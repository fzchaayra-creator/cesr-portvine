## ================================================================
##  patch_etape6.R  -  VERSION CORRIGEE (v2)
##
##  Correction du bug de la v1 : la regression poolee ne ponderait
##  pas par w_{t,m} = 1{U_{I,t} <= alpha_{I,m}}. Elle testait donc
##  E[xi_t] = 0 au lieu de E[xi_t | U_{I,t} <= alpha_{I,m}] = 0.
##
##  Ce patch v2 rend la ponderation OBLIGATOIRE si les indicateurs
##  sont fournis, et OPTIONNELLE mais avec un warning explicite si
##  aucun n'est fourni.
##
##  Convention theorique du manuscrit :
##    xi_t = e_t - q_t + (1/alpha) * 1{y_t <= q_t} * (q_t - y_t)
##    E[xi_t | U_{I,t} <= alpha_{I,m}] = 0  sous H_0
##  Estimateur pondere :
##    hat{E}[xi | U_I <= aI] = sum_t w_{t,m} xi_t / sum_t w_{t,m}
##
##  La regression poolee devient :
##    (t, m) : w_{t,m} * xi_t = w_{t,m} * beta_0 + w_{t,m} * aI_m + u_{t,m}
##  ce qui produit l'OLS pondere standard sur les seules dates
##  d'evenement. La variance est clusterisee par date.
##
##  ATTENTION D'USAGE : le xi de retour de compute_xi() garde sa
##  definition standard (non-pondere). La ponderation intervient
##  UNIQUEMENT dans les fonctions de test cesr_Z1() et cesr_pooled(),
##  qui acceptent maintenant un argument 'w' pour Z1 et une matrice
##  'W' pour Z2/Z3.
## ================================================================


## --- Residu d'identification CESR (inchange) --------------------
compute_xi <- function(var_fc, es_fc, y_real, alpha = 0.025) {
  stopifnot(length(var_fc) == length(y_real),
            length(es_fc)  == length(y_real))
  es_fc - var_fc + (1 / alpha) * (y_real <= var_fc) * (var_fc - y_real)
}


## --- Z1 : test par niveau, PONDERE ------------------------------
##  xi : vecteur des residus (longueur T_out)
##  w  : vecteur d'indicateurs 0/1 = 1{U_{I,t} <= alpha_I}
##       (si NULL, warning et retour au test non pondere de la v1)
cesr_Z1 <- function(xi, w = NULL, alpha_level = 0.05) {
  T_out <- length(xi)
  if (is.null(w)) {
    warning("cesr_Z1 : aucun indicateur d'evenement fourni. ",
            "Le test estime E[xi]=0, pas E[xi|U_I<=alpha_I]=0.",
            call. = FALSE)
    w <- rep(1, T_out)
  }
  stopifnot(length(w) == T_out, all(w %in% c(0, 1)))
  n_event <- sum(w)
  if (n_event < 5)
    return(list(spec = "Z1", K = 1, beta = NA_real_, se = NA_real_,
                T_stat = NA_real_, p_value = NA_real_, p = NA_real_,
                n_event = n_event, reject = NA))
  ## Moyenne conditionnelle sur les seules dates d'evenement
  xi_ev <- xi[w == 1]
  beta  <- mean(xi_ev)
  ## Variance Eicker-White sur les dates d'evenement uniquement
  V <- sum((xi_ev - beta)^2) / n_event^2
  if (!is.finite(V) || V <= 0)
    return(list(spec = "Z1", K = 1, beta = beta, se = NA_real_,
                T_stat = NA_real_, p_value = NA_real_, p = NA_real_,
                n_event = n_event, reject = NA))
  T_stat  <- beta^2 / V
  p_value <- pchisq(T_stat, df = 1, lower.tail = FALSE)
  list(spec = "Z1", K = 1, beta = beta, se = sqrt(V),
       T_stat = T_stat, p_value = p_value, p = p_value,
       n_event = n_event, reject = p_value < alpha_level)
}


## --- Z2 / Z3 : regression poolee PONDEREE -----------------------
##  XI  : matrice T_out x M des residus
##  aI  : vecteur des niveaux, EN POINTS DE POURCENTAGE
##  W   : matrice T_out x M des indicateurs 1{U_{I,t} <= alpha_{I,m}}
##        (si NULL, warning et retour au comportement v1 non pondere)
cesr_pooled <- function(XI, aI, W = NULL, spec = c("Z2", "Z3"),
                        alpha_level = 0.05) {
  
  spec  <- match.arg(spec)
  T_out <- nrow(XI); M <- ncol(XI)
  stopifnot(length(aI) == M)
  
  if (is.null(W)) {
    warning("cesr_pooled : aucune matrice de ponderation fournie. ",
            "Le test estime E[xi]=0 au lieu de l'estimand du manuscrit ",
            "E[xi|U_I<=alpha_I]=0.", call. = FALSE)
    W <- matrix(1, T_out, M)
  }
  stopifnot(nrow(W) == T_out, ncol(W) == M)
  
  ## Construction des regresseurs (identique a la v1)
  un <- rep(1, T_out)
  Zlist <- vector("list", M)
  for (m in seq_len(M)) {
    if (spec == "Z2") {
      Zlist[[m]] <- cbind(un, rep(aI[m], T_out))
    } else {
      xi_lag     <- c(0, XI[-T_out, m])
      Zlist[[m]] <- cbind(un, rep(aI[m], T_out),
                          rep(aI[m]^2, T_out), xi_lag)
    }
  }
  K <- ncol(Zlist[[1]])
  
  ## OLS PONDERE : chaque terme est multiplie par w_{t,m}
  ##   min sum_{t,m} w_{t,m} * (xi_{t,m} - z_{t,m} beta)^2
  ##   => A = sum_m (W_m^{1/2} Z_m)' (W_m^{1/2} Z_m) / T_out
  ##      s = sum_m (W_m^{1/2} Z_m)' (W_m^{1/2} xi_m) / T_out
  A <- matrix(0, K, K); s <- numeric(K)
  for (m in seq_len(M)) {
    Zw <- Zlist[[m]] * W[, m]                   # ponderation ligne a ligne
    A  <- A + crossprod(Zw, Zlist[[m]]) / T_out
    s  <- s + drop(crossprod(Zw, XI[, m]))  / T_out
  }
  beta <- tryCatch(drop(solve(A, s)), error = function(e) rep(NA_real_, K))
  if (anyNA(beta))
    return(list(spec = spec, K = K, beta = beta, se = rep(NA_real_, K),
                T_stat = NA_real_, p_value = NA_real_, p = NA_real_,
                reject = NA))
  
  ## Score agrege PAR DATE, avec ponderation :
  ##   g_t = sum_m w_{t,m} * z_{t,m} * (xi_{t,m} - z_{t,m}' beta)
  ## Variance clusterisee par date (chaque date est un cluster).
  g <- matrix(0, T_out, K)
  for (m in seq_len(M)) {
    u <- XI[, m] - drop(Zlist[[m]] %*% beta)
    g <- g + Zlist[[m]] * (W[, m] * u)          # ponderation dans le score
  }
  B  <- crossprod(g) / T_out
  Ai <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(Ai))
    return(list(spec = spec, K = K, beta = beta, se = rep(NA_real_, K),
                T_stat = NA_real_, p_value = NA_real_, p = NA_real_,
                reject = NA))
  V <- Ai %*% B %*% Ai
  
  T_stat  <- tryCatch(as.numeric(T_out * t(beta) %*% solve(V, beta)),
                      error = function(e) NA_real_)
  p_value <- if (is.na(T_stat)) NA_real_ else
    pchisq(T_stat, df = K, lower.tail = FALSE)
  
  list(spec = spec, K = K, beta = beta,
       se = suppressWarnings(sqrt(diag(V) / T_out)),
       T_stat = T_stat, p_value = p_value, p = p_value,
       reject = !is.na(p_value) && p_value < alpha_level)
}


## --- Application sur la grille (avec ponderation obligatoire) --
##  fc     : liste avec $VaR et $ES (matrices T_out x M)
##  y_real : rendements portefeuille out-of-sample (longueur T_out)
##  u_I    : PIT du facteur de stress (longueur T_out) - OBLIGATOIRE
##  alpha_I_grid : niveaux (proba, ex. c(0.05, 0.10, 0.25, 0.50))
run_cesr_grid <- function(fc, y_real, u_I, alpha_I_grid,
                          specs = c("Z1", "Z2", "Z3"),
                          alpha_ES = 0.025,
                          alpha_I_en_points = TRUE) {
  
  stopifnot(length(y_real) == length(u_I),
            length(y_real) == nrow(fc$VaR))
  
  M  <- length(alpha_I_grid)
  XI <- sapply(seq_len(M), function(j)
    compute_xi(fc$VaR[, j], fc$ES[, j], y_real, alpha = alpha_ES))
  XI <- matrix(XI, ncol = M)
  colnames(XI) <- paste0("aI=", alpha_I_grid)
  
  ## Matrice de ponderation par niveau
  W <- sapply(seq_len(M), function(j) as.numeric(u_I <= alpha_I_grid[j]))
  W <- matrix(W, ncol = M)
  colnames(W) <- paste0("aI=", alpha_I_grid)
  
  aI <- if (alpha_I_en_points) 100 * alpha_I_grid else alpha_I_grid
  
  results <- list()
  
  ## Z1 : un test par niveau, sur les seules dates d'evenement
  if ("Z1" %in% specs)
    for (j in seq_len(M))
      results[[paste("Z1", alpha_I_grid[j], sep = "_")]] <-
    c(cesr_Z1(XI[, j], w = W[, j]),
      list(alpha_I = alpha_I_grid[j]))
  
  ## Z2 / Z3 : un test poole, DUPLIQUE sur les cles de niveau
  ## pour retro-compatibilite avec les scripts existants
  for (sp in intersect(c("Z2", "Z3"), specs)) {
    res <- cesr_pooled(XI, aI, W = W, spec = sp)
    for (j in seq_len(M))
      results[[paste(sp, alpha_I_grid[j], sep = "_")]] <-
        c(res, list(alpha_I = alpha_I_grid[j], pooled = TRUE))
  }
  
  ## Bonferroni sur la famille Z1
  p_Z1 <- vapply(seq_len(M), function(j)
    results[[paste("Z1", alpha_I_grid[j], sep = "_")]]$p_value, numeric(1))
  p_global <- min(1, M * min(p_Z1, na.rm = TRUE))
  
  ## Diagnostic : nombre d'exceedances CONDITIONNELLES par niveau
  n_exc <- vapply(seq_len(M), function(j) {
    sum(W[, j] == 1 & y_real <= fc$VaR[, j])
  }, numeric(1))
  n_event <- colSums(W)
  diagnostics <- data.frame(
    alpha_I         = alpha_I_grid,
    n_event         = n_event,
    n_exceedances   = n_exc,
    exceed_rate_pct = round(100 * n_exc / n_event, 2),
    expected_pct    = 100 * alpha_ES,
    row.names = NULL)
  
  list(results = results, p_global = p_global,
       xi = XI, W = W,
       diagnostics = diagnostics,
       beta_Z2 = if (!is.null(results[[paste("Z2", alpha_I_grid[1], sep = "_")]]))
         results[[paste("Z2", alpha_I_grid[1], sep = "_")]]$beta else NULL)
}


## --- Execution ---------------------------------------------------
.y_oos <- NULL
for (nm in c("port_oos", "port_returns_oos", "portfolio_returns_oos",
             "port_ret_oos", "y_oos")) {
  if (exists(nm, envir = globalenv())) { .y_oos <- get(nm, envir = globalenv())
  .y_nom <- nm; break }
}
if (is.null(.y_oos))
  stop("Rendements out-of-sample introuvables.", call. = FALSE)
cat(sprintf("Rendements out-of-sample : '%s'  (%d obs, plage %.3f a %.3f)\n",
            .y_nom, length(.y_oos), min(.y_oos), max(.y_oos)))

if (max(abs(.y_oos)) < 1 && max(abs(fc_sp500$VaR)) > 1)
  warning("Les rendements semblent en DECIMALES et les previsions en POURCENTAGE. ",
          "Le residu serait faux d'un facteur 100.", call. = FALSE)

cat("\n--- Test CESR : S&P 500 (avec ponderation event) ---\n")
cesr_sp500  <- run_cesr_grid(fc_sp500,  .y_oos, u_I_sp, alpha_I_grid)
print(cesr_sp500$diagnostics)

cat("\n--- Test CESR : Euro Stoxx 50 (avec ponderation event) ---\n")
cesr_estoxx <- run_cesr_grid(fc_estoxx, .y_oos, u_I_es, alpha_I_grid)
print(cesr_estoxx$diagnostics)


## --- Recapitulatif des p-values ----------------------------------
recap <- data.frame(
  spec = c(paste0("Z1_", alpha_I_grid), "Z2 (poole)", "Z3 (poole)", "Bonferroni(Z1)"),
  SP500 = c(sapply(alpha_I_grid, function(a) cesr_sp500$results[[paste0("Z1_", a)]]$p_value),
            cesr_sp500$results[[paste0("Z2_", alpha_I_grid[1])]]$p_value,
            cesr_sp500$results[[paste0("Z3_", alpha_I_grid[1])]]$p_value,
            cesr_sp500$p_global),
  ESTOXX = c(sapply(alpha_I_grid, function(a) cesr_estoxx$results[[paste0("Z1_", a)]]$p_value),
             cesr_estoxx$results[[paste0("Z2_", alpha_I_grid[1])]]$p_value,
             cesr_estoxx$results[[paste0("Z3_", alpha_I_grid[1])]]$p_value,
             cesr_estoxx$p_global))

cat("\n=== NOUVELLES p-VALUES CESR (patch v2 pondere) ===\n")
print(recap, row.names = FALSE, digits = 4)

cat(sprintf("\nbeta_Z2  S&P 500       : b0 = %+.4f   b1 = %+.4f\n",
            cesr_sp500$beta_Z2[1],  cesr_sp500$beta_Z2[2]))
cat(sprintf("beta_Z2  Euro Stoxx 50 : b0 = %+.4f   b1 = %+.4f",
            cesr_estoxx$beta_Z2[1], cesr_estoxx$beta_Z2[2]))
cat("   [points de rendement / point de alpha_I]\n")

## Moyenne des residus CONDITIONNELS (sur les dates d'evenement seulement)
cat("\nMoyenne des residus PONDERES par niveau (positif = ES sous-estime) :\n")
res_cond <- rbind(
  `S&P 500` = sapply(seq_along(alpha_I_grid), function(m) {
    xi <- cesr_sp500$xi[, m][cesr_sp500$W[, m] == 1]
    if (length(xi) < 1) NA else mean(xi)
  }),
  `Euro Stoxx 50` = sapply(seq_along(alpha_I_grid), function(m) {
    xi <- cesr_estoxx$xi[, m][cesr_estoxx$W[, m] == 1]
    if (length(xi) < 1) NA else mean(xi)
  })
)
colnames(res_cond) <- paste0("aI=", alpha_I_grid)
print(round(res_cond, 4))
