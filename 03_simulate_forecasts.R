## =============================================================================
##  03_simulate_forecasts.R
##
##  Generates conditional VaR/ES forecasts under the event-based estimand
##  U_I <= alpha_I, using the fitted portvine objects from step 02.
##
##  Rationale : portvine's estimate_risk_roll() implements POINTWISE
##  conditioning (u_I = alpha_I fixed). For the event-based estimand of
##  the manuscript we need u_I ~ Uniform(0, alpha_I). We reuse the fitted
##  marginals and D-vines from step 02 and only re-do the simulation
##  layer here.
##
##  As a by-product, this script also computes the realised PIT of the
##  stress factor u_I,t, required by the CESR event-weight indicators.
##
##  Produces cesr_forecasts_event_2016_2025.RData containing:
##    - fc_sp500, fc_estoxx : lists with $VaR, $ES matrices (T_out x 4)
##    - u_I_sp, u_I_es       : realised PIT of each stress factor
##    - port_oos              : out-of-sample portfolio returns
##    - dates_oos             : out-of-sample dates
##    - alpha_ES, alpha_I_grid, T_out, weights, and window sizes
##
##  Runtime : ~40 minutes. Checkpoint-restartable.
## =============================================================================

library(portvine)
library(rugarch)
library(rvinecopulib)
library(xts)

set.seed(20260728)

## --- Load previous stages ---------------------------------------------------
if (!file.exists("cesr_portvine_2016_2025_fits.RData"))
  stop("Fits not found. Run 02_estimate_portvine.R first.")
load("cesr_data_2015_2025.RData")
load("cesr_portvine_2016_2025_fits.RData")

## --- Parameters -------------------------------------------------------------
noms_actifs  <- c("IBE","SAN","BBVA","ITX","CLNX","AMS","TEF","REP","FER")
weights      <- setNames(rep(1/9, 9), noms_actifs)
alpha_I_grid <- c(0.05, 0.10, 0.25, 0.50)
n_samples    <- 1000    # conditional draws per (date, level)

fichier_final <- "cesr_forecasts_event_2016_2025.RData"

dates_all <- index(ret)
T_total   <- nrow(ret)
port_ret  <- as.numeric(as.matrix(ret[, noms_actifs]) %*% weights)

## --- Extract fitted marginals and vines from portvine ----------------------
extract_marginals <- function(pv, stress_name) {
  vars <- c(noms_actifs, stress_name)
  out  <- list()
  for (v in vars) {
    roll <- fitted_marginals(pv)[[v]]
    df   <- as.data.frame(roll, which = "density")  # Mu, Sigma, Skew, Shape, Realized
    out[[v]] <- df
  }
  out
}
marg_sp <- extract_marginals(pv_sp, "SP500")
marg_es <- extract_marginals(pv_es, "ESTOXX")

## PIT of the stress factor : u_I,t = F_skewt( (r_t - mu_t)/sigma_t )
pit_of <- function(df) {
  z <- (df$Realized - df$Mu) / df$Sigma
  rugarch::pdist("sstd", q = z, mu = 0, sigma = 1,
                 skew = df$Skew, shape = df$Shape)
}
u_I_sp <- pit_of(marg_sp[["SP500"]])
u_I_es <- pit_of(marg_es[["ESTOXX"]])
T_out  <- length(u_I_sp)

cat(sprintf("T_out = %d out-of-sample dates\n", T_out))
cat(sprintf("KS uniformity of u_I : SP500 p = %.3f | ESTOXX p = %.3f\n",
            ks.test(u_I_sp, "punif")$p.value,
            ks.test(u_I_es, "punif")$p.value))

vines_sp <- fitted_vines(pv_sp)
vines_es <- fitted_vines(pv_es)
vine_index_for_t <- function(t) 1L + (t - 1L) %/% refit_vine

## --- Event-conditional simulator -------------------------------------------
##  For each date t and stress level aI :
##     u_I^(s) ~ Uniform(0, aI)                            <- the change
##     inverse Rosenblatt of the D-vine conditioned on u_I^(s)
##     return-scale mapping via each asset's marginal (skew-t)
##     portfolio aggregation and empirical VaR / ES
simulate_event <- function(vine, marg_list, stress_name, t, aI, S) {
  d  <- length(noms_actifs)
  W  <- matrix(runif(S * d), ncol = d)
  uI <- runif(S, 0, aI)
  U  <- cbind(W, uI)
  Usim <- rvinecopulib::inverse_rosenblatt(U, vine)
  r <- matrix(NA_real_, S, d)
  for (j in seq_len(d)) {
    dfj  <- marg_list[[noms_actifs[j]]]
    zsim <- rugarch::qdist("sstd", p = Usim[, j], mu = 0, sigma = 1,
                            skew = dfj$Skew[t], shape = dfj$Shape[t])
    r[, j] <- dfj$Mu[t] + dfj$Sigma[t] * zsim   # in PERCENT (fitted x100)
  }
  drop(r %*% weights) / 100                      # back to DECIMAL
}

regen_factor <- function(vines, marg_list, stress_name, tag) {
  M   <- length(alpha_I_grid)
  VaR <- ES <- matrix(NA_real_, T_out, M,
                      dimnames = list(NULL, paste0("aI_", alpha_I_grid)))
  ckpt <- sprintf("ckpt_%s.RData", tag)
  t0   <- 1L
  if (file.exists(ckpt)) {
    load(ckpt); t0 <- t_done + 1L
    cat(sprintf("  [%s] resume at t = %d\n", tag, t0))
  }
  for (t in t0:T_out) {
    vi <- vine_index_for_t(t)
    for (m in seq_len(M)) {
      draws     <- simulate_event(vines[[vi]], marg_list, stress_name,
                                  t, alpha_I_grid[m], n_samples)
      q         <- quantile(draws, probs = alpha_ES, type = 1, names = FALSE)
      VaR[t, m] <- q
      ES[t, m]  <- mean(draws[draws <= q])
    }
    if (t %% 100 == 0) {
      t_done <- t; save(VaR, ES, t_done, file = ckpt)
      cat(sprintf("  [%s] %d / %d  (%s)\n", tag, t, T_out,
                  format(Sys.time(), "%H:%M")))
    }
  }
  file.remove(ckpt)
  list(VaR = VaR, ES = ES)
}

cat("\n=== Event-based re-simulation : S&P 500 ===\n")
fc_sp500  <- regen_factor(vines_sp, marg_sp, "SP500",  "sp500")
cat("\n=== Event-based re-simulation : Euro Stoxx 50 ===\n")
fc_estoxx <- regen_factor(vines_es, marg_es, "ESTOXX", "estoxx")

## --- Align out-of-sample series and save -----------------------------------
port_oos  <- tail(port_ret,  T_out)
dates_oos <- tail(dates_all, T_out)

stopifnot(all(fc_sp500$ES  <= fc_sp500$VaR  + 1e-12),
          all(fc_estoxx$ES <= fc_estoxx$VaR + 1e-12),
          length(port_oos) == T_out)

save(fc_sp500, fc_estoxx, port_oos, dates_oos,
     u_I_sp, u_I_es, alpha_ES, alpha_I_grid,
     T_in_marg, refit_marg, T_in_vine, refit_vine, T_out, weights,
     file = fichier_final)
cat("\n=== DONE. Saved to", fichier_final, "===\n")
