## =============================================================================
##  06_figures.R
##
##  Regenerates the five figures used in the manuscript from the CESR test
##  results of step 04 :
##    fig1_residuals_sp500_HQ.pdf  : event-date CESR residuals for S&P 500
##    fig3_pvalues_HQ.pdf          : p-values across specifications (bar chart)
##    fig4_heatmap_HQ.pdf          : full p-value heatmap (used in appendix)
##    fig5_calibration_HQ.pdf      : VaR / ES calibration diagnostic
##    fig7_exceedance_HQ.pdf       : conditional exceedance rates
##
## =============================================================================

setwd("~/reel_data")

if (!file.exists("cesr_forecasts_event_2016_2025.RData"))
  stop("Forecasts not found. Run 03_simulate_forecasts.R first.")

load("cesr_forecasts_event_2016_2025.RData")
source("patch_etape6.R")  # produces cesr_sp500, cesr_estoxx

library(ggplot2)

FONT_BASE <- "sans"
theme_paper <- theme_minimal(base_size = 11, base_family = FONT_BASE) +
  theme(panel.grid.minor = element_blank(),
        axis.title       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(colour = "grey30"))

## --- Reference p-values for sanity check ------------------------------------
cat("\n=== p-value reference (should match Table 4 of the manuscript) ===\n")
cat(sprintf("  Z2 pooled S&P 500       : %.4f\n",
            cesr_sp500$results$Z2_0.05$p_value))
cat(sprintf("  Z2 pooled Euro Stoxx 50 : %.4f\n",
            cesr_estoxx$results$Z2_0.05$p_value))

## ============================================================================
## FIG 1 : CESR residuals under S&P 500 stress conditioning
## ============================================================================
xi_sp   <- cesr_sp500$xi
W_sp    <- sapply(seq_along(alpha_I_grid),
                  function(j) as.numeric(u_I_sp <= alpha_I_grid[j]))
df1 <- do.call(rbind, lapply(seq_along(alpha_I_grid), function(j) {
  idx <- which(W_sp[, j] == 1)
  data.frame(date  = seq_along(port_oos)[idx],
             xi    = xi_sp[idx, j],
             level = factor(sprintf("alpha_I = %.2f", alpha_I_grid[j]),
                            levels = sprintf("alpha_I = %.2f", alpha_I_grid)))
}))
means_sp <- data.frame(
  level = factor(sprintf("alpha_I = %.2f", alpha_I_grid),
                 levels = sprintf("alpha_I = %.2f", alpha_I_grid)),
  mean_xi = sapply(seq_along(alpha_I_grid),
                   function(j) mean(xi_sp[W_sp[, j] == 1, j])))
fig1 <- ggplot(df1, aes(date, xi)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_point(alpha = 0.4, size = 0.7, colour = "#2E5C8A") +
  geom_hline(data = means_sp, aes(yintercept = mean_xi),
             colour = "#B22222", linewidth = 0.7, linetype = "dashed") +
  facet_wrap(~ level, ncol = 2, scales = "free_y") +
  labs(title = "CESR residuals under S&P 500 stress conditioning",
       subtitle = "Blue: event-date residuals; red dashed: weighted mean",
       x = "Event date index",
       y = expression(hat(xi)[t](alpha[I]))) +
  theme_paper
ggsave("fig1_residuals_sp500_HQ.pdf", fig1,
       width = 9, height = 6, device = cairo_pdf)
cat("Fig 1 : OK\n")

## ============================================================================
## FIG 3 : p-values on log scale, by specification and stress factor
## ============================================================================
p_data <- data.frame(
  stress = rep(c("S&P 500","Euro Stoxx 50"), each = 6),
  spec   = rep(c(sprintf("Z1 (%.2f)", alpha_I_grid),
                 "Z2 pooled","Z3 pooled"), 2),
  p      = c(
    sapply(alpha_I_grid,
           function(a) cesr_sp500$results[[paste0("Z1_", a)]]$p_value),
    cesr_sp500$results$Z2_0.05$p_value,
    cesr_sp500$results$Z3_0.05$p_value,
    sapply(alpha_I_grid,
           function(a) cesr_estoxx$results[[paste0("Z1_", a)]]$p_value),
    cesr_estoxx$results$Z2_0.05$p_value,
    cesr_estoxx$results$Z3_0.05$p_value))
p_data$spec   <- factor(p_data$spec, levels = c(sprintf("Z1 (%.2f)", alpha_I_grid),
                                                 "Z2 pooled","Z3 pooled"))
p_data$stress <- factor(p_data$stress, levels = c("S&P 500","Euro Stoxx 50"))
p_data$reject <- p_data$p < 0.05
fig3 <- ggplot(p_data, aes(spec, p, fill = reject)) +
  geom_col(colour = "grey30", width = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#B22222") +
  facet_wrap(~ stress, nrow = 1) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.001)) +
  scale_fill_manual(values = c("TRUE" = "#2E5C8A","FALSE" = "grey70"),
                    labels = c("TRUE" = "Rejects at 5%",
                               "FALSE" = "Does not reject")) +
  labs(title = "CESR test p-values by specification and stress factor",
       subtitle = "Red dashed line: 5% significance threshold (log scale)",
       x = NULL, y = "p-value (log scale)", fill = NULL) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")
ggsave("fig3_pvalues_HQ.pdf", fig3,
       width = 10, height = 5.5, device = cairo_pdf)
cat("Fig 3 : OK\n")

## ============================================================================
## FIG 4 : full p-value heatmap (used in the appendix)
## ============================================================================
hm_df <- data.frame(
  ai     = factor(rep(rep(sprintf("%.2f", alpha_I_grid), each = 3), 2),
                  levels = sprintf("%.2f", alpha_I_grid)),
  spec   = factor(rep(rep(c("Z1","Z2","Z3"), 4), 2), levels = c("Z1","Z2","Z3")),
  stress = factor(rep(c("S&P 500","Euro Stoxx 50"), each = 12),
                  levels = c("S&P 500","Euro Stoxx 50")),
  pv     = c(
    unlist(lapply(alpha_I_grid, function(a)
      c(cesr_sp500$results[[paste0("Z1_", a)]]$p_value,
        cesr_sp500$results[[paste0("Z2_", a)]]$p_value,
        cesr_sp500$results[[paste0("Z3_", a)]]$p_value))),
    unlist(lapply(alpha_I_grid, function(a)
      c(cesr_estoxx$results[[paste0("Z1_", a)]]$p_value,
        cesr_estoxx$results[[paste0("Z2_", a)]]$p_value,
        cesr_estoxx$results[[paste0("Z3_", a)]]$p_value)))))
hm_df$lab <- ifelse(hm_df$pv < 0.001, "<0.001", sprintf("%.3f", hm_df$pv))
fig4 <- ggplot(hm_df, aes(ai, spec, fill = pv)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = lab, colour = pv < 0.05),
            size = 3.5, fontface = "bold", family = FONT_BASE) +
  scale_fill_gradientn(
    colours = c("#08306B","#4292C6","#F7FBFF","#FDD49E","#FDAE61"),
    values  = c(0, 0.02, 0.05, 0.20, 1),
    limits  = c(0, 1), name = "p-value") +
  scale_colour_manual(values = c("TRUE" = "white","FALSE" = "black"),
                      guide = "none") +
  facet_wrap(~ stress, nrow = 1) +
  labs(title = "CESR p-value heatmap",
       subtitle = "Rows: test specification; columns: stress level alpha_I",
       x = expression(alpha[I]), y = NULL) +
  theme_paper +
  theme(panel.border = element_rect(colour = "grey80", fill = NA,
                                    linewidth = 0.5))
ggsave("fig4_heatmap_HQ.pdf", fig4,
       width = 10, height = 4.5, device = cairo_pdf)
cat("Fig 4 : OK\n")

## ============================================================================
## FIG 5, FIG 7 : reuse the existing high-quality generators
## ============================================================================
##  Figures 5 (calibration) and 7 (exceedance rates) do not depend on the
##  p-value calculation; they depend only on port_oos, fc_sp500, fc_estoxx.
##  If your local install already contains the up-to-date PDF files, no
##  regeneration is needed. If not, run the legacy generator :
##     source("cesr_figures_event_2016_2025.R")
##  which produces fig5_calibration_HQ.pdf and fig7_exceedance_HQ.pdf.
if (!file.exists("fig5_calibration_HQ.pdf") ||
    !file.exists("fig7_exceedance_HQ.pdf")) {
  cat("\nNote : fig5_calibration and fig7_exceedance are not regenerated ",
      "by this script.\n",
      "      Use the legacy generator or open the .tex documentation.\n",
      sep = "")
}

cat("\n=== DONE. Figures regenerated ===\n")
