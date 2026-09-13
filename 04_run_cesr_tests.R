## =============================================================================
##  04_run_cesr_tests.R
##
##  Applies the CESR test (weighted, event-based) to the conditional VaR/ES
##  forecasts from step 03. Reproduces Table 4 of the manuscript.
##
##  All the test logic (Z1, Z2, Z3 with event weighting, clustered variance)
##  lives in patch_etape6.R, which this script sources.
##
##  Produces cesr_main_results.RData containing:
##    - cesr_sp500, cesr_estoxx : full test outputs for both stress factors
##    - recap                    : summary table of p-values (Table 4 content)
##
## =============================================================================

setwd("~/reel_data")

if (!file.exists("cesr_forecasts_event_2016_2025.RData"))
  stop("Forecasts not found. Run 03_simulate_forecasts.R first.")

load("cesr_forecasts_event_2016_2025.RData")
source("patch_etape6.R")   # patch_etape6.R computes cesr_sp500 and cesr_estoxx
                            # via run_cesr_grid, prints the recap table, and
                            # exports the p-value summary to the console.

## Save results for downstream scripts (05, 06)
recap_table <- data.frame(
  spec  = c(paste0("Z1_", alpha_I_grid),
            "Z2_pooled", "Z3_pooled", "Bonferroni_Z1"),
  SP500 = c(sapply(alpha_I_grid,
                   function(a) cesr_sp500$results[[paste0("Z1_", a)]]$p_value),
            cesr_sp500$results[[paste0("Z2_", alpha_I_grid[1])]]$p_value,
            cesr_sp500$results[[paste0("Z3_", alpha_I_grid[1])]]$p_value,
            cesr_sp500$p_global),
  ESTOXX = c(sapply(alpha_I_grid,
                    function(a) cesr_estoxx$results[[paste0("Z1_", a)]]$p_value),
             cesr_estoxx$results[[paste0("Z2_", alpha_I_grid[1])]]$p_value,
             cesr_estoxx$results[[paste0("Z3_", alpha_I_grid[1])]]$p_value,
             cesr_estoxx$p_global)
)

save(cesr_sp500, cesr_estoxx, recap_table,
     file = "cesr_main_results.RData")

cat("\n=== DONE. Main results saved to cesr_main_results.RData ===\n")
cat("These p-values reproduce Table 4 of the manuscript.\n")
