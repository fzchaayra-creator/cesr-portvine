## =============================================================================
##  01_download_data.R
##
##  Downloads daily equity return data from Yahoo Finance for the CESR study.
##  Produces cesr_data_2015_2025.RData containing:
##    - px  : adjusted closing prices (xts)
##    - ret : daily log-returns (xts, in DECIMAL units)
##
##
##  Reproducibility : the exact series depend on Yahoo Finance data at the
##  time of download. All symbols are publicly listed; small differences
##  (dividend adjustments, split corrections) may appear across downloads.
## =============================================================================

library(quantmod)
library(xts)

## --- Parameters --------------------------------------------------------------
tickers <- c("IBE.MC","SAN.MC","BBVA.MC","ITX.MC","CLNX.MC",
             "AMS.MC","TEF.MC","REP.MC","FER.MC")   # 9 Spanish equities
indices <- c("^GSPC","^STOXX50E")                    # 2 stress-factor indices

date_from <- "2015-01-01"
date_to   <- "2025-04-30"

fichier_donnees <- "cesr_data_2015_2025.RData"

## --- Download ---------------------------------------------------------------
if (file.exists(fichier_donnees)) {
  cat("Data file already exists:", fichier_donnees, "\n")
  cat("Delete it and re-run this script to re-download.\n")
} else {
  cat("Downloading from Yahoo Finance...\n")
  getSymbols(c(tickers, indices), src = "yahoo",
             from = date_from, to = date_to)

  ## Combine Adjusted Close from all series (Yahoo strips the '^' prefix)
  px <- do.call(merge, lapply(c(tickers, indices),
                              function(t) Ad(get(sub("\\^", "", t)))))
  colnames(px) <- c("IBE","SAN","BBVA","ITX","CLNX","AMS","TEF","REP","FER",
                    "SP500","ESTOXX")

  ## Align on common trading days and compute log-returns in DECIMAL units
  px  <- na.omit(px)
  ret <- diff(log(px))[-1, ]

  save(px, ret, file = fichier_donnees)
  cat(sprintf("Saved: %s (%d observations, %s -> %s)\n",
              fichier_donnees, nrow(ret),
              format(index(ret)[1]),
              format(index(ret)[nrow(ret)])))
}
