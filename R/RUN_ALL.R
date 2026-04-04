################################################################################
# FCI MASTER SCRIPT - Run All Analyses
################################################################################
#
# This script runs the complete FCI analysis pipeline in correct order.
#
# Usage:
#   cd R/
#   Rscript RUN_ALL.R
#
################################################################################

cat("\n")
cat("================================================================================\n")
cat("  FINANCIAL CONDITIONS INDEX - COMPLETE ANALYSIS\n")
cat("  Banco Central del Paraguay\n")
cat("================================================================================\n\n")

start_time <- Sys.time()

# Check required packages
required_packages <- c("readxl", "dplyr", "tidyr", "zoo", "FactoMineR",
                       "vars", "MARSS", "ggplot2", "gridExtra", "lmtest", "sandwich",
                       "tseries", "urca", "quantreg", "ivreg")  # ivreg for IV-LP

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# Check data file
data_file <- "../data/FCI_data_1.xlsx"
if (!file.exists(data_file)) {
  stop("ERROR: Data file not found at ", data_file)
}
cat("Data file found:", data_file, "\n\n")

# Create output directory
if (!dir.exists("../output")) {
  dir.create("../output", recursive = TRUE)
  cat("Created output directory\n\n")
}

#-------------------------------------------------------------------------------
# STEP 1: Complete FCI Calculation (Consolidated Script)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 1: Complete FCI Calculation (01_FCI_Complete.R)\n")
cat("        - Comprehensive FCI (all variables)\n")
cat("        - Level 2: Endogenous/Exogenous decomposition\n")
cat("        - Level 3: Rates/Banking/External channels\n")
cat("        - Orthogonalization\n")
cat("================================================================================\n\n")

tryCatch({
  source("01_FCI_Complete.R")
  cat("\n[OK] FCI calculation complete\n\n")
}, error = function(e) {
  cat("\n[ERROR] FCI calculation failed:", e$message, "\n")
  stop("Cannot continue without FCI calculation")
})

#-------------------------------------------------------------------------------
# STEP 1B: Stationarity Tests (Phase 1.4)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 1B: Stationarity Tests (03_FCI_Stationarity_Tests.R)\n")
cat("================================================================================\n\n")

if (file.exists("03_FCI_Stationarity_Tests.R")) {
  tryCatch({
    source("03_FCI_Stationarity_Tests.R")
    cat("\n[OK] Stationarity tests complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Stationarity tests failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Stationarity tests script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 2: Effects Analysis
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 2: Effects Analysis (02_FCI_Effects_Analysis.R)\n")
cat("================================================================================\n\n")

if (file.exists("02_FCI_Effects_Analysis.R")) {
  tryCatch({
    source("02_FCI_Effects_Analysis.R")
    cat("\n[OK] Effects analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Effects analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Effects analysis script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 3: Local Projections Analysis
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 3: Local Projections (05_FCI_Local_Projections.R)\n")
cat("================================================================================\n\n")

if (file.exists("05_FCI_Local_Projections.R")) {
  tryCatch({
    source("05_FCI_Local_Projections.R")
    cat("\n[OK] Local Projections complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Local Projections failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Local Projections script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 4: Rolling/Expanding Stability Analysis
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 4: Rolling/Expanding Stability (06_FCI_Rolling_Stability.R)\n")
cat("================================================================================\n\n")

if (file.exists("06_FCI_Rolling_Stability.R")) {
  tryCatch({
    source("06_FCI_Rolling_Stability.R")
    cat("\n[OK] Rolling Stability analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Rolling Stability analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Rolling Stability script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 5: Robustness/Endogeneity Analysis
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 5: Robustness/Endogeneity (07_FCI_Robustness_Endogeneity.R)\n")
cat("================================================================================\n\n")

if (file.exists("07_FCI_Robustness_Endogeneity.R")) {
  tryCatch({
    source("07_FCI_Robustness_Endogeneity.R")
    cat("\n[OK] Robustness analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Robustness analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Robustness script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 6: Regime and TVP Analysis
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 6: Regime and TVP Analysis (08_FCI_Regime_TVP_Analysis.R)\n")
cat("================================================================================\n\n")

if (file.exists("08_FCI_Regime_TVP_Analysis.R")) {
  tryCatch({
    source("08_FCI_Regime_TVP_Analysis.R")
    cat("\n[OK] Regime and TVP analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Regime and TVP analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Regime and TVP script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 7: Monetary Policy Interaction and Forecasting
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 7: Monetary Policy Interaction (09_FCI_Monetary_Policy_Interaction.R)\n")
cat("================================================================================\n\n")

if (file.exists("09_FCI_Monetary_Policy_Interaction.R")) {
  tryCatch({
    source("09_FCI_Monetary_Policy_Interaction.R")
    cat("\n[OK] Monetary Policy Interaction analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Monetary Policy Interaction analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Monetary Policy Interaction script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 8: Growth-at-Risk Analysis (Phase 3.1)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 8: Growth-at-Risk (10_FCI_Growth_at_Risk.R)\n")
cat("        - Quantile regressions for tail risk analysis\n")
cat("        - Tests asymmetric effects on growth distribution\n")
cat("================================================================================\n\n")

if (file.exists("10_FCI_Growth_at_Risk.R")) {
  tryCatch({
    source("10_FCI_Growth_at_Risk.R")
    cat("\n[OK] Growth-at-Risk analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Growth-at-Risk analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Growth-at-Risk script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 9: Block-Exogenous SVAR (Phase 3.2)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 9: Block-Exogenous SVAR (11_FCI_Block_SVAR.R)\n")
cat("        - External spillover analysis\n")
cat("        - Structural IRFs with bootstrap CIs\n")
cat("================================================================================\n\n")

if (file.exists("11_FCI_Block_SVAR.R")) {
  tryCatch({
    source("11_FCI_Block_SVAR.R")
    cat("\n[OK] Block-SVAR analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Block-SVAR analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Block-SVAR script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 10: Output Puzzle Investigation (Phase 3.3)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 10: Output Puzzle Investigation (12_FCI_Output_Puzzle_Investigation.R)\n")
cat("         - Local Projections with IMAEP_SANB (non-agricultural GDP)\n")
cat("         - VECM cointegration analysis (credit-output long-run)\n")
cat("         - Investment/Consumption transmission channels\n")
cat("================================================================================\n\n")

if (file.exists("12_FCI_Output_Puzzle_Investigation.R")) {
  tryCatch({
    source("12_FCI_Output_Puzzle_Investigation.R")
    cat("\n[OK] Output Puzzle Investigation complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Output Puzzle Investigation failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Output Puzzle Investigation script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 10B: Output Puzzle - Sectoral Analysis
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 10B: Output Puzzle Sectoral Analysis (14_FCI_Output_Puzzle_Sectoral.R)\n")
cat("          - Sectoral credit/output analysis\n")
cat("================================================================================\n\n")

if (file.exists("14_FCI_Output_Puzzle_Sectoral.R")) {
  tryCatch({
    source("14_FCI_Output_Puzzle_Sectoral.R")
    cat("\n[OK] Output Puzzle Sectoral Analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Output Puzzle Sectoral Analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Output Puzzle Sectoral Analysis script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 11: Post-IT Subsample Analysis
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 11: Post-IT Subsample Analysis (15_FCI_PostIT_Subsample_Analysis.R)\n")
cat("================================================================================\n\n")

if (file.exists("15_FCI_PostIT_Subsample_Analysis.R")) {
  tryCatch({
    source("15_FCI_PostIT_Subsample_Analysis.R")
    cat("\n[OK] Post-IT Subsample Analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Post-IT Subsample Analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Post-IT Subsample Analysis script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 12: TCN Reclassification Sensitivity
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 12: TCN Reclassification Sensitivity (16_FCI_TCN_Reclassification.R)\n")
cat("================================================================================\n\n")

if (file.exists("16_FCI_TCN_Reclassification.R")) {
  tryCatch({
    source("16_FCI_TCN_Reclassification.R")
    cat("\n[OK] TCN Reclassification Sensitivity complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] TCN Reclassification Sensitivity failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] TCN Reclassification Sensitivity script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 13: New External Data Construction (Phase 5)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 13: New External Data Construction (17_FCI_New_External_Data.R)\n")
cat("         - GFC index, regional proxy, disaggregated commodities\n")
cat("================================================================================\n\n")

if (file.exists("17_FCI_New_External_Data.R")) {
  tryCatch({
    source("17_FCI_New_External_Data.R")
    cat("\n[OK] New External Data Construction complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] New External Data Construction failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] New External Data script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 14: Improved IV-LP (Phase 5)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 14: Improved IV-LP (18_FCI_Improved_IV_LP.R)\n")
cat("         - Expanded instruments, weak-IV robust inference\n")
cat("================================================================================\n\n")

if (file.exists("18_FCI_Improved_IV_LP.R")) {
  tryCatch({
    source("18_FCI_Improved_IV_LP.R")
    cat("\n[OK] Improved IV-LP complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Improved IV-LP failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Improved IV-LP script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 15: Proxy-SVAR (Phase 5)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 15: Proxy-SVAR (19_FCI_Proxy_SVAR.R)\n")
cat("         - Mertens & Ravn (2013) external instruments SVAR\n")
cat("================================================================================\n\n")

if (file.exists("19_FCI_Proxy_SVAR.R")) {
  tryCatch({
    source("19_FCI_Proxy_SVAR.R")
    cat("\n[OK] Proxy-SVAR complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Proxy-SVAR failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Proxy-SVAR script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 16: Enriched Block-SVAR (Phase 5)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 16: Enriched Block-SVAR (20_FCI_Enriched_Block_SVAR.R)\n")
cat("         - 7-var external block, historical decomposition\n")
cat("================================================================================\n\n")

if (file.exists("20_FCI_Enriched_Block_SVAR.R")) {
  tryCatch({
    source("20_FCI_Enriched_Block_SVAR.R")
    cat("\n[OK] Enriched Block-SVAR complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Enriched Block-SVAR failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Enriched Block-SVAR script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 17: Commodity Puzzle Investigation (Phase 5)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 17: Commodity Puzzle Investigation (21_FCI_Commodity_Puzzle.R)\n")
cat("         - Disaggregated commodities, ToT, offsetting effects\n")
cat("================================================================================\n\n")

if (file.exists("21_FCI_Commodity_Puzzle.R")) {
  tryCatch({
    source("21_FCI_Commodity_Puzzle.R")
    cat("\n[OK] Commodity Puzzle Investigation complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Commodity Puzzle Investigation failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Commodity Puzzle script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 18: Regional Spillover Analysis (Phase 5)
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 18: Regional Spillover Analysis (22_FCI_Regional_Spillovers.R)\n")
cat("         - Brazil Selic spillovers, contagion vs common shocks\n")
cat("================================================================================\n\n")

if (file.exists("22_FCI_Regional_Spillovers.R")) {
  tryCatch({
    source("22_FCI_Regional_Spillovers.R")
    cat("\n[OK] Regional Spillover Analysis complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Regional Spillover Analysis failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Regional Spillover script not found\n\n")
}

#-------------------------------------------------------------------------------
# STEP 19: Identification Triangulation
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("STEP 19: Identification Triangulation (23_FCI_Identification_Triangulation.R)\n")
cat("         - DXY-only IV-LP, GFC-PCA IV, globally-purged FCI, triangulation figure\n")
cat("================================================================================\n\n")

if (file.exists("23_FCI_Identification_Triangulation.R")) {
  tryCatch({
    source("23_FCI_Identification_Triangulation.R")
    cat("\n[OK] Identification Triangulation complete\n\n")
  }, error = function(e) {
    cat("\n[WARNING] Identification Triangulation failed:", e$message, "\n\n")
  })
} else {
  cat("[SKIP] Identification Triangulation script not found\n\n")
}

#-------------------------------------------------------------------------------
# FINAL STEP: Organize output files by type
#-------------------------------------------------------------------------------
cat("================================================================================\n")
cat("FINAL STEP: Organizing output files by type\n")
cat("================================================================================\n\n")

out <- "../output"
dir.create(file.path(out, "png"), showWarnings = FALSE)
dir.create(file.path(out, "csv"), showWarnings = FALSE)
dir.create(file.path(out, "pdf"), showWarnings = FALSE)

move_files <- function(dir, pattern, dest) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) > 0) {
    file.rename(files, file.path(dest, basename(files)))
    cat(sprintf("[OK] Moved %d %s files to %s\n", length(files),
                toupper(sub("\\.", "", pattern)), dest))
  }
}

move_files(out, "\\.png$", file.path(out, "png"))
move_files(out, "\\.csv$", file.path(out, "csv"))
move_files(out, "\\.pdf$", file.path(out, "pdf"))

cat("\n[OK] Output organization complete\n\n")

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
end_time <- Sys.time()
duration <- round(difftime(end_time, start_time, units = "mins"), 1)

cat("================================================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Duration:", duration, "minutes\n\n")

# List output files
png_files  <- list.files("../output/png", pattern = "\\.png$", full.names = FALSE)
csv_files  <- list.files("../output/csv", pattern = "\\.csv$", full.names = FALSE)
total_files <- length(png_files) + length(csv_files)
if (total_files > 0) {
  cat("Output files generated:", total_files, "\n")
  cat("  PNG charts:", length(png_files), "(in output/png/)\n")
  cat("  CSV files: ", length(csv_files), "(in output/csv/)\n")
  cat("\nFiles saved to: ../output/png/ and ../output/csv/\n")
}

cat("\n[DONE]\n\n")
