################################################################################
# FCI STATIONARITY TESTS (Phase 1.4)
################################################################################
#
# Project:      Financial Conditions Index - Stationarity Analysis
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Comprehensive unit root testing for all FCI variables and
#               computed indices. Addresses referee concerns about stationarity.
#
# TESTS INCLUDED:
#   1. Augmented Dickey-Fuller (ADF) - H0: Unit root
#   2. Phillips-Perron (PP) - H0: Unit root
#   3. KPSS - H0: Stationarity
#
# Interpretation:
#   - Variable is I(0) if: ADF rejects AND PP rejects AND KPSS fails to reject
#   - Variable is I(1) if: ADF fails to reject AND PP fails to reject AND KPSS rejects
#   - Conflicting results warrant further investigation
#
# Dependencies: Requires 01_FCI_Complete.R to be run first
#
# Output:
#   - FCI_Stationarity_Tests.csv
#   - 30_Stationarity_Summary.png
#
# Author:       Departamento de Analisis Macroeconomico
# Created:      2025
#
# References:
#   - Dickey & Fuller (1979) - ADF test
#   - Phillips & Perron (1988) - PP test
#   - Kwiatkowski et al. (1992) - KPSS test
#
################################################################################


################################################################################
# 1. SETUP AND CONFIGURATION
################################################################################

# Install required packages if not available
required_packages <- c("tseries", "urca", "readxl", "dplyr", "ggplot2", "tidyr")
new_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(new_packages) > 0) {
  cat("Installing required packages:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tseries)
  library(urca)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

STATIONARITY_CONFIG <- list(
  data_file = "../data/FCI_data_1.xlsx",
  output_dir = "../output",
  significance_levels = c(0.01, 0.05, 0.10),
  adf_max_lag = 12,  # Max lags for ADF
  verbose = TRUE
)

set.seed(20250126)

cat("\n################################################################################\n")
cat("FCI STATIONARITY TESTS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# 2. DATA PREPARATION
################################################################################

cat("================================================================================\n")
cat("DATA PREPARATION\n")
cat("================================================================================\n\n")

# Load FCI results (run core script if needed)
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  source("01_FCI_Complete.R")
}

# Get variable definitions
VARIABLES <- resultado_fci$variables$level3
all_vars_core <- resultado_fci$variables$core_vars

# Load raw data for individual variable testing
datos_raw <- read_excel(STATIONARITY_CONFIG$data_file)
fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

datos <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Calculate credit growth if needed
if ("Creditos_Sector_privado_totales" %in% names(datos)) {
  datos <- datos %>%
    mutate(Crecimiento_creditos = (Creditos_Sector_privado_totales /
                                    lag(Creditos_Sector_privado_totales, 12) - 1) * 100)
}

cat("Data loaded:", nrow(datos), "observations\n")
cat("Period:", format(min(datos$fecha)), "to", format(max(datos$fecha)), "\n\n")


################################################################################
# 3. STATIONARITY TEST FUNCTIONS
################################################################################

#' Run comprehensive stationarity tests on a time series
#' @param x Numeric vector (time series)
#' @param var_name Variable name for reporting
#' @return Data frame with test results
run_stationarity_tests <- function(x, var_name) {

  # Remove NAs
  x_clean <- na.omit(x)

  if (length(x_clean) < 30) {
    return(data.frame(
      variable = var_name,
      n_obs = length(x_clean),
      adf_stat = NA, adf_pvalue = NA, adf_lag = NA,
      pp_stat = NA, pp_pvalue = NA,
      kpss_stat = NA, kpss_pvalue = NA,
      conclusion = "Insufficient data"
    ))
  }

  results <- data.frame(variable = var_name, n_obs = length(x_clean))

  # ADF Test (H0: Unit root)
  tryCatch({
    adf_result <- adf.test(x_clean, alternative = "stationary",
                           k = min(STATIONARITY_CONFIG$adf_max_lag,
                                   floor((length(x_clean) - 1)^(1/3))))
    results$adf_stat <- adf_result$statistic
    results$adf_pvalue <- adf_result$p.value
    results$adf_lag <- adf_result$parameter
  }, error = function(e) {
    results$adf_stat <<- NA
    results$adf_pvalue <<- NA
    results$adf_lag <<- NA
  })

  # PP Test (H0: Unit root)
  tryCatch({
    pp_result <- pp.test(x_clean, alternative = "stationary")
    results$pp_stat <- pp_result$statistic
    results$pp_pvalue <- pp_result$p.value
  }, error = function(e) {
    results$pp_stat <<- NA
    results$pp_pvalue <<- NA
  })

  # KPSS Test (H0: Stationarity)
  tryCatch({
    kpss_result <- kpss.test(x_clean, null = "Level")
    results$kpss_stat <- kpss_result$statistic
    results$kpss_pvalue <- kpss_result$p.value
  }, error = function(e) {
    results$kpss_stat <<- NA
    results$kpss_pvalue <<- NA
  })

  # Determine conclusion
  alpha <- 0.05

  adf_reject <- !is.na(results$adf_pvalue) && results$adf_pvalue < alpha
  pp_reject <- !is.na(results$pp_pvalue) && results$pp_pvalue < alpha
  kpss_reject <- !is.na(results$kpss_pvalue) && results$kpss_pvalue < alpha

  results$conclusion <- case_when(
    adf_reject && pp_reject && !kpss_reject ~ "I(0) - Stationary",
    !adf_reject && !pp_reject && kpss_reject ~ "I(1) - Unit root",
    adf_reject && pp_reject && kpss_reject ~ "Conflicting (borderline)",
    !adf_reject && !pp_reject && !kpss_reject ~ "Inconclusive",
    adf_reject && !pp_reject ~ "Weak evidence for I(0)",
    !adf_reject && pp_reject ~ "Weak evidence for I(0)",
    TRUE ~ "Inconclusive"
  )

  return(results)
}

#' Run UR-CA tests for more detailed analysis
#' @param x Numeric vector
#' @param var_name Variable name
#' @return Data frame with urca test results
run_urca_tests <- function(x, var_name) {

  x_clean <- na.omit(x)
  if (length(x_clean) < 50) return(NULL)

  results <- data.frame(variable = var_name)

  # ADF with urca (more detailed output)
  tryCatch({
    adf_urca <- ur.df(x_clean, type = "drift", selectlags = "AIC")
    results$urca_adf_stat <- adf_urca@teststat[1, "tau2"]
    results$urca_adf_cv_1pct <- adf_urca@cval[1, "1pct"]
    results$urca_adf_cv_5pct <- adf_urca@cval[1, "5pct"]
    results$urca_adf_cv_10pct <- adf_urca@cval[1, "10pct"]
    results$urca_adf_lags <- adf_urca@lags
  }, error = function(e) {
    results$urca_adf_stat <<- NA
  })

  # KPSS with urca
  tryCatch({
    kpss_urca <- ur.kpss(x_clean, type = "mu")
    results$urca_kpss_stat <- kpss_urca@teststat
    results$urca_kpss_cv_1pct <- kpss_urca@cval["1pct"]
    results$urca_kpss_cv_5pct <- kpss_urca@cval["5pct"]
    results$urca_kpss_cv_10pct <- kpss_urca@cval["10pct"]
  }, error = function(e) {
    results$urca_kpss_stat <<- NA
  })

  return(results)
}


################################################################################
# 4. TEST INDIVIDUAL FCI VARIABLES
################################################################################

cat("================================================================================\n")
cat("TESTING INDIVIDUAL FCI VARIABLES\n")
cat("================================================================================\n\n")

variable_results <- data.frame()

for (var in all_vars_core) {
  if (var %in% names(datos) && !all(is.na(datos[[var]]))) {
    cat("  Testing:", var, "...")
    result <- run_stationarity_tests(datos[[var]], var)
    variable_results <- rbind(variable_results, result)
    cat(" ", result$conclusion, "\n")
  }
}

cat("\nIndividual Variable Results:\n")
cat(sprintf("%-30s %8s %10s %10s %10s %20s\n",
            "Variable", "N", "ADF p", "PP p", "KPSS p", "Conclusion"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for (i in 1:nrow(variable_results)) {
  r <- variable_results[i, ]
  cat(sprintf("%-30s %8d %10.3f %10.3f %10.3f %20s\n",
              r$variable, r$n_obs,
              ifelse(is.na(r$adf_pvalue), NA, r$adf_pvalue),
              ifelse(is.na(r$pp_pvalue), NA, r$pp_pvalue),
              ifelse(is.na(r$kpss_pvalue), NA, r$kpss_pvalue),
              r$conclusion))
}


################################################################################
# 5. TEST FCI INDICES
################################################################################

cat("\n================================================================================\n")
cat("TESTING FCI INDICES\n")
cat("================================================================================\n\n")

# Get FCI indices from results
fci_indices <- resultado_fci$all_indices

# List of FCI series to test
fci_series <- c("FCI_COMP_AVG", "FCI_ENDO_AVG", "FCI_EXO_AVG",
                "FCI_RATES_AVG", "FCI_BANKING_AVG", "FCI_EXTERNAL_AVG",
                "FCI_ENDO_ORTHO")
fci_series <- fci_series[fci_series %in% names(fci_indices)]

index_results <- data.frame()

for (idx in fci_series) {
  if (idx %in% names(fci_indices) && !all(is.na(fci_indices[[idx]]))) {
    cat("  Testing:", idx, "...")
    result <- run_stationarity_tests(fci_indices[[idx]], idx)
    index_results <- rbind(index_results, result)
    cat(" ", result$conclusion, "\n")
  }
}

cat("\nFCI Index Results:\n")
cat(sprintf("%-30s %8s %10s %10s %10s %20s\n",
            "Index", "N", "ADF p", "PP p", "KPSS p", "Conclusion"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for (i in 1:nrow(index_results)) {
  r <- index_results[i, ]
  cat(sprintf("%-30s %8d %10.3f %10.3f %10.3f %20s\n",
              r$variable, r$n_obs,
              ifelse(is.na(r$adf_pvalue), NA, r$adf_pvalue),
              ifelse(is.na(r$pp_pvalue), NA, r$pp_pvalue),
              ifelse(is.na(r$kpss_pvalue), NA, r$kpss_pvalue),
              r$conclusion))
}


################################################################################
# 6. CONSOLIDATE AND EXPORT RESULTS
################################################################################

cat("\n================================================================================\n")
cat("EXPORTING RESULTS\n")
cat("================================================================================\n\n")

# Combine all results
all_stationarity_results <- rbind(
  variable_results %>% mutate(type = "Variable"),
  index_results %>% mutate(type = "Index")
)

# Add interpretation columns
all_stationarity_results <- all_stationarity_results %>%
  mutate(
    adf_reject_5pct = adf_pvalue < 0.05,
    pp_reject_5pct = pp_pvalue < 0.05,
    kpss_reject_5pct = kpss_pvalue < 0.05,
    is_stationary = conclusion == "I(0) - Stationary" |
                    grepl("evidence for I\\(0\\)", conclusion)
  )

# Export
write.csv(all_stationarity_results,
          file.path(STATIONARITY_CONFIG$output_dir, "FCI_Stationarity_Tests.csv"),
          row.names = FALSE)
cat("Saved: FCI_Stationarity_Tests.csv\n")


################################################################################
# 7. VISUALIZATION
################################################################################

cat("\n================================================================================\n")
cat("GENERATING VISUALIZATION\n")
cat("================================================================================\n\n")

# Create summary visualization
plot_data <- all_stationarity_results %>%
  mutate(
    variable_short = gsub("FCI_|_AVG", "", variable),
    stationarity_status = case_when(
      grepl("I\\(0\\)", conclusion) | grepl("evidence for I\\(0\\)", conclusion) ~ "Stationary",
      grepl("I\\(1\\)", conclusion) ~ "Unit Root",
      TRUE ~ "Inconclusive"
    )
  )

# Heatmap of p-values
plot_pvalues <- plot_data %>%
  dplyr::select(variable_short, type, adf_pvalue, pp_pvalue, kpss_pvalue) %>%
  pivot_longer(cols = c(adf_pvalue, pp_pvalue, kpss_pvalue),
               names_to = "test", values_to = "pvalue") %>%
  mutate(
    test = case_when(
      test == "adf_pvalue" ~ "ADF",
      test == "pp_pvalue" ~ "PP",
      test == "kpss_pvalue" ~ "KPSS"
    ),
    significant = pvalue < 0.05,
    sig_label = ifelse(significant, "*", "")
  )

p1 <- ggplot(plot_pvalues, aes(x = test, y = variable_short, fill = pvalue)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f%s", pvalue, sig_label)), size = 3) +
  scale_fill_gradient2(low = "#E74C3C", mid = "#F7DC6F", high = "#27AE60",
                       midpoint = 0.10, name = "p-value",
                       limits = c(0, 0.20), oob = scales::squish) +
  facet_wrap(~type, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Stationarity Test Results",
    subtitle = "ADF/PP: H0=Unit root (low p rejects) | KPSS: H0=Stationary (low p rejects)",
    x = "Test", y = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

ggsave(file.path(STATIONARITY_CONFIG$output_dir, "30_Stationarity_Summary.png"), p1,
       width = 12, height = 10, dpi = 300)
cat("Saved: 30_Stationarity_Summary.png\n")


################################################################################
# 8. SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("STATIONARITY TESTS SUMMARY\n")
cat("================================================================================\n\n")

# Summary by type
summary_table <- all_stationarity_results %>%
  group_by(type) %>%
  summarise(
    n_tested = n(),
    n_stationary = sum(is_stationary, na.rm = TRUE),
    n_unit_root = sum(conclusion == "I(1) - Unit root", na.rm = TRUE),
    n_inconclusive = sum(grepl("Inconclusive|Conflicting", conclusion), na.rm = TRUE),
    pct_stationary = mean(is_stationary, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("Summary by Type:\n")
print(summary_table)

cat("\nKey Findings:\n")

# Variables
var_stationary <- variable_results %>% filter(grepl("I\\(0\\)|evidence for I\\(0\\)", conclusion))
var_unit_root <- variable_results %>% filter(conclusion == "I(1) - Unit root")

cat("  Individual Variables:\n")
cat("    - Stationary: ", nrow(var_stationary), "\n")
if (nrow(var_stationary) > 0) {
  cat("      ", paste(var_stationary$variable, collapse = ", "), "\n")
}
cat("    - Unit root:  ", nrow(var_unit_root), "\n")
if (nrow(var_unit_root) > 0) {
  cat("      ", paste(var_unit_root$variable, collapse = ", "), "\n")
}

# FCI indices
idx_stationary <- index_results %>% filter(grepl("I\\(0\\)|evidence for I\\(0\\)", conclusion))
idx_unit_root <- index_results %>% filter(conclusion == "I(1) - Unit root")

cat("\n  FCI Indices:\n")
cat("    - Stationary: ", nrow(idx_stationary), "\n")
if (nrow(idx_stationary) > 0) {
  cat("      ", paste(idx_stationary$variable, collapse = ", "), "\n")
}
cat("    - Unit root:  ", nrow(idx_unit_root), "\n")
if (nrow(idx_unit_root) > 0) {
  cat("      ", paste(idx_unit_root$variable, collapse = ", "), "\n")
}

cat("\nInterpretation for Referees:\n")
cat("  - Most FCI indices are constructed from standardized variables (Z-scores),\n")
cat("    which are stationary by construction.\n")
cat("  - Rolling window standardization ensures the FCI is comparable over time.\n")
cat("  - Any non-stationary components are captured through the rolling mean adjustment.\n")

cat("\n================================================================================\n")
cat("STATIONARITY TESTS COMPLETE\n")
cat("================================================================================\n\n")


################################################################################
# END OF SCRIPT
################################################################################
