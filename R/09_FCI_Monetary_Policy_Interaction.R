################################################################################
# FCI-MONETARY POLICY INTERACTION AND SYSTEMATIC OUT-OF-SAMPLE FORECASTING
################################################################################
#
# Project:      Financial Conditions Index - Monetary Policy Interaction
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Analysis of two-way interaction between FCI and monetary policy
#               plus systematic out-of-sample RMSE comparison across models
#
#   PART 1: FCI x TPM Interaction Analysis
#     A: Does FCI amplify/dampen monetary policy transmission?
#        Model: y_{t+h} = a + b1*dTPM + b2*FCI + b3*(dTPM x FCI) + controls
#     B: Does TPM stance affect FCI transmission?
#        Model: y_{t+h} = a + b1*FCI + b2*stance + b3*(FCI x stance) + controls
#
#   PART 2: Systematic Out-of-Sample RMSE Comparison
#     - 7 models: AR, AR+FCI_COMP, AR+FCI_ENDO, AR+FCI_EXO, VAR, VAR+FCI, Factor-AR
#     - Horizons: h = 1, 3, 6, 12 months
#     - Windows: Expanding (T0=80) and Rolling (60 months)
#     - Diebold-Mariano tests for equal predictive accuracy
#
# Key Research Questions:
#   1. Does FCI amplify or dampen monetary policy transmission?
#   2. Does the monetary policy stance affect FCI transmission to the economy?
#   3. Which FCI specification provides the best forecasting performance?
#
# References:
#   - Jorda (2005) - Local Projections
#   - Diebold & Mariano (1995) - Comparing Predictive Accuracy
#   - Newey & West (1987) - HAC Standard Errors
#
# Dependencies: Requires 01_FCI_Complete.R to be run first
#
# Output:
#   - 90-99_*.png (visualizations)
#   - FCI_TPM_Interaction_*.csv (interaction results)
#   - Forecast_*.csv (forecasting results)
#
# Author:       Departamento de Analisis Macroeconomico
# Created:      2025-01-23
#
################################################################################


################################################################################
# SECTION 1: SETUP AND CONFIGURATION
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(ggplot2)
  library(gridExtra)
  library(sandwich)
  library(lmtest)
  library(vars)
})

# Configuration parameters
MP_CONFIG <- list(
  # Data
  data_file = "../data/FCI_data_1.xlsx",
  macro_sheet = "Datos_macro",

  # LP parameters
  horizons_interaction = c(3, 6, 12, 18, 24),
  max_horizon = 24,
  n_lags = 2,
  confidence_level = 0.90,

  # Forecasting parameters
  horizons_forecast = c(1, 3, 6, 12),
  min_train = 80,           # Minimum training observations for expanding window
  rolling_window = 60,      # Rolling window size

  # HP filter parameter for trend extraction
  hp_lambda = 14400,        # Monthly data standard

  # FCI types for robustness
  fci_types = c("FCI_COMP_AVG", "FCI_ENDO_AVG", "FCI_BANKING_AVG"),

  # Target variables
  target_vars = c("IMAEP_yoy", "IPC_yoy", "Credit_yoy"),

  # Output
  output_dir = "../output",
  verbose = TRUE
)

set.seed(20250123)

cat("\n################################################################################\n")
cat("FCI-MONETARY POLICY INTERACTION AND SYSTEMATIC FORECASTING\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# SECTION 2: DATA PREPARATION (TPM TRANSFORMATIONS)
################################################################################

cat("================================================================================\n")
cat("DATA PREPARATION\n")
cat("================================================================================\n\n")

# Load FCI results (run core script if needed)
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  source("01_FCI_Complete.R")
}

# Extract FCI data
fci_data <- resultado_fci$all_indices %>%
  dplyr::select(fecha,
                FCI_COMP = FCI_COMP_AVG,
                FCI_ENDO = FCI_ENDO_AVG,
                FCI_EXO = FCI_EXO_AVG,
                FCI_exCredit = FCI_exCredit_AVG)

# Check for FCI_BANKING_AVG
if ("FCI_BANKING_AVG" %in% names(resultado_fci$all_indices)) {
  fci_data <- fci_data %>%
    left_join(resultado_fci$all_indices %>%
                dplyr::select(fecha, FCI_BANKING = FCI_BANKING_AVG),
              by = "fecha")
}

# Load raw data for TPM
datos_raw <- read_excel(MP_CONFIG$data_file)
fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

datos <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Load macro data
macro_raw <- tryCatch({
  read_excel(MP_CONFIG$data_file, sheet = MP_CONFIG$macro_sheet)
}, error = function(e) {
  cat("Note: Macro sheet not found, using main sheet for macro data\n")
  datos_raw
})

fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

macro_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

#' Calculate HP filter trend
#' @param x Numeric vector
#' @param lambda Smoothing parameter
#' @return Trend component
hp_filter <- function(x, lambda = 14400) {
  n <- length(x)
  if (n < 10) return(rep(NA, n))

  # Build penalty matrix
  I <- diag(n)
  D <- matrix(0, n - 2, n)
  for (i in 1:(n - 2)) {
    D[i, i] <- 1
    D[i, i + 1] <- -2
    D[i, i + 2] <- 1
  }

  # HP filter: (I + lambda * D'D)^{-1} * x
  trend <- tryCatch({
    solve(I + lambda * t(D) %*% D) %*% x
  }, error = function(e) {
    # Fallback to rolling mean
    zoo::rollmean(x, k = 24, fill = NA, align = "center")
  })

  return(as.vector(trend))
}

#' Calculate TPM stance measures
#' @param data Data frame with TPM
#' @return Data frame with stance measures
calculate_tpm_stance <- function(data) {
  data %>%
    mutate(
      # Raw TPM
      TPM = TPM,

      # Monthly change
      dTPM_m = TPM - lag(TPM, 1),

      # YoY change
      dTPM_yoy = TPM - lag(TPM, 12),

      # Standardized level stance
      TPM_std = (TPM - mean(TPM, na.rm = TRUE)) / sd(TPM, na.rm = TRUE),

      # HP deviation stance
      TPM_trend = hp_filter(TPM, MP_CONFIG$hp_lambda),
      TPM_hp_dev = TPM - TPM_trend,
      TPM_hp_std = (TPM_hp_dev - mean(TPM_hp_dev, na.rm = TRUE)) / sd(TPM_hp_dev, na.rm = TRUE),

      # Binary stance: 1 if TPM > median
      TPM_binary = as.integer(TPM > median(TPM, na.rm = TRUE))
    )
}

# Prepare TPM data
if ("TPM" %in% names(datos)) {
  tpm_data <- datos %>%
    dplyr::select(fecha, TPM) %>%
    calculate_tpm_stance()
} else {
  stop("TPM not found in data file")
}

# Calculate credit growth
if ("Creditos_Sector_privado_totales" %in% names(datos)) {
  credit_data <- datos %>%
    dplyr::select(fecha, Creditos_Sector_privado_totales) %>%
    mutate(Credit_yoy = (Creditos_Sector_privado_totales /
                           lag(Creditos_Sector_privado_totales, 12) - 1) * 100) %>%
    dplyr::select(fecha, Credit_yoy)
} else {
  credit_data <- data.frame(fecha = datos$fecha, Credit_yoy = NA)
}

# Calculate macro YoY growth
for (var in c("IMAEP", "IPC")) {
  if (var %in% names(macro_data)) {
    macro_data <- macro_data %>%
      mutate(!!paste0(var, "_yoy") := (!!sym(var) / lag(!!sym(var), 12) - 1) * 100)
  }
}

# Merge all data
analysis_data <- fci_data %>%
  inner_join(tpm_data, by = "fecha") %>%
  inner_join(credit_data, by = "fecha") %>%
  left_join(macro_data %>% dplyr::select(fecha, any_of(c("IMAEP_yoy", "IPC_yoy"))),
            by = "fecha") %>%
  arrange(fecha) %>%
  na.omit()

cat("Analysis data prepared:\n")
cat("  Observations:", nrow(analysis_data), "\n")
cat("  Period:", format(min(analysis_data$fecha)), "to",
    format(max(analysis_data$fecha)), "\n\n")

# Summary statistics
cat("TPM Statistics:\n")
cat(sprintf("  Level: Mean=%.2f%%, SD=%.2f%%\n",
            mean(analysis_data$TPM), sd(analysis_data$TPM)))
cat(sprintf("  Monthly change: Mean=%.3f, SD=%.3f\n",
            mean(analysis_data$dTPM_m, na.rm = TRUE),
            sd(analysis_data$dTPM_m, na.rm = TRUE)))
cat(sprintf("  Tightening months: %.1f%%\n",
            mean(analysis_data$TPM_binary, na.rm = TRUE) * 100))


################################################################################
# SECTION 3: PART 1A - FCI AMPLIFICATION OF MONETARY POLICY
################################################################################

cat("\n================================================================================\n")
cat("PART 1A: FCI AMPLIFICATION OF MONETARY POLICY\n")
cat("================================================================================\n\n")

cat("Model: y_{t+h} = a + b1*dTPM + b2*FCI + b3*(dTPM x FCI) + controls + e\n")
cat("Interpretation:\n")
cat("  b3 > 0 with dTPM > 0: FCI AMPLIFIES contractionary policy\n")
cat("  b3 < 0: FCI DAMPENS policy transmission\n\n")

#' Run LP with TPM x FCI interaction (amplification model)
#' @param data Data frame
#' @param y_var Dependent variable name
#' @param fci_var FCI variable name
#' @param tpm_change_var TPM change variable (dTPM_m or dTPM_yoy)
#' @param max_h Maximum horizon
#' @param n_lags Number of lags
#' @return Data frame with results
run_lp_tpm_fci_interaction <- function(data, y_var, fci_var, tpm_change_var = "dTPM_yoy",
                                        max_h, n_lags = 2) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - MP_CONFIG$confidence_level) / 2)

  for (h in 1:max_h) {
    # Create variables
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        fci_lag1 = lag(!!sym(fci_var), 1),
        tpm_x_fci = !!sym(tpm_change_var) * !!sym(fci_var)
      )

    # Add more lags
    for (i in 2:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", i) := lag(!!sym(y_var), i))
    }

    # Build formula
    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")
    formula_str <- paste("y_fwd ~", tpm_change_var, "+", fci_var, "+ tpm_x_fci +",
                         paste(lag_vars, collapse = " + "))

    # Estimate
    reg_data <- data_h %>%
      dplyr::select(y_fwd, !!sym(tpm_change_var), !!sym(fci_var), tpm_x_fci,
                    all_of(lag_vars)) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    # Extract coefficients
    idx_tpm <- which(rownames(coef_test) == tpm_change_var)
    idx_fci <- which(rownames(coef_test) == fci_var)
    idx_interact <- which(rownames(coef_test) == "tpm_x_fci")

    results <- rbind(results, data.frame(
      horizon = h,
      # TPM effect (when FCI = 0)
      coef_tpm = coef_test[idx_tpm, 1],
      se_tpm = coef_test[idx_tpm, 2],
      p_tpm = coef_test[idx_tpm, 4],
      # FCI effect (when dTPM = 0)
      coef_fci = coef_test[idx_fci, 1],
      se_fci = coef_test[idx_fci, 2],
      p_fci = coef_test[idx_fci, 4],
      # Interaction effect
      coef_interact = coef_test[idx_interact, 1],
      se_interact = coef_test[idx_interact, 2],
      p_interact = coef_test[idx_interact, 4],
      ci_interact_lo = coef_test[idx_interact, 1] - z_crit * coef_test[idx_interact, 2],
      ci_interact_hi = coef_test[idx_interact, 1] + z_crit * coef_test[idx_interact, 2],
      n_obs = nrow(reg_data),
      r_squared = summary(model)$r.squared
    ))
  }

  return(results)
}

# Run amplification analysis for each target and FCI type
amplification_results <- list()

# Primary analysis: FCI_COMP with YoY TPM change
# Use FCI_exCredit when target is credit (endogeneity correction)
for (target in c("IMAEP_yoy", "IPC_yoy", "Credit_yoy")) {
  if (!target %in% names(analysis_data)) next

  fci_use <- if (target == "Credit_yoy" && "FCI_exCredit" %in% names(analysis_data)) "FCI_exCredit" else "FCI_COMP"
  cat("  ", target, "with", fci_use, "... ")

  result <- run_lp_tpm_fci_interaction(
    analysis_data, target, fci_use, "dTPM_yoy",
    MP_CONFIG$max_horizon, MP_CONFIG$n_lags
  )

  if (nrow(result) > 0) {
    result <- result %>%
      mutate(target_var = target, fci_type = "FCI_COMP", tpm_change = "YoY")
    amplification_results[[paste(target, "FCI_COMP", sep = "_")]] <- result
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

# Robustness: FCI_BANKING (excludes TPM)
for (target in c("IMAEP_yoy", "IPC_yoy", "Credit_yoy")) {
  if (!target %in% names(analysis_data)) next
  if (!"FCI_BANKING" %in% names(analysis_data)) next

  cat("  ", target, "with FCI_BANKING (robustness) ... ")

  result <- run_lp_tpm_fci_interaction(
    analysis_data, target, "FCI_BANKING", "dTPM_yoy",
    MP_CONFIG$max_horizon, MP_CONFIG$n_lags
  )

  if (nrow(result) > 0) {
    result <- result %>%
      mutate(target_var = target, fci_type = "FCI_BANKING", tpm_change = "YoY")
    amplification_results[[paste(target, "FCI_BANKING", sep = "_")]] <- result
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

# Combine results
all_amplification <- bind_rows(amplification_results)

# Results summary
cat("\n================================================================================\n")
cat("AMPLIFICATION RESULTS SUMMARY (at key horizons)\n")
cat("================================================================================\n\n")

if (nrow(all_amplification) > 0) {
  summary_amp <- all_amplification %>%
    filter(horizon %in% MP_CONFIG$horizons_interaction, fci_type == "FCI_COMP") %>%
    mutate(
      stars = case_when(
        p_interact < 0.01 ~ "***",
        p_interact < 0.05 ~ "**",
        p_interact < 0.10 ~ "*",
        TRUE ~ ""
      )
    )

  cat("dTPM x FCI Interaction Coefficient (b3)\n")
  cat("*** p<0.01, ** p<0.05, * p<0.10\n\n")
  cat(sprintf("%-15s %8s %8s %8s %8s %8s\n",
              "Target", "h=3", "h=6", "h=9", "h=12", "h=18"))
  cat(paste(rep("-", 65), collapse = ""), "\n")

  for (target in unique(summary_amp$target_var)) {
    vals <- summary_amp %>% filter(target_var == target)
    row_str <- sprintf("%-15s", gsub("_yoy", "", target))

    for (h in c(3, 6, 9, 12, 18)) {
      v <- vals %>% filter(horizon == h)
      if (nrow(v) > 0) {
        row_str <- paste0(row_str, sprintf(" %+7.3f%s", v$coef_interact, v$stars))
      } else {
        row_str <- paste0(row_str, "       NA")
      }
    }
    cat(row_str, "\n")
  }

  # Interpretation
  cat("\nINTERPRETATION:\n")

  for (target in unique(summary_amp$target_var)) {
    v12 <- summary_amp %>% filter(target_var == target, horizon == 12)
    if (nrow(v12) == 0) next

    target_label <- gsub("_yoy", "", target)

    if (v12$coef_interact > 0 && v12$p_interact < 0.10) {
      cat(sprintf("  %s: FCI AMPLIFIES monetary policy (b3=%+.3f, p=%.3f)\n",
                  target_label, v12$coef_interact, v12$p_interact))
    } else if (v12$coef_interact < 0 && v12$p_interact < 0.10) {
      cat(sprintf("  %s: FCI DAMPENS monetary policy (b3=%+.3f, p=%.3f)\n",
                  target_label, v12$coef_interact, v12$p_interact))
    } else {
      cat(sprintf("  %s: No significant interaction (b3=%+.3f, p=%.3f)\n",
                  target_label, v12$coef_interact, v12$p_interact))
    }
  }
}


################################################################################
# SECTION 4: PART 1B - TPM STANCE EFFECT ON FCI TRANSMISSION
################################################################################

cat("\n================================================================================\n")
cat("PART 1B: TPM STANCE EFFECT ON FCI TRANSMISSION\n")
cat("================================================================================\n\n")

cat("Model: y_{t+h} = a + b1*FCI + b2*stance + b3*(FCI x stance) + controls + e\n")
cat("Interpretation:\n")
cat("  b1 = FCI effect when stance is low/neutral\n")
cat("  b1 + b3 = FCI effect when stance is high/tight\n")
cat("  b3 != 0 = FCI effects depend on policy stance\n\n")

#' Run LP with FCI x TPM stance interaction
#' @param data Data frame
#' @param y_var Dependent variable name
#' @param fci_var FCI variable name
#' @param stance_var Stance variable name
#' @param max_h Maximum horizon
#' @param n_lags Number of lags
#' @return Data frame with results
run_lp_fci_stance_interaction <- function(data, y_var, fci_var, stance_var,
                                           max_h, n_lags = 2) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - MP_CONFIG$confidence_level) / 2)

  for (h in 1:max_h) {
    # Create variables
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        fci_lag1 = lag(!!sym(fci_var), 1),
        fci_x_stance = !!sym(fci_var) * !!sym(stance_var)
      )

    # Add more lags
    for (i in 2:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", i) := lag(!!sym(y_var), i))
    }

    # Build formula
    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")
    formula_str <- paste("y_fwd ~", fci_var, "+", stance_var, "+ fci_x_stance +",
                         paste(lag_vars, collapse = " + "))

    # Estimate
    reg_data <- data_h %>%
      dplyr::select(y_fwd, !!sym(fci_var), !!sym(stance_var), fci_x_stance,
                    all_of(lag_vars)) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    # Extract coefficients
    idx_fci <- which(rownames(coef_test) == fci_var)
    idx_stance <- which(rownames(coef_test) == stance_var)
    idx_interact <- which(rownames(coef_test) == "fci_x_stance")

    # Effect in low stance (stance = 0 for binary, or at mean for continuous)
    effect_low <- coef_test[idx_fci, 1]
    se_low <- coef_test[idx_fci, 2]

    # Effect in high stance (for binary: stance = 1)
    effect_high <- coef_test[idx_fci, 1] + coef_test[idx_interact, 1]
    se_high <- sqrt(vcov_hac[idx_fci, idx_fci] + vcov_hac[idx_interact, idx_interact] +
                      2 * vcov_hac[idx_fci, idx_interact])

    results <- rbind(results, data.frame(
      horizon = h,
      # FCI effect in low stance
      coef_fci_low = effect_low,
      se_fci_low = se_low,
      p_fci_low = coef_test[idx_fci, 4],
      ci_fci_low_lo = effect_low - z_crit * se_low,
      ci_fci_low_hi = effect_low + z_crit * se_low,
      # FCI effect in high stance
      coef_fci_high = effect_high,
      se_fci_high = se_high,
      ci_fci_high_lo = effect_high - z_crit * se_high,
      ci_fci_high_hi = effect_high + z_crit * se_high,
      # Interaction (differential)
      coef_interact = coef_test[idx_interact, 1],
      se_interact = coef_test[idx_interact, 2],
      p_interact = coef_test[idx_interact, 4],
      n_obs = nrow(reg_data)
    ))
  }

  return(results)
}

# Run stance analysis with three definitions
stance_results <- list()

stance_definitions <- c(
  "TPM_std" = "Level (Standardized)",
  "TPM_hp_std" = "HP Deviation",
  "TPM_binary" = "Binary (> median)"
)

for (target in c("IMAEP_yoy", "IPC_yoy", "Credit_yoy")) {
  if (!target %in% names(analysis_data)) next

  # Endogeneity correction: use FCI_exCredit when target is credit
  fci_use <- if (target == "Credit_yoy" && "FCI_exCredit" %in% names(analysis_data)) "FCI_exCredit" else "FCI_COMP"

  for (stance_var in names(stance_definitions)) {
    if (!stance_var %in% names(analysis_data)) next

    cat("  ", target, "with stance:", stance_definitions[stance_var], "... ")

    result <- run_lp_fci_stance_interaction(
      analysis_data, target, fci_use, stance_var,
      MP_CONFIG$max_horizon, MP_CONFIG$n_lags
    )

    if (nrow(result) > 0) {
      result <- result %>%
        mutate(target_var = target,
               stance_type = stance_var,
               stance_label = stance_definitions[stance_var])
      stance_results[[paste(target, stance_var, sep = "_")]] <- result
      cat("done\n")
    } else {
      cat("insufficient data\n")
    }
  }
}

# Combine results
all_stance <- bind_rows(stance_results)

# Results summary
cat("\n================================================================================\n")
cat("STANCE INTERACTION RESULTS (at h=12)\n")
cat("================================================================================\n\n")

if (nrow(all_stance) > 0) {
  stance_12 <- all_stance %>%
    filter(horizon == 12) %>%
    mutate(
      stars = case_when(
        p_interact < 0.01 ~ "***",
        p_interact < 0.05 ~ "**",
        p_interact < 0.10 ~ "*",
        TRUE ~ ""
      )
    )

  cat("FCI Effect in Low vs High TPM Stance (h=12)\n")
  cat("*** p<0.01, ** p<0.05, * p<0.10 for interaction\n\n")
  cat(sprintf("%-12s %-20s %10s %10s %10s\n",
              "Target", "Stance Definition", "Low", "High", "Interact."))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in 1:nrow(stance_12)) {
    r <- stance_12[i, ]
    cat(sprintf("%-12s %-20s %+10.2f %+10.2f %+10.3f%s\n",
                gsub("_yoy", "", r$target_var),
                r$stance_label,
                r$coef_fci_low,
                r$coef_fci_high,
                r$coef_interact,
                r$stars))
  }

  # Key finding
  cat("\nKEY FINDING:\n")
  binary_credit <- stance_12 %>%
    filter(target_var == "Credit_yoy", stance_type == "TPM_binary")

  if (nrow(binary_credit) > 0) {
    if (binary_credit$p_interact < 0.10) {
      direction <- ifelse(abs(binary_credit$coef_fci_high) > abs(binary_credit$coef_fci_low),
                          "STRONGER", "WEAKER")
      cat(sprintf("  Credit: FCI effects are %s when TPM is above median\n", direction))
    } else {
      cat("  Credit: FCI effects do not significantly depend on TPM stance\n")
    }
  }
}


################################################################################
# SECTION 5: PART 2A - FORECASTING MODEL SPECIFICATIONS
################################################################################

cat("\n================================================================================\n")
cat("PART 2: SYSTEMATIC OUT-OF-SAMPLE FORECASTING\n")
cat("================================================================================\n\n")

cat("Models to compare:\n")
cat("  1. AR (baseline): y_t = a + sum(b_i * y_{t-i}) + e\n")
cat("  2. AR + FCI_COMP: y_t = a + sum(b_i * y_{t-i}) + sum(g_j * FCI_{t-j}) + e\n")
cat("  3. AR + FCI_ENDO: Same with FCI_ENDO\n")
cat("  4. AR + FCI_EXO: Same with FCI_EXO\n")
cat("  5. VAR without FCI: VAR(p) on [y, macro vars]\n")
cat("  6. VAR with FCI: VAR(p) on [y, macro vars, FCI]\n")
cat("  7. Factor-AR: AR + PC of all FCI types\n\n")

#' Forecast with AR model
#' @param train_data Training data
#' @param target Target variable name
#' @param h Forecast horizon
#' @param n_lags Number of AR lags
#' @return Forecast value
forecast_ar <- function(train_data, target, h, n_lags = 2) {
  y <- train_data[[target]]
  n <- length(y)

  if (n < n_lags + 10) return(NA)

  # Create lagged data
  df <- data.frame(y = y)
  for (i in 1:n_lags) {
    df[[paste0("y_lag", i)]] <- c(rep(NA, i), y[1:(n - i)])
  }
  df <- na.omit(df)

  if (nrow(df) < 20) return(NA)

  # Estimate AR
  formula_str <- paste("y ~", paste(paste0("y_lag", 1:n_lags), collapse = " + "))
  model <- tryCatch({
    lm(as.formula(formula_str), data = df)
  }, error = function(e) NULL)

  if (is.null(model)) return(NA)

  # Iterated forecast
  forecast_vals <- tail(y, n_lags)

  for (step in 1:h) {
    newdata <- data.frame(t(forecast_vals[length(forecast_vals):(length(forecast_vals) - n_lags + 1)]))
    names(newdata) <- paste0("y_lag", 1:n_lags)
    pred <- predict(model, newdata = newdata)
    forecast_vals <- c(forecast_vals, pred)
  }

  return(tail(forecast_vals, 1))
}

#' Forecast with AR + FCI (ADL model)
#' @param train_data Training data
#' @param target Target variable name
#' @param fci_var FCI variable name
#' @param h Forecast horizon
#' @param n_lags Number of lags
#' @return Forecast value
forecast_ar_fci <- function(train_data, target, fci_var, h, n_lags = 2) {
  y <- train_data[[target]]
  fci <- train_data[[fci_var]]
  n <- length(y)

  if (n < n_lags + 10 || length(fci) != n) return(NA)

  # Create lagged data
  df <- data.frame(y = y, fci = fci)
  for (i in 1:n_lags) {
    df[[paste0("y_lag", i)]] <- c(rep(NA, i), y[1:(n - i)])
    df[[paste0("fci_lag", i)]] <- c(rep(NA, i), fci[1:(n - i)])
  }
  df <- na.omit(df)

  if (nrow(df) < 20) return(NA)

  # Estimate ADL
  y_lags <- paste(paste0("y_lag", 1:n_lags), collapse = " + ")
  fci_lags <- paste(paste0("fci_lag", 1:n_lags), collapse = " + ")
  formula_str <- paste("y ~ fci +", y_lags, "+", fci_lags)

  model <- tryCatch({
    lm(as.formula(formula_str), data = df)
  }, error = function(e) NULL)

  if (is.null(model)) return(NA)

  # Direct forecast (use current FCI for h-step ahead)
  last_y <- tail(y, n_lags)
  last_fci <- tail(fci, n_lags)

  newdata <- data.frame(
    fci = tail(fci, 1)
  )
  for (i in 1:n_lags) {
    newdata[[paste0("y_lag", i)]] <- last_y[n_lags - i + 1]
    newdata[[paste0("fci_lag", i)]] <- last_fci[n_lags - i + 1]
  }

  pred <- tryCatch({
    predict(model, newdata = newdata)
  }, error = function(e) NA)

  return(pred)
}

#' Forecast with VAR
#' @param train_data Training data
#' @param target Target variable name
#' @param vars_include Variables to include
#' @param h Forecast horizon
#' @param p VAR lag order
#' @return Forecast value
forecast_var <- function(train_data, target, vars_include, h, p = 2) {
  var_data <- train_data[, vars_include, drop = FALSE]
  var_data <- var_data[complete.cases(var_data), , drop = FALSE]

  if (nrow(var_data) < p * length(vars_include) + 20) return(NA)

  model <- tryCatch({
    VAR(var_data, p = p, type = "const")
  }, error = function(e) NULL)

  if (is.null(model)) return(NA)

  forecast <- tryCatch({
    predict(model, n.ahead = h)$fcst[[target]][h, "fcst"]
  }, error = function(e) NA)

  return(forecast)
}

#' Calculate Factor-AR forecast (AR + PC of FCI types)
#' @param train_data Training data
#' @param target Target variable name
#' @param fci_vars Vector of FCI variable names
#' @param h Forecast horizon
#' @param n_lags Number of lags
#' @return Forecast value
forecast_factor_ar <- function(train_data, target, fci_vars, h, n_lags = 2) {
  y <- train_data[[target]]
  n <- length(y)

  # Extract FCI data
  fci_data <- train_data[, fci_vars, drop = FALSE]
  fci_data <- fci_data[complete.cases(fci_data), , drop = FALSE]

  if (nrow(fci_data) < 30 || nrow(fci_data) != n) return(NA)

  # Calculate first PC
  pca <- tryCatch({
    prcomp(fci_data, center = TRUE, scale. = TRUE)
  }, error = function(e) NULL)

  if (is.null(pca)) return(NA)

  fci_pc1 <- pca$x[, 1]

  # Create lagged data
  df <- data.frame(y = y, fci_pc1 = fci_pc1)
  for (i in 1:n_lags) {
    df[[paste0("y_lag", i)]] <- c(rep(NA, i), y[1:(n - i)])
    df[[paste0("pc1_lag", i)]] <- c(rep(NA, i), fci_pc1[1:(n - i)])
  }
  df <- na.omit(df)

  if (nrow(df) < 20) return(NA)

  # Estimate
  y_lags <- paste(paste0("y_lag", 1:n_lags), collapse = " + ")
  pc_lags <- paste(paste0("pc1_lag", 1:n_lags), collapse = " + ")
  formula_str <- paste("y ~ fci_pc1 +", y_lags, "+", pc_lags)

  model <- tryCatch({
    lm(as.formula(formula_str), data = df)
  }, error = function(e) NULL)

  if (is.null(model)) return(NA)

  # Forecast
  last_y <- tail(y, n_lags)
  last_pc1 <- tail(fci_pc1, n_lags)

  newdata <- data.frame(fci_pc1 = tail(fci_pc1, 1))
  for (i in 1:n_lags) {
    newdata[[paste0("y_lag", i)]] <- last_y[n_lags - i + 1]
    newdata[[paste0("pc1_lag", i)]] <- last_pc1[n_lags - i + 1]
  }

  pred <- tryCatch({
    predict(model, newdata = newdata)
  }, error = function(e) NA)

  return(pred)
}


################################################################################
# SECTION 6: PART 2B - ROLLING/EXPANDING FORECASTING
################################################################################

cat("================================================================================\n")
cat("RUNNING OUT-OF-SAMPLE FORECASTING\n")
cat("================================================================================\n\n")

#' Run out-of-sample forecasting comparison
#' @param data Full dataset
#' @param target Target variable name
#' @param horizons Forecast horizons
#' @param window_type "expanding" or "rolling"
#' @param min_train Minimum training size
#' @param roll_window Rolling window size (if applicable)
#' @return Data frame with forecasts and errors
run_oos_forecasting <- function(data, target, horizons, window_type = "expanding",
                                 min_train = 80, roll_window = 60) {

  n_total <- nrow(data)
  results <- data.frame()

  # Determine available FCI variables
  fci_vars <- intersect(c("FCI_COMP", "FCI_ENDO", "FCI_EXO"), names(data))
  macro_vars_base <- intersect(c("IMAEP_yoy", "IPC_yoy", "Credit_yoy"), names(data))
  macro_vars_base <- setdiff(macro_vars_base, target)

  for (h in horizons) {
    cat(sprintf("  Horizon %d months (%s)... ", h, window_type))

    # Determine OOS period
    if (window_type == "expanding") {
      start_idx <- min_train + 1
    } else {
      start_idx <- roll_window + 1
    }

    n_forecasts <- 0

    for (t in start_idx:(n_total - h)) {
      # Training window
      if (window_type == "expanding") {
        train_start <- 1
      } else {
        train_start <- t - roll_window
      }
      train_end <- t
      test_idx <- t + h

      if (test_idx > n_total) break
      if (train_end - train_start + 1 < 40) next

      train_data <- data[train_start:train_end, ]
      actual <- data[[target]][test_idx]

      if (is.na(actual)) next

      # Model 1: AR
      pred_ar <- forecast_ar(train_data, target, h, n_lags = 2)

      # Model 2: AR + FCI_COMP
      pred_ar_comp <- forecast_ar_fci(train_data, target, "FCI_COMP", h, n_lags = 2)

      # Model 3: AR + FCI_ENDO
      pred_ar_endo <- if ("FCI_ENDO" %in% names(train_data)) {
        forecast_ar_fci(train_data, target, "FCI_ENDO", h, n_lags = 2)
      } else NA

      # Model 4: AR + FCI_EXO
      pred_ar_exo <- if ("FCI_EXO" %in% names(train_data)) {
        forecast_ar_fci(train_data, target, "FCI_EXO", h, n_lags = 2)
      } else NA

      # Model 5: VAR without FCI
      var_vars_wo <- c(target, macro_vars_base)
      var_vars_wo <- intersect(var_vars_wo, names(train_data))
      pred_var_wo <- if (length(var_vars_wo) >= 2) {
        forecast_var(train_data, target, var_vars_wo, h, p = 2)
      } else NA

      # Model 6: VAR with FCI
      var_vars_w <- c(target, macro_vars_base, "FCI_COMP")
      var_vars_w <- intersect(var_vars_w, names(train_data))
      pred_var_w <- if (length(var_vars_w) >= 3) {
        forecast_var(train_data, target, var_vars_w, h, p = 2)
      } else NA

      # Model 7: Factor-AR
      pred_factor <- if (length(fci_vars) >= 2) {
        forecast_factor_ar(train_data, target, fci_vars, h, n_lags = 2)
      } else NA

      results <- rbind(results, data.frame(
        fecha = data$fecha[test_idx],
        target = target,
        horizon = h,
        window_type = window_type,
        actual = actual,
        pred_AR = pred_ar,
        pred_AR_COMP = pred_ar_comp,
        pred_AR_ENDO = pred_ar_endo,
        pred_AR_EXO = pred_ar_exo,
        pred_VAR_wo = pred_var_wo,
        pred_VAR_w = pred_var_w,
        pred_Factor = pred_factor
      ))

      n_forecasts <- n_forecasts + 1
    }

    cat(n_forecasts, "forecasts\n")
  }

  return(results)
}

# Prepare forecast data (align FCI names)
forecast_data <- analysis_data %>%
  rename(FCI_COMP = FCI_COMP, FCI_ENDO = FCI_ENDO, FCI_EXO = FCI_EXO)

# Run forecasting for each target and window type
forecast_results <- list()

for (target in c("IMAEP_yoy", "IPC_yoy", "Credit_yoy")) {
  if (!target %in% names(forecast_data)) {
    cat("  Skipping", target, "(not available)\n")
    next
  }

  cat("\nTarget:", target, "\n")

  # Endogeneity correction: when forecasting credit, swap FCI_COMP with FCI_exCredit
  fd <- forecast_data
  if (target == "Credit_yoy" && "FCI_exCredit" %in% names(analysis_data)) {
    fd$FCI_COMP <- analysis_data$FCI_exCredit
    cat("  (Using FCI_exCredit for endogeneity correction)\n")
  }

  # Expanding window
  exp_results <- run_oos_forecasting(
    fd, target, MP_CONFIG$horizons_forecast,
    window_type = "expanding", min_train = MP_CONFIG$min_train
  )
  if (nrow(exp_results) > 0) {
    forecast_results[[paste(target, "expanding", sep = "_")]] <- exp_results
  }

  # Rolling window
  roll_results <- run_oos_forecasting(
    fd, target, MP_CONFIG$horizons_forecast,
    window_type = "rolling", roll_window = MP_CONFIG$rolling_window
  )
  if (nrow(roll_results) > 0) {
    forecast_results[[paste(target, "rolling", sep = "_")]] <- roll_results
  }
}

# Combine all forecast results
all_forecasts <- bind_rows(forecast_results)


################################################################################
# SECTION 7: PART 2C - STATISTICAL COMPARISON (DIEBOLD-MARIANO)
################################################################################

cat("\n================================================================================\n")
cat("CALCULATING RMSE AND DIEBOLD-MARIANO TESTS\n")
cat("================================================================================\n\n")

#' Diebold-Mariano test for equal predictive accuracy
#' @param e1 Forecast errors from model 1
#' @param e2 Forecast errors from model 2
#' @param h Forecast horizon
#' @return List with test statistic and p-value
dm_test <- function(e1, e2, h = 1) {
  # Loss differential
  d <- e1^2 - e2^2
  d <- d[!is.na(d)]

  if (length(d) < 20) {
    return(list(statistic = NA, p_value = NA))
  }

  # HAC variance (Newey-West with h-1 lags)
  n <- length(d)
  mean_d <- mean(d)

  # Autocovariance
  gamma <- function(k) {
    if (k >= n) return(0)
    cov(d[1:(n-k)], d[(k+1):n])
  }

  # HAC variance
  var_hac <- gamma(0)
  for (k in 1:min(h-1, n-1)) {
    var_hac <- var_hac + 2 * (1 - k/h) * gamma(k)
  }

  if (var_hac <= 0) {
    var_hac <- var(d)
  }

  # DM statistic
  dm_stat <- mean_d / sqrt(var_hac / n)
  p_value <- 2 * (1 - pnorm(abs(dm_stat)))

  return(list(statistic = dm_stat, p_value = p_value))
}

# Calculate RMSE by model
if (nrow(all_forecasts) > 0) {

  # Calculate errors
  all_forecasts <- all_forecasts %>%
    mutate(
      e_AR = actual - pred_AR,
      e_AR_COMP = actual - pred_AR_COMP,
      e_AR_ENDO = actual - pred_AR_ENDO,
      e_AR_EXO = actual - pred_AR_EXO,
      e_VAR_wo = actual - pred_VAR_wo,
      e_VAR_w = actual - pred_VAR_w,
      e_Factor = actual - pred_Factor
    )

  # RMSE summary
  rmse_summary <- all_forecasts %>%
    group_by(target, horizon, window_type) %>%
    summarise(
      n_forecasts = n(),
      RMSE_AR = sqrt(mean(e_AR^2, na.rm = TRUE)),
      RMSE_AR_COMP = sqrt(mean(e_AR_COMP^2, na.rm = TRUE)),
      RMSE_AR_ENDO = sqrt(mean(e_AR_ENDO^2, na.rm = TRUE)),
      RMSE_AR_EXO = sqrt(mean(e_AR_EXO^2, na.rm = TRUE)),
      RMSE_VAR_wo = sqrt(mean(e_VAR_wo^2, na.rm = TRUE)),
      RMSE_VAR_w = sqrt(mean(e_VAR_w^2, na.rm = TRUE)),
      RMSE_Factor = sqrt(mean(e_Factor^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      # Relative RMSE (vs AR baseline)
      rel_AR_COMP = RMSE_AR_COMP / RMSE_AR,
      rel_AR_ENDO = RMSE_AR_ENDO / RMSE_AR,
      rel_AR_EXO = RMSE_AR_EXO / RMSE_AR,
      rel_VAR_wo = RMSE_VAR_wo / RMSE_AR,
      rel_VAR_w = RMSE_VAR_w / RMSE_AR,
      rel_Factor = RMSE_Factor / RMSE_AR
    )

  cat("RMSE Summary (Expanding Window):\n\n")

  rmse_exp <- rmse_summary %>% filter(window_type == "expanding")

  cat(sprintf("%-12s %3s %8s %8s %8s %8s %8s %8s %8s\n",
              "Target", "h", "AR", "AR+COMP", "AR+ENDO", "AR+EXO", "VAR", "VAR+FCI", "Factor"))
  cat(paste(rep("-", 85), collapse = ""), "\n")

  for (i in 1:nrow(rmse_exp)) {
    r <- rmse_exp[i, ]
    cat(sprintf("%-12s %3d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
                gsub("_yoy", "", r$target), r$horizon,
                r$RMSE_AR, r$RMSE_AR_COMP, r$RMSE_AR_ENDO, r$RMSE_AR_EXO,
                r$RMSE_VAR_wo, r$RMSE_VAR_w, r$RMSE_Factor))
  }

  # Diebold-Mariano tests (AR+COMP vs AR)
  cat("\n\nDiebold-Mariano Tests: AR+FCI_COMP vs AR (Expanding Window)\n")
  cat("H0: Equal predictive accuracy\n\n")

  dm_results <- data.frame()

  for (target in unique(all_forecasts$target)) {
    for (h in MP_CONFIG$horizons_forecast) {
      for (wtype in c("expanding", "rolling")) {
        subset_data <- all_forecasts %>%
          filter(target == !!target, horizon == h, window_type == wtype)

        if (nrow(subset_data) < 20) next

        dm <- dm_test(subset_data$e_AR, subset_data$e_AR_COMP, h)

        dm_results <- rbind(dm_results, data.frame(
          target = target,
          horizon = h,
          window_type = wtype,
          comparison = "AR_COMP vs AR",
          dm_stat = dm$statistic,
          p_value = dm$p_value,
          n_forecasts = nrow(subset_data)
        ))
      }
    }
  }

  if (nrow(dm_results) > 0) {
    dm_results$significant <- dm_results$p_value < 0.10

    cat(sprintf("%-12s %3s %-10s %10s %10s\n",
                "Target", "h", "Window", "DM stat", "p-value"))
    cat(paste(rep("-", 55), collapse = ""), "\n")

    dm_exp <- dm_results %>% filter(window_type == "expanding")
    for (i in 1:nrow(dm_exp)) {
      r <- dm_exp[i, ]
      stars <- ifelse(r$p_value < 0.01, "***",
                      ifelse(r$p_value < 0.05, "**",
                             ifelse(r$p_value < 0.10, "*", "")))
      cat(sprintf("%-12s %3d %-10s %+10.2f %10.3f%s\n",
                  gsub("_yoy", "", r$target), r$horizon, r$window_type,
                  r$dm_stat, r$p_value, stars))
    }
  }

  # Model rankings
  cat("\n\nMODEL RANKINGS (by Relative RMSE, lower is better):\n\n")

  model_rankings <- rmse_summary %>%
    filter(window_type == "expanding") %>%
    dplyr::select(target, horizon, rel_AR_COMP, rel_AR_ENDO, rel_AR_EXO,
                  rel_VAR_wo, rel_VAR_w, rel_Factor) %>%
    pivot_longer(cols = starts_with("rel_"),
                 names_to = "model", values_to = "rel_rmse") %>%
    mutate(model = gsub("rel_", "", model)) %>%
    group_by(target, horizon) %>%
    arrange(rel_rmse) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  # Best model for each target/horizon
  best_models <- model_rankings %>%
    filter(rank == 1) %>%
    dplyr::select(target, horizon, model, rel_rmse)

  cat(sprintf("%-12s %3s %-12s %10s\n", "Target", "h", "Best Model", "Rel.RMSE"))
  cat(paste(rep("-", 45), collapse = ""), "\n")

  for (i in 1:nrow(best_models)) {
    r <- best_models[i, ]
    cat(sprintf("%-12s %3d %-12s %10.3f\n",
                gsub("_yoy", "", r$target), r$horizon, r$model, r$rel_rmse))
  }
}


################################################################################
# SECTION 8: VISUALIZATIONS
################################################################################

cat("\n================================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

if (!dir.exists(MP_CONFIG$output_dir)) {
  dir.create(MP_CONFIG$output_dir, recursive = TRUE)
}

# Color palette
colors_models <- c(
  "AR" = "#7F8C8D",
  "AR_COMP" = "#E74C3C",
  "AR_ENDO" = "#3498DB",
  "AR_EXO" = "#27AE60",
  "VAR_wo" = "#9B59B6",
  "VAR_w" = "#F39C12",
  "Factor" = "#1ABC9C"
)

colors_stance <- c(
  "Low TPM" = "#3498DB",
  "High TPM" = "#E74C3C"
)

# Plot 90: TPM x FCI Interaction IRFs
if (nrow(all_amplification) > 0) {

  amp_plot_data <- all_amplification %>%
    filter(fci_type == "FCI_COMP") %>%
    mutate(target_label = gsub("_yoy", "", target_var))

  p90 <- ggplot(amp_plot_data, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_interact_lo, ymax = ci_interact_hi),
                fill = "#E74C3C", alpha = 0.2) +
    geom_line(aes(y = coef_interact), color = "#E74C3C", linewidth = 1) +
    geom_point(aes(y = coef_interact), color = "#E74C3C", size = 2) +
    facet_wrap(~target_label, scales = "free_y", ncol = 3) +
    theme_minimal(base_size = 11) +
    labs(title = "FCI x TPM Interaction: Amplification Effect",
         subtitle = "Coefficient on dTPM x FCI interaction | 90% CI",
         x = "Months", y = "Interaction coefficient (b3)") +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(MP_CONFIG$output_dir, "90_TPM_FCI_Interaction_IRF.png"), p90,
         width = 12, height = 4, dpi = 300)
  cat("Saved: 90_TPM_FCI_Interaction_IRF.png\n")
}

# Plot 91: Amplification Heatmap
if (nrow(all_amplification) > 0) {

  heatmap_amp <- all_amplification %>%
    filter(fci_type == "FCI_COMP", horizon %in% MP_CONFIG$horizons_interaction) %>%
    mutate(
      target_label = gsub("_yoy", "", target_var),
      stars = case_when(
        p_interact < 0.01 ~ "***",
        p_interact < 0.05 ~ "**",
        p_interact < 0.10 ~ "*",
        TRUE ~ ""
      ),
      label = sprintf("%+.3f%s", coef_interact, stars)
    )

  p91 <- ggplot(heatmap_amp, aes(x = factor(horizon), y = target_label, fill = coef_interact)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                         midpoint = 0, name = "Coefficient") +
    theme_minimal(base_size = 12) +
    labs(title = "FCI Amplification of Monetary Policy",
         subtitle = "Interaction coefficient (dTPM x FCI) | *** p<0.01, ** p<0.05, * p<0.10",
         x = "Horizon (months)", y = NULL) +
    theme(plot.title = element_text(face = "bold", size = 14))

  ggsave(file.path(MP_CONFIG$output_dir, "91_TPM_FCI_Amplification_Heatmap.png"), p91,
         width = 10, height = 5, dpi = 300)
  cat("Saved: 91_TPM_FCI_Amplification_Heatmap.png\n")
}

# Plot 92: FCI Stance Transmission
if (nrow(all_stance) > 0) {

  # Focus on binary stance for clearest interpretation
  stance_binary <- all_stance %>%
    filter(stance_type == "TPM_binary") %>%
    dplyr::select(horizon, target_var, coef_fci_low, ci_fci_low_lo, ci_fci_low_hi,
                  coef_fci_high, ci_fci_high_lo, ci_fci_high_hi) %>%
    pivot_longer(cols = c(coef_fci_low, coef_fci_high),
                 names_to = "stance", values_to = "coef") %>%
    mutate(
      stance_label = ifelse(stance == "coef_fci_low", "Low TPM", "High TPM"),
      ci_lo = ifelse(stance == "coef_fci_low", ci_fci_low_lo, ci_fci_high_lo),
      ci_hi = ifelse(stance == "coef_fci_low", ci_fci_low_hi, ci_fci_high_hi),
      target_label = gsub("_yoy", "", target_var)
    )

  p92 <- ggplot(stance_binary, aes(x = horizon, y = coef, color = stance_label, fill = stance_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~target_label, scales = "free_y", ncol = 3) +
    scale_color_manual(values = colors_stance) +
    scale_fill_manual(values = colors_stance) +
    theme_minimal(base_size = 11) +
    labs(title = "FCI Effect by Monetary Policy Stance",
         subtitle = "FCI effect in high vs low TPM periods (binary stance) | 90% CI",
         x = "Months", y = "FCI effect (pp)",
         color = "TPM Stance", fill = "TPM Stance") +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(MP_CONFIG$output_dir, "92_FCI_Stance_Transmission.png"), p92,
         width = 12, height = 4, dpi = 300)
  cat("Saved: 92_FCI_Stance_Transmission.png\n")
}

# Plot 93: Stance Comparison (all 3 definitions)
if (nrow(all_stance) > 0) {

  stance_compare_12 <- all_stance %>%
    filter(horizon == 12) %>%
    mutate(target_label = gsub("_yoy", "", target_var))

  p93 <- ggplot(stance_compare_12, aes(x = stance_label, y = coef_interact, fill = target_label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, color = "gray30") +
    geom_errorbar(aes(ymin = coef_interact - 1.645 * se_interact,
                      ymax = coef_interact + 1.645 * se_interact),
                  position = position_dodge(width = 0.8), width = 0.2) +
    scale_fill_manual(values = c("IMAEP" = "#3498DB", "IPC" = "#E74C3C", "Credit" = "#27AE60")) +
    theme_minimal(base_size = 11) +
    labs(title = "Stance Interaction Comparison (h=12)",
         subtitle = "Comparing three TPM stance definitions | 90% CI",
         x = "Stance Definition", y = "Interaction coefficient",
         fill = "Target") +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 15, hjust = 1),
          legend.position = "bottom")

  ggsave(file.path(MP_CONFIG$output_dir, "93_Stance_Comparison.png"), p93,
         width = 10, height = 6, dpi = 300)
  cat("Saved: 93_Stance_Comparison.png\n")
}

# Plot 94: RMSE Comparison (Expanding)
if (exists("rmse_summary") && nrow(rmse_summary) > 0) {

  rmse_long <- rmse_summary %>%
    filter(window_type == "expanding") %>%
    dplyr::select(target, horizon, RMSE_AR, RMSE_AR_COMP, RMSE_AR_ENDO,
                  RMSE_AR_EXO, RMSE_VAR_wo, RMSE_VAR_w, RMSE_Factor) %>%
    pivot_longer(cols = starts_with("RMSE_"),
                 names_to = "model", values_to = "RMSE") %>%
    mutate(
      model = gsub("RMSE_", "", model),
      target_label = gsub("_yoy", "", target)
    )

  p94 <- ggplot(rmse_long, aes(x = factor(horizon), y = RMSE, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~target_label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = colors_models) +
    theme_minimal(base_size = 11) +
    labs(title = "Out-of-Sample RMSE Comparison (Expanding Window)",
         subtitle = "Lower is better",
         x = "Horizon (months)", y = "RMSE",
         fill = "Model") +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(MP_CONFIG$output_dir, "94_RMSE_Comparison_Expanding.png"), p94,
         width = 14, height = 5, dpi = 300)
  cat("Saved: 94_RMSE_Comparison_Expanding.png\n")
}

# Plot 95: RMSE Comparison (Rolling)
if (exists("rmse_summary") && nrow(rmse_summary) > 0) {

  rmse_long_roll <- rmse_summary %>%
    filter(window_type == "rolling") %>%
    dplyr::select(target, horizon, RMSE_AR, RMSE_AR_COMP, RMSE_AR_ENDO,
                  RMSE_AR_EXO, RMSE_VAR_wo, RMSE_VAR_w, RMSE_Factor) %>%
    pivot_longer(cols = starts_with("RMSE_"),
                 names_to = "model", values_to = "RMSE") %>%
    mutate(
      model = gsub("RMSE_", "", model),
      target_label = gsub("_yoy", "", target)
    )

  p95 <- ggplot(rmse_long_roll, aes(x = factor(horizon), y = RMSE, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~target_label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = colors_models) +
    theme_minimal(base_size = 11) +
    labs(title = "Out-of-Sample RMSE Comparison (Rolling Window, 60 months)",
         subtitle = "Lower is better",
         x = "Horizon (months)", y = "RMSE",
         fill = "Model") +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(MP_CONFIG$output_dir, "95_RMSE_Comparison_Rolling.png"), p95,
         width = 14, height = 5, dpi = 300)
  cat("Saved: 95_RMSE_Comparison_Rolling.png\n")
}

# Plot 96: DM Test Heatmap
if (exists("dm_results") && nrow(dm_results) > 0) {

  dm_heatmap <- dm_results %>%
    filter(window_type == "expanding") %>%
    mutate(
      target_label = gsub("_yoy", "", target),
      stars = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.10 ~ "*",
        TRUE ~ ""
      ),
      label = sprintf("%+.2f%s", dm_stat, stars)
    )

  p96 <- ggplot(dm_heatmap, aes(x = factor(horizon), y = target_label, fill = -dm_stat)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "DM stat\n(neg = FCI better)") +
    theme_minimal(base_size = 12) +
    labs(title = "Diebold-Mariano Test: AR+FCI_COMP vs AR",
         subtitle = "Negative DM = FCI model better | *** p<0.01, ** p<0.05, * p<0.10",
         x = "Horizon (months)", y = NULL) +
    theme(plot.title = element_text(face = "bold", size = 14))

  ggsave(file.path(MP_CONFIG$output_dir, "96_DM_Test_Heatmap.png"), p96,
         width = 10, height = 5, dpi = 300)
  cat("Saved: 96_DM_Test_Heatmap.png\n")
}

# Plot 97: Cumulative Squared Errors
if (nrow(all_forecasts) > 0) {

  cse_data <- all_forecasts %>%
    filter(window_type == "expanding", horizon == 12) %>%
    arrange(fecha) %>%
    group_by(target) %>%
    mutate(
      CSE_AR = cumsum(e_AR^2),
      CSE_AR_COMP = cumsum(e_AR_COMP^2),
      target_label = gsub("_yoy", "", target)
    ) %>%
    ungroup() %>%
    dplyr::select(fecha, target_label, CSE_AR, CSE_AR_COMP) %>%
    pivot_longer(cols = c(CSE_AR, CSE_AR_COMP),
                 names_to = "model", values_to = "CSE") %>%
    mutate(model = ifelse(model == "CSE_AR", "AR", "AR + FCI"))

  p97 <- ggplot(cse_data, aes(x = fecha, y = CSE, color = model)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~target_label, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("AR" = "#7F8C8D", "AR + FCI" = "#E74C3C")) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    theme_minimal(base_size = 11) +
    labs(title = "Cumulative Squared Forecast Errors (h=12, Expanding)",
         subtitle = "Lower line = better forecasting performance",
         x = NULL, y = "Cumulative Squared Error",
         color = "Model") +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(MP_CONFIG$output_dir, "97_Forecast_Cumulative_Error.png"), p97,
         width = 12, height = 4, dpi = 300)
  cat("Saved: 97_Forecast_Cumulative_Error.png\n")
}

# Plot 98: Model Ranking Summary
if (exists("model_rankings") && nrow(model_rankings) > 0) {

  rank_summary <- model_rankings %>%
    group_by(model) %>%
    summarise(
      avg_rank = mean(rank),
      avg_rel_rmse = mean(rel_rmse, na.rm = TRUE),
      n_best = sum(rank == 1),
      .groups = "drop"
    ) %>%
    arrange(avg_rank)

  p98 <- ggplot(rank_summary, aes(x = reorder(model, avg_rank), y = avg_rank, fill = model)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", avg_rank)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = colors_models) +
    theme_minimal(base_size = 12) +
    labs(title = "Model Rankings by Average Rank",
         subtitle = "Lower is better (1 = best)",
         x = "Model", y = "Average Rank") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "none") +
    coord_cartesian(ylim = c(0, max(rank_summary$avg_rank) * 1.2))

  ggsave(file.path(MP_CONFIG$output_dir, "98_Model_Ranking_Summary.png"), p98,
         width = 10, height = 6, dpi = 300)
  cat("Saved: 98_Model_Ranking_Summary.png\n")
}

# Plot 99: Summary Dashboard
p99_a <- if (exists("p91")) {
  p91 + theme(legend.position = "none") + labs(title = "A: Amplification")
} else {
  ggplot() + theme_void() + labs(title = "A: Amplification (N/A)")
}

p99_b <- if (exists("p92")) {
  p92 + theme(legend.position = "none") + labs(title = "B: Stance Effect")
} else {
  ggplot() + theme_void() + labs(title = "B: Stance Effect (N/A)")
}

p99_c <- if (exists("p94")) {
  p94 + theme(legend.position = "none") + labs(title = "C: RMSE Comparison")
} else {
  ggplot() + theme_void() + labs(title = "C: RMSE (N/A)")
}

p99_d <- if (exists("p96")) {
  p96 + theme(legend.position = "none") + labs(title = "D: DM Tests")
} else {
  ggplot() + theme_void() + labs(title = "D: DM Tests (N/A)")
}

p99 <- grid.arrange(
  p99_a, p99_b, p99_c, p99_d,
  ncol = 2,
  top = grid::textGrob("Monetary Policy Interaction and Forecasting Summary",
                       gp = grid::gpar(fontsize = 16, fontface = "bold"))
)

ggsave(file.path(MP_CONFIG$output_dir, "99_Monetary_Policy_Dashboard.png"), p99,
       width = 14, height = 12, dpi = 300)
cat("Saved: 99_Monetary_Policy_Dashboard.png\n")


################################################################################
# SECTION 9: EXPORT RESULTS
################################################################################

cat("\n================================================================================\n")
cat("EXPORTING RESULTS\n")
cat("================================================================================\n\n")

# Export amplification results
if (nrow(all_amplification) > 0) {
  write.csv(all_amplification,
            file.path(MP_CONFIG$output_dir, "FCI_TPM_Interaction_Amplification.csv"),
            row.names = FALSE)
  cat("Saved: FCI_TPM_Interaction_Amplification.csv\n")
}

# Export stance results
if (nrow(all_stance) > 0) {
  write.csv(all_stance,
            file.path(MP_CONFIG$output_dir, "FCI_TPM_Interaction_Stance.csv"),
            row.names = FALSE)
  cat("Saved: FCI_TPM_Interaction_Stance.csv\n")
}

# Export RMSE comparison
if (exists("rmse_summary") && nrow(rmse_summary) > 0) {
  write.csv(rmse_summary,
            file.path(MP_CONFIG$output_dir, "Forecast_RMSE_Comparison.csv"),
            row.names = FALSE)
  cat("Saved: Forecast_RMSE_Comparison.csv\n")
}

# Export DM test results
if (exists("dm_results") && nrow(dm_results) > 0) {
  write.csv(dm_results,
            file.path(MP_CONFIG$output_dir, "Forecast_DM_Tests.csv"),
            row.names = FALSE)
  cat("Saved: Forecast_DM_Tests.csv\n")
}

# Export model rankings
if (exists("model_rankings") && nrow(model_rankings) > 0) {
  write.csv(model_rankings,
            file.path(MP_CONFIG$output_dir, "Forecast_Model_Rankings.csv"),
            row.names = FALSE)
  cat("Saved: Forecast_Model_Rankings.csv\n")
}


################################################################################
# SECTION 10: SUMMARY AND KEY FINDINGS
################################################################################

cat("\n================================================================================\n")
cat("KEY FINDINGS SUMMARY\n")
cat("================================================================================\n\n")

cat("PART 1: FCI x MONETARY POLICY INTERACTION\n\n")

cat("1A. FCI AMPLIFICATION OF MONETARY POLICY:\n")
if (nrow(all_amplification) > 0) {
  amp_12 <- all_amplification %>%
    filter(fci_type == "FCI_COMP", horizon == 12)

  for (i in 1:nrow(amp_12)) {
    r <- amp_12[i, ]
    target_label <- gsub("_yoy", "", r$target_var)

    if (r$p_interact < 0.10) {
      direction <- ifelse(r$coef_interact > 0, "AMPLIFIES", "DAMPENS")
      cat(sprintf("   %s: FCI %s monetary policy (b3=%+.3f, p=%.3f)\n",
                  target_label, direction, r$coef_interact, r$p_interact))
    } else {
      cat(sprintf("   %s: No significant interaction (b3=%+.3f, p=%.3f)\n",
                  target_label, r$coef_interact, r$p_interact))
    }
  }
}

cat("\n1B. TPM STANCE EFFECT ON FCI TRANSMISSION:\n")
if (nrow(all_stance) > 0) {
  stance_12 <- all_stance %>%
    filter(stance_type == "TPM_binary", horizon == 12)

  for (i in 1:nrow(stance_12)) {
    r <- stance_12[i, ]
    target_label <- gsub("_yoy", "", r$target_var)

    if (r$p_interact < 0.10) {
      direction <- ifelse(abs(r$coef_fci_high) > abs(r$coef_fci_low), "STRONGER", "WEAKER")
      cat(sprintf("   %s: FCI effect %s in high TPM periods\n", target_label, direction))
    } else {
      cat(sprintf("   %s: FCI effect similar across TPM stances\n", target_label))
    }
  }
}

cat("\nPART 2: FORECASTING COMPARISON\n\n")
if (exists("best_models") && nrow(best_models) > 0) {
  cat("Best performing models (by Relative RMSE):\n")
  for (i in 1:nrow(best_models)) {
    r <- best_models[i, ]
    cat(sprintf("   %s (h=%d): %s (Rel.RMSE=%.3f)\n",
                gsub("_yoy", "", r$target), r$horizon, r$model, r$rel_rmse))
  }
}

cat("\nENDOGENEITY NOTE:\n")
cat("   TPM is a component of FCI_RATES and FCI_ENDO (1 of 12 in FCI_COMP)\n")
cat("   Robustness check with FCI_BANKING (excludes TPM) available in results\n")

cat("\n================================================================================\n")
cat("MONETARY POLICY INTERACTION ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
cat("Plots: 90-99\n")
cat("Output:", MP_CONFIG$output_dir, "\n\n")


################################################################################
# END OF SCRIPT
################################################################################
