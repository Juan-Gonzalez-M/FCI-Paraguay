################################################################################
# FCI LOCAL PROJECTIONS ANALYSIS - COMPREHENSIVE
################################################################################
#
# Project:      Financial Conditions Index - Local Projections (Jordà, 2005)
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Comprehensive Local Projections analysis including:
#
#   PART A: MACRO EFFECTS
#     - FCI effects on IMAEP, IPC, Total Credit
#     - Standard and asymmetric (non-linear) LP
#
#   PART B: CREDIT CHANNEL DECOMPOSITION
#     - 3 FCI types: Comprehensive, Endogenous, Exogenous
#     - 3 Credit types: Total, MN (local currency), USD (dollarized)
#     - Tests differential effects of domestic vs external conditions
#
# Methodology:
#   - Standard LP: y_{t+h} = α + β×FCI_t + γ×X_t + ε
#   - Asymmetric LP: y_{t+h} = α + β⁺×FCI⁺_t + β⁻×|FCI⁻|_t + γ×X_t + ε
#   - Newey-West HAC standard errors
#
# References:
#   - Jordà (2005) AER
#   - Ramey & Zubairy (2018) JPE
#
################################################################################


################################################################################
# 1. SETUP AND CONFIGURATION
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(sandwich)
  library(lmtest)
})

LP_CONFIG <- list(
  # Data
  data_file = "../data/FCI_data_1.xlsx",
  macro_sheet = "Datos_macro",

  # LP parameters
  max_horizon = 24,
  horizons_report = c(3, 6, 12, 18, 24),
  n_lags = 2,
  confidence_level = 0.90,

  # Analysis settings
  macro_variables = c("IMAEP", "IPC"),
  fci_types = c("FCI_COMP", "FCI_ENDO", "FCI_EXO"),
  credit_types = c("Total", "USD", "Real_MN"),

  # Output
  output_dir = "../output"
)

# Labels for display
FCI_LABELS <- c(
  "FCI_COMP" = "Comprehensive",
  "FCI_ENDO" = "Endogenous (Domestic)",
  "FCI_EXO" = "Exogenous (External)"
)

CREDIT_LABELS <- c(
  "Total" = "Real Total Credit",
  "USD" = "USD (Dollarized)",
  "Real_MN" = "Real MN Credit"
)

set.seed(20250116)

cat("\n################################################################################\n")
cat("FCI LOCAL PROJECTIONS - COMPREHENSIVE ANALYSIS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# 2. DATA PREPARATION
################################################################################

cat("================================================================================\n")
cat("DATA PREPARATION\n")
cat("================================================================================\n\n")

# Load FCI
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  source("01_FCI_Complete.R")
}

# Extract all FCI types (including endogeneity-corrected versions)
fci_data <- resultado_fci$all_indices %>%
  dplyr::select(fecha,
                FCI_COMP = FCI_COMP_AVG,
                FCI_ENDO = FCI_ENDO_AVG,
                FCI_EXO = FCI_EXO_AVG,
                FCI_exCredit = FCI_exCredit_AVG,
                FCI_ENDO_exCredit = FCI_ENDO_exCredit_AVG,
                FCI_RATES = FCI_RATES_AVG,
                FCI_ENDO_exCredit_exTCN = FCI_ENDO_exCredit_exTCN_AVG)

# Also extract individual method FCIs for method-by-method LP
fci_methods <- resultado_fci$all_indices %>%
  dplyr::select(fecha,
                matches("^FCI_exCredit_(ZS|PCA|VAR|DFM)_norm$"))

# Load raw data for credit and macro variables
datos_raw <- read_excel(LP_CONFIG$data_file)
fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

# Load macro data (needed for IPC deflation of credit)
macro_raw <- tryCatch({
  read_excel(LP_CONFIG$data_file, sheet = LP_CONFIG$macro_sheet)
}, error = function(e) {
  cat("Note: Macro sheet not found, using main sheet\n")
  datos_raw
})

fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

# Prepare credit data (YoY growth from stock)
# Use deflated credit for Total, add Real MN
credit_data <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Get IPC from macro sheet for deflation
ipc_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  dplyr::select(fecha, IPC, Creditos_deflactados) %>%
  arrange(fecha)

credit_data <- credit_data %>%
  left_join(ipc_data, by = "fecha") %>%
  mutate(
    Cred_Total = (Creditos_deflactados /
                  lag(Creditos_deflactados, 12) - 1) * 100,
    Cred_MN = (Creditos_Sector_privado_MN /
               lag(Creditos_Sector_privado_MN, 12) - 1) * 100,
    Cred_USD = (Creditos_Sector_privado_USD_equivalente /
                lag(Creditos_Sector_privado_USD_equivalente, 12) - 1) * 100,
    Cred_Real_MN = ((Creditos_Sector_privado_MN / IPC) /
                     lag(Creditos_Sector_privado_MN / IPC, 12) - 1) * 100
  ) %>%
  dplyr::select(fecha, Cred_Total, Cred_MN, Cred_USD, Cred_Real_MN)

macro_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Calculate YoY growth for macro variables
for (var in LP_CONFIG$macro_variables) {
  if (var %in% names(macro_data)) {
    macro_data <- macro_data %>%
      mutate(!!paste0(var, "_yoy") := (!!sym(var) / lag(!!sym(var), 12) - 1) * 100)
  }
}

# Also compute Cred_Real_yoy (real total credit growth) for use as LP control
macro_data <- macro_data %>%
  left_join(credit_data %>% dplyr::select(fecha, Cred_Real_yoy = Cred_Total), by = "fecha") %>%
  dplyr::select(fecha, ends_with("_yoy"))

# Extract external variables (instruments for IV-LP)
external_instruments <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  dplyr::select(fecha, VIX, FFER, Commodities)

# Merge all data
analysis_data <- fci_data %>%
  left_join(fci_methods, by = "fecha") %>%
  inner_join(credit_data, by = "fecha") %>%
  inner_join(macro_data, by = "fecha") %>%
  inner_join(external_instruments, by = "fecha") %>%
  arrange(fecha) %>%
  mutate(
    # Create asymmetric FCI components for each type
    FCI_COMP_pos = pmax(FCI_COMP, 0),
    FCI_COMP_neg = abs(pmin(FCI_COMP, 0)),
    FCI_ENDO_pos = pmax(FCI_ENDO, 0),
    FCI_ENDO_neg = abs(pmin(FCI_ENDO, 0)),
    FCI_EXO_pos = pmax(FCI_EXO, 0),
    FCI_EXO_neg = abs(pmin(FCI_EXO, 0)),
    FCI_exCredit_pos = pmax(FCI_exCredit, 0),
    FCI_exCredit_neg = abs(pmin(FCI_exCredit, 0)),
    FCI_ENDO_exCredit_pos = pmax(FCI_ENDO_exCredit, 0),
    FCI_ENDO_exCredit_neg = abs(pmin(FCI_ENDO_exCredit, 0))
  ) %>%
  na.omit()

cat("Data prepared:\n")
cat("  Observations:", nrow(analysis_data), "\n")
cat("  Period:", format(min(analysis_data$fecha)), "to",
    format(max(analysis_data$fecha)), "\n\n")

cat("Credit Growth Statistics (YoY %):\n")
cat(sprintf("  Real Total: Mean=%.1f, SD=%.1f\n",
            mean(analysis_data$Cred_Total), sd(analysis_data$Cred_Total)))
cat(sprintf("  USD:        Mean=%.1f, SD=%.1f\n",
            mean(analysis_data$Cred_USD), sd(analysis_data$Cred_USD)))
cat(sprintf("  Real MN:    Mean=%.1f, SD=%.1f\n\n",
            mean(analysis_data$Cred_Real_MN), sd(analysis_data$Cred_Real_MN)))


################################################################################
# 3. LOCAL PROJECTION FUNCTIONS
################################################################################

#' Standard Local Projection
run_lp_standard <- function(data, y_var, fci_var, max_h, n_lags = 2, control_vars = NULL) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - LP_CONFIG$confidence_level) / 2)

  for (h in 1:max_h) {
    # Create variables
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        fci_lag1 = lag(!!sym(fci_var), 1)
      )

    # Add more lags if needed
    for (i in 2:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", i) := lag(!!sym(y_var), i))
    }

    # Build formula
    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")

    # Add control variables and their lags
    control_cols <- c()
    if (!is.null(control_vars)) {
      for (cv in control_vars) {
        if (cv %in% names(data)) {
          data_h <- data_h %>%
            mutate(!!cv := !!sym(cv))
          control_cols <- c(control_cols, cv)
          for (j in 1:n_lags) {
            cv_lag_name <- paste0(cv, "_lag", j)
            data_h <- data_h %>%
              mutate(!!cv_lag_name := lag(!!sym(cv), j))
            control_cols <- c(control_cols, cv_lag_name)
          }
        }
      }
    }

    formula_str <- paste("y_fwd ~", fci_var, "+", paste(c(lag_vars, control_cols), collapse = " + "))

    # Estimate
    reg_data <- data_h %>% dplyr::select(y_fwd, !!sym(fci_var), all_of(lag_vars), all_of(control_cols)) %>% na.omit()

    if (nrow(reg_data) < 30) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    idx <- which(rownames(coef_test) == fci_var)

    results <- rbind(results, data.frame(
      horizon = h,
      coef = coef_test[idx, 1],
      se = coef_test[idx, 2],
      p_value = coef_test[idx, 4],
      ci_lower = coef_test[idx, 1] - z_crit * coef_test[idx, 2],
      ci_upper = coef_test[idx, 1] + z_crit * coef_test[idx, 2],
      n_obs = nrow(reg_data)
    ))
  }

  return(results)
}

#' Asymmetric (Non-linear) Local Projection
run_lp_asymmetric <- function(data, y_var, fci_var, max_h, n_lags = 2, control_vars = NULL) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - LP_CONFIG$confidence_level) / 2)

  fci_pos <- paste0(fci_var, "_pos")
  fci_neg <- paste0(fci_var, "_neg")

  for (h in 1:max_h) {
    # Create variables
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        fci_lag1 = lag(!!sym(fci_var), 1)
      )

    for (i in 2:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", i) := lag(!!sym(y_var), i))
    }

    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")

    # Add control variables and their lags
    control_cols <- c()
    if (!is.null(control_vars)) {
      for (cv in control_vars) {
        if (cv %in% names(data)) {
          data_h <- data_h %>%
            mutate(!!cv := !!sym(cv))
          control_cols <- c(control_cols, cv)
          for (j in 1:n_lags) {
            cv_lag_name <- paste0(cv, "_lag", j)
            data_h <- data_h %>%
              mutate(!!cv_lag_name := lag(!!sym(cv), j))
            control_cols <- c(control_cols, cv_lag_name)
          }
        }
      }
    }

    formula_str <- paste("y_fwd ~", fci_pos, "+", fci_neg, "+",
                         paste(c(lag_vars, control_cols), collapse = " + "))

    reg_data <- data_h %>%
      dplyr::select(y_fwd, !!sym(fci_pos), !!sym(fci_neg), all_of(lag_vars), all_of(control_cols)) %>%
      na.omit()

    if (nrow(reg_data) < 30) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    idx_pos <- which(rownames(coef_test) == fci_pos)
    idx_neg <- which(rownames(coef_test) == fci_neg)

    # Wald test for asymmetry
    coef_diff <- coef_test[idx_pos, 1] - coef_test[idx_neg, 1]
    var_diff <- vcov_hac[idx_pos, idx_pos] + vcov_hac[idx_neg, idx_neg] -
                2 * vcov_hac[idx_pos, idx_neg]
    wald_p <- 1 - pchisq(coef_diff^2 / var_diff, df = 1)

    results <- rbind(results, data.frame(
      horizon = h,
      coef_tight = coef_test[idx_pos, 1],
      se_tight = coef_test[idx_pos, 2],
      p_tight = coef_test[idx_pos, 4],
      ci_tight_lo = coef_test[idx_pos, 1] - z_crit * coef_test[idx_pos, 2],
      ci_tight_hi = coef_test[idx_pos, 1] + z_crit * coef_test[idx_pos, 2],
      coef_ease = coef_test[idx_neg, 1],
      se_ease = coef_test[idx_neg, 2],
      p_ease = coef_test[idx_neg, 4],
      ci_ease_lo = coef_test[idx_neg, 1] - z_crit * coef_test[idx_neg, 2],
      ci_ease_hi = coef_test[idx_neg, 1] + z_crit * coef_test[idx_neg, 2],
      asym_p = wald_p,
      n_obs = nrow(reg_data)
    ))
  }

  return(results)
}


#' IV Local Projection (2SLS with external instruments)
run_lp_iv <- function(data, y_var, endo_var, instruments, max_h, n_lags = 2, control_vars = NULL) {

  if (!requireNamespace("ivreg", quietly = TRUE)) {
    stop("Package 'ivreg' required for IV-LP. Install with install.packages('ivreg')")
  }

  coefs <- data.frame()
  diag <- data.frame()
  z_crit <- qnorm(1 - (1 - LP_CONFIG$confidence_level) / 2)

  for (h in 1:max_h) {
    # Create variables
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        fci_lag1 = lag(!!sym(endo_var), 1)
      )

    for (i in 2:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", i) := lag(!!sym(y_var), i))
    }

    # Exogenous regressors (lags + controls)
    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")

    control_cols <- c()
    if (!is.null(control_vars)) {
      for (cv in control_vars) {
        if (cv %in% names(data)) {
          data_h <- data_h %>%
            mutate(!!cv := !!sym(cv))
          control_cols <- c(control_cols, cv)
          for (j in 1:n_lags) {
            cv_lag_name <- paste0(cv, "_lag", j)
            data_h <- data_h %>%
              mutate(!!cv_lag_name := lag(!!sym(cv), j))
            control_cols <- c(control_cols, cv_lag_name)
          }
        }
      }
    }

    exog_vars <- c(lag_vars, control_cols)

    # Select complete cases
    all_vars <- c("y_fwd", endo_var, instruments, exog_vars)
    reg_data <- data_h %>% dplyr::select(all_of(all_vars)) %>% na.omit()

    if (nrow(reg_data) < 30) next

    # IV formula: y_fwd ~ endo_var + exog | instruments + exog
    exog_str <- paste(exog_vars, collapse = " + ")
    iv_formula <- as.formula(paste(
      "y_fwd ~", endo_var, "+", exog_str,
      "|", paste(instruments, collapse = " + "), "+", exog_str
    ))

    # Estimate 2SLS
    iv_model <- tryCatch(
      ivreg::ivreg(iv_formula, data = reg_data),
      error = function(e) NULL
    )
    if (is.null(iv_model)) next

    # HAC standard errors
    vcov_hac <- tryCatch(
      sandwich::NeweyWest(iv_model, lag = h + 1, prewhite = FALSE),
      error = function(e) sandwich::vcovHC(iv_model, type = "HC1")
    )

    coef_test <- lmtest::coeftest(iv_model, vcov = vcov_hac)
    idx <- which(rownames(coef_test) == endo_var)
    if (length(idx) == 0) next

    coefs <- rbind(coefs, data.frame(
      horizon = h,
      coef = coef_test[idx, 1],
      se = coef_test[idx, 2],
      p_value = coef_test[idx, 4],
      ci_lower = coef_test[idx, 1] - z_crit * coef_test[idx, 2],
      ci_upper = coef_test[idx, 1] + z_crit * coef_test[idx, 2],
      n_obs = nrow(reg_data)
    ))

    # Diagnostics from ivreg summary
    iv_summ <- tryCatch(
      summary(iv_model, vcov = vcov_hac, diagnostics = TRUE),
      error = function(e) NULL
    )

    first_F <- first_p <- wu_h <- wu_p <- sargan_s <- sargan_p_val <- NA

    if (!is.null(iv_summ) && !is.null(iv_summ$diagnostics)) {
      d <- iv_summ$diagnostics
      if ("Weak instruments" %in% rownames(d)) {
        first_F <- d["Weak instruments", "statistic"]
        first_p <- d["Weak instruments", "p-value"]
      }
      if ("Wu-Hausman" %in% rownames(d)) {
        wu_h <- d["Wu-Hausman", "statistic"]
        wu_p <- d["Wu-Hausman", "p-value"]
      }
      if ("Sargan" %in% rownames(d)) {
        sargan_s <- d["Sargan", "statistic"]
        sargan_p_val <- d["Sargan", "p-value"]
      }
    }

    # Partial R² from first stage
    partial_r2 <- NA
    tryCatch({
      fs_full_formula <- as.formula(paste(endo_var, "~",
        paste(c(instruments, exog_vars), collapse = " + ")))
      fs_restricted_formula <- as.formula(paste(endo_var, "~",
        paste(exog_vars, collapse = " + ")))
      fs_full <- lm(fs_full_formula, data = reg_data)
      fs_restricted <- lm(fs_restricted_formula, data = reg_data)
      partial_r2 <- summary(fs_full)$r.squared - summary(fs_restricted)$r.squared
    }, error = function(e) NULL)

    diag <- rbind(diag, data.frame(
      horizon = h,
      first_stage_F = first_F,
      first_stage_p = first_p,
      partial_R2 = partial_r2,
      wu_hausman = wu_h,
      wu_hausman_p = wu_p,
      sargan_stat = sargan_s,
      sargan_p = sargan_p_val,
      n_obs = nrow(reg_data)
    ))
  }

  return(list(coefs = coefs, diag = diag))
}


################################################################################
# 4. PART A: MACRO EFFECTS (Standard Analysis)
################################################################################

cat("================================================================================\n")
cat("PART A: MACRO EFFECTS\n")
cat("================================================================================\n\n")

# Run LP for macro variables using comprehensive FCI
macro_results_std <- list()
macro_results_asym <- list()

for (var in LP_CONFIG$macro_variables) {
  y_var <- paste0(var, "_yoy")
  if (!y_var %in% names(analysis_data)) next

  cat("Processing:", var, "... ")

  # Determine controls based on LHS variable
  if (y_var == "IMAEP_yoy") {
    ctrl <- c("IPC_yoy", "Cred_Real_yoy")
  } else if (y_var == "IPC_yoy") {
    ctrl <- c("IMAEP_yoy", "Cred_Real_yoy")
  } else {
    ctrl <- NULL
  }

  macro_results_std[[var]] <- run_lp_standard(
    analysis_data, y_var, "FCI_COMP", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
    control_vars = ctrl
  ) %>% mutate(variable = var)

  macro_results_asym[[var]] <- run_lp_asymmetric(
    analysis_data, y_var, "FCI_COMP", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
    control_vars = ctrl
  ) %>% mutate(variable = var)

  cat("done\n")
}

# Add credit (total) to macro results
# NOTE: Uses FCI_exCredit to address mechanical correlation
#       (Crecimiento_creditos is a component of FCI_COMP)
cat("Processing: Credit (Total) [using FCI_exCredit] ... ")
macro_results_std[["Creditos"]] <- run_lp_standard(
  analysis_data, "Cred_Total", "FCI_exCredit", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
  control_vars = c("IMAEP_yoy", "IPC_yoy")
) %>% mutate(variable = "Creditos")

macro_results_asym[["Creditos"]] <- run_lp_asymmetric(
  analysis_data, "Cred_Total", "FCI_exCredit", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
  control_vars = c("IMAEP_yoy", "IPC_yoy")
) %>% mutate(variable = "Creditos")
cat("done\n\n")

all_macro_std <- bind_rows(macro_results_std)
all_macro_asym <- bind_rows(macro_results_asym)


################################################################################
# 5. PART B: CREDIT CHANNEL DECOMPOSITION
################################################################################

cat("================================================================================\n")
cat("PART B: CREDIT CHANNEL DECOMPOSITION (3 FCI × 3 Credit types)\n")
cat("================================================================================\n\n")

credit_results_std <- list()
credit_results_asym <- list()

# Endogeneity correction: when credit is on LHS, swap FCI variants
# FCI_COMP -> FCI_exCredit (credit growth excluded)
# FCI_ENDO -> FCI_ENDO_exCredit (credit growth excluded from domestic)
# FCI_EXO  -> unchanged (no credit growth in exogenous)
fci_credit_swap <- c("FCI_COMP" = "FCI_exCredit",
                     "FCI_ENDO" = "FCI_ENDO_exCredit",
                     "FCI_EXO" = "FCI_EXO")

for (fci in LP_CONFIG$fci_types) {
  for (cred in LP_CONFIG$credit_types) {

    y_var <- paste0("Cred_", cred)
    key <- paste(fci, cred, sep = "_")

    # Use endogeneity-corrected FCI for credit regressions
    fci_actual <- fci_credit_swap[fci]

    cat(sprintf("  %s [%s] → %s ... ", FCI_LABELS[fci], fci_actual, CREDIT_LABELS[cred]))

    # Standard LP
    credit_results_std[[key]] <- run_lp_standard(
      analysis_data, y_var, fci_actual, LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    ) %>% mutate(fci_type = fci, credit_type = cred,
                 fci_label = FCI_LABELS[fci], credit_label = CREDIT_LABELS[cred])

    # Asymmetric LP
    credit_results_asym[[key]] <- run_lp_asymmetric(
      analysis_data, y_var, fci_actual, LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    ) %>% mutate(fci_type = fci, credit_type = cred,
                 fci_label = FCI_LABELS[fci], credit_label = CREDIT_LABELS[cred])

    cat("done\n")
  }
}

all_credit_std <- bind_rows(credit_results_std)
all_credit_asym <- bind_rows(credit_results_asym)


################################################################################
# 6. RESULTS SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("RESULTS SUMMARY\n")
cat("================================================================================\n\n")

# Part A: Macro effects
cat("PART A: MACRO EFFECTS (Comprehensive FCI)\n")
cat("Response to 1 SD FCI increase | *** p<0.01, ** p<0.05, * p<0.10\n\n")

for (var in c(LP_CONFIG$macro_variables, "Creditos")) {
  var_res <- all_macro_std %>% filter(variable == var, horizon %in% LP_CONFIG$horizons_report)
  if (nrow(var_res) == 0) next

  cat(var, ":\n")
  cat(sprintf("  %3s  %8s  %8s  %6s\n", "h", "Coef", "SE", ""))

  for (i in 1:nrow(var_res)) {
    r <- var_res[i,]
    stars <- ifelse(r$p_value < 0.01, "***", ifelse(r$p_value < 0.05, "**",
                                                     ifelse(r$p_value < 0.10, "*", "")))
    cat(sprintf("  %3d  %8.2f  %8.2f  %s\n", r$horizon, r$coef, r$se, stars))
  }
  cat("\n")
}

# Part B: Credit channel (12-month horizon)
cat("\nPART B: CREDIT CHANNEL at 12-MONTH HORIZON\n")
cat("Effect (pp) of 1 SD FCI increase on credit growth\n\n")

credit_12m <- all_credit_std %>% filter(horizon == 12)

cat(sprintf("%-25s %12s %12s %12s\n", "FCI Type", "Total", "USD", "Real_MN"))
cat(paste(rep("-", 63), collapse = ""), "\n")

for (fci in LP_CONFIG$fci_types) {
  vals <- credit_12m %>% filter(fci_type == fci)

  total_v <- vals %>% filter(credit_type == "Total")
  usd_v <- vals %>% filter(credit_type == "USD")
  real_mn_v <- vals %>% filter(credit_type == "Real_MN")

  fmt_val <- function(v) {
    if (nrow(v) == 0) return("NA")
    stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                                                     ifelse(v$p_value < 0.10, "*", "")))
    sprintf("%+.2f%s", v$coef, stars)
  }

  cat(sprintf("%-25s %12s %12s %12s\n",
              FCI_LABELS[fci], fmt_val(total_v), fmt_val(usd_v), fmt_val(real_mn_v)))
}

# Asymmetric effects summary
cat("\n\nASYMMETRIC EFFECTS: CREDIT CHANNEL at 12 MONTHS\n")
cat("Comparing Tightening (FCI>0) vs Easing (FCI<0) effects\n\n")

asym_12m <- all_credit_asym %>% filter(horizon == 12)

cat(sprintf("%-20s %-8s %10s %10s %8s\n", "FCI", "Credit", "Tight", "Ease", "Asym.p"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (fci in LP_CONFIG$fci_types) {
  for (cred in LP_CONFIG$credit_types) {
    v <- asym_12m %>% filter(fci_type == fci, credit_type == cred)
    if (nrow(v) == 0) next

    asym_star <- ifelse(v$asym_p < 0.05, "**", ifelse(v$asym_p < 0.10, "*", ""))

    cat(sprintf("%-20s %-8s %+10.2f %+10.2f %8.3f%s\n",
                FCI_LABELS[fci], cred, v$coef_tight, v$coef_ease, v$asym_p, asym_star))
  }
}


################################################################################
# 7. KEY ECONOMIC FINDINGS
################################################################################

cat("\n================================================================================\n")
cat("KEY ECONOMIC FINDINGS\n")
cat("================================================================================\n\n")

# Compare Endogenous vs Exogenous effects
cat("1. DOMESTIC vs EXTERNAL CONDITIONS (at 12 months):\n\n")

for (cred in LP_CONFIG$credit_types) {
  endo <- all_credit_std %>% filter(horizon == 12, fci_type == "FCI_ENDO", credit_type == cred)
  exo <- all_credit_std %>% filter(horizon == 12, fci_type == "FCI_EXO", credit_type == cred)

  if (nrow(endo) > 0 && nrow(exo) > 0) {
    ratio <- abs(endo$coef / exo$coef)
    cat(sprintf("   %s: Endogenous=%.2f pp, Exogenous=%.2f pp, Ratio=%.1fx\n",
                CREDIT_LABELS[cred], endo$coef, exo$coef, ratio))
  }
}

# Compare MN vs USD sensitivity
cat("\n2. MN vs USD SENSITIVITY:\n\n")

for (fci in LP_CONFIG$fci_types) {
  mn <- all_credit_std %>% filter(horizon == 12, fci_type == fci, credit_type == "MN")
  usd <- all_credit_std %>% filter(horizon == 12, fci_type == fci, credit_type == "USD")

  if (nrow(mn) > 0 && nrow(usd) > 0) {
    more_sens <- ifelse(abs(mn$coef) > abs(usd$coef), "MN", "USD")
    cat(sprintf("   %s: MN=%.2f, USD=%.2f → %s more sensitive\n",
                FCI_LABELS[fci], mn$coef, usd$coef, more_sens))
  }
}

# Asymmetry findings
cat("\n3. SIGNIFICANT ASYMMETRIES (p < 0.10):\n\n")

sig_asym <- all_credit_asym %>%
  filter(asym_p < 0.10, horizon == 12) %>%
  arrange(asym_p)

if (nrow(sig_asym) > 0) {
  for (i in 1:nrow(sig_asym)) {
    r <- sig_asym[i,]
    direction <- ifelse(abs(r$coef_tight) > abs(r$coef_ease),
                        "Tightening stronger", "Easing stronger")
    cat(sprintf("   %s → %s: %s (p=%.3f)\n",
                r$fci_label, r$credit_label, direction, r$asym_p))
  }
} else {
  cat("   No significant asymmetries at 12-month horizon\n")
}


################################################################################
# 8. VISUALIZATIONS
################################################################################

cat("\n================================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

if (!dir.exists(LP_CONFIG$output_dir)) dir.create(LP_CONFIG$output_dir, recursive = TRUE)

# Colors (colorblind-safe palette)
colors_fci <- c("Comprehensive" = "#2C3E50",
                "Endogenous (Domestic)" = "#D55E00",
                "Exogenous (External)" = "#0072B2")

colors_credit <- c("Real Total Credit" = "#2C3E50",
                   "USD (Dollarized)" = "#CC79A7",
                   "Real MN Credit" = "#009E73")

colors_asym <- c("Tightening" = "#D55E00", "Easing" = "#0072B2")

# Plot 1: Macro IRFs
p1 <- all_macro_std %>%
  ggplot(aes(x = horizon)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#377EB8", alpha = 0.2) +
  geom_line(aes(y = coef), color = "#377EB8", linewidth = 1) +
  geom_point(aes(y = coef), color = "#377EB8", size = 2) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  theme_minimal(base_size = 11) +
  labs(title = "Local Projections: FCI Effect on Macro Variables",
       subtitle = "Response to 1 SD FCI shock (tightening) | 90% CI",
       x = "Months", y = "Effect (pp)") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(LP_CONFIG$output_dir, "50_LP_Macro_Standard.png"), p1,
       width = 12, height = 4, dpi = 300)
cat("Saved: 50_LP_Macro_Standard.png\n")

# Plot 2: Credit channel by credit type
p2 <- all_credit_std %>%
  ggplot(aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~credit_label, ncol = 3) +
  scale_color_manual(values = colors_fci) +
  scale_fill_manual(values = colors_fci) +
  theme_minimal(base_size = 11) +
  labs(title = "Credit Channel: Effect of Different FCI Types",
       subtitle = "Response of credit growth (YoY) to 1 SD FCI shock | 90% CI",
       x = "Months", y = "Effect on credit growth (pp)",
       color = "FCI Type", fill = "FCI Type") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(LP_CONFIG$output_dir, "51_LP_Credit_byType.png"), p2,
       width = 14, height = 5, dpi = 300)
cat("Saved: 51_LP_Credit_byType.png\n")

# Plot 3: Credit channel by FCI type
p3 <- all_credit_std %>%
  ggplot(aes(x = horizon, y = coef, color = credit_label, fill = credit_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~fci_label, ncol = 3) +
  scale_color_manual(values = colors_credit) +
  scale_fill_manual(values = colors_credit) +
  theme_minimal(base_size = 11) +
  labs(title = "Credit Channel: Different Credit Types' Response",
       subtitle = "Comparing Real Total, USD, and Real MN credit | 90% CI",
       x = "Months", y = "Effect on credit growth (pp)",
       color = "Credit Type", fill = "Credit Type") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(LP_CONFIG$output_dir, "52_LP_Credit_byFCI.png"), p3,
       width = 14, height = 5, dpi = 300)
cat("Saved: 52_LP_Credit_byFCI.png\n")

# Plot 4: Heatmap at 12 months
heatmap_data <- all_credit_std %>%
  filter(horizon == 12) %>%
  mutate(label = sprintf("%.1f%s", coef,
                         ifelse(p_value < 0.01, "***",
                                ifelse(p_value < 0.05, "**",
                                       ifelse(p_value < 0.10, "*", "")))))

p4 <- ggplot(heatmap_data, aes(x = credit_label, y = fci_label, fill = coef)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = label), size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                       midpoint = 0, name = "Effect (pp)") +
  theme_minimal(base_size = 12) +
  labs(title = "FCI Effects on Credit Growth (12-month horizon)",
       subtitle = "Effect in percentage points | *** p<0.01, ** p<0.05, * p<0.10",
       x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_text(size = 11))

ggsave(file.path(LP_CONFIG$output_dir, "53_LP_Credit_Heatmap.png"), p4,
       width = 10, height = 6, dpi = 300)
cat("Saved: 53_LP_Credit_Heatmap.png\n")

# Plot 5: Asymmetric effects - Endogenous FCI on different credits
asym_endo <- all_credit_asym %>%
  filter(fci_type == "FCI_ENDO") %>%
  pivot_longer(cols = c(coef_tight, coef_ease),
               names_to = "effect_type", values_to = "coef") %>%
  mutate(
    effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
    ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
    ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
  )

p5 <- ggplot(asym_endo, aes(x = horizon, y = coef, color = effect_label, fill = effect_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~credit_label, ncol = 3) +
  scale_color_manual(values = colors_asym) +
  scale_fill_manual(values = colors_asym) +
  theme_minimal(base_size = 11) +
  labs(title = "Asymmetric Effects: Endogenous FCI on Credit",
       subtitle = "Comparing tightening (FCI>0) vs easing (FCI<0) episodes | 90% CI",
       x = "Months", y = "Effect on credit growth (pp)",
       color = "FCI State", fill = "FCI State") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(LP_CONFIG$output_dir, "54_LP_Asymmetric_Endo.png"), p5,
       width = 14, height = 5, dpi = 300)
cat("Saved: 54_LP_Asymmetric_Endo.png\n")

# Plot 6: Asymmetric effects - Exogenous FCI
asym_exo <- all_credit_asym %>%
  filter(fci_type == "FCI_EXO") %>%
  pivot_longer(cols = c(coef_tight, coef_ease),
               names_to = "effect_type", values_to = "coef") %>%
  mutate(
    effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
    ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
    ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
  )

p6 <- ggplot(asym_exo, aes(x = horizon, y = coef, color = effect_label, fill = effect_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~credit_label, ncol = 3) +
  scale_color_manual(values = colors_asym) +
  scale_fill_manual(values = colors_asym) +
  theme_minimal(base_size = 11) +
  labs(title = "Asymmetric Effects: Exogenous FCI on Credit",
       subtitle = "Comparing tightening vs easing of external conditions | 90% CI",
       x = "Months", y = "Effect on credit growth (pp)",
       color = "FCI State", fill = "FCI State") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(LP_CONFIG$output_dir, "55_LP_Asymmetric_Exo.png"), p6,
       width = 14, height = 5, dpi = 300)
cat("Saved: 55_LP_Asymmetric_Exo.png\n")


################################################################################
# 9. EXPORT RESULTS
################################################################################

cat("\n================================================================================\n")
cat("EXPORTING RESULTS\n")
cat("================================================================================\n\n")

write.csv(all_macro_std, file.path(LP_CONFIG$output_dir, "LP_Macro_Standard.csv"), row.names = FALSE)
write.csv(all_macro_asym, file.path(LP_CONFIG$output_dir, "LP_Macro_Asymmetric.csv"), row.names = FALSE)
write.csv(all_credit_std, file.path(LP_CONFIG$output_dir, "LP_Credit_Standard.csv"), row.names = FALSE)
write.csv(all_credit_asym, file.path(LP_CONFIG$output_dir, "LP_Credit_Asymmetric.csv"), row.names = FALSE)

cat("Saved: LP_Macro_Standard.csv\n")
cat("Saved: LP_Macro_Asymmetric.csv\n")
cat("Saved: LP_Credit_Standard.csv\n")
cat("Saved: LP_Credit_Asymmetric.csv\n")

################################################################################
# 10. PART C: NPL (MOROSIDAD) LOCAL PROJECTIONS
################################################################################
#
# Financial stability angle - policy-relevant
# NPL is used in FCI but never tested as outcome
# Creates "feedback loop" analysis unique to this paper
#
################################################################################

cat("\n================================================================================\n")
cat("PART C: NPL (MOROSIDAD) LOCAL PROJECTIONS\n")
cat("================================================================================\n\n")

# Load NPL data from original dataset
datos_npl <- tryCatch({
  read_excel(LP_CONFIG$data_file)
}, error = function(e) {
  cat("Error loading NPL data:", e$message, "\n")
  NULL
})

if (!is.null(datos_npl)) {
  fecha_col_npl <- names(datos_npl)[grepl("fecha|date", names(datos_npl), ignore.case = TRUE)][1]

  npl_data <- datos_npl %>%
    rename(fecha = !!sym(fecha_col_npl)) %>%
    mutate(fecha = as.Date(fecha)) %>%
    arrange(fecha) %>%
    dplyr::select(fecha, Morosidad) %>%
    mutate(
      # NPL level
      NPL_level = Morosidad,
      # NPL YoY change (difference, not growth rate since it's already a ratio)
      NPL_change = Morosidad - lag(Morosidad, 12)
    ) %>%
    na.omit()

  cat("NPL Data:\n")
  cat("  Observations:", nrow(npl_data), "\n")
  cat("  NPL Level: Mean=", round(mean(npl_data$NPL_level), 2), "%, SD=",
      round(sd(npl_data$NPL_level), 2), "%\n")
  cat("  NPL Change (YoY): Mean=", round(mean(npl_data$NPL_change, na.rm = TRUE), 3),
      ", SD=", round(sd(npl_data$NPL_change, na.rm = TRUE), 3), "\n\n")

  # Merge with FCI data
  npl_analysis <- fci_data %>%
    inner_join(npl_data, by = "fecha") %>%
    inner_join(macro_data, by = "fecha") %>%
    mutate(
      # Create asymmetric FCI components
      FCI_COMP_pos = pmax(FCI_COMP, 0),
      FCI_COMP_neg = abs(pmin(FCI_COMP, 0)),
      FCI_ENDO_pos = pmax(FCI_ENDO, 0),
      FCI_ENDO_neg = abs(pmin(FCI_ENDO, 0)),
      FCI_EXO_pos = pmax(FCI_EXO, 0),
      FCI_EXO_neg = abs(pmin(FCI_EXO, 0))
    ) %>%
    na.omit()

  cat("NPL Analysis Data: ", nrow(npl_analysis), " observations\n\n")

  # Run LP for NPL
  npl_results_std <- list()
  npl_results_asym <- list()

  for (fci in LP_CONFIG$fci_types) {
    cat(sprintf("  %s → NPL ... ", FCI_LABELS[fci]))

    # NPL level
    npl_results_std[[paste0(fci, "_level")]] <- run_lp_standard(
      npl_analysis, "NPL_level", fci, LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    ) %>% mutate(fci_type = fci, fci_label = FCI_LABELS[fci], npl_type = "Level")

    # NPL change
    npl_results_std[[paste0(fci, "_change")]] <- run_lp_standard(
      npl_analysis, "NPL_change", fci, LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    ) %>% mutate(fci_type = fci, fci_label = FCI_LABELS[fci], npl_type = "YoY Change")

    # Asymmetric for NPL level
    npl_results_asym[[fci]] <- run_lp_asymmetric(
      npl_analysis, "NPL_level", fci, LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    ) %>% mutate(fci_type = fci, fci_label = FCI_LABELS[fci])

    cat("done\n")
  }

  all_npl_std <- bind_rows(npl_results_std)
  all_npl_asym <- bind_rows(npl_results_asym)

  #---------------------------------------------------------------------------
  # NPL Results Summary
  #---------------------------------------------------------------------------
  cat("\n================================================================================\n")
  cat("NPL LOCAL PROJECTIONS RESULTS\n")
  cat("================================================================================\n\n")

  cat("Effect of 1 SD FCI shock on NPL (percentage points)\n")
  cat("*** p<0.01, ** p<0.05, * p<0.10\n\n")

  # Table for NPL level at key horizons
  npl_summary <- all_npl_std %>%
    filter(npl_type == "Level", horizon %in% c(6, 12, 18, 24)) %>%
    mutate(
      stars = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.10 ~ "*",
        TRUE ~ ""
      ),
      value = sprintf("%+.3f%s", coef, stars)
    )

  cat("NPL LEVEL Response at h=6, 12, 18, 24 months:\n")
  cat(sprintf("%-25s %12s %12s %12s %12s\n", "FCI Type", "h=6", "h=12", "h=18", "h=24"))
  cat(paste(rep("-", 72), collapse = ""), "\n")

  for (fci in LP_CONFIG$fci_types) {
    vals <- npl_summary %>% filter(fci_type == fci)
    h6 <- vals %>% filter(horizon == 6)
    h12 <- vals %>% filter(horizon == 12)
    h18 <- vals %>% filter(horizon == 18)
    h24 <- vals %>% filter(horizon == 24)

    cat(sprintf("%-25s %12s %12s %12s %12s\n",
                FCI_LABELS[fci],
                ifelse(nrow(h6) > 0, h6$value, "NA"),
                ifelse(nrow(h12) > 0, h12$value, "NA"),
                ifelse(nrow(h18) > 0, h18$value, "NA"),
                ifelse(nrow(h24) > 0, h24$value, "NA")))
  }

  # Asymmetric effects
  cat("\n\nASYMMETRIC EFFECTS ON NPL at 12 MONTHS:\n")
  cat(sprintf("%-25s %12s %12s %12s\n", "FCI Type", "Tightening", "Easing", "Asym.p"))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  asym_npl_12 <- all_npl_asym %>% filter(horizon == 12)

  for (fci in LP_CONFIG$fci_types) {
    v <- asym_npl_12 %>% filter(fci_type == fci)
    if (nrow(v) == 0) next

    asym_star <- ifelse(v$asym_p < 0.05, "**", ifelse(v$asym_p < 0.10, "*", ""))
    cat(sprintf("%-25s %+12.3f %+12.3f %12.3f%s\n",
                FCI_LABELS[fci], v$coef_tight, v$coef_ease, v$asym_p, asym_star))
  }

  #---------------------------------------------------------------------------
  # NPL Key Economic Findings
  #---------------------------------------------------------------------------
  cat("\n\nKEY FINDINGS - FINANCIAL STABILITY FEEDBACK:\n\n")

  # Check if tightening leads to higher NPL
  npl_12 <- all_npl_std %>% filter(npl_type == "Level", horizon == 12)

  for (fci in LP_CONFIG$fci_types) {
    v <- npl_12 %>% filter(fci_type == fci)
    if (nrow(v) == 0) next

    if (v$coef > 0 && v$p_value < 0.10) {
      cat(sprintf("  - %s: Tightening INCREASES NPL by %.3f pp (p=%.3f)\n",
                  FCI_LABELS[fci], v$coef, v$p_value))
    } else if (v$coef < 0 && v$p_value < 0.10) {
      cat(sprintf("  - %s: Tightening DECREASES NPL by %.3f pp (p=%.3f)\n",
                  FCI_LABELS[fci], abs(v$coef), v$p_value))
    } else {
      cat(sprintf("  - %s: No significant effect (coef=%.3f, p=%.3f)\n",
                  FCI_LABELS[fci], v$coef, v$p_value))
    }
  }

  #---------------------------------------------------------------------------
  # NPL Visualizations
  #---------------------------------------------------------------------------
  cat("\n\nGenerating NPL visualizations...\n")

  # Plot: NPL IRFs by FCI type
  p_npl_1 <- all_npl_std %>%
    filter(npl_type == "Level") %>%
    ggplot(aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = colors_fci) +
    scale_fill_manual(values = colors_fci) +
    theme_minimal(base_size = 11) +
    labs(title = "Financial Stability Feedback: FCI Effect on NPL",
         subtitle = "Response of NPL ratio (level) to 1 SD FCI shock | 90% CI",
         x = "Months", y = "Effect on NPL (pp)",
         color = "FCI Type", fill = "FCI Type") +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(LP_CONFIG$output_dir, "56_LP_NPL_byFCI.png"), p_npl_1,
         width = 12, height = 6, dpi = 300)
  cat("Saved: 56_LP_NPL_byFCI.png\n")

  # Plot: Asymmetric effects on NPL
  asym_npl_plot <- all_npl_asym %>%
    pivot_longer(cols = c(coef_tight, coef_ease),
                 names_to = "effect_type", values_to = "coef") %>%
    mutate(
      effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
      ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
      ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
    )

  p_npl_2 <- ggplot(asym_npl_plot, aes(x = horizon, y = coef,
                                        color = effect_label, fill = effect_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~fci_label, ncol = 3) +
    scale_color_manual(values = colors_asym) +
    scale_fill_manual(values = colors_asym) +
    theme_minimal(base_size = 11) +
    labs(title = "Asymmetric FCI Effects on NPL",
         subtitle = "Comparing tightening (FCI>0) vs easing (FCI<0) effects | 90% CI",
         x = "Months", y = "Effect on NPL (pp)",
         color = "FCI State", fill = "FCI State") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(LP_CONFIG$output_dir, "57_LP_NPL_Asymmetric.png"), p_npl_2,
         width = 14, height = 5, dpi = 300)
  cat("Saved: 57_LP_NPL_Asymmetric.png\n")

  # Export NPL results
  write.csv(all_npl_std, file.path(LP_CONFIG$output_dir, "LP_NPL_Standard.csv"), row.names = FALSE)
  write.csv(all_npl_asym, file.path(LP_CONFIG$output_dir, "LP_NPL_Asymmetric.csv"), row.names = FALSE)
  cat("Saved: LP_NPL_Standard.csv\n")
  cat("Saved: LP_NPL_Asymmetric.csv\n")

  #---------------------------------------------------------------------------
  # ROBUSTNESS CHECK: FCI Excluding NPL (Addresses Endogeneity Concern)
  #---------------------------------------------------------------------------
  cat("\n\n================================================================================\n")
  cat("ROBUSTNESS: FCI EXCLUDING NPL (Endogeneity Test)\n")
  cat("================================================================================\n\n")

  cat("ISSUE: NPL (Morosidad) is a component of FCI with sign +1.\n")
  cat("       When FCI increases, it may be partly because NPL increased.\n")
  cat("       This creates mechanical correlation in the LP: FCI_t -> NPL_{t+h}\n\n")

  cat("SOLUTION: Create FCI_exNPL = FCI calculated WITHOUT Morosidad\n\n")

  # Calculate FCI excluding NPL
  # We need to recalculate from the raw standardized data
  tryCatch({
    # Load raw data again
    datos_raw_npl <- read_excel(LP_CONFIG$data_file)
    fecha_col_raw <- names(datos_raw_npl)[grepl("fecha|date", names(datos_raw_npl), ignore.case = TRUE)][1]

    # Variable definitions (from 01_FCI_Complete.R)
    banking_vars_exNPL <- c("Crecimiento_creditos", "Ratio_Cred_Depo", "Rentabilidad", "Liquidez")
    banking_signs_exNPL <- c(-1, -1, -1, -1)  # Removed Morosidad

    rates_vars <- c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM")
    rates_signs <- c(+1, +1, +1)

    external_vars <- c("TCN", "Commodities", "FFER", "VIX")
    external_signs <- c(+1, -1, +1, +1)

    # All vars excluding NPL
    all_vars_exNPL <- c(rates_vars, banking_vars_exNPL, external_vars)
    all_signs_exNPL <- c(rates_signs, banking_signs_exNPL, external_signs)

    # Prepare data
    datos_exNPL <- datos_raw_npl %>%
      rename(fecha = !!sym(fecha_col_raw)) %>%
      mutate(fecha = as.Date(fecha)) %>%
      arrange(fecha)

    # Calculate credit growth
    if ("Creditos_Sector_privado_totales" %in% names(datos_exNPL)) {
      datos_exNPL <- datos_exNPL %>%
        mutate(Crecimiento_creditos = (Creditos_Sector_privado_totales /
                                        lag(Creditos_Sector_privado_totales, 12) - 1) * 100)
    }

    # Apply signs
    for (i in seq_along(all_vars_exNPL)) {
      if (all_vars_exNPL[i] %in% names(datos_exNPL)) {
        datos_exNPL[[all_vars_exNPL[i]]] <- datos_exNPL[[all_vars_exNPL[i]]] * all_signs_exNPL[i]
      }
    }

    # Standardize (rolling 60-month window)
    datos_exNPL <- datos_exNPL %>%
      dplyr::select(fecha, any_of(all_vars_exNPL)) %>%
      mutate(across(where(is.numeric), ~ {
        rm <- zoo::rollapply(., width = 60, FUN = mean, na.rm = TRUE,
                             fill = NA, align = "right", partial = TRUE)
        rsd <- zoo::rollapply(., width = 60, FUN = sd, na.rm = TRUE,
                              fill = NA, align = "right", partial = TRUE)
        (. - rm) / rsd
      }))

    # Calculate FCI_exNPL as simple average (Z-Score method - most stable)
    vars_available <- intersect(all_vars_exNPL, names(datos_exNPL))
    datos_exNPL$FCI_exNPL <- rowMeans(datos_exNPL[, vars_available, drop = FALSE], na.rm = TRUE)

    cat("FCI_exNPL calculated using", length(vars_available), "variables (excluding Morosidad)\n\n")

    # Merge with NPL data
    npl_robustness <- datos_exNPL %>%
      dplyr::select(fecha, FCI_exNPL) %>%
      inner_join(npl_data, by = "fecha") %>%
      inner_join(macro_data, by = "fecha") %>%
      mutate(
        FCI_exNPL_pos = pmax(FCI_exNPL, 0),
        FCI_exNPL_neg = abs(pmin(FCI_exNPL, 0))
      ) %>%
      na.omit()

    cat("Robustness data:", nrow(npl_robustness), "observations\n\n")

    # Run LP with FCI_exNPL
    cat("Running LP: FCI_exNPL → NPL...\n")

    npl_exNPL_std <- run_lp_standard(
      npl_robustness, "NPL_level", "FCI_exNPL", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    ) %>% mutate(fci_type = "FCI_exNPL", fci_label = "Excl. NPL", npl_type = "Level")

    npl_exNPL_asym <- run_lp_asymmetric(
      npl_robustness, "NPL_level", "FCI_exNPL", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    ) %>% mutate(fci_type = "FCI_exNPL", fci_label = "Excl. NPL")

    # Compare results
    cat("\n================================================================================\n")
    cat("COMPARISON: FCI with NPL vs FCI excluding NPL\n")
    cat("================================================================================\n\n")

    # Get ENDO results for comparison (since ENDO includes Morosidad)
    endo_12 <- all_npl_std %>%
      filter(fci_type == "FCI_ENDO", npl_type == "Level", horizon == 12)

    exNPL_12 <- npl_exNPL_std %>%
      filter(horizon == 12)

    cat("Effect on NPL at h=12 months:\n\n")
    cat(sprintf("  FCI_ENDO (includes NPL): coef = %+.4f (p = %.3f)\n",
                endo_12$coef, endo_12$p_value))
    cat(sprintf("  FCI_exNPL (excludes NPL): coef = %+.4f (p = %.3f)\n",
                exNPL_12$coef, exNPL_12$p_value))

    # Interpretation
    cat("\nINTERPRETATION:\n")
    if (exNPL_12$p_value < 0.10 && exNPL_12$coef > 0) {
      cat("  >>> EFFECT IS ROBUST: Even after excluding NPL from FCI,\n")
      cat("      tighter financial conditions still predict higher NPL.\n")
      cat("      This is NOT just mechanical correlation.\n")
    } else if (exNPL_12$p_value >= 0.10) {
      cat("  >>> CAUTION: Effect becomes insignificant when excluding NPL.\n")
      cat("      The original result may be partially driven by mechanical correlation.\n")
      cat("      The FCI-NPL relationship should be interpreted with care.\n")
    } else {
      cat("  >>> MIXED: Effect changes sign or magnitude significantly.\n")
      cat("      Further investigation needed.\n")
    }

    # Correlation between FCI_ENDO and FCI_exNPL
    fci_comparison <- npl_analysis %>%
      inner_join(npl_robustness %>% dplyr::select(fecha, FCI_exNPL), by = "fecha")

    cat(sprintf("\nCorrelation FCI_ENDO vs FCI_exNPL: %.3f\n",
                cor(fci_comparison$FCI_ENDO, fci_comparison$FCI_exNPL, use = "complete")))

    # Plot comparison
    exNPL_plot <- npl_exNPL_std %>%
      mutate(fci_label = "Excluding NPL")

    endo_plot <- all_npl_std %>%
      filter(fci_type == "FCI_ENDO", npl_type == "Level") %>%
      mutate(fci_label = "Including NPL")

    comparison_plot <- bind_rows(exNPL_plot, endo_plot)

    p_robustness <- ggplot(comparison_plot, aes(x = horizon, y = coef,
                                                 color = fci_label, fill = fci_label)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      scale_color_manual(values = c("Including NPL" = "#E74C3C", "Excluding NPL" = "#3498DB")) +
      scale_fill_manual(values = c("Including NPL" = "#E74C3C", "Excluding NPL" = "#3498DB")) +
      theme_minimal(base_size = 11) +
      labs(title = "Robustness Check: FCI Effect on NPL",
           subtitle = "Comparing FCI with and without NPL as a component | 90% CI",
           x = "Months", y = "Effect on NPL (pp)",
           color = "FCI Specification", fill = "FCI Specification") +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom")

    ggsave(file.path(LP_CONFIG$output_dir, "58_LP_NPL_Robustness_exNPL.png"), p_robustness,
           width = 10, height = 6, dpi = 300)
    cat("\nSaved: 58_LP_NPL_Robustness_exNPL.png\n")

    # Export
    write.csv(npl_exNPL_std, file.path(LP_CONFIG$output_dir, "LP_NPL_exNPL.csv"), row.names = FALSE)
    cat("Saved: LP_NPL_exNPL.csv\n")

  }, error = function(e) {
    cat("Error in robustness check:", e$message, "\n")
  })

} else {
  cat("NPL data not available - skipping Part C\n")
}


################################################################################
# PART E: IV-LP EXTERNAL INSTRUMENT IDENTIFICATION
################################################################################
#
# Strategy: Use external variables (VIX, FFER, Commodities) as excluded
# instruments for FCI_ENDO_exCredit under the small-open-economy exclusion
# restriction: these variables affect Paraguayan credit only through their
# effect on domestic financial conditions.
#
# 3 instruments - 1 endogenous = 2 overidentification degrees of freedom
#
################################################################################

cat("\n================================================================================\n")
cat("PART E: IV-LP EXTERNAL INSTRUMENT IDENTIFICATION\n")
cat("================================================================================\n\n")

if (!require(ivreg, quietly = TRUE)) {
  install.packages("ivreg", repos = "https://cloud.r-project.org")
  library(ivreg)
}

instruments <- c("VIX", "FFER", "Commodities")
iv_endo_var <- "FCI_ENDO_exCredit"
iv_controls <- c("IMAEP_yoy", "IPC_yoy")

cat("Endogenous variable:", iv_endo_var, "\n")
cat("Instruments:", paste(instruments, collapse = ", "), "\n")
cat("Controls:", paste(iv_controls, collapse = ", "), "\n\n")

# --- Run IV-LP for Total Credit ---
cat("Running IV-LP: FCI_ENDO_exCredit → Cred_Total ... ")
iv_total <- run_lp_iv(
  analysis_data, "Cred_Total", iv_endo_var, instruments,
  LP_CONFIG$max_horizon, LP_CONFIG$n_lags, control_vars = iv_controls
)
cat("done\n")

# --- Run IV-LP for Real MN Credit ---
cat("Running IV-LP: FCI_ENDO_exCredit → Cred_Real_MN ... ")
iv_real_mn <- run_lp_iv(
  analysis_data, "Cred_Real_MN", iv_endo_var, instruments,
  LP_CONFIG$max_horizon, LP_CONFIG$n_lags, control_vars = iv_controls
)
cat("done\n")

# --- Run OLS baselines for apples-to-apples comparison ---
cat("Running OLS baselines (same spec) ... ")
ols_total <- run_lp_standard(
  analysis_data, "Cred_Total", iv_endo_var,
  LP_CONFIG$max_horizon, LP_CONFIG$n_lags, control_vars = iv_controls
)
ols_real_mn <- run_lp_standard(
  analysis_data, "Cred_Real_MN", iv_endo_var,
  LP_CONFIG$max_horizon, LP_CONFIG$n_lags, control_vars = iv_controls
)
cat("done\n\n")

# --- Diagnostic Summary ---
cat("================================================================================\n")
cat("IV DIAGNOSTICS at KEY HORIZONS\n")
cat("================================================================================\n\n")

cat(sprintf("%-8s %-10s %10s %10s %10s %10s %10s\n",
            "Credit", "Horizon", "1st-F", "Partial-R2", "Sargan-p", "Wu-Haus-p", "N"))
cat(paste(rep("-", 72), collapse = ""), "\n")

for (ctype in c("Total", "Real_MN")) {
  d <- if (ctype == "Total") iv_total$diag else iv_real_mn$diag
  for (hh in c(6, 12, 18)) {
    row <- d %>% filter(horizon == hh)
    if (nrow(row) == 0) next
    cat(sprintf("%-8s h=%-8d %10.2f %10.3f %10.3f %10.3f %10d\n",
                ctype, hh,
                row$first_stage_F, row$partial_R2,
                row$sargan_p, row$wu_hausman_p, row$n_obs))
  }
}

# --- OLS vs IV Comparison ---
cat("\n================================================================================\n")
cat("OLS vs IV COEFFICIENT COMPARISON\n")
cat("================================================================================\n\n")

cat(sprintf("%-8s %-8s %10s %10s %10s %10s\n",
            "Credit", "Horizon", "OLS-coef", "IV-coef", "OLS-p", "IV-p"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (ctype in c("Total", "Real_MN")) {
  iv_c <- if (ctype == "Total") iv_total$coefs else iv_real_mn$coefs
  ols_c <- if (ctype == "Total") ols_total else ols_real_mn
  for (hh in c(6, 12, 18)) {
    iv_row <- iv_c %>% filter(horizon == hh)
    ols_row <- ols_c %>% filter(horizon == hh)
    if (nrow(iv_row) == 0 || nrow(ols_row) == 0) next
    cat(sprintf("%-8s h=%-8d %+10.2f %+10.2f %10.3f %10.3f\n",
                ctype, hh,
                ols_row$coef, iv_row$coef,
                ols_row$p_value, iv_row$p_value))
  }
}

# --- Visualization: 4-panel figure ---
cat("\nGenerating IV-LP visualization...\n")

# Combine OLS and IV results
iv_total_plot <- iv_total$coefs %>% mutate(method = "IV (2SLS)", credit_type = "Real Total Credit")
ols_total_plot <- ols_total %>% mutate(method = "OLS", credit_type = "Real Total Credit")
iv_mn_plot <- iv_real_mn$coefs %>% mutate(method = "IV (2SLS)", credit_type = "Real MN Credit")
ols_mn_plot <- ols_real_mn %>% mutate(method = "OLS", credit_type = "Real MN Credit")

all_iv_ols <- bind_rows(iv_total_plot, ols_total_plot, iv_mn_plot, ols_mn_plot)

colors_method <- c("OLS" = "#2C3E50", "IV (2SLS)" = "#E74C3C")

# Panel A: OLS vs IV for Total Credit
pA <- all_iv_ols %>%
  filter(credit_type == "Real Total Credit") %>%
  ggplot(aes(x = horizon, y = coef, color = method, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = colors_method) +
  scale_fill_manual(values = colors_method) +
  theme_minimal(base_size = 10) +
  labs(title = "A. Real Total Credit", x = "Months", y = "Effect (pp)",
       color = NULL, fill = NULL) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# Panel B: OLS vs IV for Real MN Credit
pB <- all_iv_ols %>%
  filter(credit_type == "Real MN Credit") %>%
  ggplot(aes(x = horizon, y = coef, color = method, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = colors_method) +
  scale_fill_manual(values = colors_method) +
  theme_minimal(base_size = 10) +
  labs(title = "B. Real MN Credit", x = "Months", y = "Effect (pp)",
       color = NULL, fill = NULL) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# Panel C: First-stage F across horizons
diag_all <- bind_rows(
  iv_total$diag %>% mutate(credit_type = "Real Total Credit"),
  iv_real_mn$diag %>% mutate(credit_type = "Real MN Credit")
)

pC <- ggplot(diag_all, aes(x = horizon, y = first_stage_F, color = credit_type)) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  annotate("text", x = 1, y = 11, label = "Stock-Yogo F=10", hjust = 0,
           size = 3, color = "#E74C3C") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Real Total Credit" = "#2C3E50", "Real MN Credit" = "#27AE60")) +
  theme_minimal(base_size = 10) +
  labs(title = "C. First-Stage F-Statistic", x = "Months", y = "F-statistic",
       color = NULL) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# Panel D: Sargan J p-value across horizons
pD <- ggplot(diag_all, aes(x = horizon, y = sargan_p, color = credit_type)) +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  annotate("text", x = 1, y = 0.13, label = "10% rejection", hjust = 0,
           size = 3, color = "#E74C3C") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Real Total Credit" = "#2C3E50", "Real MN Credit" = "#27AE60")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 10) +
  labs(title = "D. Sargan J-Test p-value", x = "Months", y = "p-value",
       color = NULL) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

p_iv <- gridExtra::grid.arrange(
  pA, pB, pC, pD, ncol = 2,
  top = grid::textGrob(
    "IV-LP Identification: External Instruments for Domestic FCI",
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  ),
  bottom = grid::textGrob(
    "Instruments: VIX, FFER, Commodities | Endogenous: FCI_ENDO_exCredit | 90% CI",
    gp = grid::gpar(fontsize = 9, col = "gray40")
  )
)

ggsave(file.path(LP_CONFIG$output_dir, "59_LP_IV_Identification.png"), p_iv,
       width = 12, height = 10, dpi = 300)
cat("Saved: 59_LP_IV_Identification.png\n")

# --- Export CSVs ---
iv_results_export <- bind_rows(
  iv_total$coefs %>% mutate(credit_type = "Total", method = "IV"),
  iv_real_mn$coefs %>% mutate(credit_type = "Real_MN", method = "IV"),
  ols_total %>% mutate(credit_type = "Total", method = "OLS"),
  ols_real_mn %>% mutate(credit_type = "Real_MN", method = "OLS")
)

iv_diag_export <- bind_rows(
  iv_total$diag %>% mutate(credit_type = "Total"),
  iv_real_mn$diag %>% mutate(credit_type = "Real_MN")
)

write.csv(iv_results_export, file.path(LP_CONFIG$output_dir, "LP_IV_Results.csv"), row.names = FALSE)
write.csv(iv_diag_export, file.path(LP_CONFIG$output_dir, "LP_IV_Diagnostics.csv"), row.names = FALSE)
cat("Saved: LP_IV_Results.csv\n")
cat("Saved: LP_IV_Diagnostics.csv\n")

# --- Summary interpretation ---
cat("\n================================================================================\n")
cat("IV-LP INTERPRETATION SUMMARY\n")
cat("================================================================================\n\n")

# Check first-stage strength
median_F_total <- median(iv_total$diag$first_stage_F, na.rm = TRUE)
median_F_mn <- median(iv_real_mn$diag$first_stage_F, na.rm = TRUE)
cat(sprintf("Median first-stage F: Total=%.1f, Real_MN=%.1f\n", median_F_total, median_F_mn))

if (median_F_total > 10) {
  cat("  -> Instruments are STRONG (F > Stock-Yogo threshold of 10)\n")
} else {
  cat("  -> WARNING: Instruments may be WEAK (F < 10)\n")
}

# Check Sargan
median_sargan_total <- median(iv_total$diag$sargan_p, na.rm = TRUE)
cat(sprintf("Median Sargan p-value: Total=%.3f\n", median_sargan_total))

if (median_sargan_total > 0.10) {
  cat("  -> Exclusion restriction NOT rejected (Sargan p > 0.10)\n")
} else {
  cat("  -> CAUTION: Exclusion restriction rejected at some horizons\n")
}

# Check Wu-Hausman
median_wu_total <- median(iv_total$diag$wu_hausman_p, na.rm = TRUE)
cat(sprintf("Median Wu-Hausman p-value: Total=%.3f\n", median_wu_total))

if (median_wu_total > 0.10) {
  cat("  -> OLS is consistent (Wu-Hausman p > 0.10); IV confirms robustness\n")
} else {
  cat("  -> OLS may be biased (Wu-Hausman p < 0.10); IV estimates preferred\n")
}

# Compare IV vs OLS at h=12
iv_12_total <- iv_total$coefs %>% filter(horizon == 12)
ols_12_total <- ols_total %>% filter(horizon == 12)
if (nrow(iv_12_total) > 0 && nrow(ols_12_total) > 0) {
  cat(sprintf("\nAt h=12 (Total Credit): OLS=%.2f pp, IV=%.2f pp (ratio=%.2f)\n",
              ols_12_total$coef, iv_12_total$coef,
              iv_12_total$coef / ols_12_total$coef))
}


################################################################################
# PART F: PLACEBO TEST (BACKWARD LOCAL PROJECTIONS)
################################################################################
#
# If FCI genuinely predicts *future* credit, it should NOT predict *past* credit.
# We regress lag(credit, h) on FCI_t — reversing the temporal direction.
# Significant backward coefficients would indicate omitted-variable bias or
# mechanical contamination rather than genuine forward-looking content.
#
################################################################################

cat("\n================================================================================\n")
cat("PART F: PLACEBO TEST (BACKWARD LOCAL PROJECTIONS)\n")
cat("================================================================================\n\n")

z_crit <- qnorm(1 - (1 - LP_CONFIG$confidence_level) / 2)

placebo_results <- data.frame()

for (h in 1:LP_CONFIG$max_horizon) {
  data_h <- analysis_data %>%
    mutate(
      y_back = lag(Cred_Total, h),
      y_lag1 = lag(Cred_Total, 1),
      y_lag2 = lag(Cred_Total, 2),
      fci_lag1 = lag(FCI_exCredit, 1),
      IMAEP_yoy_c = IMAEP_yoy,
      IMAEP_yoy_lag1 = lag(IMAEP_yoy, 1),
      IMAEP_yoy_lag2 = lag(IMAEP_yoy, 2),
      IPC_yoy_c = IPC_yoy,
      IPC_yoy_lag1 = lag(IPC_yoy, 1),
      IPC_yoy_lag2 = lag(IPC_yoy, 2)
    )

  formula_str <- "y_back ~ FCI_exCredit + y_lag1 + y_lag2 + fci_lag1 + IMAEP_yoy_c + IMAEP_yoy_lag1 + IMAEP_yoy_lag2 + IPC_yoy_c + IPC_yoy_lag1 + IPC_yoy_lag2"

  reg_data <- data_h %>%
    dplyr::select(y_back, FCI_exCredit, y_lag1, y_lag2, fci_lag1,
                  IMAEP_yoy_c, IMAEP_yoy_lag1, IMAEP_yoy_lag2,
                  IPC_yoy_c, IPC_yoy_lag1, IPC_yoy_lag2) %>%
    na.omit()

  if (nrow(reg_data) < 30) next

  model <- lm(as.formula(formula_str), data = reg_data)
  vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
  coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

  idx <- which(rownames(coef_test) == "FCI_exCredit")

  placebo_results <- rbind(placebo_results, data.frame(
    horizon = h,
    direction = "Placebo (backward)",
    coef = coef_test[idx, 1],
    se = coef_test[idx, 2],
    p_value = coef_test[idx, 4],
    ci_lower = coef_test[idx, 1] - z_crit * coef_test[idx, 2],
    ci_upper = coef_test[idx, 1] + z_crit * coef_test[idx, 2],
    n_obs = nrow(reg_data)
  ))
}

# Combine with forward LP results
forward_results <- macro_results_std[["Creditos"]] %>%
  mutate(direction = "Forward (standard)")

placebo_combined <- bind_rows(
  forward_results %>% dplyr::select(horizon, direction, coef, se, p_value, ci_lower, ci_upper, n_obs),
  placebo_results
)

# Console summary
cat("Placebo Test Results (Backward LP) at key horizons:\n\n")
cat(sprintf("%-5s %12s %12s %12s %12s\n", "h", "Forward", "p-fwd", "Placebo", "p-placebo"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (hh in c(3, 6, 12)) {
  fwd <- forward_results %>% filter(horizon == hh)
  plac <- placebo_results %>% filter(horizon == hh)
  if (nrow(fwd) > 0 && nrow(plac) > 0) {
    cat(sprintf("h=%-3d %+12.2f %12.3f %+12.2f %12.3f\n",
                hh, fwd$coef, fwd$p_value, plac$coef, plac$p_value))
  }
}

cat("\nInterpretation: Forward coefficients should be large and significant;\n")
cat("placebo (backward) coefficients should be near zero and insignificant.\n")

# Visualization: single panel with forward vs placebo
p_placebo <- ggplot(placebo_combined, aes(x = horizon, y = coef,
                                           color = direction, fill = direction,
                                           linetype = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Forward (standard)" = "#2C3E50", "Placebo (backward)" = "#95A5A6")) +
  scale_fill_manual(values = c("Forward (standard)" = "#2C3E50", "Placebo (backward)" = "#95A5A6")) +
  scale_linetype_manual(values = c("Forward (standard)" = "solid", "Placebo (backward)" = "dashed")) +
  theme_minimal(base_size = 11) +
  labs(title = "Placebo Test: Forward vs Backward Local Projections",
       subtitle = "FCI_exCredit effect on real credit growth | 90% CI | Backward should be zero",
       x = "Horizon (months)", y = "Effect (pp)",
       color = "Direction", fill = "Direction", linetype = "Direction") +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(LP_CONFIG$output_dir, "59b_LP_Placebo_Test.png"), p_placebo,
       width = 10, height = 6, dpi = 300)
cat("\nSaved: 59b_LP_Placebo_Test.png\n")

write.csv(placebo_combined, file.path(LP_CONFIG$output_dir, "LP_Placebo_Test.csv"), row.names = FALSE)
cat("Saved: LP_Placebo_Test.csv\n")


################################################################################
# PART G: PURGED FCI (MACRO-ORTHOGONALIZED)
################################################################################
#
# Conservative lower-bound test: purge FCI_exCredit of all macro information
# by regressing it on IMAEP_yoy and IPC_yoy (+ 2 lags each). The residuals
# ("FCI_purged") contain only the financial conditions signal that is
# orthogonal to current and recent macroeconomic conditions.
#
# This is deliberately over-conservative — it removes both genuine FCI
# content AND any causal effect that macro → FCI, potentially biasing
# the estimate toward zero. A significant purged-FCI coefficient provides
# a strong lower bound on the true predictive relationship.
#
################################################################################

cat("\n================================================================================\n")
cat("PART G: PURGED FCI (MACRO-ORTHOGONALIZED)\n")
cat("================================================================================\n\n")

# Purging regression: remove macro information from FCI
purge_formula <- FCI_exCredit ~ IMAEP_yoy + lag(IMAEP_yoy, 1) + lag(IMAEP_yoy, 2) +
                                IPC_yoy + lag(IPC_yoy, 1) + lag(IPC_yoy, 2)

purge_data <- analysis_data %>%
  mutate(
    IMAEP_yoy_L1 = lag(IMAEP_yoy, 1),
    IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
    IPC_yoy_L1 = lag(IPC_yoy, 1),
    IPC_yoy_L2 = lag(IPC_yoy, 2)
  ) %>%
  na.omit()

purge_model <- lm(FCI_exCredit ~ IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +
                                  IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2,
                   data = purge_data)

purge_r2 <- summary(purge_model)$r.squared
purge_data$FCI_purged <- residuals(purge_model)

# Add asymmetric components for purged FCI
purge_data$FCI_purged_pos <- pmax(purge_data$FCI_purged, 0)
purge_data$FCI_purged_neg <- abs(pmin(purge_data$FCI_purged, 0))

cat(sprintf("Purging regression R-squared: %.3f\n", purge_r2))
cat(sprintf("  -> Macro variables explain %.1f%% of FCI_exCredit variation\n", purge_r2 * 100))
cat(sprintf("  -> FCI_purged retains %.1f%% of original variation\n\n", (1 - purge_r2) * 100))

# Run LP with purged FCI
cat("Running LP: FCI_purged -> Cred_Total ... ")
purged_total <- run_lp_standard(
  purge_data, "Cred_Total", "FCI_purged", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
  control_vars = c("IMAEP_yoy", "IPC_yoy")
) %>% mutate(fci_type = "Purged", credit_type = "Total")
cat("done\n")

cat("Running LP: FCI_purged -> Cred_Real_MN ... ")
purged_mn <- run_lp_standard(
  purge_data, "Cred_Real_MN", "FCI_purged", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
  control_vars = c("IMAEP_yoy", "IPC_yoy")
) %>% mutate(fci_type = "Purged", credit_type = "Real_MN")
cat("done\n")

# Get baselines for comparison (from credit_results_std)
baseline_total <- credit_results_std[["FCI_COMP_Total"]] %>%
  mutate(fci_type = "Baseline (FCI_exCredit)", credit_type = "Total")
baseline_mn <- credit_results_std[["FCI_COMP_Real_MN"]] %>%
  mutate(fci_type = "Baseline (FCI_exCredit)", credit_type = "Real_MN")

# Console summary
cat("\nPurged FCI vs Baseline at key horizons:\n\n")
cat(sprintf("%-8s %-10s %12s %12s %8s\n", "Credit", "Horizon", "Baseline", "Purged", "Ratio"))
cat(paste(rep("-", 52), collapse = ""), "\n")

for (ctype in c("Total", "Real_MN")) {
  base <- if (ctype == "Total") baseline_total else baseline_mn
  purg <- if (ctype == "Total") purged_total else purged_mn
  for (hh in c(6, 12)) {
    b <- base %>% filter(horizon == hh)
    p <- purg %>% filter(horizon == hh)
    if (nrow(b) > 0 && nrow(p) > 0) {
      ratio <- ifelse(b$coef != 0, p$coef / b$coef, NA)
      cat(sprintf("%-8s h=%-8d %+12.2f %+12.2f %8.1f%%\n",
                  ctype, hh, b$coef, p$coef, ratio * 100))
    }
  }
}

cat("\nCAVEAT: Purging is deliberately over-conservative. It removes both:\n")
cat("  (a) macro → FCI feedback (appropriate to remove)\n")
cat("  (b) FCI → macro causal content (inappropriate to remove)\n")
cat("Significant purged-FCI effects represent a conservative lower bound.\n")

# Visualization: 2-panel faceted plot
purged_plot_data <- bind_rows(
  baseline_total %>% dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, fci_type, credit_type),
  baseline_mn %>% dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, fci_type, credit_type),
  purged_total %>% dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, fci_type, credit_type),
  purged_mn %>% dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, fci_type, credit_type)
) %>%
  mutate(credit_label = ifelse(credit_type == "Total", "Real Total Credit", "Real MN Credit"))

p_purged <- ggplot(purged_plot_data, aes(x = horizon, y = coef,
                                          color = fci_type, fill = fci_type,
                                          linetype = fci_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~credit_label, ncol = 2) +
  scale_color_manual(values = c("Baseline (FCI_exCredit)" = "#2C3E50", "Purged" = "#E67E22")) +
  scale_fill_manual(values = c("Baseline (FCI_exCredit)" = "#2C3E50", "Purged" = "#E67E22")) +
  scale_linetype_manual(values = c("Baseline (FCI_exCredit)" = "solid", "Purged" = "dashed")) +
  theme_minimal(base_size = 11) +
  labs(title = "Purged FCI: Macro-Orthogonalized Predictor (Conservative Lower Bound)",
       subtitle = sprintf("Purging R² = %.3f | FCI_purged = residuals after removing IMAEP, IPC + 2 lags | 90%% CI", purge_r2),
       x = "Horizon (months)", y = "Effect on credit growth (pp)",
       color = "FCI Specification", fill = "FCI Specification",
       linetype = "FCI Specification") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(LP_CONFIG$output_dir, "59c_LP_Purged_FCI.png"), p_purged,
       width = 12, height = 6, dpi = 300)
cat("\nSaved: 59c_LP_Purged_FCI.png\n")

purged_export <- bind_rows(
  purged_total %>% dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, n_obs, fci_type, credit_type),
  purged_mn %>% dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, n_obs, fci_type, credit_type)
)
write.csv(purged_export, file.path(LP_CONFIG$output_dir, "LP_Purged_FCI.csv"), row.names = FALSE)
cat("Saved: LP_Purged_FCI.csv\n")


################################################################################
# PART H: RATES-ONLY FCI, exTCN, LAG SENSITIVITY, METHOD-BY-METHOD LP
################################################################################

cat("\n================================================================================\n")
cat("PART H: ADDITIONAL ROBUSTNESS — RATES-ONLY, exTCN, LAGS, METHODS\n")
cat("================================================================================\n\n")

# --- H1: Rates-only FCI LP (3 price variables only: TPM, 2 spreads) ---
cat("H1: Rates-only FCI credit LP...\n")

rates_lp_results <- data.frame()
for (cred in c("Cred_Total", "Cred_Real_MN", "Cred_USD")) {
  res <- run_lp_standard(
    analysis_data, cred, "FCI_RATES", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  )
  if (nrow(res) > 0) {
    res$credit_type <- cred
    res$fci_type <- "FCI_RATES"
    rates_lp_results <- rbind(rates_lp_results, res)
  }
}

cat("  Rates-only LP at h=12:\n")
for (ct in unique(rates_lp_results$credit_type)) {
  r <- rates_lp_results[rates_lp_results$credit_type == ct & rates_lp_results$horizon == 12, ]
  if (nrow(r) > 0) {
    stars <- ifelse(r$p_value < 0.01, "***", ifelse(r$p_value < 0.05, "**",
                                                     ifelse(r$p_value < 0.10, "*", "")))
    cat(sprintf("    %s: %.2f pp (p=%.3f) %s\n", ct, r$coef, r$p_value, stars))
  }
}

# --- H2: FCI_ENDO_exCredit_exTCN LP (7 vars, no TCN, no credit) ---
cat("\nH2: FCI_ENDO_exCredit_exTCN credit LP...\n")

exTCN_lp_results <- data.frame()
if ("FCI_ENDO_exCredit_exTCN" %in% names(analysis_data)) {
  for (cred in c("Cred_Total", "Cred_Real_MN", "Cred_USD")) {
    res <- run_lp_standard(
      analysis_data, cred, "FCI_ENDO_exCredit_exTCN", LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    )
    if (nrow(res) > 0) {
      res$credit_type <- cred
      res$fci_type <- "FCI_ENDO_exCredit_exTCN"
      exTCN_lp_results <- rbind(exTCN_lp_results, res)
    }
  }

  cat("  exTCN LP at h=12:\n")
  for (ct in unique(exTCN_lp_results$credit_type)) {
    r <- exTCN_lp_results[exTCN_lp_results$credit_type == ct & exTCN_lp_results$horizon == 12, ]
    if (nrow(r) > 0) {
      stars <- ifelse(r$p_value < 0.01, "***", ifelse(r$p_value < 0.05, "**",
                                                       ifelse(r$p_value < 0.10, "*", "")))
      cat(sprintf("    %s: %.2f pp (p=%.3f) %s\n", ct, r$coef, r$p_value, stars))
    }
  }
} else {
  cat("  FCI_ENDO_exCredit_exTCN not available in data.\n")
}

# --- H3: Lag sensitivity (n_lags = 1 and 4) ---
cat("\nH3: Lag sensitivity for FCI_exCredit -> Cred_Total...\n")

lag_sensitivity <- data.frame()
for (nl in c(1, 2, 4)) {
  res <- run_lp_standard(
    analysis_data, "Cred_Total", "FCI_exCredit", 18, n_lags = nl,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  )
  if (nrow(res) > 0) {
    res$n_lags_used <- nl
    lag_sensitivity <- rbind(lag_sensitivity, res)
  }
}

cat("  Lag sensitivity at h=12 (FCI_exCredit -> Cred_Total):\n")
for (nl in c(1, 2, 4)) {
  r <- lag_sensitivity[lag_sensitivity$n_lags_used == nl & lag_sensitivity$horizon == 12, ]
  if (nrow(r) > 0) {
    cat(sprintf("    n_lags=%d: %.2f pp (p=%.3f)\n", nl, r$coef, r$p_value))
  }
}

# --- H4: Method-by-method LP (ZS, PCA, VAR, DFM separately) ---
cat("\nH4: Method-by-method credit LP...\n")

method_cols <- grep("^FCI_exCredit_(ZS|PCA|VAR|DFM)_norm$", names(analysis_data), value = TRUE)
method_lp_results <- data.frame()

for (mcol in method_cols) {
  method_name <- gsub("FCI_exCredit_|_norm", "", mcol)
  res <- run_lp_standard(
    analysis_data, "Cred_Total", mcol, LP_CONFIG$max_horizon, LP_CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  )
  if (nrow(res) > 0) {
    res$method <- method_name
    res$fci_var <- mcol
    method_lp_results <- rbind(method_lp_results, res)
  }
}

if (nrow(method_lp_results) > 0) {
  cat("  Method-by-method at h=12 (-> Cred_Total):\n")
  for (m in unique(method_lp_results$method)) {
    r <- method_lp_results[method_lp_results$method == m & method_lp_results$horizon == 12, ]
    if (nrow(r) > 0) {
      stars <- ifelse(r$p_value < 0.01, "***", ifelse(r$p_value < 0.05, "**",
                                                       ifelse(r$p_value < 0.10, "*", "")))
      cat(sprintf("    %s: %.2f pp (p=%.3f) %s\n", m, r$coef, r$p_value, stars))
    }
  }
} else {
  cat("  No method-level FCI columns available.\n")
}

# --- Save Part H results ---
partH_results <- rbind(
  if (nrow(rates_lp_results) > 0) cbind(rates_lp_results, spec = "Rates-only") else NULL,
  if (nrow(exTCN_lp_results) > 0) cbind(exTCN_lp_results, spec = "ENDO_exCredit_exTCN") else NULL
)
if (nrow(partH_results) > 0) {
  write.csv(partH_results, file.path(LP_CONFIG$output_dir, "LP_Robustness_Additional.csv"), row.names = FALSE)
  cat("\nSaved: LP_Robustness_Additional.csv\n")
}

write.csv(lag_sensitivity, file.path(LP_CONFIG$output_dir, "LP_Lag_Sensitivity.csv"), row.names = FALSE)
cat("Saved: LP_Lag_Sensitivity.csv\n")

if (nrow(method_lp_results) > 0) {
  write.csv(method_lp_results, file.path(LP_CONFIG$output_dir, "LP_Method_by_Method.csv"), row.names = FALSE)
  cat("Saved: LP_Method_by_Method.csv\n")
}

cat("\n================================================================================\n")
cat("LOCAL PROJECTIONS ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
cat("Plots: 50-59c\n")
cat("Output:", LP_CONFIG$output_dir, "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
