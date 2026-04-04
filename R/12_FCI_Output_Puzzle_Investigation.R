################################################################################
# FCI OUTPUT PUZZLE INVESTIGATION - EXTENDED ANALYSIS
################################################################################
#
# Project:      Financial Conditions Index - Output Puzzle Deep Dive
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Comprehensive investigation of why FCI predicts credit but not
#               aggregate output (IMAEP). Tests multiple hypotheses:
#
#   HYPOTHESIS 1: Agricultural Contamination
#     - Aggregate IMAEP includes crops (~7% GDP at constant 2014 prices);
#       manufacturing (~20% GDP) also shows no FCI response, partly because
#       a substantial agro-processing component depends on international financing
#     - TEST: Use IMAEP_SANB (excluding agriculture and binational entities)
#
#   HYPOTHESIS 2: Long-Run Cointegration
#     - Credit-output relationship may be long-term, requiring VECM
#     - Standard VAR in levels or differences may miss this
#     - TEST: Johansen cointegration test + VECM estimation
#
#   HYPOTHESIS 3: Transmission Channels
#     - FCI may affect output via investment/consumption, not directly
#     - TEST: LP for FBKf (investment) and Consumo separately
#
#   HYPOTHESIS 4: Sector Heterogeneity
#     - Financial conditions may matter more for specific sectors
#     - TEST: Compare effects across IMAEP, IMAEP_SANB, FBKf, Consumo
#
# References:
#   - Jordà (2005) - Local Projections
#   - Johansen (1988, 1991) - Cointegration Tests
#   - Lütkepohl (2005) - VECM Methodology
#   - Adrian et al. (2019) - Financial Conditions and Growth
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
  library(urca)        # For Johansen cointegration and VECM
  library(vars)        # For VAR estimation
  library(zoo)
  library(lubridate)   # For date manipulation
})

PUZZLE_CONFIG <- list(
  # Data
  data_file = "../data/FCI_data_1.xlsx",
  macro_sheet = "Datos_macro",

  # LP parameters
  max_horizon = 24,
  horizons_report = c(3, 6, 12, 18, 24),
  n_lags = 2,
  confidence_level = 0.90,

  # VECM parameters
  vecm_lag_max = 8,

  # Output variables to test
  output_vars = c("IMAEP", "IMAEP_SANB", "FBKf", "Consumo"),

  # Output
  output_dir = "../output"
)

# Labels for display
OUTPUT_LABELS <- c(
  "IMAEP" = "Aggregate GDP (IMAE-P)",
  "IMAEP_SANB" = "Non-Agro GDP (IMAE-P ex Agro/Binational)",
  "FBKf" = "Investment (FBKf)",
  "Consumo" = "Consumption"
)

set.seed(20260203)

cat("\n################################################################################\n")
cat("FCI OUTPUT PUZZLE INVESTIGATION\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# 2. DATA PREPARATION
################################################################################

cat("================================================================================\n")
cat("DATA PREPARATION\n")
cat("================================================================================\n\n")

# Load FCI from main script
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  source("01_FCI_Complete.R")
}

# Extract FCI indices (including endogeneity-corrected versions)
fci_data <- resultado_fci$all_indices %>%
  dplyr::select(fecha,
                FCI_COMP = FCI_COMP_AVG,
                FCI_ENDO = FCI_ENDO_AVG,
                FCI_EXO = FCI_EXO_AVG,
                FCI_exCredit = FCI_exCredit_AVG)

# Load macro data
cat("Loading macro data...\n")
macro_raw <- tryCatch({
  read_excel(PUZZLE_CONFIG$data_file, sheet = PUZZLE_CONFIG$macro_sheet)
}, error = function(e) {
  cat("Error loading macro sheet:", e$message, "\n")
  NULL
})

if (is.null(macro_raw)) {
  stop("Could not load macro data. Please check data file.")
}

fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

# Check which output variables are available
available_vars <- intersect(PUZZLE_CONFIG$output_vars, names(macro_raw))
missing_vars <- setdiff(PUZZLE_CONFIG$output_vars, names(macro_raw))

cat("\nAvailable output variables:\n")
for (v in available_vars) {
  n_obs <- sum(!is.na(macro_raw[[v]]))
  cat(sprintf("  %s: %d observations\n", OUTPUT_LABELS[v], n_obs))
}

if (length(missing_vars) > 0) {
  cat("\nMISSING VARIABLES (please add to Datos_macro sheet):\n")
  for (v in missing_vars) {
    cat(sprintf("  - %s (%s)\n", v, OUTPUT_LABELS[v]))
  }
  cat("\n*** Analysis will proceed with available variables ***\n")
}

# Prepare macro data with YoY growth rates
macro_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Calculate YoY growth for all available variables
for (var in available_vars) {
  macro_data <- macro_data %>%
    mutate(!!paste0(var, "_yoy") := (!!sym(var) / lag(!!sym(var), 12) - 1) * 100)
}

# Also calculate credit growth from main sheet for VECM
datos_raw <- read_excel(PUZZLE_CONFIG$data_file)
fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

credit_data <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  mutate(
    Credit_yoy = (Creditos_Sector_privado_totales /
                   lag(Creditos_Sector_privado_totales, 12) - 1) * 100
  ) %>%
  dplyr::select(fecha, Credit_yoy)

# Compute IPC_yoy and Cred_Real_yoy for use as LP controls (consistent with Script 05)
ipc_credit_data <- tryCatch({
  ipc_raw <- macro_raw %>%
    rename(fecha = !!sym(fecha_col_macro)) %>%
    mutate(fecha = as.Date(fecha)) %>%
    arrange(fecha)

  # IPC YoY inflation
  ipc_yoy_df <- ipc_raw %>%
    mutate(IPC_yoy = (IPC / lag(IPC, 12) - 1) * 100) %>%
    dplyr::select(fecha, IPC_yoy)

  # Real (deflated) credit YoY growth
  cred_real_df <- ipc_raw %>%
    mutate(Cred_Real_yoy = (Creditos_deflactados /
                             lag(Creditos_deflactados, 12) - 1) * 100) %>%
    dplyr::select(fecha, Cred_Real_yoy)

  ipc_yoy_df %>% inner_join(cred_real_df, by = "fecha")
}, error = function(e) {
  cat("Warning: Could not compute IPC_yoy / Cred_Real_yoy:", e$message, "\n")
  NULL
})

# Merge all data
analysis_data <- fci_data %>%
  inner_join(macro_data, by = "fecha") %>%
  inner_join(credit_data, by = "fecha")

if (!is.null(ipc_credit_data)) {
  analysis_data <- analysis_data %>%
    left_join(ipc_credit_data, by = "fecha")
  cat("  Macro controls added: IPC_yoy, Cred_Real_yoy\n")
}

analysis_data <- analysis_data %>%
  arrange(fecha) %>%
  na.omit()

cat("\nAnalysis dataset:\n")
cat("  Observations:", nrow(analysis_data), "\n")
cat("  Period:", format(min(analysis_data$fecha)), "to",
    format(max(analysis_data$fecha)), "\n\n")


################################################################################
# 3. LOCAL PROJECTION FUNCTIONS
################################################################################

#' Standard Local Projection with HAC standard errors
run_lp_standard <- function(data, y_var, fci_var, max_h, n_lags = 2,
                            conf_level = 0.90, control_vars = NULL) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - conf_level) / 2)

  for (h in 1:max_h) {
    # Create forward and lagged variables
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
    reg_data <- data_h %>%
      dplyr::select(y_fwd, !!sym(fci_var), all_of(lag_vars), all_of(control_cols)) %>%
      na.omit()

    if (nrow(reg_data) < 30) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    idx <- which(rownames(coef_test) == fci_var)

    results <- rbind(results, data.frame(
      horizon = h,
      coef = coef_test[idx, 1],
      se = coef_test[idx, 2],
      t_stat = coef_test[idx, 3],
      p_value = coef_test[idx, 4],
      ci_lower = coef_test[idx, 1] - z_crit * coef_test[idx, 2],
      ci_upper = coef_test[idx, 1] + z_crit * coef_test[idx, 2],
      n_obs = nrow(reg_data),
      r_squared = summary(model)$r.squared
    ))
  }

  return(results)
}


################################################################################
# 4. HYPOTHESIS 1: AGRICULTURAL CONTAMINATION
################################################################################

cat("================================================================================\n")
cat("HYPOTHESIS 1: AGRICULTURAL CONTAMINATION\n")
cat("================================================================================\n\n")

cat("RATIONALE:\n")
cat("  Aggregate IMAEP includes crops (~7% GDP at constant 2014 prices)\n")
cat("  which depend on rainfall, global prices, and international financing.\n")
cat("  Manufacturing (~20% GDP) also shows no significant FCI response,\n")
cat("  partly because a substantial agro-processing component (soybean oil/flour,\n")
cat("  alcohol, food products) inherits the primary sector's financing structure.\n\n")
cat("  Together with binational electricity (~9%) and taxes on products (~8%),\n")
cat("  nearly half of GDP is structurally insensitive to domestic financial conditions.\n\n")
cat("  TEST: Compare LP results for IMAEP vs IMAEP_SANB\n\n")

# Run LP for all available output variables
output_lp_results <- list()

for (var in available_vars) {
  y_var <- paste0(var, "_yoy")
  if (!y_var %in% names(analysis_data)) next

  cat(sprintf("Running LP: FCI_COMP -> %s ... ", OUTPUT_LABELS[var]))

  # Use macro controls consistent with Script 05: IPC_yoy + Cred_Real_yoy for output regressions
  output_ctrl <- if (all(c("IPC_yoy", "Cred_Real_yoy") %in% names(analysis_data))) {
    c("IPC_yoy", "Cred_Real_yoy")
  } else {
    NULL
  }

  output_lp_results[[var]] <- run_lp_standard(
    analysis_data, y_var, "FCI_COMP",
    PUZZLE_CONFIG$max_horizon, PUZZLE_CONFIG$n_lags,
    control_vars = output_ctrl
  ) %>%
    mutate(
      variable = var,
      variable_label = OUTPUT_LABELS[var]
    )

  cat("done\n")
}

all_output_lp <- bind_rows(output_lp_results)

# Summary table
cat("\n================================================================================\n")
cat("LOCAL PROJECTIONS RESULTS: FCI EFFECT ON OUTPUT MEASURES\n")
cat("================================================================================\n\n")

cat("Response (pp) to 1 SD FCI tightening | *** p<0.01, ** p<0.05, * p<0.10\n\n")

for (var in available_vars) {
  var_res <- all_output_lp %>%
    filter(variable == var, horizon %in% PUZZLE_CONFIG$horizons_report)

  if (nrow(var_res) == 0) next

  cat(sprintf("%s:\n", OUTPUT_LABELS[var]))
  cat(sprintf("  %3s  %10s  %8s  %6s  %6s\n", "h", "Coef", "SE", "p-val", ""))

  for (i in 1:nrow(var_res)) {
    r <- var_res[i,]
    stars <- case_when(
      r$p_value < 0.01 ~ "***",
      r$p_value < 0.05 ~ "**",
      r$p_value < 0.10 ~ "*",
      TRUE ~ ""
    )
    cat(sprintf("  %3d  %+10.3f  %8.3f  %6.3f  %s\n",
                r$horizon, r$coef, r$se, r$p_value, stars))
  }
  cat("\n")
}

# Key comparison: IMAEP vs IMAEP_SANB at horizon 12
if ("IMAEP" %in% available_vars && "IMAEP_SANB" %in% available_vars) {
  cat("\n================================================================================\n")
  cat("KEY COMPARISON: IMAEP vs IMAEP_SANB (Non-Agricultural GDP)\n")
  cat("================================================================================\n\n")

  imaep_12 <- all_output_lp %>%
    filter(variable == "IMAEP", horizon == 12)
  sanb_12 <- all_output_lp %>%
    filter(variable == "IMAEP_SANB", horizon == 12)

  if (nrow(imaep_12) > 0 && nrow(sanb_12) > 0) {
    cat(sprintf("At 12-month horizon:\n"))
    cat(sprintf("  IMAEP (aggregate):     coef = %+.3f, p = %.3f\n",
                imaep_12$coef, imaep_12$p_value))
    cat(sprintf("  IMAEP_SANB (non-agro): coef = %+.3f, p = %.3f\n",
                sanb_12$coef, sanb_12$p_value))

    cat("\nINTERPRETATION:\n")
    if (sanb_12$p_value < 0.10 && imaep_12$p_value >= 0.10) {
      cat("  >>> HYPOTHESIS CONFIRMED: FCI significantly affects non-agricultural GDP\n")
      cat("      but NOT aggregate GDP. Agricultural sector masks the credit channel.\n")
    } else if (sanb_12$p_value < 0.10 && imaep_12$p_value < 0.10) {
      cat("  >>> FCI affects BOTH aggregate and non-agricultural GDP.\n")
      cat("      Effect on non-agro is stronger (coef: %.3f vs %.3f)\n",
          abs(sanb_12$coef), abs(imaep_12$coef))
    } else if (sanb_12$p_value >= 0.10 && imaep_12$p_value >= 0.10) {
      cat("  >>> Neither effect is significant. Puzzle may have other explanations.\n")
    } else {
      cat("  >>> Mixed results - requires further investigation.\n")
    }
  }
}


################################################################################
# 5. HYPOTHESIS 2: LONG-RUN COINTEGRATION (VECM)
################################################################################

cat("\n================================================================================\n")
cat("HYPOTHESIS 2: LONG-RUN COINTEGRATION (VECM ANALYSIS)\n")
cat("================================================================================\n\n")

cat("RATIONALE:\n")
cat("  - Credit-output relationship may be long-term\n")
cat("  - Standard VAR in levels/differences may miss cointegration\n")
cat("  - VECM captures both short-run dynamics and long-run equilibrium\n\n")

# Prepare data for VECM
# Use levels for cointegration test
vecm_vars <- c("Credit_yoy", "IMAEP_yoy")

# Check if IMAEP_SANB is available and add it
if ("IMAEP_SANB_yoy" %in% names(analysis_data)) {
  vecm_vars <- c(vecm_vars, "IMAEP_SANB_yoy")
}

# Add FCI
vecm_data <- analysis_data %>%
  dplyr::select(fecha, FCI_COMP, all_of(vecm_vars)) %>%
  na.omit()

cat("VECM Variables:\n")
for (v in c("FCI_COMP", vecm_vars)) {
  cat(sprintf("  - %s\n", v))
}
cat(sprintf("\nObservations: %d\n\n", nrow(vecm_data)))

# Convert to time series matrix
vecm_ts <- as.matrix(vecm_data[, c("FCI_COMP", vecm_vars)])

#-----------------------------------------------------------------------------
# Step 1: Determine optimal lag length
#-----------------------------------------------------------------------------
cat("Step 1: Optimal Lag Selection\n")

tryCatch({
  # Use VARselect for lag selection
  lag_select <- vars::VARselect(vecm_ts, lag.max = PUZZLE_CONFIG$vecm_lag_max, type = "const")
  cat("  Information Criteria:\n")
  print(lag_select$selection)

  bic_lag <- lag_select$selection["SC(n)"]
  optimal_lag <- max(bic_lag, 2)  # Johansen test requires K >= 2
  cat(sprintf("\n  Selected lag (BIC): %d", bic_lag))
  if (bic_lag < 2) cat(" (adjusted to 2 for Johansen test)")
  cat("\n\n")
}, error = function(e) {
  cat("  Error in lag selection:", e$message, "\n")
  optimal_lag <<- 2
})

#-----------------------------------------------------------------------------
# Step 2: Johansen Cointegration Test
#-----------------------------------------------------------------------------
cat("Step 2: Johansen Cointegration Test\n\n")

# Ensure K is valid for ca.jo (must be >= 2)
K_johansen <- max(min(optimal_lag, 4), 2)
cat(sprintf("Using K = %d lags for Johansen test\n\n", K_johansen))

tryCatch({
  # Test with trace statistic
  johansen_trace <- urca::ca.jo(vecm_ts, type = "trace", ecdet = "const",
                                 K = K_johansen)

  # Test with eigenvalue statistic
  johansen_eigen <- urca::ca.jo(vecm_ts, type = "eigen", ecdet = "const",
                                 K = K_johansen)

  cat("Trace Test Results:\n")
  print(summary(johansen_trace))

  cat("\n\nEigenvalue Test Results:\n")
  print(summary(johansen_eigen))

  # Determine cointegration rank
  trace_stats <- johansen_trace@teststat
  trace_crit <- johansen_trace@cval[, "5pct"]

  # Count how many trace stats exceed critical values
  coint_rank <- sum(trace_stats > trace_crit)

  n_vars <- ncol(vecm_ts)
  cat(sprintf("\n\nCointegration Rank (5%% level): r = %d (out of %d variables)\n", coint_rank, n_vars))

  # Handle edge case: rank = n_vars means variables are stationary, no cointegration
  if (coint_rank >= n_vars) {
    cat("\n>>> FULL RANK: All variables appear stationary\n")
    cat("    This indicates no cointegrating relationships exist.\n")
    cat("    VECM not applicable - standard VAR in levels is appropriate.\n")
    cat("    Interpretation: Variables move independently (no long-run equilibrium).\n\n")
    vecm_results <- list(rank = coint_rank, full_rank = TRUE, error = "Full rank - no cointegration")
  } else if (coint_rank > 0) {
    cat("\n>>> COINTEGRATION DETECTED\n")
    cat("    There exists a long-run equilibrium relationship between\n")
    cat("    FCI, credit, and output.\n")

    #---------------------------------------------------------------------------
    # Step 3: Estimate VECM
    #---------------------------------------------------------------------------
    cat("\n================================================================================\n")
    cat("Step 3: VECM Estimation\n")
    cat("================================================================================\n\n")

    # Convert Johansen to VECM
    vecm_model <- urca::cajorls(johansen_trace, r = coint_rank)

    cat("Cointegrating Vector(s):\n")
    print(vecm_model$beta)

    cat("\nError Correction Coefficients (Alpha - adjustment speeds):\n")
    print(vecm_model$rlm$coefficients[1:coint_rank, ])

    # Extract adjustment coefficients
    alpha <- vecm_model$rlm$coefficients[1:coint_rank, ]

    cat("\nINTERPRETATION:\n")
    cat("  The alpha coefficients show how each variable adjusts to\n")
    cat("  deviations from long-run equilibrium:\n\n")

    var_names <- colnames(vecm_ts)
    for (i in 1:ncol(alpha)) {
      cat(sprintf("  %s: alpha = %.4f\n", var_names[i], alpha[1, i]))
      if (abs(alpha[1, i]) > 0.1) {
        if (alpha[1, i] < 0) {
          cat("    -> Adjusts TOWARD equilibrium (stabilizing)\n")
        } else {
          cat("    -> Adjusts AWAY FROM equilibrium (destabilizing)\n")
        }
      } else {
        cat("    -> Weak adjustment (nearly weakly exogenous)\n")
      }
    }

    #---------------------------------------------------------------------------
    # Step 4: Analyze Long-Run Relationship
    #---------------------------------------------------------------------------
    cat("\n================================================================================\n")
    cat("LONG-RUN EQUILIBRIUM RELATIONSHIP\n")
    cat("================================================================================\n\n")

    beta <- vecm_model$beta
    cat("Normalized cointegrating vector (first variable = 1):\n")
    print(beta / beta[1, 1])

    cat("\nThis implies the long-run equilibrium:\n")
    cat("  FCI + ")
    for (i in 2:nrow(beta)) {
      coef <- beta[i, 1] / beta[1, 1]
      cat(sprintf("%.3f*%s ", coef, rownames(beta)[i]))
      if (i < nrow(beta)) cat("+ ")
    }
    cat("= constant\n")

    # Save VECM results
    vecm_results <- list(
      rank = coint_rank,
      beta = beta,
      alpha = alpha,
      johansen = johansen_trace
    )

  } else {
    cat("\n>>> NO COINTEGRATION DETECTED\n")
    cat("    The variables do not share a long-run equilibrium.\n")
    cat("    Standard VAR in differences may be appropriate.\n")

    vecm_results <- list(rank = 0)
  }

}, error = function(e) {
  cat("Error in VECM analysis:", e$message, "\n")
  vecm_results <- list(rank = NA, error = e$message)
})


################################################################################
# 6. HYPOTHESIS 3: TRANSMISSION CHANNELS (Investment & Consumo)
################################################################################

if ("FBKf" %in% available_vars || "Consumo" %in% available_vars) {

  cat("\n================================================================================\n")
  cat("HYPOTHESIS 3: TRANSMISSION CHANNELS\n")
  cat("================================================================================\n\n")

  cat("RATIONALE:\n")
  cat("  FCI may affect output indirectly through:\n")
  cat("    - Investment (FBKf): Credit-constrained firms reduce investment\n")
  cat("    - Consumo: Households reduce spending when credit tightens\n\n")

  # Already computed in Section 4, just highlight results here
  cat("Results (from LP analysis above):\n\n")

  channel_vars <- intersect(c("FBKf", "Consumo"), available_vars)

  for (var in channel_vars) {
    var_12 <- all_output_lp %>%
      filter(variable == var, horizon == 12)

    if (nrow(var_12) > 0) {
      sig <- ifelse(var_12$p_value < 0.10, "SIGNIFICANT", "not significant")
      cat(sprintf("  %s: coef = %+.3f, p = %.3f (%s)\n",
                  OUTPUT_LABELS[var], var_12$coef, var_12$p_value, sig))
    }
  }

  cat("\nINTERPRETATION:\n")

  # Check which channels are significant
  sig_channels <- all_output_lp %>%
    filter(variable %in% channel_vars, horizon == 12, p_value < 0.10)

  if (nrow(sig_channels) > 0) {
    cat("  >>> FCI affects output through:\n")
    for (i in 1:nrow(sig_channels)) {
      cat(sprintf("      - %s\n", OUTPUT_LABELS[sig_channels$variable[i]]))
    }
  } else {
    cat("  >>> No significant transmission through investment or consumption.\n")
  }
}


################################################################################
# 6B. HYPOTHESIS 4: CREDIT CHANNEL MEDIATION
################################################################################
#
# The user correctly identifies that transmission should be sequential:
#   FCI tightening (t) → Credit contraction (t+k) → Output reduction (t+k+j)
#
# This section tests:
#   1. Credit → Output LP (does credit growth predict output growth?)
#   2. FCI → Credit LP (confirm the first link)
#   3. Mediation test: Does including credit reduce the FCI→Output coefficient?
#   4. Granger causality: Credit ↔ Output
#   5. Sequential timing analysis

cat("\n================================================================================\n")
cat("HYPOTHESIS 4: CREDIT CHANNEL MEDIATION (TRANSMISSION CHAIN)\n")
cat("================================================================================\n\n")

cat("RATIONALE:\n")
cat("  The transmission mechanism should be SEQUENTIAL with lags:\n\n")
cat("    FCI ↑ (t=0) → Credit ↓ (t+3 to t+6) → Output ↓ (t+6 to t+12)\n\n")
cat("  If FCI affects output only through credit, we should see:\n")
cat("    1. FCI strongly predicts credit (already established)\n")
cat("    2. Credit strongly predicts output (TEST THIS)\n")
cat("    3. Adding credit to FCI→Output regression reduces FCI coefficient\n\n")

#-----------------------------------------------------------------------------
# 6B.1: Credit → Output LP
#-----------------------------------------------------------------------------
cat("--------------------------------------------------------------------------------\n")
cat("TEST 1: Does Credit Growth Predict Output Growth?\n")
cat("--------------------------------------------------------------------------------\n\n")

# Check if Credit_yoy is in analysis data
if (!"Credit_yoy" %in% names(analysis_data)) {
  # Merge credit data if not already present
  analysis_data <- analysis_data %>%
    left_join(credit_data, by = "fecha")
}

# Run LP: Credit → Output for each output variable
credit_output_lp <- data.frame()

for (var in available_vars) {
  cat(sprintf("Running LP: Credit_yoy -> %s ... ", OUTPUT_LABELS[var]))

  for (h in 0:PUZZLE_CONFIG$max_horizon) {
    # Create lead of dependent variable
    dep_var <- paste0(var, "_yoy")
    if (!dep_var %in% names(analysis_data)) next

    temp_data <- analysis_data %>%
      mutate(y_lead = lead(!!sym(dep_var), h)) %>%
      filter(!is.na(y_lead), !is.na(Credit_yoy))

    if (nrow(temp_data) < 50) next

    # LP regression: y_{t+h} = α + β*Credit_yoy_t + γ*y_t + controls + ε
    formula_str <- paste0("y_lead ~ Credit_yoy + ", dep_var)
    for (lag in 1:PUZZLE_CONFIG$n_lags) {
      formula_str <- paste0(formula_str, " + lag(Credit_yoy, ", lag, ")")
      formula_str <- paste0(formula_str, " + lag(", dep_var, ", ", lag, ")")
    }

    tryCatch({
      model <- lm(as.formula(formula_str), data = temp_data)

      # Newey-West standard errors
      nw_vcov <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
      coef_test <- lmtest::coeftest(model, vcov = nw_vcov)

      # Extract Credit_yoy coefficient
      credit_row <- which(rownames(coef_test) == "Credit_yoy")
      if (length(credit_row) > 0) {
        coef_val <- coef_test[credit_row, 1]
        se_val <- coef_test[credit_row, 2]
        p_val <- coef_test[credit_row, 4]

        credit_output_lp <- rbind(credit_output_lp, data.frame(
          output_var = var,
          output_label = OUTPUT_LABELS[var],
          horizon = h,
          coef = coef_val,
          se = se_val,
          p_value = p_val,
          ci_lower = coef_val - 1.645 * se_val,
          ci_upper = coef_val + 1.645 * se_val,
          n_obs = nrow(temp_data),
          r_squared = summary(model)$r.squared
        ))
      }
    }, error = function(e) NULL)
  }
  cat("done\n")
}

# Report results
cat("\nRESULTS: Effect of 1 pp Credit Growth on Output (YoY)\n")
cat("*** p<0.01, ** p<0.05, * p<0.10\n\n")

for (var in available_vars) {
  var_results <- credit_output_lp %>%
    filter(output_var == var, horizon %in% PUZZLE_CONFIG$horizons_report)

  if (nrow(var_results) > 0) {
    cat(sprintf("%s:\n", OUTPUT_LABELS[var]))
    for (i in 1:nrow(var_results)) {
      r <- var_results[i, ]
      stars <- case_when(
        r$p_value < 0.01 ~ "***",
        r$p_value < 0.05 ~ "**",
        r$p_value < 0.10 ~ "*",
        TRUE ~ ""
      )
      cat(sprintf("    h=%2d: %+.4f (SE=%.4f, p=%.3f) %s\n",
                  r$horizon, r$coef, r$se, r$p_value, stars))
    }
    cat("\n")
  }
}

#-----------------------------------------------------------------------------
# 6B.2: Granger Causality - Credit ↔ Output
#-----------------------------------------------------------------------------
cat("--------------------------------------------------------------------------------\n")
cat("TEST 2: Granger Causality Tests (Credit ↔ Output)\n")
cat("--------------------------------------------------------------------------------\n\n")

granger_results <- data.frame()

for (var in available_vars) {
  dep_var <- paste0(var, "_yoy")
  if (!dep_var %in% names(analysis_data)) next

  # Prepare data for Granger test
  granger_data <- analysis_data %>%
    dplyr::select(all_of(c("Credit_yoy", dep_var))) %>%
    na.omit()

  if (nrow(granger_data) < 50) next

  tryCatch({
    # VAR model
    var_model <- vars::VAR(granger_data, p = 4, type = "const")

    # Granger test: Credit → Output
    granger_credit_to_output <- vars::causality(var_model, cause = "Credit_yoy")
    p_credit_to_output <- granger_credit_to_output$Granger$p.value

    # Granger test: Output → Credit
    granger_output_to_credit <- vars::causality(var_model, cause = dep_var)
    p_output_to_credit <- granger_output_to_credit$Granger$p.value

    granger_results <- rbind(granger_results, data.frame(
      output_var = var,
      output_label = OUTPUT_LABELS[var],
      p_credit_to_output = p_credit_to_output,
      p_output_to_credit = p_output_to_credit,
      credit_causes_output = ifelse(p_credit_to_output < 0.10, "YES", "NO"),
      output_causes_credit = ifelse(p_output_to_credit < 0.10, "YES", "NO")
    ))

  }, error = function(e) NULL)
}

if (nrow(granger_results) > 0) {
  cat("Granger Causality Results (4 lags, 10% significance):\n\n")
  cat(sprintf("%-40s %15s %15s\n", "Output Variable", "Credit→Output", "Output→Credit"))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (i in 1:nrow(granger_results)) {
    r <- granger_results[i, ]
    cat(sprintf("%-40s %12s (p=%.3f) %12s (p=%.3f)\n",
                r$output_label, r$credit_causes_output, r$p_credit_to_output,
                r$output_causes_credit, r$p_output_to_credit))
  }
  cat("\n")
}

#-----------------------------------------------------------------------------
# 6B.3: Mediation Test - Does Credit Mediate FCI → Output?
#-----------------------------------------------------------------------------
cat("--------------------------------------------------------------------------------\n")
cat("TEST 3: Mediation Analysis (Does Credit Mediate FCI → Output?)\n")
cat("--------------------------------------------------------------------------------\n\n")

cat("LOGIC:\n")
cat("  If FCI affects output ONLY through credit:\n")
cat("    - Direct effect (FCI → Output) should be zero when controlling for credit\n")
cat("    - Indirect effect (FCI → Credit → Output) should be significant\n\n")

mediation_results <- data.frame()

for (var in available_vars) {
  dep_var <- paste0(var, "_yoy")
  if (!dep_var %in% names(analysis_data)) next

  # Test at h=12 (where we expect to see output effects)
  h <- 12

  # Use FCI_exCredit for mediation (credit is mediator, so FCI containing credit biases the test)
  fci_med <- if ("FCI_exCredit" %in% names(analysis_data)) "FCI_exCredit" else "FCI_COMP"

  temp_data <- analysis_data %>%
    mutate(
      y_lead = lead(!!sym(dep_var), h),
      credit_lag6 = lag(Credit_yoy, 6)  # Credit lagged 6 months (FCI→Credit timing)
    ) %>%
    filter(!is.na(y_lead), !is.na(!!sym(fci_med)), !is.na(credit_lag6))

  if (nrow(temp_data) < 50) next

  tryCatch({
    # Model 1: FCI only (total effect)
    model_total <- lm(as.formula(paste0("y_lead ~ ", fci_med, " + ", dep_var)), data = temp_data)
    nw_total <- sandwich::NeweyWest(model_total, lag = h + 1, prewhite = FALSE)
    coef_total <- lmtest::coeftest(model_total, vcov = nw_total)
    fci_total <- coef_total[fci_med, 1]
    se_total <- coef_total[fci_med, 2]
    p_total <- coef_total[fci_med, 4]

    # Model 2: FCI + Credit (direct effect, controlling for mediator)
    model_direct <- lm(as.formula(paste0("y_lead ~ ", fci_med, " + credit_lag6 + ", dep_var)), data = temp_data)
    nw_direct <- sandwich::NeweyWest(model_direct, lag = h + 1, prewhite = FALSE)
    coef_direct <- lmtest::coeftest(model_direct, vcov = nw_direct)
    fci_direct <- coef_direct[fci_med, 1]
    se_direct <- coef_direct[fci_med, 2]
    p_direct <- coef_direct[fci_med, 4]
    credit_effect <- coef_direct["credit_lag6", 1]
    credit_p <- coef_direct["credit_lag6", 4]

    # Percentage reduction in FCI coefficient
    pct_reduction <- ifelse(abs(fci_total) > 0.001,
                            (1 - abs(fci_direct) / abs(fci_total)) * 100,
                            NA)

    mediation_results <- rbind(mediation_results, data.frame(
      output_var = var,
      output_label = OUTPUT_LABELS[var],
      horizon = h,
      fci_total = fci_total,
      fci_total_se = se_total,
      fci_total_p = p_total,
      fci_direct = fci_direct,
      fci_direct_se = se_direct,
      fci_direct_p = p_direct,
      credit_effect = credit_effect,
      credit_p = credit_p,
      pct_reduction = pct_reduction,
      n_obs = nrow(temp_data)
    ))

  }, error = function(e) NULL)
}

if (nrow(mediation_results) > 0) {
  cat("Mediation Results at h=12:\n\n")
  cat(sprintf("%-35s %12s %12s %12s %12s\n",
              "Output Variable", "FCI Total", "FCI Direct", "Credit(t-6)", "Mediation %"))
  cat(paste(rep("-", 90), collapse = ""), "\n")

  for (i in 1:nrow(mediation_results)) {
    r <- mediation_results[i, ]
    fci_total_sig <- ifelse(r$fci_total_p < 0.10, "*", "")
    fci_direct_sig <- ifelse(r$fci_direct_p < 0.10, "*", "")
    credit_sig <- ifelse(r$credit_p < 0.10, "*", "")

    cat(sprintf("%-35s %+10.3f%s %+10.3f%s %+10.4f%s %10.1f%%\n",
                r$output_label, r$fci_total, fci_total_sig,
                r$fci_direct, fci_direct_sig,
                r$credit_effect, credit_sig,
                ifelse(is.na(r$pct_reduction), 0, r$pct_reduction)))
  }
  cat("\nNote: * indicates p < 0.10\n")
  cat("Mediation % = reduction in FCI coefficient when controlling for credit\n\n")

  # Interpretation
  cat("INTERPRETATION:\n")
  avg_mediation <- mean(mediation_results$pct_reduction, na.rm = TRUE)
  any_credit_sig <- any(mediation_results$credit_p < 0.10)

  if (any_credit_sig && avg_mediation > 30) {
    cat("  >>> CREDIT CHANNEL CONFIRMED: Credit mediates FCI → Output relationship.\n")
    cat(sprintf("      Average mediation: %.1f%% of FCI effect works through credit.\n", avg_mediation))
  } else if (any_credit_sig) {
    cat("  >>> PARTIAL MEDIATION: Credit has some mediating role.\n")
    cat(sprintf("      Average mediation: %.1f%%\n", avg_mediation))
  } else {
    cat("  >>> NO MEDIATION DETECTED: Credit does not significantly affect output.\n")
    cat("      FCI may affect output through other channels (expectations, risk premia).\n")
  }
}

#-----------------------------------------------------------------------------
# 6B.4: Sequential Timing Analysis
#-----------------------------------------------------------------------------
cat("\n--------------------------------------------------------------------------------\n")
cat("TEST 4: Sequential Timing Analysis (Lag Structure)\n")
cat("--------------------------------------------------------------------------------\n\n")

cat("Testing the transmission chain timing:\n")
cat("  Step 1: FCI(t) → Credit(t+k) - At what lag k does FCI affect credit?\n")
cat("  Step 2: Credit(t) → Output(t+j) - At what lag j does credit affect output?\n")
cat("  Expected: Total lag = k + j\n\n")

# FCI → Credit peak effect (using FCI_exCredit for endogeneity correction)
fci_cr_col <- if ("FCI_exCredit" %in% names(analysis_data)) "FCI_exCredit" else "FCI_COMP"
fci_credit_lp <- data.frame()
for (h in 0:PUZZLE_CONFIG$max_horizon) {
  temp_data <- analysis_data %>%
    mutate(credit_lead = lead(Credit_yoy, h)) %>%
    filter(!is.na(credit_lead), !is.na(!!sym(fci_cr_col)))

  tryCatch({
    fml <- as.formula(paste0("credit_lead ~ ", fci_cr_col,
                             " + Credit_yoy + lag(", fci_cr_col, ", 1) + lag(Credit_yoy, 1)"))
    model <- lm(fml, data = temp_data)
    nw_vcov <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = nw_vcov)

    fci_credit_lp <- rbind(fci_credit_lp, data.frame(
      horizon = h,
      coef = coef_test[fci_cr_col, 1],
      se = coef_test[fci_cr_col, 2],
      p_value = coef_test[fci_cr_col, 4]
    ))
  }, error = function(e) NULL)
}

# Find peak FCI→Credit effect
if (nrow(fci_credit_lp) > 0) {
  fci_credit_sig <- fci_credit_lp %>% filter(p_value < 0.10)
  if (nrow(fci_credit_sig) > 0) {
    peak_fci_credit <- fci_credit_sig[which.min(fci_credit_sig$coef), ]
    cat(sprintf("FCI → Credit: Peak effect at h=%d (coef=%.2f, p=%.3f)\n",
                peak_fci_credit$horizon, peak_fci_credit$coef, peak_fci_credit$p_value))
  } else {
    cat("FCI → Credit: No significant effect found (unexpected)\n")
  }
}

# Credit → IMAEP_SANB peak effect (use non-agro for cleaner signal)
if ("IMAEP_SANB_yoy" %in% names(analysis_data)) {
  credit_output_timing <- credit_output_lp %>%
    filter(output_var == "IMAEP_SANB", p_value < 0.10)

  if (nrow(credit_output_timing) > 0) {
    peak_credit_output <- credit_output_timing[which.max(abs(credit_output_timing$coef)), ]
    cat(sprintf("Credit → IMAEP_SANB: Peak effect at h=%d (coef=%.4f, p=%.3f)\n",
                peak_credit_output$horizon, peak_credit_output$coef, peak_credit_output$p_value))

    if (exists("peak_fci_credit")) {
      total_lag <- peak_fci_credit$horizon + peak_credit_output$horizon
      cat(sprintf("\n>>> EXPECTED TOTAL LAG (FCI → Output): %d months\n", total_lag))
      cat("    This explains why direct FCI → Output effects appear weak at standard horizons.\n")
    }
  } else {
    cat("Credit → IMAEP_SANB: No significant effect found at any horizon.\n")
    cat("    This suggests credit growth may not be the primary transmission channel.\n")
  }
}

# Save credit channel results
write.csv(credit_output_lp,
          file.path(PUZZLE_CONFIG$output_dir, "Credit_Output_LP_Results.csv"),
          row.names = FALSE)
cat("\nSaved: Credit_Output_LP_Results.csv\n")

if (nrow(granger_results) > 0) {
  write.csv(granger_results,
            file.path(PUZZLE_CONFIG$output_dir, "Credit_Output_Granger.csv"),
            row.names = FALSE)
  cat("Saved: Credit_Output_Granger.csv\n")
}

if (nrow(mediation_results) > 0) {
  write.csv(mediation_results,
            file.path(PUZZLE_CONFIG$output_dir, "Credit_Output_Mediation.csv"),
            row.names = FALSE)
  cat("Saved: Credit_Output_Mediation.csv\n")
}


################################################################################
# 6C. LP vs SVAR COMPARISON: CREDIT-OUTPUT TRANSMISSION
################################################################################
#
# Compare Local Projections with Structural VAR to verify robustness
# of the credit-output relationship. Both methods should give similar
# impulse responses if the relationship is well-identified.
#
# Identification: Cholesky (recursive) ordering
#   1. FCI (financial conditions - most exogenous)
#   2. Credit_yoy (responds to FCI, affects output)
#   3. Output_yoy (most endogenous)

cat("\n================================================================================\n")
cat("LP vs SVAR COMPARISON: CREDIT-OUTPUT TRANSMISSION\n")
cat("================================================================================\n\n")

cat("METHODOLOGY:\n")
cat("  LOCAL PROJECTIONS (Jordà 2005):\n")
cat("    - Direct estimation of h-step ahead responses\n")
cat("    - Robust to misspecification of VAR dynamics\n")
cat("    - Allows for non-linearities and state-dependence\n\n")

cat("  STRUCTURAL VAR (Sims 1980):\n")
cat("    - Joint estimation of system dynamics\n")
cat("    - More efficient if VAR is correctly specified\n")
cat("    - Cholesky identification: FCI → Credit → Output\n\n")

cat("  If both methods agree, the credit-output relationship is robust.\n\n")

#-----------------------------------------------------------------------------
# 6C.1: SVAR Estimation - FCI → Credit → Output
#-----------------------------------------------------------------------------
cat("--------------------------------------------------------------------------------\n")
cat("SVAR ESTIMATION: FCI → Credit → Output\n")
cat("--------------------------------------------------------------------------------\n\n")

# Store SVAR results for comparison
svar_results <- list()
lp_vs_svar_comparison <- data.frame()

# Test for each output variable
for (output_var in c("IMAEP", "IMAEP_SANB")) {
  dep_var <- paste0(output_var, "_yoy")
  if (!dep_var %in% names(analysis_data)) next

  cat(sprintf("=== %s ===\n\n", OUTPUT_LABELS[output_var]))

  # Prepare data for SVAR
  # Ordering: FCI → Credit_yoy → Output_yoy (Cholesky)
  # Use FCI_exCredit to address mechanical correlation (credit is in the SVAR system)
  svar_fci_col <- if ("FCI_exCredit" %in% names(analysis_data)) "FCI_exCredit" else "FCI_COMP"
  svar_data <- analysis_data %>%
    dplyr::select(fecha, FCI_COMP = !!sym(svar_fci_col), Credit_yoy, !!sym(dep_var)) %>%
    na.omit()

  if (nrow(svar_data) < 60) {
    cat("  Insufficient observations for SVAR. Skipping.\n\n")
    next
  }

  cat(sprintf("  Sample: %d observations (%s to %s)\n",
              nrow(svar_data),
              format(min(svar_data$fecha), "%Y-%m"),
              format(max(svar_data$fecha), "%Y-%m")))

  # Create time series matrix
  svar_ts <- ts(svar_data[, c("FCI_COMP", "Credit_yoy", dep_var)],
                start = c(year(min(svar_data$fecha)), month(min(svar_data$fecha))),
                frequency = 12)

  tryCatch({
    #---------------------------------------------------------------------------
    # Step 1: Estimate reduced-form VAR
    #---------------------------------------------------------------------------

    # Select optimal lag using BIC
    lag_select <- vars::VARselect(svar_ts, lag.max = 12, type = "const")
    optimal_p <- max(lag_select$selection["SC(n)"], 2)  # Min 2 lags
    cat(sprintf("  Optimal lag (BIC): %d\n", optimal_p))

    # Estimate VAR
    var_model <- vars::VAR(svar_ts, p = optimal_p, type = "const")

    #---------------------------------------------------------------------------
    # Step 2: Structural identification (Cholesky)
    #---------------------------------------------------------------------------
    # Ordering: FCI → Credit → Output
    # This assumes:
    #   - FCI is most exogenous (financial conditions)
    #   - Credit responds to FCI contemporaneously
    #   - Output responds to both with a lag (most endogenous)

    # Compute structural IRFs using Cholesky decomposition
    irf_horizon <- 24

    # IRF: FCI shock → Credit
    irf_fci_credit <- vars::irf(var_model, impulse = "FCI_COMP",
                                 response = "Credit_yoy",
                                 n.ahead = irf_horizon, ortho = TRUE,
                                 boot = TRUE, runs = 500, ci = 0.90)

    # IRF: FCI shock → Output
    irf_fci_output <- vars::irf(var_model, impulse = "FCI_COMP",
                                 response = dep_var,
                                 n.ahead = irf_horizon, ortho = TRUE,
                                 boot = TRUE, runs = 500, ci = 0.90)

    # IRF: Credit shock → Output
    irf_credit_output <- vars::irf(var_model, impulse = "Credit_yoy",
                                    response = dep_var,
                                    n.ahead = irf_horizon, ortho = TRUE,
                                    boot = TRUE, runs = 500, ci = 0.90)

    #---------------------------------------------------------------------------
    # Step 3: Extract and store IRF results
    #---------------------------------------------------------------------------

    # FCI → Credit IRF
    svar_fci_credit <- data.frame(
      horizon = 0:irf_horizon,
      impulse = "FCI_COMP",
      response = "Credit_yoy",
      irf = irf_fci_credit$irf$FCI_COMP[, "Credit_yoy"],
      lower = irf_fci_credit$Lower$FCI_COMP[, "Credit_yoy"],
      upper = irf_fci_credit$Upper$FCI_COMP[, "Credit_yoy"],
      output_var = output_var,
      method = "SVAR"
    )

    # FCI → Output IRF
    svar_fci_output <- data.frame(
      horizon = 0:irf_horizon,
      impulse = "FCI_COMP",
      response = dep_var,
      irf = irf_fci_output$irf$FCI_COMP[, dep_var],
      lower = irf_fci_output$Lower$FCI_COMP[, dep_var],
      upper = irf_fci_output$Upper$FCI_COMP[, dep_var],
      output_var = output_var,
      method = "SVAR"
    )

    # Credit → Output IRF
    svar_credit_output <- data.frame(
      horizon = 0:irf_horizon,
      impulse = "Credit_yoy",
      response = dep_var,
      irf = irf_credit_output$irf$Credit_yoy[, dep_var],
      lower = irf_credit_output$Lower$Credit_yoy[, dep_var],
      upper = irf_credit_output$Upper$Credit_yoy[, dep_var],
      output_var = output_var,
      method = "SVAR"
    )

    # Store results
    svar_results[[output_var]] <- list(
      var_model = var_model,
      fci_credit = svar_fci_credit,
      fci_output = svar_fci_output,
      credit_output = svar_credit_output
    )

    #---------------------------------------------------------------------------
    # Step 4: Report SVAR results
    #---------------------------------------------------------------------------

    cat("\n  SVAR Impulse Response Functions (Cholesky identification):\n\n")

    # FCI → Credit
    cat("  FCI → Credit (1 SD shock):\n")
    for (h in c(0, 3, 6, 12, 18, 24)) {
      if (h <= irf_horizon) {
        row <- svar_fci_credit[svar_fci_credit$horizon == h, ]
        sig <- ifelse(row$lower > 0 | row$upper < 0, "*", "")
        cat(sprintf("    h=%2d: %+6.2f [%+6.2f, %+6.2f] %s\n",
                    h, row$irf, row$lower, row$upper, sig))
      }
    }

    # Credit → Output
    cat(sprintf("\n  Credit → %s (1 SD shock):\n", output_var))
    for (h in c(0, 3, 6, 12, 18, 24)) {
      if (h <= irf_horizon) {
        row <- svar_credit_output[svar_credit_output$horizon == h, ]
        sig <- ifelse(row$lower > 0 | row$upper < 0, "*", "")
        cat(sprintf("    h=%2d: %+6.3f [%+6.3f, %+6.3f] %s\n",
                    h, row$irf, row$lower, row$upper, sig))
      }
    }

    # FCI → Output (total effect)
    cat(sprintf("\n  FCI → %s (1 SD shock, total effect):\n", output_var))
    for (h in c(0, 3, 6, 12, 18, 24)) {
      if (h <= irf_horizon) {
        row <- svar_fci_output[svar_fci_output$horizon == h, ]
        sig <- ifelse(row$lower > 0 | row$upper < 0, "*", "")
        cat(sprintf("    h=%2d: %+6.3f [%+6.3f, %+6.3f] %s\n",
                    h, row$irf, row$lower, row$upper, sig))
      }
    }

    #---------------------------------------------------------------------------
    # Step 5: Compare with LP estimates
    #---------------------------------------------------------------------------

    cat("\n  LP vs SVAR Comparison (FCI → Output):\n\n")
    cat(sprintf("  %5s %12s %12s %12s\n", "h", "LP", "SVAR", "Difference"))
    cat(paste(rep("-", 50), collapse = ""), "\n")

    for (h in c(3, 6, 12, 18, 24)) {
      # Get LP estimate
      lp_row <- all_output_lp %>%
        filter(variable == output_var, horizon == h)

      # Get SVAR estimate
      svar_row <- svar_fci_output[svar_fci_output$horizon == h, ]

      if (nrow(lp_row) > 0 && nrow(svar_row) > 0) {
        diff <- lp_row$coef - svar_row$irf
        cat(sprintf("  %5d %+12.3f %+12.3f %+12.3f\n",
                    h, lp_row$coef, svar_row$irf, diff))

        # Store comparison
        lp_vs_svar_comparison <- rbind(lp_vs_svar_comparison, data.frame(
          output_var = output_var,
          horizon = h,
          lp_coef = lp_row$coef,
          lp_se = lp_row$se,
          svar_irf = svar_row$irf,
          svar_lower = svar_row$lower,
          svar_upper = svar_row$upper,
          difference = diff
        ))
      }
    }
    cat("\n")

  }, error = function(e) {
    cat(sprintf("  Error in SVAR estimation: %s\n\n", e$message))
  })
}

#-----------------------------------------------------------------------------
# 6C.2: Forecast Error Variance Decomposition
#-----------------------------------------------------------------------------
cat("--------------------------------------------------------------------------------\n")
cat("FORECAST ERROR VARIANCE DECOMPOSITION (FEVD)\n")
cat("--------------------------------------------------------------------------------\n\n")

cat("How much of output variance is explained by FCI and Credit shocks?\n\n")

fevd_results <- data.frame()

for (output_var in names(svar_results)) {
  if (is.null(svar_results[[output_var]]$var_model)) next

  dep_var <- paste0(output_var, "_yoy")
  var_model <- svar_results[[output_var]]$var_model

  tryCatch({
    # Compute FEVD
    fevd_obj <- vars::fevd(var_model, n.ahead = 24)

    # Extract FEVD for output variable
    fevd_output <- as.data.frame(fevd_obj[[dep_var]])
    fevd_output$horizon <- 0:(nrow(fevd_output) - 1)
    fevd_output$output_var <- output_var

    cat(sprintf("%s - Variance Decomposition:\n", OUTPUT_LABELS[output_var]))
    cat(sprintf("  %5s %12s %12s %12s\n", "h", "FCI", "Credit", "Own"))
    cat(paste(rep("-", 50), collapse = ""), "\n")

    for (h in c(1, 6, 12, 24)) {
      if (h <= nrow(fevd_output)) {
        row <- fevd_output[h, ]
        cat(sprintf("  %5d %11.1f%% %11.1f%% %11.1f%%\n",
                    h - 1,
                    row$FCI_COMP * 100,
                    row$Credit_yoy * 100,
                    row[[dep_var]] * 100))
      }
    }
    cat("\n")

    # Store FEVD results
    fevd_results <- rbind(fevd_results, fevd_output)

  }, error = function(e) {
    cat(sprintf("  Error computing FEVD for %s: %s\n\n", output_var, e$message))
  })
}

#-----------------------------------------------------------------------------
# 6C.3: Historical Decomposition
#-----------------------------------------------------------------------------
cat("--------------------------------------------------------------------------------\n")
cat("CUMULATIVE EFFECTS: Credit Channel Contribution\n")
cat("--------------------------------------------------------------------------------\n\n")

for (output_var in names(svar_results)) {
  if (is.null(svar_results[[output_var]])) next

  # Calculate cumulative IRF (total effect over time)
  credit_output_irf <- svar_results[[output_var]]$credit_output

  # Sum of IRF coefficients (cumulative multiplier)
  cum_effect_12 <- sum(credit_output_irf$irf[1:13])  # h=0 to h=12
  cum_effect_24 <- sum(credit_output_irf$irf[1:25])  # h=0 to h=24

  cat(sprintf("%s:\n", OUTPUT_LABELS[output_var]))
  cat(sprintf("  Cumulative Credit → Output effect (h=0-12): %+.3f\n", cum_effect_12))
  cat(sprintf("  Cumulative Credit → Output effect (h=0-24): %+.3f\n\n", cum_effect_24))
}

#-----------------------------------------------------------------------------
# 6C.4: Summary and Interpretation
#-----------------------------------------------------------------------------
cat("--------------------------------------------------------------------------------\n")
cat("LP vs SVAR: SUMMARY AND INTERPRETATION\n")
cat("--------------------------------------------------------------------------------\n\n")

if (nrow(lp_vs_svar_comparison) > 0) {
  # Calculate average absolute difference
  avg_diff <- mean(abs(lp_vs_svar_comparison$difference), na.rm = TRUE)

  # Calculate correlation
  cor_lp_svar <- cor(lp_vs_svar_comparison$lp_coef,
                     lp_vs_svar_comparison$svar_irf,
                     use = "complete.obs")

  cat(sprintf("Correlation between LP and SVAR estimates: %.3f\n", cor_lp_svar))
  cat(sprintf("Average absolute difference: %.3f pp\n\n", avg_diff))

  if (cor_lp_svar > 0.7) {
    cat(">>> STRONG AGREEMENT: LP and SVAR give similar impulse responses.\n")
    cat("    The credit-output transmission is robust to methodology.\n\n")
  } else if (cor_lp_svar > 0.4) {
    cat(">>> MODERATE AGREEMENT: LP and SVAR give broadly similar results.\n")
    cat("    Some differences may reflect non-linearities (captured by LP) or\n")
    cat("    VAR dynamics (captured by SVAR).\n\n")
  } else {
    cat(">>> WEAK AGREEMENT: LP and SVAR give different results.\n")
    cat("    This may indicate model misspecification or non-linear dynamics.\n")
    cat("    LP estimates may be more robust in this case.\n\n")
  }

  # Key insight about transmission
  cat("KEY INSIGHTS FROM SVAR:\n")

  if (nrow(fevd_results) > 0) {
    for (output_var in unique(fevd_results$output_var)) {
      fevd_24 <- fevd_results %>%
        filter(output_var == !!output_var, horizon == 23)

      if (nrow(fevd_24) > 0) {
        fci_share <- fevd_24$FCI_COMP * 100
        credit_share <- fevd_24$Credit_yoy * 100
        total_financial <- fci_share + credit_share

        cat(sprintf("\n  %s (at h=24):\n", OUTPUT_LABELS[output_var]))
        cat(sprintf("    FCI explains:    %5.1f%% of variance\n", fci_share))
        cat(sprintf("    Credit explains: %5.1f%% of variance\n", credit_share))
        cat(sprintf("    TOTAL FINANCIAL: %5.1f%% of output variance\n", total_financial))
      }
    }
  }
}

# Save comparison results
if (nrow(lp_vs_svar_comparison) > 0) {
  write.csv(lp_vs_svar_comparison,
            file.path(PUZZLE_CONFIG$output_dir, "LP_vs_SVAR_Comparison.csv"),
            row.names = FALSE)
  cat("\nSaved: LP_vs_SVAR_Comparison.csv\n")
}

if (nrow(fevd_results) > 0) {
  write.csv(fevd_results,
            file.path(PUZZLE_CONFIG$output_dir, "SVAR_FEVD_Results.csv"),
            row.names = FALSE)
  cat("Saved: SVAR_FEVD_Results.csv\n")
}

# Save SVAR IRF results
svar_irf_all <- data.frame()
for (output_var in names(svar_results)) {
  if (!is.null(svar_results[[output_var]])) {
    svar_irf_all <- rbind(svar_irf_all,
                          svar_results[[output_var]]$fci_credit,
                          svar_results[[output_var]]$fci_output,
                          svar_results[[output_var]]$credit_output)
  }
}
if (nrow(svar_irf_all) > 0) {
  write.csv(svar_irf_all,
            file.path(PUZZLE_CONFIG$output_dir, "SVAR_IRF_Results.csv"),
            row.names = FALSE)
  cat("Saved: SVAR_IRF_Results.csv\n")
}


################################################################################
# 7. COMPREHENSIVE COMPARISON: ALL OUTPUT MEASURES
################################################################################

cat("\n================================================================================\n")
cat("COMPREHENSIVE COMPARISON: ALL OUTPUT MEASURES\n")
cat("================================================================================\n\n")

# Compare effects at 12-month horizon
comparison_12m <- all_output_lp %>%
  filter(horizon == 12) %>%
  arrange(p_value)

cat("Effect of 1 SD FCI tightening at 12-month horizon:\n\n")
cat(sprintf("%-35s %10s %10s %10s\n", "Variable", "Coef", "p-value", "Signif"))
cat(paste(rep("-", 70), collapse = ""), "\n")

for (i in 1:nrow(comparison_12m)) {
  r <- comparison_12m[i, ]
  stars <- case_when(
    r$p_value < 0.01 ~ "***",
    r$p_value < 0.05 ~ "**",
    r$p_value < 0.10 ~ "*",
    TRUE ~ ""
  )
  cat(sprintf("%-35s %+10.3f %10.3f %10s\n",
              r$variable_label, r$coef, r$p_value, stars))
}


################################################################################
# 8. VISUALIZATIONS
################################################################################

cat("\n================================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

if (!dir.exists(PUZZLE_CONFIG$output_dir)) {
  dir.create(PUZZLE_CONFIG$output_dir, recursive = TRUE)
}

# Color palette
colors_output <- c(
  "Aggregate GDP (IMAE-P)" = "#2C3E50",
  "Non-Agro GDP (IMAE-P ex Agro/Binational)" = "#E74C3C",
  "Investment (FBKf)" = "#3498DB",
  "Consumption" = "#27AE60"
)

# Plot 1: Comparison of LP IRFs for all output measures
p1 <- all_output_lp %>%
  ggplot(aes(x = horizon, y = coef, color = variable_label, fill = variable_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = colors_output) +
  scale_fill_manual(values = colors_output) +
  theme_minimal(base_size = 12) +
  labs(title = "Output Puzzle Investigation: FCI Effect on Different Output Measures",
       subtitle = "Response to 1 SD FCI tightening | 90% CI | Testing agricultural contamination hypothesis",
       x = "Months", y = "Effect on YoY growth (pp)",
       color = "Output Measure", fill = "Output Measure") +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2))

ggsave(file.path(PUZZLE_CONFIG$output_dir, "120_Output_Puzzle_LP_Comparison.png"), p1,
       width = 14, height = 7, dpi = 300)
cat("Saved: 120_Output_Puzzle_LP_Comparison.png\n")

# Plot 2: Faceted view by output measure
p2 <- all_output_lp %>%
  ggplot(aes(x = horizon)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#377EB8", alpha = 0.2) +
  geom_line(aes(y = coef), color = "#377EB8", linewidth = 1) +
  geom_point(aes(y = coef), color = "#377EB8", size = 2) +
  facet_wrap(~variable_label, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 11) +
  labs(title = "FCI Effect on Output: By Measure",
       subtitle = "Response to 1 SD FCI shock (tightening) | 90% CI",
       x = "Months", y = "Effect (pp)") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(PUZZLE_CONFIG$output_dir, "121_Output_Puzzle_LP_Faceted.png"), p2,
       width = 12, height = 8, dpi = 300)
cat("Saved: 121_Output_Puzzle_LP_Faceted.png\n")

# Plot 3: Bar chart comparison at key horizons
bar_data <- all_output_lp %>%
  filter(horizon %in% c(6, 12, 18)) %>%
  mutate(horizon_label = paste0("h=", horizon, "m"),
         signif = ifelse(p_value < 0.10, "Significant", "Not significant"))

p3 <- ggplot(bar_data, aes(x = variable_label, y = coef, fill = signif)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~horizon_label, ncol = 3) +
  scale_fill_manual(values = c("Significant" = "#E74C3C", "Not significant" = "#BDC3C7")) +
  theme_minimal(base_size = 11) +
  labs(title = "FCI Effect on Output Measures at Key Horizons",
       subtitle = "Error bars show 90% CI",
       x = NULL, y = "Effect (pp)",
       fill = NULL) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(PUZZLE_CONFIG$output_dir, "122_Output_Puzzle_BarChart.png"), p3,
       width = 14, height = 7, dpi = 300)
cat("Saved: 122_Output_Puzzle_BarChart.png\n")

# Plot 4: IMAEP vs IMAEP_SANB direct comparison (if both available)
if ("IMAEP" %in% available_vars && "IMAEP_SANB" %in% available_vars) {

  comparison_data <- all_output_lp %>%
    filter(variable %in% c("IMAEP", "IMAEP_SANB"))

  p4 <- ggplot(comparison_data, aes(x = horizon, y = coef,
                                     color = variable_label, fill = variable_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("Aggregate GDP (IMAE-P)" = "#2C3E50",
                                   "Non-Agro GDP (IMAE-P ex Agro/Binational)" = "#E74C3C")) +
    scale_fill_manual(values = c("Aggregate GDP (IMAE-P)" = "#2C3E50",
                                  "Non-Agro GDP (IMAE-P ex Agro/Binational)" = "#E74C3C")) +
    annotate("rect", xmin = 10, xmax = 14, ymin = -Inf, ymax = Inf,
             alpha = 0.1, fill = "gray50") +
    theme_minimal(base_size = 12) +
    labs(title = "Agricultural Contamination Hypothesis",
         subtitle = "If FCI affects non-agro GDP but not aggregate, agriculture masks the credit channel",
         x = "Months", y = "Effect on YoY growth (pp)",
         color = NULL, fill = NULL) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")

  ggsave(file.path(PUZZLE_CONFIG$output_dir, "123_IMAEP_vs_SANB_Comparison.png"), p4,
         width = 12, height = 7, dpi = 300)
  cat("Saved: 123_IMAEP_vs_SANB_Comparison.png\n")
}

# Plot 5: LP vs SVAR Comparison
if (exists("lp_vs_svar_comparison") && nrow(lp_vs_svar_comparison) > 0 &&
    exists("svar_results") && length(svar_results) > 0) {

  # Prepare LP data for plotting
  lp_plot_data <- all_output_lp %>%
    filter(variable %in% c("IMAEP", "IMAEP_SANB")) %>%
    mutate(method = "Local Projection",
           irf = coef,
           lower = ci_lower,
           upper = ci_upper) %>%
    dplyr::select(output_var = variable, horizon, irf, lower, upper, method)

  # Prepare SVAR data for plotting
  svar_plot_data <- data.frame()
  for (output_var in names(svar_results)) {
    if (!is.null(svar_results[[output_var]]$fci_output)) {
      temp <- svar_results[[output_var]]$fci_output %>%
        mutate(method = "SVAR (Cholesky)") %>%
        dplyr::select(output_var, horizon, irf, lower, upper, method)
      svar_plot_data <- rbind(svar_plot_data, temp)
    }
  }

  # Combine data
  comparison_plot_data <- rbind(lp_plot_data, svar_plot_data) %>%
    mutate(output_label = ifelse(output_var == "IMAEP",
                                 "Aggregate GDP (IMAE-P)",
                                 "Non-Agro GDP (IMAE-P ex Agro/Binational)"))

  # Create comparison plot
  p5 <- ggplot(comparison_plot_data,
               aes(x = horizon, y = irf, color = method, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~output_label, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("Local Projection" = "#E74C3C",
                                   "SVAR (Cholesky)" = "#3498DB")) +
    scale_fill_manual(values = c("Local Projection" = "#E74C3C",
                                  "SVAR (Cholesky)" = "#3498DB")) +
    theme_minimal(base_size = 12) +
    labs(title = "LP vs SVAR Comparison: FCI Effect on Output",
         subtitle = "Both methods estimate impulse response to 1 SD FCI shock | 90% CI",
         x = "Months", y = "Effect on YoY growth (pp)",
         color = "Method", fill = "Method") +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          strip.text = element_text(face = "bold"))

  ggsave(file.path(PUZZLE_CONFIG$output_dir, "124_LP_vs_SVAR_Comparison.png"), p5,
         width = 12, height = 10, dpi = 300)
  cat("Saved: 124_LP_vs_SVAR_Comparison.png\n")

  # Plot 6: Credit → Output IRF from SVAR
  credit_output_plot_data <- data.frame()
  for (output_var in names(svar_results)) {
    if (!is.null(svar_results[[output_var]]$credit_output)) {
      temp <- svar_results[[output_var]]$credit_output %>%
        mutate(output_label = ifelse(output_var == "IMAEP",
                                     "Aggregate GDP (IMAE-P)",
                                     "Non-Agro GDP (IMAE-P ex Agro/Binational)"))
      credit_output_plot_data <- rbind(credit_output_plot_data, temp)
    }
  }

  if (nrow(credit_output_plot_data) > 0) {
    p6 <- ggplot(credit_output_plot_data,
                 aes(x = horizon, y = irf)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#27AE60", alpha = 0.2) +
      geom_line(color = "#27AE60", linewidth = 1) +
      geom_point(color = "#27AE60", size = 2) +
      facet_wrap(~output_label, scales = "free_y", ncol = 1) +
      theme_minimal(base_size = 12) +
      labs(title = "SVAR: Credit Shock Effect on Output",
           subtitle = "Response to 1 SD credit growth shock | Cholesky identification | 90% CI",
           x = "Months", y = "Effect on YoY output growth (pp)") +
      theme(plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold"))

    ggsave(file.path(PUZZLE_CONFIG$output_dir, "125_SVAR_Credit_Output.png"), p6,
           width = 12, height = 8, dpi = 300)
    cat("Saved: 125_SVAR_Credit_Output.png\n")
  }

  # Plot 7: FCI → Credit IRF from SVAR
  fci_credit_plot_data <- data.frame()
  for (output_var in names(svar_results)) {
    if (!is.null(svar_results[[output_var]]$fci_credit)) {
      temp <- svar_results[[output_var]]$fci_credit %>%
        mutate(output_label = ifelse(output_var == "IMAEP",
                                     "Using IMAEP model",
                                     "Using IMAEP_SANB model"))
      fci_credit_plot_data <- rbind(fci_credit_plot_data, temp)
    }
  }

  if (nrow(fci_credit_plot_data) > 0) {
    # Use only one model (they should be very similar)
    fci_credit_single <- fci_credit_plot_data %>%
      filter(output_label == "Using IMAEP_SANB model")

    p7 <- ggplot(fci_credit_single, aes(x = horizon, y = irf)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#9B59B6", alpha = 0.2) +
      geom_line(color = "#9B59B6", linewidth = 1.2) +
      geom_point(color = "#9B59B6", size = 2.5) +
      theme_minimal(base_size = 12) +
      labs(title = "SVAR: FCI Shock Effect on Credit Growth",
           subtitle = "Response of credit growth (YoY) to 1 SD FCI shock | Cholesky | 90% CI",
           x = "Months", y = "Effect on credit growth (pp)") +
      theme(plot.title = element_text(face = "bold"))

    ggsave(file.path(PUZZLE_CONFIG$output_dir, "126_SVAR_FCI_Credit.png"), p7,
           width = 10, height = 6, dpi = 300)
    cat("Saved: 126_SVAR_FCI_Credit.png\n")
  }
}


################################################################################
# 8B. QUARTERLY TRANSMISSION CHANNEL ANALYSIS
################################################################################
#
# Consumo and FBKf are no longer in the Datos_macro sheet (monthly).
# They ARE available at quarterly frequency in the Quarterly_SA sheet.
# This section loads quarterly national accounts, aggregates the monthly FCI
# to quarterly, and runs LP to test the transmission channel hypothesis
# at quarterly frequency (max 8 quarters = 24 months equivalent).
#

cat("\n================================================================================\n")
cat("QUARTERLY TRANSMISSION CHANNEL ANALYSIS\n")
cat("================================================================================\n\n")

cat("RATIONALE:\n")
cat("  Consumo (private consumption) and FBKf (gross fixed capital formation)\n")
cat("  are available only at quarterly frequency (Quarterly_SA sheet).\n")
cat("  We aggregate the monthly FCI to quarterly and run LP at quarterly horizons\n")
cat("  to test whether FCI transmits to output via investment or consumption.\n\n")

tryCatch({

  # ── 1. Load Quarterly_SA sheet ──────────────────────────────────────────────
  cat("Loading Quarterly_SA sheet...\n")
  quarterly_raw <- read_excel(PUZZLE_CONFIG$data_file, sheet = "Quarterly_SA")

  # Rename columns
  quarterly_raw <- quarterly_raw %>%
    rename(
      fecha            = Date,
      Consumo          = `Consumo Privado`,
      FBKf             = `Formación bruta de capital fijo`
    ) %>%
    mutate(fecha = as.Date(fecha)) %>%
    arrange(fecha)

  cat(sprintf("  Quarterly observations: %d (%s to %s)\n",
              nrow(quarterly_raw),
              format(min(quarterly_raw$fecha)),
              format(max(quarterly_raw$fecha))))

  # ── 2. Compute YoY growth rates ────────────────────────────────────────────
  quarterly_na <- quarterly_raw %>%
    mutate(
      Consumo_yoy = (Consumo / lag(Consumo, 4) - 1) * 100,
      FBKf_yoy    = (FBKf    / lag(FBKf, 4)    - 1) * 100
    ) %>%
    dplyr::select(fecha, Consumo_yoy, FBKf_yoy)

  # ── 3. Aggregate monthly FCI to quarterly (end-of-quarter months) ──────────
  cat("Aggregating monthly FCI to quarterly (end-of-quarter months)...\n")

  fci_quarterly <- fci_data %>%
    mutate(month = as.integer(format(fecha, "%m"))) %>%
    filter(month %in% c(3, 6, 9, 12)) %>%
    dplyr::select(fecha, FCI_COMP, FCI_exCredit)

  cat(sprintf("  Monthly FCI observations retained: %d\n", nrow(fci_quarterly)))

  # ── 4. Aggregate IPC_yoy control to quarterly ──────────────────────────────
  has_ipc_control <- FALSE
  if (!is.null(ipc_credit_data) && "IPC_yoy" %in% names(ipc_credit_data)) {
    ipc_quarterly <- ipc_credit_data %>%
      mutate(month = as.integer(format(fecha, "%m"))) %>%
      filter(month %in% c(3, 6, 9, 12)) %>%
      dplyr::select(fecha, IPC_yoy)
    has_ipc_control <- TRUE
    cat("  IPC_yoy control aggregated to quarterly.\n")
  }

  # ── 5. Merge quarterly datasets ────────────────────────────────────────────
  quarterly_data <- fci_quarterly %>%
    inner_join(quarterly_na, by = "fecha")

  if (has_ipc_control) {
    quarterly_data <- quarterly_data %>%
      left_join(ipc_quarterly, by = "fecha")
  }

  quarterly_data <- quarterly_data %>%
    arrange(fecha) %>%
    na.omit()

  cat(sprintf("  Merged quarterly dataset: %d observations (%s to %s)\n\n",
              nrow(quarterly_data),
              format(min(quarterly_data$fecha)),
              format(max(quarterly_data$fecha))))

  # ── 6. Quarterly LP function (reuses run_lp_standard logic) ────────────────

  run_lp_quarterly <- function(data, y_var, fci_var, max_h, n_lags = 2,
                               conf_level = 0.90, control_vars = NULL) {
    results <- data.frame()
    z_crit <- qnorm(1 - (1 - conf_level) / 2)

    for (h in 1:max_h) {
      data_h <- data %>%
        mutate(
          y_fwd    = lead(!!sym(y_var), h),
          y_lag1   = lag(!!sym(y_var), 1),
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

      formula_str <- paste("y_fwd ~", fci_var, "+",
                           paste(c(lag_vars, control_cols), collapse = " + "))

      reg_data <- data_h %>%
        dplyr::select(y_fwd, !!sym(fci_var), all_of(lag_vars), all_of(control_cols)) %>%
        na.omit()

      if (nrow(reg_data) < 20) next

      model <- lm(as.formula(formula_str), data = reg_data)
      vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
      coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

      idx <- which(rownames(coef_test) == fci_var)

      results <- rbind(results, data.frame(
        horizon_q = h,
        horizon_m = h * 3,
        coef      = coef_test[idx, 1],
        se        = coef_test[idx, 2],
        t_stat    = coef_test[idx, 3],
        p_value   = coef_test[idx, 4],
        ci_lower  = coef_test[idx, 1] - z_crit * coef_test[idx, 2],
        ci_upper  = coef_test[idx, 1] + z_crit * coef_test[idx, 2],
        n_obs     = nrow(reg_data),
        r_squared = summary(model)$r.squared
      ))
    }
    return(results)
  }

  # ── 7. Run quarterly LP for Consumo and FBKf ──────────────────────────────
  max_h_q <- 8   # 8 quarters = 24 months

  quarterly_lp_results <- data.frame()

  for (var in c("Consumo", "FBKf")) {
    y_var <- paste0(var, "_yoy")
    cat(sprintf("Running quarterly LP: FCI_COMP -> %s ... ", OUTPUT_LABELS[var]))

    q_ctrl <- if (has_ipc_control) "IPC_yoy" else NULL

    lp_res <- run_lp_quarterly(
      quarterly_data, y_var, "FCI_COMP",
      max_h = max_h_q, n_lags = 2,
      conf_level = 0.90,
      control_vars = q_ctrl
    ) %>%
      mutate(
        variable       = var,
        variable_label = OUTPUT_LABELS[var],
        frequency      = "quarterly"
      )

    quarterly_lp_results <- rbind(quarterly_lp_results, lp_res)
    cat("done\n")
  }

  # ── 8. Print results ──────────────────────────────────────────────────────
  cat("\n================================================================================\n")
  cat("QUARTERLY LP RESULTS: FCI EFFECT ON TRANSMISSION CHANNELS\n")
  cat("================================================================================\n\n")

  cat("Response (pp) to 1 SD FCI tightening | *** p<0.01, ** p<0.05, * p<0.10\n\n")

  for (var in c("Consumo", "FBKf")) {
    var_res <- quarterly_lp_results %>%
      filter(variable == var)

    if (nrow(var_res) == 0) next

    cat(sprintf("%s (quarterly LP, max h = %d quarters = %d months):\n",
                OUTPUT_LABELS[var], max_h_q, max_h_q * 3))
    cat(sprintf("  %3s  %5s  %10s  %8s  %6s  %6s\n",
                "h_q", "h_m", "Coef", "SE", "p-val", ""))

    for (i in 1:nrow(var_res)) {
      r <- var_res[i, ]
      stars <- case_when(
        r$p_value < 0.01 ~ "***",
        r$p_value < 0.05 ~ "**",
        r$p_value < 0.10 ~ "*",
        TRUE ~ ""
      )
      cat(sprintf("  %3d  %5d  %+10.3f  %8.3f  %6.3f  %s\n",
                  r$horizon_q, r$horizon_m, r$coef, r$se, r$p_value, stars))
    }
    cat("\n")
  }

  # Key comparison at h=4 quarters (= 12 months)
  cat("KEY RESULTS (h = 4 quarters = 12 months equivalent):\n")
  for (var in c("Consumo", "FBKf")) {
    r4 <- quarterly_lp_results %>% filter(variable == var, horizon_q == 4)
    if (nrow(r4) > 0) {
      sig <- ifelse(r4$p_value < 0.10, "SIGNIFICANT", "not significant")
      cat(sprintf("  %s: coef = %+.3f, p = %.3f (%s)\n",
                  OUTPUT_LABELS[var], r4$coef, r4$p_value, sig))
    }
  }

  # ── 9. Export CSV ──────────────────────────────────────────────────────────
  write.csv(quarterly_lp_results,
            file.path(PUZZLE_CONFIG$output_dir, "Output_Puzzle_Quarterly_Transmission.csv"),
            row.names = FALSE)
  cat("\nSaved: Output_Puzzle_Quarterly_Transmission.csv\n")

  # ── 10. 2-panel plot: Consumo LP + FBKf LP ────────────────────────────────
  cat("Generating quarterly transmission plot...\n")

  pq <- quarterly_lp_results %>%
    ggplot(aes(x = horizon_q, y = coef)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#377EB8", alpha = 0.20) +
    geom_line(color = "#377EB8", linewidth = 1) +
    geom_point(aes(shape = ifelse(p_value < 0.10, "Significant", "Not significant")),
               color = "#377EB8", size = 3) +
    scale_shape_manual(values = c("Significant" = 16, "Not significant" = 1)) +
    facet_wrap(~variable_label, scales = "free_y", ncol = 2) +
    scale_x_continuous(breaks = 1:max_h_q,
                       labels = paste0("Q", 1:max_h_q, "\n(", (1:max_h_q)*3, "m)")) +
    theme_minimal(base_size = 12) +
    labs(
      title    = "Quarterly Transmission Channels: FCI Effect on Consumption and Investment",
      subtitle = "Response to 1 SD FCI tightening | Quarterly LP | 90% CI",
      x        = "Horizon (quarters / months)",
      y        = "Effect on YoY growth (pp)",
      shape    = NULL
    ) +
    theme(
      plot.title    = element_text(face = "bold"),
      strip.text    = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(file.path(PUZZLE_CONFIG$output_dir, "127_Output_Puzzle_Quarterly_Transmission.png"),
         pq, width = 14, height = 7, dpi = 300)
  cat("Saved: 127_Output_Puzzle_Quarterly_Transmission.png\n")

}, error = function(e) {
  cat(sprintf("WARNING: Quarterly transmission analysis failed: %s\n", e$message))
  cat("Continuing with remaining script...\n\n")
})


################################################################################
# 9. EXPORT RESULTS
################################################################################

cat("\n================================================================================\n")
cat("EXPORTING RESULTS\n")
cat("================================================================================\n\n")

# LP Results
write.csv(all_output_lp,
          file.path(PUZZLE_CONFIG$output_dir, "Output_Puzzle_LP_Results.csv"),
          row.names = FALSE)
cat("Saved: Output_Puzzle_LP_Results.csv\n")

# Summary comparison
comparison_summary <- all_output_lp %>%
  filter(horizon %in% c(6, 12, 18, 24)) %>%
  dplyr::select(variable, variable_label, horizon, coef, se, p_value,
                ci_lower, ci_upper, n_obs, r_squared) %>%
  mutate(significant = ifelse(p_value < 0.10, "*", ""))

write.csv(comparison_summary,
          file.path(PUZZLE_CONFIG$output_dir, "Output_Puzzle_Summary.csv"),
          row.names = FALSE)
cat("Saved: Output_Puzzle_Summary.csv\n")

# VECM results (if available)
if (exists("vecm_results") && !is.null(vecm_results$rank) && !is.na(vecm_results$rank)) {
  vecm_summary <- data.frame(
    test = "Johansen Cointegration",
    rank = vecm_results$rank,
    interpretation = ifelse(vecm_results$rank > 0,
                            "Long-run equilibrium exists",
                            "No cointegration")
  )
  write.csv(vecm_summary,
            file.path(PUZZLE_CONFIG$output_dir, "Output_Puzzle_VECM_Summary.csv"),
            row.names = FALSE)
  cat("Saved: Output_Puzzle_VECM_Summary.csv\n")
}


################################################################################
# 10. FINAL SUMMARY AND CONCLUSIONS
################################################################################

cat("\n################################################################################\n")
cat("OUTPUT PUZZLE INVESTIGATION - CONCLUSIONS\n")
cat("################################################################################\n\n")

# Hypothesis 1 conclusion
cat("HYPOTHESIS 1: AGRICULTURAL CONTAMINATION\n")
if ("IMAEP_SANB" %in% available_vars) {
  sanb_sig <- any(all_output_lp$variable == "IMAEP_SANB" &
                    all_output_lp$horizon %in% 6:18 &
                    all_output_lp$p_value < 0.10)
  imaep_sig <- any(all_output_lp$variable == "IMAEP" &
                     all_output_lp$horizon %in% 6:18 &
                     all_output_lp$p_value < 0.10)

  if (sanb_sig && !imaep_sig) {
    cat("  CONFIRMED: FCI affects non-agricultural GDP but not aggregate.\n")
    cat("  The agricultural sector masks the credit channel in Paraguay.\n\n")
  } else if (sanb_sig && imaep_sig) {
    cat("  PARTIAL: FCI affects both, but non-agro effect is likely stronger.\n\n")
  } else if (!sanb_sig && !imaep_sig) {
    cat("  NOT CONFIRMED: FCI does not significantly affect either measure.\n")
    cat("  The puzzle may have other explanations.\n\n")
  }
} else {
  cat("  UNTESTED: IMAEP_SANB variable not available.\n")
  cat("  Please add IMAEP_SANB (non-agricultural GDP) to the data.\n\n")
}

# Hypothesis 2 conclusion
cat("HYPOTHESIS 2: LONG-RUN COINTEGRATION\n")
if (exists("vecm_results") && !is.null(vecm_results$rank) && !is.na(vecm_results$rank)) {
  if (!is.null(vecm_results$full_rank) && vecm_results$full_rank) {
    cat("  NOT APPLICABLE: Full rank indicates variables are stationary.\n")
    cat("  No cointegrating relationships exist - variables move independently.\n")
    cat("  Standard VAR/LP approach is appropriate.\n\n")
  } else if (vecm_results$rank > 0) {
    cat(sprintf("  CONFIRMED: %d cointegrating relationship(s) detected.\n", vecm_results$rank))
    cat("  Credit, output, and FCI share a long-run equilibrium.\n")
    cat("  Short-run LP may underestimate the true credit-output link.\n\n")
  } else {
    cat("  NOT CONFIRMED: No cointegration detected.\n")
    cat("  Standard VAR/LP approach is appropriate.\n\n")
  }
} else {
  cat("  INCONCLUSIVE: VECM analysis could not be completed.\n\n")
}

# Hypothesis 3 conclusion
cat("HYPOTHESIS 3: TRANSMISSION CHANNELS\n")
channel_vars <- intersect(c("FBKf", "Consumo"), available_vars)
if (length(channel_vars) > 0) {
  sig_channels <- all_output_lp %>%
    filter(variable %in% channel_vars, horizon == 12, p_value < 0.10)

  if (nrow(sig_channels) > 0) {
    cat("  CONFIRMED (monthly): FCI affects output through:\n")
    for (v in sig_channels$variable) {
      cat(sprintf("    - %s\n", OUTPUT_LABELS[v]))
    }
    cat("\n")
  } else {
    cat("  NOT CONFIRMED (monthly): No significant effect on investment or consumption.\n\n")
  }
} else {
  cat("  Monthly data not available in Datos_macro.\n")
}

# Check quarterly transmission results
if (exists("quarterly_lp_results") && nrow(quarterly_lp_results) > 0) {
  q_sig <- quarterly_lp_results %>%
    filter(horizon_q == 4, p_value < 0.10)
  if (nrow(q_sig) > 0) {
    cat("  CONFIRMED (quarterly): Significant FCI effect at h=4Q through:\n")
    for (i in 1:nrow(q_sig)) {
      cat(sprintf("    - %s (coef=%+.3f, p=%.3f)\n",
                  q_sig$variable_label[i], q_sig$coef[i], q_sig$p_value[i]))
    }
    cat("\n")
  } else {
    cat("  NOT CONFIRMED (quarterly): No significant effect at h=4Q.\n\n")
  }
} else {
  cat("  Quarterly analysis: not available.\n\n")
}

# Hypothesis 4 conclusion (Credit Channel Mediation)
cat("HYPOTHESIS 4: CREDIT CHANNEL MEDIATION\n")
if (exists("credit_output_lp") && nrow(credit_output_lp) > 0) {
  # Check if credit Granger-causes any output measure
  credit_causes_any <- FALSE
  if (exists("granger_results") && nrow(granger_results) > 0) {
    credit_causes_any <- any(granger_results$p_credit_to_output < 0.10)
  }

  # Check LP significance
  credit_lp_sig <- credit_output_lp %>%
    filter(horizon %in% 6:18, p_value < 0.10)

  if (credit_causes_any || nrow(credit_lp_sig) > 0) {
    cat("  CONFIRMED: Credit growth predicts output growth.\n")
    if (exists("granger_results") && nrow(granger_results) > 0) {
      sig_granger <- granger_results %>% filter(p_credit_to_output < 0.10)
      if (nrow(sig_granger) > 0) {
        cat("  Granger causality (Credit → Output):\n")
        for (i in 1:nrow(sig_granger)) {
          cat(sprintf("    - %s (p=%.3f)\n",
                      sig_granger$output_label[i], sig_granger$p_credit_to_output[i]))
        }
      }
    }
    cat("\n  TRANSMISSION CHAIN: FCI → Credit → Output is supported.\n\n")
  } else {
    cat("  NOT CONFIRMED: Credit growth does not significantly predict output.\n")
    cat("  The credit channel may be weak in Paraguay, or transmission is non-linear.\n\n")
  }

  # Report mediation results
  if (exists("mediation_results") && nrow(mediation_results) > 0) {
    avg_med <- mean(mediation_results$pct_reduction, na.rm = TRUE)
    if (!is.na(avg_med) && avg_med > 20) {
      cat(sprintf("  MEDIATION: Credit mediates %.0f%% of FCI → Output effect.\n\n", avg_med))
    }
  }
} else {
  cat("  UNTESTED: Credit channel analysis could not be completed.\n\n")
}

# Overall conclusion
cat("================================================================================\n")
cat("OVERALL ASSESSMENT\n")
cat("================================================================================\n\n")

n_tested <- length(available_vars)
n_significant <- sum(all_output_lp %>%
                       filter(horizon == 12, p_value < 0.10) %>%
                       pull(variable) %in% available_vars)

cat(sprintf("  Variables tested: %d\n", n_tested))
cat(sprintf("  Significant effects at h=12: %d\n\n", n_significant))

if (n_tested < 4) {
  cat("  NOTE: Monthly analysis covers", n_tested, "of 4 output variables.\n")
  cat("    Missing from Datos_macro (monthly):", paste(missing_vars, collapse = ", "), "\n")
  if (exists("quarterly_lp_results") && nrow(quarterly_lp_results) > 0) {
    cat("    However, Consumo and FBKf were tested via quarterly LP (section 8B).\n")
  }
}

cat("\n################################################################################\n")
cat("OUTPUT PUZZLE INVESTIGATION COMPLETE\n")
cat("################################################################################\n\n")
cat("Plots: 120-127\n")
cat("Output:", PUZZLE_CONFIG$output_dir, "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
