################################################################################
# FCI REGIME-SWITCHING AND TIME-VARYING PARAMETER ANALYSIS
################################################################################
#
# Project:      Financial Conditions Index - Regime and TVP Analysis
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Comprehensive analysis of regime-dependent FCI behavior and
#               time-varying parameter estimation including:
#
#   PART A: Regime-Split PCA Analysis
#     - Pre-IT vs Post-IT loading comparison
#     - Bootstrap tests for loading differences
#
#   PART B: Regime Dummies Storage
#     - IT regime indicator (May 2011+)
#     - COVID period indicator (Mar 2020 - Dec 2021)
#
#   PART C: Regime-Dependent Local Projections
#     - Interaction models: FCI × regime
#     - Split-sample LP comparisons
#
#   PART D: Time-Varying Parameter PCA
#     - Exponential weighting approach
#     - State-space MARSS approach (if convergent)
#
#   PART E: Threshold Effects Testing
#     - Crisis (FCI > 1 SD) vs normal times
#     - Asymmetric threshold LP
#
# Key Dates for Paraguay:
#   - Inflation Targeting (IT): May 2011
#   - COVID-19 Pandemic: March 2020 - December 2021
#
# Dependencies: Requires 01_FCI_Complete.R to be run first
#
# Output:
#   - 80-89_*.png (visualizations)
#   - FCI_Regime_*.csv (results)
#   - FCI_TVP_*.csv (time-varying loadings)
#
# Author:       Departamento de Analisis Macroeconomico
# Created:      2025-01-23
#
# References:
#   - Hamilton (1989) - Regime Switching Models
#   - Stock & Watson (2002) - Time-Varying Parameters
#   - Jordà (2005) - Local Projections
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
})

# Configuration parameters
REGIME_CONFIG <- list(
  # Key regime dates for Paraguay
  IT_START = as.Date("2011-05-01"),        # Inflation Targeting adoption
  COVID_START = as.Date("2020-03-01"),     # COVID pandemic start
  COVID_END = as.Date("2021-12-31"),       # COVID pandemic end

  # Analysis parameters
  FCI_THRESHOLD = 1.0,                     # Standard deviations for crisis
  BOOTSTRAP_REPS = 1000,                   # Bootstrap replications
  TVP_LAMBDA = 0.97,                       # Exponential decay (33-month eff. window)

  # LP parameters
  max_horizon = 24,
  n_lags = 2,
  confidence_level = 0.90,

  # Data
  data_file = "../data/FCI_data_1.xlsx",

  # Output
  output_dir = "../output",
  verbose = TRUE
)

set.seed(20250123)

cat("\n################################################################################\n")
cat("FCI REGIME-SWITCHING AND TIME-VARYING PARAMETER ANALYSIS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Key dates:\n")
cat("  Inflation Targeting: ", format(REGIME_CONFIG$IT_START, "%B %Y"), "\n")
cat("  COVID period: ", format(REGIME_CONFIG$COVID_START, "%B %Y"), " - ",
    format(REGIME_CONFIG$COVID_END, "%B %Y"), "\n\n")


################################################################################
# SECTION 2: DATA PREPARATION AND REGIME DUMMIES
################################################################################

cat("================================================================================\n")
cat("DATA PREPARATION AND REGIME DUMMIES\n")
cat("================================================================================\n\n")

# Load FCI results (run core script if needed)
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  source("01_FCI_Complete.R")
}

# Get variable definitions
VARIABLES <- resultado_fci$variables$level3
all_vars_core <- resultado_fci$variables$core_vars
all_signs_core <- c(VARIABLES$rates$signs, VARIABLES$banking$signs,
                    VARIABLES$external$signs)

# Load raw data for regime analysis
datos_raw <- read_excel(REGIME_CONFIG$data_file)
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

#' Create regime dummy variables
#' @param data Data frame with fecha column
#' @return Data frame with added regime indicators
create_regime_dummies <- function(data) {
  data %>%
    mutate(
      # IT regime: 0 before May 2011, 1 from May 2011+
      regime_IT = as.integer(fecha >= REGIME_CONFIG$IT_START),

      # COVID period: 1 during Mar 2020 - Dec 2021, 0 otherwise
      regime_COVID = as.integer(fecha >= REGIME_CONFIG$COVID_START &
                                  fecha <= REGIME_CONFIG$COVID_END),

      # IT period excluding COVID (for clean comparison)
      regime_IT_exCOVID = as.integer(fecha >= REGIME_CONFIG$IT_START &
                                       !(fecha >= REGIME_CONFIG$COVID_START &
                                           fecha <= REGIME_CONFIG$COVID_END)),

      # Period labels for analysis
      period_label = case_when(
        fecha < REGIME_CONFIG$IT_START ~ "Pre-IT",
        fecha >= REGIME_CONFIG$COVID_START & fecha <= REGIME_CONFIG$COVID_END ~ "COVID",
        TRUE ~ "Post-IT"
      )
    )
}

# Apply regime dummies to FCI data
fci_data <- resultado_fci$all_indices %>%
  create_regime_dummies()

# Apply signs to raw data for PCA analysis
apply_signs <- function(data, vars, signs) {
  data_adj <- data
  for (i in seq_along(vars)) {
    if (vars[i] %in% names(data_adj)) {
      data_adj[[vars[i]]] <- data_adj[[vars[i]]] * signs[i]
    }
  }
  return(data_adj)
}

datos_signed <- apply_signs(datos, all_vars_core, all_signs_core) %>%
  create_regime_dummies()

# Report regime splits
regime_summary <- fci_data %>%
  group_by(period_label) %>%
  summarise(
    n_obs = n(),
    start_date = min(fecha),
    end_date = max(fecha),
    .groups = "drop"
  )

cat("Regime Summary:\n")
print(regime_summary)

cat("\nDetailed counts:\n")
cat("  Pre-IT (before May 2011):    ", sum(fci_data$regime_IT == 0), " observations\n")
cat("  Post-IT (from May 2011):     ", sum(fci_data$regime_IT == 1), " observations\n")
cat("  COVID (Mar 2020 - Dec 2021): ", sum(fci_data$regime_COVID == 1), " observations\n")
cat("  IT excluding COVID:          ", sum(fci_data$regime_IT_exCOVID == 1), " observations\n\n")


################################################################################
# SECTION 3: PART A - REGIME-SPLIT PCA ANALYSIS
################################################################################

cat("================================================================================\n")
cat("PART A: REGIME-SPLIT PCA ANALYSIS\n")
cat("================================================================================\n\n")

#' Calculate PCA loadings for a given period
#' @param data Data frame with variables
#' @param vars Vector of variable names
#' @param start_date Start date filter
#' @param end_date End date filter
#' @return List with loadings, variance explained, n_obs
calc_regime_pca <- function(data, vars, start_date = NULL, end_date = NULL) {

  # Filter by date if specified
  data_filtered <- data
  if (!is.null(start_date)) {
    data_filtered <- data_filtered %>% filter(fecha >= start_date)
  }
  if (!is.null(end_date)) {
    data_filtered <- data_filtered %>% filter(fecha <= end_date)
  }

  vars_available <- intersect(vars, names(data_filtered))
  if (length(vars_available) < 2) return(NULL)

  datos_subset <- data_filtered %>%
    dplyr::select(all_of(vars_available))

  # Standardize within period
  datos_std <- scale(datos_subset)

  # Remove rows with NAs
  complete_rows <- complete.cases(datos_std)
  if (sum(complete_rows) < 20) return(NULL)

  datos_clean <- datos_std[complete_rows, , drop = FALSE]

  tryCatch({
    pca_model <- prcomp(datos_clean, center = FALSE, scale. = FALSE)
    loadings <- pca_model$rotation[, 1]
    var_explained <- summary(pca_model)$importance[2, 1]

    # Sign correction: ensure positive correlation with row means
    avg_data <- rowMeans(datos_clean)
    scores <- datos_clean %*% loadings
    if (cor(scores, avg_data) < 0) {
      loadings <- -loadings
    }

    return(list(
      loadings = loadings,
      var_explained = var_explained,
      n_obs = nrow(datos_clean),
      period_start = min(data_filtered$fecha),
      period_end = max(data_filtered$fecha)
    ))
  }, error = function(e) {
    warning("PCA failed: ", e$message)
    return(NULL)
  })
}

#' Compare loadings between two regimes
#' @param loadings_pre Loadings from pre-regime
#' @param loadings_post Loadings from post-regime
#' @return Comparison data frame
compare_regime_loadings <- function(loadings_pre, loadings_post) {
  vars <- intersect(names(loadings_pre), names(loadings_post))

  comparison <- data.frame(
    variable = vars,
    loading_pre = loadings_pre[vars],
    loading_post = loadings_post[vars],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      abs_diff = abs(loading_post - loading_pre),
      pct_change = (loading_post - loading_pre) / abs(loading_pre) * 100,
      sign_change = sign(loading_pre) != sign(loading_post),
      direction = case_when(
        loading_post > loading_pre ~ "Increased",
        loading_post < loading_pre ~ "Decreased",
        TRUE ~ "No change"
      )
    ) %>%
    arrange(desc(abs_diff))

  return(comparison)
}

#' Bootstrap test for loading differences
#' @param data_pre Pre-regime data
#' @param data_post Post-regime data
#' @param vars Variables to include
#' @param n_boot Number of bootstrap replications
#' @return Data frame with p-values for each variable
bootstrap_loading_test <- function(data_pre, data_post, vars, n_boot = 1000) {

  vars_available <- intersect(vars, names(data_pre))
  vars_available <- intersect(vars_available, names(data_post))

  # Original loadings
  pca_pre <- calc_regime_pca(data_pre, vars_available)
  pca_post <- calc_regime_pca(data_post, vars_available)

  if (is.null(pca_pre) || is.null(pca_post)) {
    return(NULL)
  }

  original_diff <- pca_post$loadings[vars_available] - pca_pre$loadings[vars_available]

  # Bootstrap under null (pooled data)
  data_pooled <- rbind(
    data_pre %>% dplyr::select(all_of(vars_available)),
    data_post %>% dplyr::select(all_of(vars_available))
  )

  n_pre <- nrow(data_pre)
  n_post <- nrow(data_post)
  n_total <- n_pre + n_post

  boot_diffs <- matrix(NA, nrow = n_boot, ncol = length(vars_available))
  colnames(boot_diffs) <- vars_available

  for (b in 1:n_boot) {
    # Resample indices
    idx_pre <- sample(1:n_total, n_pre, replace = TRUE)
    idx_post <- sample(1:n_total, n_post, replace = TRUE)

    boot_pre <- data_pooled[idx_pre, , drop = FALSE]
    boot_post <- data_pooled[idx_post, , drop = FALSE]

    # Calculate loadings
    pca_boot_pre <- tryCatch({
      std_pre <- scale(boot_pre)
      std_pre[is.na(std_pre)] <- 0
      pca_m <- prcomp(std_pre, center = FALSE, scale. = FALSE)
      loadings <- pca_m$rotation[, 1]
      if (cor(std_pre %*% loadings, rowMeans(std_pre)) < 0) loadings <- -loadings
      loadings
    }, error = function(e) NULL)

    pca_boot_post <- tryCatch({
      std_post <- scale(boot_post)
      std_post[is.na(std_post)] <- 0
      pca_m <- prcomp(std_post, center = FALSE, scale. = FALSE)
      loadings <- pca_m$rotation[, 1]
      if (cor(std_post %*% loadings, rowMeans(std_post)) < 0) loadings <- -loadings
      loadings
    }, error = function(e) NULL)

    if (!is.null(pca_boot_pre) && !is.null(pca_boot_post)) {
      boot_diffs[b, ] <- pca_boot_post[vars_available] - pca_boot_pre[vars_available]
    }
  }

  # Calculate p-values (two-sided)
  p_values <- sapply(vars_available, function(v) {
    boot_vals <- boot_diffs[, v]
    boot_vals <- boot_vals[!is.na(boot_vals)]
    if (length(boot_vals) < 100) return(NA)
    mean(abs(boot_vals) >= abs(original_diff[v]))
  })

  return(data.frame(
    variable = vars_available,
    diff = original_diff,
    p_value = p_values,
    significant_10 = p_values < 0.10,
    significant_05 = p_values < 0.05
  ))
}

# Calculate regime-split PCA
cat("Calculating Pre-IT PCA loadings...\n")
pca_pre_IT <- calc_regime_pca(
  datos_signed,
  all_vars_core,
  end_date = REGIME_CONFIG$IT_START - 1
)

cat("Calculating Post-IT PCA loadings...\n")
pca_post_IT <- calc_regime_pca(
  datos_signed,
  all_vars_core,
  start_date = REGIME_CONFIG$IT_START
)

# Store results
regime_pca_results <- list(
  pre_IT = pca_pre_IT,
  post_IT = pca_post_IT
)

if (!is.null(pca_pre_IT) && !is.null(pca_post_IT)) {

  cat("\nPre-IT Period:\n")
  cat("  Observations:", pca_pre_IT$n_obs, "\n")
  cat("  Variance explained (PC1):", round(pca_pre_IT$var_explained * 100, 1), "%\n")

  cat("\nPost-IT Period:\n")
  cat("  Observations:", pca_post_IT$n_obs, "\n")
  cat("  Variance explained (PC1):", round(pca_post_IT$var_explained * 100, 1), "%\n")

  # Compare loadings
  loading_comparison <- compare_regime_loadings(
    pca_pre_IT$loadings,
    pca_post_IT$loadings
  )

  cat("\nLoading Comparison (Pre-IT vs Post-IT):\n")
  cat(sprintf("%-25s %10s %10s %10s %8s\n",
              "Variable", "Pre-IT", "Post-IT", "Diff", "Sign"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in 1:nrow(loading_comparison)) {
    r <- loading_comparison[i, ]
    sign_flag <- ifelse(r$sign_change, "CHANGE", "")
    cat(sprintf("%-25s %+10.3f %+10.3f %+10.3f %8s\n",
                r$variable, r$loading_pre, r$loading_post,
                r$loading_post - r$loading_pre, sign_flag))
  }

  # Correlation between loading vectors
  loading_corr <- cor(loading_comparison$loading_pre, loading_comparison$loading_post)
  cat(sprintf("\nLoading vector correlation: %.3f\n", loading_corr))

  if (loading_corr > 0.9) {
    cat(">>> HIGH STABILITY: Loading structure similar across regimes\n")
  } else if (loading_corr > 0.7) {
    cat(">>> MODERATE STABILITY: Some structural changes\n")
  } else {
    cat(">>> LOW STABILITY: Significant structural break\n")
  }

  # Bootstrap test
  cat("\nRunning bootstrap test for loading differences...\n")
  cat("(", REGIME_CONFIG$BOOTSTRAP_REPS, " replications)\n")

  data_pre_subset <- datos_signed %>%
    filter(fecha < REGIME_CONFIG$IT_START)
  data_post_subset <- datos_signed %>%
    filter(fecha >= REGIME_CONFIG$IT_START)

  bootstrap_results <- bootstrap_loading_test(
    data_pre_subset,
    data_post_subset,
    all_vars_core,
    n_boot = REGIME_CONFIG$BOOTSTRAP_REPS
  )

  if (!is.null(bootstrap_results)) {
    cat("\nBootstrap Test Results (H0: loading_pre = loading_post):\n")
    cat(sprintf("%-25s %10s %10s\n", "Variable", "Diff", "p-value"))
    cat(paste(rep("-", 50), collapse = ""), "\n")

    for (i in 1:nrow(bootstrap_results)) {
      r <- bootstrap_results[i, ]
      stars <- ifelse(r$p_value < 0.01, "***",
                      ifelse(r$p_value < 0.05, "**",
                             ifelse(r$p_value < 0.10, "*", "")))
      cat(sprintf("%-25s %+10.3f %10.3f%s\n",
                  r$variable, r$diff, r$p_value, stars))
    }

    n_sig <- sum(bootstrap_results$significant_05, na.rm = TRUE)
    cat(sprintf("\nSignificant changes (p<0.05): %d/%d variables\n",
                n_sig, nrow(bootstrap_results)))

    regime_pca_results$bootstrap <- bootstrap_results
    regime_pca_results$loading_comparison <- loading_comparison
  }
}


################################################################################
# SECTION 4: PART B - REGIME DUMMIES STORAGE
################################################################################

cat("\n================================================================================\n")
cat("PART B: REGIME DUMMIES STORAGE\n")
cat("================================================================================\n\n")

# Add regime columns to resultado_fci$all_indices
resultado_fci$all_indices <- resultado_fci$all_indices %>%
  create_regime_dummies()

# Export regime indicators
regime_indicators <- fci_data %>%
  dplyr::select(fecha, regime_IT, regime_COVID, regime_IT_exCOVID, period_label)

write.csv(regime_indicators,
          file.path(REGIME_CONFIG$output_dir, "FCI_Regime_Indicators.csv"),
          row.names = FALSE)
cat("Saved: FCI_Regime_Indicators.csv\n")

# Print regime period statistics
cat("\nFCI Statistics by Regime:\n\n")

regime_stats <- fci_data %>%
  group_by(period_label) %>%
  summarise(
    n = n(),
    FCI_mean = mean(FCI_COMP_AVG, na.rm = TRUE),
    FCI_sd = sd(FCI_COMP_AVG, na.rm = TRUE),
    FCI_min = min(FCI_COMP_AVG, na.rm = TRUE),
    FCI_max = max(FCI_COMP_AVG, na.rm = TRUE),
    pct_tightening = mean(FCI_COMP_AVG > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

print(regime_stats)


################################################################################
# SECTION 5: PART C - REGIME-DEPENDENT LOCAL PROJECTIONS
################################################################################

cat("\n================================================================================\n")
cat("PART C: REGIME-DEPENDENT LOCAL PROJECTIONS\n")
cat("================================================================================\n\n")

# Load macro data for LP
macro_raw <- tryCatch({
  read_excel(REGIME_CONFIG$data_file, sheet = "Datos_macro")
}, error = function(e) {
  read_excel(REGIME_CONFIG$data_file)
})

fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

macro_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Calculate YoY growth for macro variables
if ("IMAEP" %in% names(macro_data)) {
  macro_data <- macro_data %>%
    mutate(IMAEP_yoy = (IMAEP / lag(IMAEP, 12) - 1) * 100)
}

# Prepare credit growth
credit_data <- datos %>%
  dplyr::select(fecha, Creditos_Sector_privado_totales) %>%
  mutate(Credit_yoy = (Creditos_Sector_privado_totales /
                         lag(Creditos_Sector_privado_totales, 12) - 1) * 100) %>%
  dplyr::select(fecha, Credit_yoy)

# Prepare NPL data
npl_data <- datos %>%
  dplyr::select(fecha, Morosidad) %>%
  rename(NPL = Morosidad)

# Merge all data for LP analysis
lp_data <- fci_data %>%
  dplyr::select(fecha, FCI_COMP_AVG, FCI_ENDO_AVG, FCI_EXO_AVG,
                regime_IT, regime_COVID, regime_IT_exCOVID) %>%
  left_join(credit_data, by = "fecha") %>%
  left_join(npl_data, by = "fecha") %>%
  left_join(macro_data %>% dplyr::select(fecha, any_of("IMAEP_yoy")), by = "fecha") %>%
  na.omit()

cat("LP analysis data prepared:", nrow(lp_data), "observations\n\n")

#' Run Local Projection with regime interaction
#' Model: y_{t+h} = α + β₁·FCI_t + β₂·regime_t + β₃·FCI_t×regime_t + γ·controls + ε
#' @param data Data frame
#' @param y_var Dependent variable name
#' @param fci_var FCI variable name
#' @param regime_var Regime indicator name
#' @param max_h Maximum horizon
#' @param n_lags Number of lags
#' @return Data frame with results
run_lp_regime_interaction <- function(data, y_var, fci_var, regime_var, max_h, n_lags = 2) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - REGIME_CONFIG$confidence_level) / 2)

  for (h in 1:max_h) {
    # Create variables
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        fci_lag1 = lag(!!sym(fci_var), 1),
        fci_x_regime = !!sym(fci_var) * !!sym(regime_var)
      )

    # Add more lags
    for (i in 2:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", i) := lag(!!sym(y_var), i))
    }

    # Build formula
    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")
    formula_str <- paste("y_fwd ~", fci_var, "+", regime_var, "+ fci_x_regime +",
                         paste(lag_vars, collapse = " + "))

    # Estimate
    reg_data <- data_h %>%
      dplyr::select(y_fwd, !!sym(fci_var), !!sym(regime_var), fci_x_regime,
                    all_of(lag_vars)) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    # Extract coefficients
    idx_fci <- which(rownames(coef_test) == fci_var)
    idx_regime <- which(rownames(coef_test) == regime_var)
    idx_interact <- which(rownames(coef_test) == "fci_x_regime")

    # Effect in regime=0 is β₁
    # Effect in regime=1 is β₁ + β₃
    effect_base <- coef_test[idx_fci, 1]
    effect_interact <- coef_test[idx_interact, 1]
    effect_regime1 <- effect_base + effect_interact

    # SE for regime=1 effect: sqrt(var(β₁) + var(β₃) + 2*cov(β₁,β₃))
    se_regime1 <- sqrt(vcov_hac[idx_fci, idx_fci] + vcov_hac[idx_interact, idx_interact] +
                         2 * vcov_hac[idx_fci, idx_interact])

    results <- rbind(results, data.frame(
      horizon = h,
      # Base effect (regime=0)
      coef_base = effect_base,
      se_base = coef_test[idx_fci, 2],
      p_base = coef_test[idx_fci, 4],
      ci_base_lo = effect_base - z_crit * coef_test[idx_fci, 2],
      ci_base_hi = effect_base + z_crit * coef_test[idx_fci, 2],
      # Regime=1 effect (base + interaction)
      coef_regime1 = effect_regime1,
      se_regime1 = se_regime1,
      ci_regime1_lo = effect_regime1 - z_crit * se_regime1,
      ci_regime1_hi = effect_regime1 + z_crit * se_regime1,
      # Interaction (differential effect)
      coef_interact = effect_interact,
      se_interact = coef_test[idx_interact, 2],
      p_interact = coef_test[idx_interact, 4],
      n_obs = nrow(reg_data)
    ))
  }

  return(results)
}

#' Run split-sample Local Projections
#' @param data Data frame
#' @param y_var Dependent variable name
#' @param fci_var FCI variable name
#' @param regime_var Regime indicator name
#' @param max_h Maximum horizon
#' @param n_lags Number of lags
#' @return Data frame with results for both subsamples
run_lp_split_sample <- function(data, y_var, fci_var, regime_var, max_h, n_lags = 2) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - REGIME_CONFIG$confidence_level) / 2)

  # Split data
  data_0 <- data %>% filter(!!sym(regime_var) == 0)
  data_1 <- data %>% filter(!!sym(regime_var) == 1)

  for (h in 1:max_h) {
    for (regime_val in c(0, 1)) {
      data_subset <- if (regime_val == 0) data_0 else data_1

      if (nrow(data_subset) < 40) next

      # Create variables
      data_h <- data_subset %>%
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
      formula_str <- paste("y_fwd ~", fci_var, "+", paste(lag_vars, collapse = " + "))

      reg_data <- data_h %>%
        dplyr::select(y_fwd, !!sym(fci_var), all_of(lag_vars)) %>%
        na.omit()

      if (nrow(reg_data) < 30) next

      model <- lm(as.formula(formula_str), data = reg_data)
      vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
      coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

      idx_fci <- which(rownames(coef_test) == fci_var)

      results <- rbind(results, data.frame(
        horizon = h,
        regime = regime_val,
        coef = coef_test[idx_fci, 1],
        se = coef_test[idx_fci, 2],
        p_value = coef_test[idx_fci, 4],
        ci_lower = coef_test[idx_fci, 1] - z_crit * coef_test[idx_fci, 2],
        ci_upper = coef_test[idx_fci, 1] + z_crit * coef_test[idx_fci, 2],
        n_obs = nrow(reg_data)
      ))
    }
  }

  return(results)
}

# Run regime-dependent LP analyses
lp_regime_results <- list()

# Analysis 1: IT regime effect on credit
cat("Running LP: FCI_COMP → Credit by IT regime...\n")
lp_regime_results$IT_credit_interact <- run_lp_regime_interaction(
  lp_data, "Credit_yoy", "FCI_COMP_AVG", "regime_IT",
  REGIME_CONFIG$max_horizon, REGIME_CONFIG$n_lags
) %>% mutate(analysis = "IT_Credit", type = "interaction")

lp_regime_results$IT_credit_split <- run_lp_split_sample(
  lp_data, "Credit_yoy", "FCI_COMP_AVG", "regime_IT",
  REGIME_CONFIG$max_horizon, REGIME_CONFIG$n_lags
) %>% mutate(analysis = "IT_Credit", type = "split")

# Analysis 2: IT regime effect using endogenous FCI
cat("Running LP: FCI_ENDO → Credit by IT regime...\n")
lp_regime_results$IT_endo_interact <- run_lp_regime_interaction(
  lp_data, "Credit_yoy", "FCI_ENDO_AVG", "regime_IT",
  REGIME_CONFIG$max_horizon, REGIME_CONFIG$n_lags
) %>% mutate(analysis = "IT_Endo_Credit", type = "interaction")

# Analysis 3: COVID period effect
cat("Running LP: FCI_COMP → Credit by COVID regime...\n")
lp_regime_results$COVID_credit_interact <- run_lp_regime_interaction(
  lp_data, "Credit_yoy", "FCI_COMP_AVG", "regime_COVID",
  REGIME_CONFIG$max_horizon, REGIME_CONFIG$n_lags
) %>% mutate(analysis = "COVID_Credit", type = "interaction")

# Analysis 4: Clean IT comparison (excluding COVID)
cat("Running LP: FCI_COMP → Credit by IT (excl. COVID)...\n")
lp_data_exCOVID <- lp_data %>% filter(regime_COVID == 0)

lp_regime_results$IT_exCOVID_split <- run_lp_split_sample(
  lp_data_exCOVID, "Credit_yoy", "FCI_COMP_AVG", "regime_IT",
  REGIME_CONFIG$max_horizon, REGIME_CONFIG$n_lags
) %>% mutate(analysis = "IT_exCOVID_Credit", type = "split")

# Analysis 5: NPL effects by regime
if ("NPL" %in% names(lp_data)) {
  cat("Running LP: FCI_COMP → NPL by IT regime...\n")
  lp_regime_results$IT_NPL_interact <- run_lp_regime_interaction(
    lp_data, "NPL", "FCI_COMP_AVG", "regime_IT",
    REGIME_CONFIG$max_horizon, REGIME_CONFIG$n_lags
  ) %>% mutate(analysis = "IT_NPL", type = "interaction")
}

# Combine and export results
all_lp_regime <- bind_rows(lp_regime_results)

cat("\n================================================================================\n")
cat("REGIME LP RESULTS SUMMARY (at h=12)\n")
cat("================================================================================\n\n")

# IT regime results
it_12 <- lp_regime_results$IT_credit_interact %>% filter(horizon == 12)
if (nrow(it_12) > 0) {
  cat("IT REGIME EFFECT ON CREDIT (FCI_COMP → Credit_yoy, h=12):\n")
  cat(sprintf("  Pre-IT effect:  %+.2f pp (p=%.3f)\n", it_12$coef_base, it_12$p_base))
  cat(sprintf("  Post-IT effect: %+.2f pp\n", it_12$coef_regime1))
  cat(sprintf("  Differential:   %+.2f pp (p=%.3f)\n", it_12$coef_interact, it_12$p_interact))

  if (it_12$p_interact < 0.10) {
    cat("  >>> SIGNIFICANT regime difference detected\n")
  } else {
    cat("  >>> No significant regime difference\n")
  }
}

# COVID effect
covid_12 <- lp_regime_results$COVID_credit_interact %>% filter(horizon == 12)
if (nrow(covid_12) > 0) {
  cat("\nCOVID PERIOD EFFECT ON CREDIT (h=12):\n")
  cat(sprintf("  Non-COVID effect: %+.2f pp (p=%.3f)\n", covid_12$coef_base, covid_12$p_base))
  cat(sprintf("  COVID effect:     %+.2f pp\n", covid_12$coef_regime1))
  cat(sprintf("  Differential:     %+.2f pp (p=%.3f)\n", covid_12$coef_interact, covid_12$p_interact))
}

# Export LP regime results
write.csv(all_lp_regime %>% filter(type == "interaction"),
          file.path(REGIME_CONFIG$output_dir, "LP_Regime_IT_Results.csv"),
          row.names = FALSE)
cat("\nSaved: LP_Regime_IT_Results.csv\n")

write.csv(lp_regime_results$COVID_credit_interact,
          file.path(REGIME_CONFIG$output_dir, "LP_Regime_COVID_Results.csv"),
          row.names = FALSE)
cat("Saved: LP_Regime_COVID_Results.csv\n")


################################################################################
# SECTION 6: PART D - TIME-VARYING PARAMETER PCA
################################################################################

cat("\n================================================================================\n")
cat("PART D: TIME-VARYING PARAMETER PCA\n")
cat("================================================================================\n\n")

#' Calculate TVP-PCA using exponential weighting
#' @param data Data frame with fecha and variables
#' @param vars Variables to include
#' @param lambda Decay parameter (default 0.97)
#' @return Data frame with time-varying loadings
tvp_pca_exponential <- function(data, vars, lambda = 0.97) {

  vars_available <- intersect(vars, names(data))
  if (length(vars_available) < 2) return(NULL)

  data_sorted <- data %>% arrange(fecha)
  n_obs <- nrow(data_sorted)

  # Minimum observations to start
  min_start <- 36

  results <- list()

  for (t in min_start:n_obs) {
    # Calculate weights: w_s = λ^(t-s) for s ≤ t
    weights <- lambda^((t-1):0)
    weights <- weights / sum(weights)  # Normalize

    # Get data up to time t
    data_t <- data_sorted[1:t, vars_available, drop = FALSE]

    # Weighted covariance matrix
    data_centered <- scale(data_t, center = TRUE, scale = FALSE)

    # Handle NAs
    complete_rows <- complete.cases(data_centered)
    if (sum(complete_rows) < 20) next

    data_clean <- data_centered[complete_rows, , drop = FALSE]
    weights_clean <- weights[complete_rows]
    weights_clean <- weights_clean / sum(weights_clean)

    # Weighted covariance
    weighted_cov <- t(data_clean) %*% diag(weights_clean) %*% data_clean

    # Eigen decomposition
    eig <- tryCatch({
      eigen(weighted_cov)
    }, error = function(e) NULL)

    if (is.null(eig)) next

    # First eigenvector (loadings)
    loadings <- eig$vectors[, 1]
    names(loadings) <- vars_available
    var_explained <- eig$values[1] / sum(eig$values)

    # Sign correction
    scores <- data_clean %*% loadings
    if (cor(scores, rowMeans(data_clean)) < 0) {
      loadings <- -loadings
    }

    for (v in vars_available) {
      results[[length(results) + 1]] <- data.frame(
        fecha = data_sorted$fecha[t],
        variable = v,
        loading = loadings[v],
        var_explained = var_explained,
        n_effective = sum(weights_clean > 0.001)
      )
    }
  }

  return(bind_rows(results))
}

cat("Calculating TVP-PCA with exponential weighting (lambda =",
    REGIME_CONFIG$TVP_LAMBDA, ")...\n")

tvp_loadings <- tvp_pca_exponential(
  datos_signed,
  all_vars_core,
  lambda = REGIME_CONFIG$TVP_LAMBDA
)

if (!is.null(tvp_loadings) && nrow(tvp_loadings) > 0) {
  cat("TVP loadings calculated:", nrow(tvp_loadings), "observations\n")
  cat("Date range:", format(min(tvp_loadings$fecha)), "to",
      format(max(tvp_loadings$fecha)), "\n\n")

  # Summary statistics
  tvp_summary <- tvp_loadings %>%
    group_by(variable) %>%
    summarise(
      mean_loading = mean(loading, na.rm = TRUE),
      sd_loading = sd(loading, na.rm = TRUE),
      min_loading = min(loading, na.rm = TRUE),
      max_loading = max(loading, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(mean_loading)))

  cat("TVP Loading Summary:\n")
  cat(sprintf("%-25s %10s %10s %10s\n", "Variable", "Mean", "SD", "Range"))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  for (i in 1:nrow(tvp_summary)) {
    r <- tvp_summary[i, ]
    cat(sprintf("%-25s %+10.3f %10.3f %10.3f\n",
                r$variable, r$mean_loading, r$sd_loading,
                r$max_loading - r$min_loading))
  }

  # Export TVP loadings
  write.csv(tvp_loadings,
            file.path(REGIME_CONFIG$output_dir, "FCI_TVP_Loadings.csv"),
            row.names = FALSE)
  cat("\nSaved: FCI_TVP_Loadings.csv\n")

  # Compare TVP with rolling loadings if available
  rolling_file <- file.path(REGIME_CONFIG$output_dir, "FCI_Rolling_Loadings.csv")
  if (file.exists(rolling_file)) {
    cat("\nComparing TVP with rolling loadings...\n")

    rolling_loadings <- read.csv(rolling_file) %>%
      mutate(fecha = as.Date(fecha)) %>%
      filter(window_size == 60)

    # Merge and compare
    comparison <- tvp_loadings %>%
      inner_join(rolling_loadings %>% dplyr::select(fecha, variable, loading),
                 by = c("fecha", "variable"),
                 suffix = c("_tvp", "_rolling"))

    if (nrow(comparison) > 100) {
      corr_by_var <- comparison %>%
        group_by(variable) %>%
        summarise(
          correlation = cor(loading_tvp, loading_rolling, use = "complete.obs"),
          .groups = "drop"
        )

      cat("\nTVP vs Rolling (60m) Correlation by Variable:\n")
      print(corr_by_var)

      write.csv(corr_by_var,
                file.path(REGIME_CONFIG$output_dir, "FCI_TVP_vs_Rolling_Comparison.csv"),
                row.names = FALSE)
      cat("Saved: FCI_TVP_vs_Rolling_Comparison.csv\n")
    }
  }
}


################################################################################
# SECTION 6B: SIGN REVERSAL INVESTIGATION (Phase 1.2: Critical Fix)
################################################################################

cat("\n================================================================================\n")
cat("SIGN REVERSAL INVESTIGATION\n")
cat("================================================================================\n\n")

# Expected theoretical signs based on CLAUDE.md documentation
expected_signs <- data.frame(
  variable = c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM",
               "Crecimiento_creditos", "Ratio_Cred_Depo", "Morosidad",
               "Rentabilidad", "Liquidez", "TCN", "Commodities", "FFER", "VIX"),
  expected_sign = c(+1, +1, +1, -1, -1, +1, -1, -1, +1, -1, +1, +1),
  rationale = c(
    "Higher rates = tighter",
    "Higher spread = tighter",
    "Higher market spread = tighter",
    "Higher growth = looser",
    "Higher intermediation = looser",
    "Higher NPL = tighter",
    "Higher profitability = looser",
    "Higher liquidity = looser",
    "Depreciation = tighter (net importer)",
    "Higher prices = looser (soy/beef exporter)",
    "Higher US rates = tighter",
    "Higher volatility = tighter"
  ),
  stringsAsFactors = FALSE
)

# Analyze sign flips in regime PCA comparison
sign_flip_analysis <- NULL
if (!is.null(pca_pre_IT) && !is.null(pca_post_IT) && exists("loading_comparison")) {

  sign_flip_analysis <- loading_comparison %>%
    left_join(expected_signs, by = "variable") %>%
    mutate(
      # Determine actual signs in each period
      sign_pre = sign(loading_pre),
      sign_post = sign(loading_post),

      # Check for sign flips
      sign_flipped = sign_pre != sign_post,

      # Check consistency with expected sign
      consistent_pre = (sign_pre == expected_sign) | is.na(expected_sign),
      consistent_post = (sign_post == expected_sign) | is.na(expected_sign),

      # Categorize the issue
      issue_category = case_when(
        !sign_flipped & consistent_pre & consistent_post ~ "No issue",
        !sign_flipped & !consistent_pre & !consistent_post ~ "Both periods: unexpected sign",
        sign_flipped & consistent_pre & !consistent_post ~ "SIGN FLIP: Post-IT deviates",
        sign_flipped & !consistent_pre & consistent_post ~ "SIGN FLIP: Pre-IT deviates",
        sign_flipped ~ "SIGN FLIP: Both deviate from expected",
        !consistent_pre ~ "Pre-IT: unexpected sign",
        !consistent_post ~ "Post-IT: unexpected sign",
        TRUE ~ "Unknown"
      )
    )

  cat("Sign Flip Analysis Results:\n")
  cat(sprintf("%-25s %8s %8s %8s %8s %s\n",
              "Variable", "Pre-IT", "Post-IT", "Exp.", "Flip?", "Issue"))
  cat(paste(rep("-", 90), collapse = ""), "\n")

  for (i in 1:nrow(sign_flip_analysis)) {
    r <- sign_flip_analysis[i, ]
    flip_flag <- ifelse(r$sign_flipped, "YES", "")
    exp_sign <- ifelse(!is.na(r$expected_sign),
                       ifelse(r$expected_sign > 0, "+", "-"), "?")
    pre_sign <- ifelse(r$sign_pre > 0, "+", "-")
    post_sign <- ifelse(r$sign_post > 0, "+", "-")

    cat(sprintf("%-25s %8s %8s %8s %8s %s\n",
                r$variable, pre_sign, post_sign, exp_sign, flip_flag,
                ifelse(r$issue_category != "No issue", r$issue_category, "")))
  }

  # Count sign flips
  n_flips <- sum(sign_flip_analysis$sign_flipped, na.rm = TRUE)
  cat(sprintf("\nTotal sign flips: %d/%d variables\n", n_flips, nrow(sign_flip_analysis)))

  # Identify problematic variables
  problematic <- sign_flip_analysis %>%
    filter(issue_category != "No issue")

  if (nrow(problematic) > 0) {
    cat("\n*** PROBLEMATIC VARIABLES REQUIRING ATTENTION ***\n")
    for (i in 1:nrow(problematic)) {
      r <- problematic[i, ]
      cat(sprintf("  - %s: %s\n    Rationale: %s\n",
                  r$variable, r$issue_category,
                  ifelse(!is.na(r$rationale), r$rationale, "N/A")))
    }
  }

  # Economic interpretation
  cat("\n================================================================================\n")
  cat("ECONOMIC INTERPRETATION OF SIGN FLIPS\n")
  cat("================================================================================\n\n")

  # Check specific variables mentioned in the plan
  spread_flip <- sign_flip_analysis %>% filter(variable == "Spread_activas_pasivas")
  commodities_flip <- sign_flip_analysis %>% filter(variable == "Commodities")

  if (nrow(spread_flip) > 0 && spread_flip$sign_flipped) {
    cat("1. Spread_activas_pasivas sign flip:\n")
    cat("   Pre-IT:  Loading = ", round(spread_flip$loading_pre, 3), "\n")
    cat("   Post-IT: Loading = ", round(spread_flip$loading_post, 3), "\n")
    cat("   HYPOTHESIS: Post-IT, high spreads may reflect strong credit DEMAND\n")
    cat("   (endogenous to growth) rather than supply constraints.\n")
    cat("   Under IT, banks may widen spreads during expansion phases.\n\n")
  }

  if (nrow(commodities_flip) > 0 && commodities_flip$sign_flipped) {
    cat("2. Commodities sign flip:\n")
    cat("   Pre-IT:  Loading = ", round(commodities_flip$loading_pre, 3), "\n")
    cat("   Post-IT: Loading = ", round(commodities_flip$loading_post, 3), "\n")
    cat("   HYPOTHESIS: Pre-IT, BCP may have raised rates aggressively during\n")
    cat("   commodity booms (policy reaction function artifact).\n")
    cat("   Under IT, commodity booms translate to looser conditions (exports).\n\n")
  }

  cat("3. Structural interpretation:\n")
  cat("   The May 2011 IT adoption changed Paraguay's monetary policy framework.\n")
  cat("   - Pre-IT: Exchange rate as de facto nominal anchor\n")
  cat("   - Post-IT: Interest rate as primary policy instrument\n")
  cat("   This explains many loading changes - they reflect genuine\n")
  cat("   structural shifts, not data errors.\n\n")

  # Export sign flip analysis
  write.csv(sign_flip_analysis,
            file.path(REGIME_CONFIG$output_dir, "FCI_Sign_Flip_Diagnostic.csv"),
            row.names = FALSE)
  cat("Saved: FCI_Sign_Flip_Diagnostic.csv\n")

  # Export problematic variables separately
  if (nrow(problematic) > 0) {
    write.csv(problematic %>% dplyr::select(variable, loading_pre, loading_post,
                                             expected_sign, issue_category, rationale),
              file.path(REGIME_CONFIG$output_dir, "FCI_Inconsistent_Loadings.csv"),
              row.names = FALSE)
    cat("Saved: FCI_Inconsistent_Loadings.csv\n")
  }
}

# Additional TVP sign consistency analysis
if (!is.null(tvp_loadings) && nrow(tvp_loadings) > 0) {
  cat("\n================================================================================\n")
  cat("TVP SIGN CONSISTENCY OVER TIME\n")
  cat("================================================================================\n\n")

  tvp_sign_analysis <- tvp_loadings %>%
    left_join(expected_signs, by = "variable") %>%
    group_by(variable) %>%
    summarise(
      pct_positive = mean(loading > 0, na.rm = TRUE) * 100,
      pct_expected_sign = mean(sign(loading) == expected_sign, na.rm = TRUE) * 100,
      n_sign_changes = sum(diff(sign(loading)) != 0, na.rm = TRUE),
      expected = first(expected_sign),
      .groups = "drop"
    ) %>%
    mutate(
      expected_str = ifelse(expected > 0, "+", "-"),
      stability = case_when(
        pct_expected_sign >= 90 ~ "Highly stable",
        pct_expected_sign >= 70 ~ "Moderately stable",
        pct_expected_sign >= 50 ~ "Unstable",
        TRUE ~ "Inverted"
      )
    ) %>%
    arrange(pct_expected_sign)

  cat("Sign Consistency with Theory (TVP Analysis):\n")
  cat(sprintf("%-25s %8s %12s %12s %15s\n",
              "Variable", "Expected", "% Correct", "# Changes", "Stability"))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (i in 1:nrow(tvp_sign_analysis)) {
    r <- tvp_sign_analysis[i, ]
    cat(sprintf("%-25s %8s %12.1f%% %12d %15s\n",
                r$variable, r$expected_str, r$pct_expected_sign,
                r$n_sign_changes, r$stability))
  }

  # Save TVP sign analysis
  write.csv(tvp_sign_analysis,
            file.path(REGIME_CONFIG$output_dir, "FCI_TVP_Sign_Consistency.csv"),
            row.names = FALSE)
  cat("\nSaved: FCI_TVP_Sign_Consistency.csv\n")
}


################################################################################
# SECTION 7: PART E - THRESHOLD EFFECTS TESTING
################################################################################

cat("\n================================================================================\n")
cat("PART E: THRESHOLD EFFECTS TESTING\n")
cat("================================================================================\n\n")

#' Create threshold indicator
#' @param fci FCI vector
#' @param threshold_sd Threshold in standard deviations
#' @return List with indicator and decomposed FCI
create_threshold_indicator <- function(fci, threshold_sd = 1.0) {
  fci_std <- (fci - mean(fci, na.rm = TRUE)) / sd(fci, na.rm = TRUE)

  high <- as.integer(fci_std > threshold_sd)
  low <- as.integer(fci_std <= threshold_sd)

  fci_high <- fci * high
  fci_low <- fci * low

  return(list(
    indicator = high,
    fci_high = fci_high,
    fci_low = fci_low,
    pct_high = mean(high, na.rm = TRUE) * 100
  ))
}

# Create threshold variables
threshold_result <- create_threshold_indicator(
  lp_data$FCI_COMP_AVG,
  REGIME_CONFIG$FCI_THRESHOLD
)

lp_data <- lp_data %>%
  mutate(
    FCI_high = threshold_result$fci_high,
    FCI_low = threshold_result$fci_low,
    crisis_indicator = threshold_result$indicator
  )

cat("Threshold Analysis (FCI > ", REGIME_CONFIG$FCI_THRESHOLD, " SD):\n")
cat("  Crisis periods:", round(threshold_result$pct_high, 1), "%\n")
cat("  Normal periods:", round(100 - threshold_result$pct_high, 1), "%\n\n")

#' Run threshold LP
#' Model: y_{t+h} = α + β_high·FCI_high + β_low·FCI_low + controls + ε
run_lp_threshold <- function(data, y_var, max_h, n_lags = 2) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - REGIME_CONFIG$confidence_level) / 2)

  for (h in 1:max_h) {
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        fci_lag1 = lag(FCI_COMP_AVG, 1)
      )

    for (i in 2:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", i) := lag(!!sym(y_var), i))
    }

    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")
    formula_str <- paste("y_fwd ~ FCI_high + FCI_low +", paste(lag_vars, collapse = " + "))

    reg_data <- data_h %>%
      dplyr::select(y_fwd, FCI_high, FCI_low, all_of(lag_vars)) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    idx_high <- which(rownames(coef_test) == "FCI_high")
    idx_low <- which(rownames(coef_test) == "FCI_low")

    # Test H0: β_high = β_low
    coef_diff <- coef_test[idx_high, 1] - coef_test[idx_low, 1]
    var_diff <- vcov_hac[idx_high, idx_high] + vcov_hac[idx_low, idx_low] -
      2 * vcov_hac[idx_high, idx_low]
    wald_stat <- coef_diff^2 / var_diff
    wald_p <- 1 - pchisq(wald_stat, df = 1)

    results <- rbind(results, data.frame(
      horizon = h,
      coef_high = coef_test[idx_high, 1],
      se_high = coef_test[idx_high, 2],
      p_high = coef_test[idx_high, 4],
      ci_high_lo = coef_test[idx_high, 1] - z_crit * coef_test[idx_high, 2],
      ci_high_hi = coef_test[idx_high, 1] + z_crit * coef_test[idx_high, 2],
      coef_low = coef_test[idx_low, 1],
      se_low = coef_test[idx_low, 2],
      p_low = coef_test[idx_low, 4],
      ci_low_lo = coef_test[idx_low, 1] - z_crit * coef_test[idx_low, 2],
      ci_low_hi = coef_test[idx_low, 1] + z_crit * coef_test[idx_low, 2],
      diff_p = wald_p,
      n_obs = nrow(reg_data)
    ))
  }

  return(results)
}

# Run threshold LP for credit
cat("Running threshold LP: FCI → Credit...\n")
threshold_credit <- run_lp_threshold(lp_data, "Credit_yoy",
                                      REGIME_CONFIG$max_horizon, REGIME_CONFIG$n_lags)

# Run threshold LP for NPL if available
if ("NPL" %in% names(lp_data)) {
  cat("Running threshold LP: FCI → NPL...\n")
  threshold_npl <- run_lp_threshold(lp_data, "NPL",
                                     REGIME_CONFIG$max_horizon, REGIME_CONFIG$n_lags)
}

# Results summary
cat("\nTHRESHOLD LP RESULTS (at h=12):\n\n")

credit_12 <- threshold_credit %>% filter(horizon == 12)
if (nrow(credit_12) > 0) {
  cat("Credit Growth Response:\n")
  cat(sprintf("  Crisis (FCI>1SD):  %+.2f pp (p=%.3f)\n",
              credit_12$coef_high, credit_12$p_high))
  cat(sprintf("  Normal (FCI≤1SD):  %+.2f pp (p=%.3f)\n",
              credit_12$coef_low, credit_12$p_low))
  cat(sprintf("  Equality test p-value: %.3f\n", credit_12$diff_p))

  if (credit_12$diff_p < 0.10) {
    if (abs(credit_12$coef_high) > abs(credit_12$coef_low)) {
      cat("  >>> FCI effects STRONGER in crisis periods\n")
    } else {
      cat("  >>> FCI effects WEAKER in crisis periods\n")
    }
  } else {
    cat("  >>> No significant difference between crisis and normal periods\n")
  }
}

# Export threshold results
threshold_indicator_df <- lp_data %>%
  dplyr::select(fecha, FCI_COMP_AVG, crisis_indicator, FCI_high, FCI_low)

write.csv(threshold_indicator_df,
          file.path(REGIME_CONFIG$output_dir, "FCI_Threshold_Indicator.csv"),
          row.names = FALSE)
cat("\nSaved: FCI_Threshold_Indicator.csv\n")

write.csv(threshold_credit,
          file.path(REGIME_CONFIG$output_dir, "LP_Threshold_Effects.csv"),
          row.names = FALSE)
cat("Saved: LP_Threshold_Effects.csv\n")

# =============================================================================
# ALTERNATIVE THRESHOLD SENSITIVITY (Phase 2.2)
# =============================================================================

cat("\n================================================================================\n")
cat("ALTERNATIVE THRESHOLD SENSITIVITY (1.0, 1.5, 2.0 SD)\n")
cat("================================================================================\n\n")

threshold_levels <- c(1.0, 1.5, 2.0)
threshold_sensitivity_results <- list()

for (thresh in threshold_levels) {
  cat(sprintf("Testing threshold: %.1f SD...\n", thresh))

  # Create threshold variables for this level
  thresh_result <- create_threshold_indicator(lp_data$FCI_COMP_AVG, thresh)

  cat(sprintf("  Crisis periods (FCI > %.1f SD): %.1f%%\n", thresh, thresh_result$pct_high))

  # Create temporary data with this threshold
  lp_data_temp <- lp_data %>%
    mutate(
      FCI_high = thresh_result$fci_high,
      FCI_low = thresh_result$fci_low
    )

  # Run LP for credit
  thresh_lp <- tryCatch({
    run_lp_threshold(lp_data_temp, "Credit_yoy",
                     REGIME_CONFIG$max_horizon, REGIME_CONFIG$n_lags)
  }, error = function(e) NULL)

  if (!is.null(thresh_lp) && nrow(thresh_lp) > 0) {
    thresh_lp$threshold_sd <- thresh
    thresh_lp$pct_crisis <- thresh_result$pct_high
    threshold_sensitivity_results[[as.character(thresh)]] <- thresh_lp
  }
}

# Combine and export
if (length(threshold_sensitivity_results) > 0) {
  all_threshold_sensitivity <- bind_rows(threshold_sensitivity_results)

  # Summary at h=12
  cat("\nThreshold Sensitivity Results at h=12:\n")
  cat(sprintf("%-12s %12s %12s %12s %12s %12s\n",
              "Threshold", "% Crisis", "Coef High", "Coef Low", "Diff", "p(Diff)"))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (thresh in threshold_levels) {
    result_12 <- all_threshold_sensitivity %>%
      filter(threshold_sd == thresh, horizon == 12)
    if (nrow(result_12) > 0) {
      diff_effect <- result_12$coef_high - result_12$coef_low
      cat(sprintf("%-12.1f %12.1f%% %+12.2f %+12.2f %+12.2f %12.3f\n",
                  thresh, result_12$pct_crisis,
                  result_12$coef_high, result_12$coef_low,
                  diff_effect, result_12$diff_p))
    }
  }

  # Export sensitivity results
  write.csv(all_threshold_sensitivity,
            file.path(REGIME_CONFIG$output_dir, "LP_Threshold_Sensitivity.csv"),
            row.names = FALSE)
  cat("\nSaved: LP_Threshold_Sensitivity.csv\n")

  # Interpretation
  cat("\nInterpretation:\n")
  cat("  - Higher thresholds isolate more extreme crisis episodes\n")
  cat("  - Coefficient magnitudes should increase if FCI effects are non-linear\n")
  cat("  - Consistent sign across thresholds supports robustness\n")
}


################################################################################
# SECTION 8: VISUALIZATIONS
################################################################################

cat("\n================================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

if (!dir.exists(REGIME_CONFIG$output_dir)) {
  dir.create(REGIME_CONFIG$output_dir, recursive = TRUE)
}

# Color palettes
colors_regime <- c("Pre-IT" = "#E74C3C", "Post-IT" = "#3498DB", "COVID" = "#9B59B6")
colors_threshold <- c("Crisis" = "#E74C3C", "Normal" = "#27AE60")

# Plot 80: Regime-Split PCA Loadings
if (!is.null(pca_pre_IT) && !is.null(pca_post_IT)) {

  loading_plot_data <- loading_comparison %>%
    pivot_longer(cols = c(loading_pre, loading_post),
                 names_to = "period", values_to = "loading") %>%
    mutate(
      period = ifelse(period == "loading_pre", "Pre-IT", "Post-IT"),
      variable = factor(variable, levels = loading_comparison$variable)
    )

  p80 <- ggplot(loading_plot_data, aes(x = variable, y = loading, fill = period)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, color = "gray30") +
    scale_fill_manual(values = c("Pre-IT" = "#E74C3C", "Post-IT" = "#3498DB")) +
    coord_flip() +
    theme_minimal(base_size = 11) +
    labs(title = "PCA Loadings by Monetary Policy Regime",
         subtitle = paste0("Pre-IT (n=", pca_pre_IT$n_obs, ") vs Post-IT (n=", pca_post_IT$n_obs, ")"),
         x = NULL, y = "PC1 Loading", fill = "Period") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")

  ggsave(file.path(REGIME_CONFIG$output_dir, "80_Regime_PCA_Loadings.png"), p80,
         width = 10, height = 8, dpi = 300)
  cat("Saved: 80_Regime_PCA_Loadings.png\n")
}

# Plot 81: Regime Timeline with FCI
regime_plot_data <- fci_data %>%
  dplyr::select(fecha, FCI_COMP_AVG, regime_IT, regime_COVID)

p81 <- ggplot(regime_plot_data, aes(x = fecha)) +
  # Shade IT period
  annotate("rect", xmin = REGIME_CONFIG$IT_START, xmax = max(regime_plot_data$fecha),
           ymin = -Inf, ymax = Inf, fill = "#3498DB", alpha = 0.1) +
  # Shade COVID period
  annotate("rect", xmin = REGIME_CONFIG$COVID_START, xmax = REGIME_CONFIG$COVID_END,
           ymin = -Inf, ymax = Inf, fill = "#9B59B6", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "gray70") +
  geom_line(aes(y = FCI_COMP_AVG), color = "#2C3E50", linewidth = 0.8) +
  geom_vline(xintercept = REGIME_CONFIG$IT_START, linetype = "dashed", color = "#3498DB") +
  annotate("text", x = REGIME_CONFIG$IT_START, y = max(regime_plot_data$FCI_COMP_AVG, na.rm = TRUE),
           label = "IT Start", hjust = -0.1, vjust = 1, color = "#3498DB", size = 3) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  theme_minimal(base_size = 11) +
  labs(title = "FCI Timeline with Regime Periods",
       subtitle = "Blue shading = IT period | Purple shading = COVID period",
       x = NULL, y = "FCI (standard deviations)") +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(REGIME_CONFIG$output_dir, "81_Regime_Timeline_FCI.png"), p81,
       width = 12, height = 6, dpi = 300)
cat("Saved: 81_Regime_Timeline_FCI.png\n")

# Plot 82: LP IT Regime Effects
if (nrow(lp_regime_results$IT_credit_split) > 0) {

  split_data <- lp_regime_results$IT_credit_split %>%
    mutate(regime_label = ifelse(regime == 0, "Pre-IT", "Post-IT"))

  p82 <- ggplot(split_data, aes(x = horizon, y = coef, color = regime_label, fill = regime_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Pre-IT" = "#E74C3C", "Post-IT" = "#3498DB")) +
    scale_fill_manual(values = c("Pre-IT" = "#E74C3C", "Post-IT" = "#3498DB")) +
    theme_minimal(base_size = 11) +
    labs(title = "FCI Effect on Credit: Pre-IT vs Post-IT",
         subtitle = "Split-sample Local Projections | 90% CI",
         x = "Months", y = "Effect on credit growth (pp)",
         color = "Period", fill = "Period") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")

  ggsave(file.path(REGIME_CONFIG$output_dir, "82_LP_IT_Regime_Effects.png"), p82,
         width = 10, height = 6, dpi = 300)
  cat("Saved: 82_LP_IT_Regime_Effects.png\n")
}

# Plot 83: LP COVID Regime Effects
if (nrow(lp_regime_results$COVID_credit_interact) > 0) {

  covid_plot_data <- lp_regime_results$COVID_credit_interact %>%
    dplyr::select(horizon, coef_base, ci_base_lo, ci_base_hi, coef_regime1, ci_regime1_lo, ci_regime1_hi) %>%
    pivot_longer(cols = -horizon,
                 names_to = c("metric", "type"),
                 names_sep = "_(?=[^_]+$)",
                 values_to = "value") %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(period = ifelse(type %in% c("base", "lo", "hi"), "Non-COVID", "COVID"))

  # Simpler approach
  covid_plot_non <- lp_regime_results$COVID_credit_interact %>%
    dplyr::select(horizon, coef = coef_base, ci_lo = ci_base_lo, ci_hi = ci_base_hi) %>%
    mutate(period = "Non-COVID")

  covid_plot_covid <- lp_regime_results$COVID_credit_interact %>%
    dplyr::select(horizon, coef = coef_regime1, ci_lo = ci_regime1_lo, ci_hi = ci_regime1_hi) %>%
    mutate(period = "COVID")

  covid_plot_data <- bind_rows(covid_plot_non, covid_plot_covid)

  p83 <- ggplot(covid_plot_data, aes(x = horizon, y = coef, color = period, fill = period)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Non-COVID" = "#27AE60", "COVID" = "#9B59B6")) +
    scale_fill_manual(values = c("Non-COVID" = "#27AE60", "COVID" = "#9B59B6")) +
    theme_minimal(base_size = 11) +
    labs(title = "FCI Effect on Credit: COVID vs Non-COVID Periods",
         subtitle = "Interaction model Local Projections | 90% CI",
         x = "Months", y = "Effect on credit growth (pp)",
         color = "Period", fill = "Period") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")

  ggsave(file.path(REGIME_CONFIG$output_dir, "83_LP_COVID_Regime_Effects.png"), p83,
         width = 10, height = 6, dpi = 300)
  cat("Saved: 83_LP_COVID_Regime_Effects.png\n")
}

# Plot 84: TVP-PCA Loadings for Key Variables
if (!is.null(tvp_loadings) && nrow(tvp_loadings) > 0) {

  key_vars <- c("TPM", "TCN", "VIX", "Crecimiento_creditos", "Morosidad")
  key_vars <- key_vars[key_vars %in% unique(tvp_loadings$variable)]

  tvp_plot_data <- tvp_loadings %>%
    filter(variable %in% key_vars)

  p84 <- ggplot(tvp_plot_data, aes(x = fecha, y = loading, color = variable)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = REGIME_CONFIG$IT_START, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, color = "gray70") +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    theme_minimal(base_size = 11) +
    labs(title = "Time-Varying PCA Loadings (Exponential Weighting)",
         subtitle = paste0("Lambda = ", REGIME_CONFIG$TVP_LAMBDA,
                           " | Vertical line = IT adoption (May 2011)"),
         x = NULL, y = "PC1 Loading", color = "Variable") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(REGIME_CONFIG$output_dir, "84_TVP_PCA_Loadings.png"), p84,
         width = 12, height = 6, dpi = 300)
  cat("Saved: 84_TVP_PCA_Loadings.png\n")
}

# Plot 85: TVP vs Rolling Comparison (if both available)
rolling_file <- file.path(REGIME_CONFIG$output_dir, "FCI_Rolling_Loadings.csv")
if (!is.null(tvp_loadings) && file.exists(rolling_file)) {

  rolling_loadings <- read.csv(rolling_file) %>%
    mutate(fecha = as.Date(fecha)) %>%
    filter(window_size == 60)

  # Select one variable for comparison
  comparison_var <- "TPM"
  if (comparison_var %in% unique(tvp_loadings$variable)) {

    tvp_comp <- tvp_loadings %>%
      filter(variable == comparison_var) %>%
      dplyr::select(fecha, loading) %>%
      mutate(method = "TVP (Exponential)")

    roll_comp <- rolling_loadings %>%
      filter(variable == comparison_var) %>%
      dplyr::select(fecha, loading) %>%
      mutate(method = "Rolling (60m)")

    comp_data <- bind_rows(tvp_comp, roll_comp)

    p85 <- ggplot(comp_data, aes(x = fecha, y = loading, color = method)) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = REGIME_CONFIG$IT_START, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = 0, color = "gray70") +
      scale_color_manual(values = c("TVP (Exponential)" = "#E74C3C", "Rolling (60m)" = "#3498DB")) +
      scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
      theme_minimal(base_size = 11) +
      labs(title = paste("TVP vs Rolling Loadings:", comparison_var),
           subtitle = "Comparing exponential weighting with 60-month rolling window",
           x = NULL, y = "PC1 Loading", color = "Method") +
      theme(plot.title = element_text(face = "bold", size = 14),
            legend.position = "bottom")

    ggsave(file.path(REGIME_CONFIG$output_dir, "85_TVP_vs_Rolling_Comparison.png"), p85,
           width = 12, height = 6, dpi = 300)
    cat("Saved: 85_TVP_vs_Rolling_Comparison.png\n")
  }
}

# Plot 86: FCI with Threshold
p86 <- ggplot(lp_data, aes(x = fecha)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = REGIME_CONFIG$FCI_THRESHOLD, linetype = "dotted", color = "#E74C3C") +
  geom_hline(yintercept = -REGIME_CONFIG$FCI_THRESHOLD, linetype = "dotted", color = "#E74C3C") +
  geom_area(data = lp_data %>% filter(crisis_indicator == 1),
            aes(y = FCI_COMP_AVG), fill = "#E74C3C", alpha = 0.3) +
  geom_line(aes(y = FCI_COMP_AVG), color = "#2C3E50", linewidth = 0.8) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  theme_minimal(base_size = 11) +
  labs(title = "FCI with Crisis Threshold",
       subtitle = paste0("Shaded = Crisis periods (FCI > ", REGIME_CONFIG$FCI_THRESHOLD, " SD)"),
       x = NULL, y = "FCI (standard deviations)") +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(REGIME_CONFIG$output_dir, "86_FCI_Threshold.png"), p86,
       width = 12, height = 6, dpi = 300)
cat("Saved: 86_FCI_Threshold.png\n")

# Plot 87: Threshold LP Credit
if (nrow(threshold_credit) > 0) {

  threshold_plot_data <- threshold_credit %>%
    dplyr::select(horizon,
                  Crisis = coef_high, Crisis_lo = ci_high_lo, Crisis_hi = ci_high_hi,
                  Normal = coef_low, Normal_lo = ci_low_lo, Normal_hi = ci_low_hi) %>%
    pivot_longer(cols = c(Crisis, Normal),
                 names_to = "state", values_to = "coef") %>%
    mutate(
      ci_lo = ifelse(state == "Crisis", Crisis_lo, Normal_lo),
      ci_hi = ifelse(state == "Crisis", Crisis_hi, Normal_hi)
    ) %>%
    dplyr::select(horizon, state, coef, ci_lo, ci_hi)

  # Fix: properly reshape
  crisis_data <- threshold_credit %>%
    dplyr::select(horizon, coef = coef_high, ci_lo = ci_high_lo, ci_hi = ci_high_hi) %>%
    mutate(state = "Crisis (FCI > 1 SD)")

  normal_data <- threshold_credit %>%
    dplyr::select(horizon, coef = coef_low, ci_lo = ci_low_lo, ci_hi = ci_low_hi) %>%
    mutate(state = "Normal (FCI <= 1 SD)")

  threshold_plot_data <- bind_rows(crisis_data, normal_data)

  p87 <- ggplot(threshold_plot_data, aes(x = horizon, y = coef, color = state, fill = state)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Crisis (FCI > 1 SD)" = "#E74C3C",
                                   "Normal (FCI <= 1 SD)" = "#27AE60")) +
    scale_fill_manual(values = c("Crisis (FCI > 1 SD)" = "#E74C3C",
                                  "Normal (FCI <= 1 SD)" = "#27AE60")) +
    theme_minimal(base_size = 11) +
    labs(title = "Threshold LP: FCI Effect on Credit",
         subtitle = "Comparing crisis vs normal periods | 90% CI",
         x = "Months", y = "Effect on credit growth (pp)",
         color = "State", fill = "State") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")

  ggsave(file.path(REGIME_CONFIG$output_dir, "87_Threshold_LP_Credit.png"), p87,
         width = 10, height = 6, dpi = 300)
  cat("Saved: 87_Threshold_LP_Credit.png\n")
}

# Plot 88: Loading Stability Heatmap by Regime
if (!is.null(pca_pre_IT) && !is.null(pca_post_IT)) {

  heatmap_data <- loading_comparison %>%
    dplyr::select(variable, `Pre-IT` = loading_pre, `Post-IT` = loading_post) %>%
    pivot_longer(cols = c(`Pre-IT`, `Post-IT`),
                 names_to = "regime", values_to = "loading") %>%
    mutate(variable = factor(variable, levels = rev(loading_comparison$variable)))

  p88 <- ggplot(heatmap_data, aes(x = regime, y = variable, fill = loading)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%+.2f", loading)), size = 3.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Loading") +
    theme_minimal(base_size = 11) +
    labs(title = "PCA Loading Stability by Regime",
         subtitle = "First principal component loadings",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold", size = 14))

  ggsave(file.path(REGIME_CONFIG$output_dir, "88_Loading_Stability_Heatmap.png"), p88,
         width = 8, height = 10, dpi = 300)
  cat("Saved: 88_Loading_Stability_Heatmap.png\n")
}

# Plot 89: Summary Dashboard
p89_a <- p81 + theme(legend.position = "none") +
  labs(title = "FCI with Regime Periods")

p89_b <- if (exists("p82")) {
  p82 + theme(legend.position = "none") + labs(title = "IT Regime LP")
} else {
  ggplot() + theme_void() + labs(title = "IT Regime LP (N/A)")
}

p89_c <- if (exists("p87")) {
  p87 + theme(legend.position = "none") + labs(title = "Threshold LP")
} else {
  ggplot() + theme_void() + labs(title = "Threshold LP (N/A)")
}

p89_d <- if (exists("p80")) {
  p80 + theme(legend.position = "none") + labs(title = "Regime PCA Loadings")
} else {
  ggplot() + theme_void() + labs(title = "PCA Loadings (N/A)")
}

p89 <- grid.arrange(
  p89_a, p89_b, p89_c, p89_d,
  ncol = 2,
  top = grid::textGrob("Regime and TVP Analysis Summary",
                       gp = grid::gpar(fontsize = 16, fontface = "bold"))
)

ggsave(file.path(REGIME_CONFIG$output_dir, "89_Summary_Dashboard.png"), p89,
       width = 14, height = 12, dpi = 300)
cat("Saved: 89_Summary_Dashboard.png\n")


################################################################################
# SECTION 9: EXPORT CONSOLIDATED RESULTS
################################################################################

cat("\n================================================================================\n")
cat("EXPORTING CONSOLIDATED RESULTS\n")
cat("================================================================================\n\n")

# Export regime PCA comparison
if (exists("loading_comparison")) {
  regime_pca_export <- loading_comparison
  if (!is.null(regime_pca_results$bootstrap)) {
    regime_pca_export <- regime_pca_export %>%
      left_join(regime_pca_results$bootstrap %>% dplyr::select(variable, p_value, significant_05),
                by = "variable")
  }

  write.csv(regime_pca_export,
            file.path(REGIME_CONFIG$output_dir, "FCI_Regime_PCA_Comparison.csv"),
            row.names = FALSE)
  cat("Saved: FCI_Regime_PCA_Comparison.csv\n")
}

# List all generated outputs
cat("\n================================================================================\n")
cat("OUTPUT FILES SUMMARY\n")
cat("================================================================================\n\n")

cat("CSV Files:\n")
csv_files <- c(
  "FCI_Regime_Indicators.csv",
  "FCI_Regime_PCA_Comparison.csv",
  "LP_Regime_IT_Results.csv",
  "LP_Regime_COVID_Results.csv",
  "FCI_TVP_Loadings.csv",
  "FCI_TVP_vs_Rolling_Comparison.csv",
  "FCI_Threshold_Indicator.csv",
  "LP_Threshold_Effects.csv"
)
for (f in csv_files) {
  if (file.exists(file.path(REGIME_CONFIG$output_dir, f))) {
    cat("  [OK]", f, "\n")
  }
}

cat("\nPNG Visualizations:\n")
png_files <- paste0(80:89, "_*.png")
actual_pngs <- list.files(REGIME_CONFIG$output_dir, pattern = "^8[0-9]_.*\\.png$")
for (f in actual_pngs) {
  cat("  [OK]", f, "\n")
}


################################################################################
# SECTION 10: SUMMARY AND KEY FINDINGS
################################################################################

cat("\n================================================================================\n")
cat("KEY FINDINGS SUMMARY\n")
cat("================================================================================\n\n")

cat("1. REGIME-SPLIT PCA:\n")
if (!is.null(pca_pre_IT) && !is.null(pca_post_IT)) {
  cat(sprintf("   Loading correlation (Pre vs Post-IT): %.3f\n", loading_corr))
  if (exists("bootstrap_results") && !is.null(bootstrap_results)) {
    n_sig <- sum(bootstrap_results$significant_05, na.rm = TRUE)
    cat(sprintf("   Significant loading changes: %d/%d variables\n",
                n_sig, nrow(bootstrap_results)))
  }
}

cat("\n2. REGIME-DEPENDENT LP EFFECTS (at h=12):\n")
if (exists("it_12") && nrow(it_12) > 0) {
  cat(sprintf("   IT regime differential: %+.2f pp (p=%.3f)\n",
              it_12$coef_interact, it_12$p_interact))
}
if (exists("covid_12") && nrow(covid_12) > 0) {
  cat(sprintf("   COVID differential: %+.2f pp (p=%.3f)\n",
              covid_12$coef_interact, covid_12$p_interact))
}

cat("\n3. THRESHOLD EFFECTS (at h=12):\n")
if (exists("credit_12") && nrow(credit_12) > 0) {
  cat(sprintf("   Crisis effect: %+.2f pp\n", credit_12$coef_high))
  cat(sprintf("   Normal effect: %+.2f pp\n", credit_12$coef_low))
  cat(sprintf("   Equality test p: %.3f\n", credit_12$diff_p))
}

cat("\n4. INTERPRETATION:\n")
cat("   - FCI uses 12 variables for full-sample consistency\n")
cat("   - IT adoption (May 2011) marks potential structural change in transmission\n")
cat("   - COVID period (Mar 2020 - Dec 2021) treated separately due to unusual dynamics\n")
cat("   - Threshold analysis identifies potential non-linearities in FCI effects\n")

cat("\n================================================================================\n")
cat("REGIME AND TVP ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
cat("Plots: 80-89\n")
cat("Output:", REGIME_CONFIG$output_dir, "\n\n")


################################################################################
# END OF SCRIPT
################################################################################
