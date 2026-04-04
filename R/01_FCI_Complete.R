################################################################################
# FINANCIAL CONDITIONS INDEX (FCI) FOR PARAGUAY - CONSOLIDATED VERSION
################################################################################
#
# Project:      Financial Conditions Index - Complete Analysis
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Consolidated script that calculates the FCI at multiple
#               aggregation levels using four econometric methodologies:
#               1. Z-Score (equal-weighted average)
#               2. Principal Component Analysis (PCA)
#               3. Vector Autoregression (VAR)
#               4. Dynamic Factor Model (DFM)
#
# OUTPUT STRUCTURE:
#   - Level 1: Comprehensive FCI (all variables)
#   - Level 2: Endogenous vs Exogenous decomposition
#   - Level 3: By economic channel (Rates, Banking, External)
#   - Orthogonalized indices (domestic net of external spillovers)
#
# Author:       Departamento de Analisis Macroeconomico
# Created:      2025
# Last Updated: 2025-01-16 (Consolidated from multiple scripts)
#
# References:
#   - Brave & Butters (2011) - Chicago Fed NFCI
#   - Hatzius et al. (2010) - Goldman Sachs FCI
#   - Arrigoni, Bobasu & Venditti (2022) - FCIs for Latin America
#
################################################################################


################################################################################
# 1. SETUP AND CONFIGURATION
################################################################################

# Clear workspace (optional)
# rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(readxl)      # Excel file reading
  library(dplyr)       # Data manipulation
  library(tidyr)       # Data reshaping
  library(zoo)         # Time series functions
  library(FactoMineR)  # PCA
  library(vars)        # VAR models
  library(MARSS)       # Dynamic Factor Models
  library(ggplot2)     # Visualization
  library(gridExtra)   # Multiple plots
})

# Configuration parameters
CONFIG <- list(
  # Data parameters
  data_file = "../data/FCI_data_1.xlsx",
  data_sheet = "Main_variables",

  # Sample truncation: paper reports January 1996 - December 2025 (360 obs)
  sample_end = as.Date("2025-12-01"),

  # Analysis parameters
  rolling_window = 60,        # 5 years for standardization
  var_max_lags = 12,          # Max lags for VAR selection
  dfm_max_iter = 1000,        # DFM convergence iterations

  # Output parameters
  output_dir = "../output",
  export_results = TRUE,
  verbose = TRUE
)

# Set seed for reproducibility
set.seed(20250116)

# Print session info
if (CONFIG$verbose) {
  cat("\n################################################################################\n")
  cat("FINANCIAL CONDITIONS INDEX - COMPLETE ANALYSIS\n")
  cat("################################################################################\n\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("R version:", R.version.string, "\n")
  cat("Working directory:", getwd(), "\n\n")
}


################################################################################
# 2. VARIABLE DEFINITIONS
################################################################################

# Sign conventions:
#   +1 = Higher value indicates TIGHTENING of financial conditions
#   -1 = Higher value indicates LOOSENING of financial conditions

# LEVEL 3: By Economic Channel
VARIABLES <- list(

  # Interest Rate Channel
  rates = list(
    vars = c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM"),
    signs = c(+1, +1, +1),
    description = "Interest rate channel: cost of credit"
  ),

  # Banking/Credit Channel
  banking = list(
    vars = c("Crecimiento_creditos", "Ratio_Cred_Depo", "Morosidad",
             "Rentabilidad", "Liquidez"),
    signs = c(-1, -1, +1, -1, -1),
    description = "Banking channel: credit supply and bank health"
  ),

  # External Channel (4 variables)
  external = list(
    vars = c("TCN", "Commodities", "FFER", "VIX"),
    signs = c(+1, -1, +1, +1),
    description = "External channel: global conditions and FX"
  )
)

# LEVEL 2: Endogenous vs Exogenous
LEVEL2 <- list(

  endogenous = list(
    vars = c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM",
             "Crecimiento_creditos", "Ratio_Cred_Depo", "Morosidad",
             "Rentabilidad", "Liquidez", "TCN"),
    signs = c(+1, +1, +1, -1, -1, +1, -1, -1, +1),
    description = "Domestic conditions (policy-influenced)"
  ),

  # Exogenous (3 variables, exogenous to Paraguay)
  exogenous = list(
    vars = c("FFER", "VIX", "Commodities"),
    signs = c(+1, +1, -1),
    description = "External conditions (exogenous to Paraguay)"
  )
)

# Print variable summary
if (CONFIG$verbose) {
  cat("Variable Classification:\n")
  cat("  Level 3 - Rates:      ", paste(VARIABLES$rates$vars, collapse = ", "), "\n")
  cat("  Level 3 - Banking:    ", paste(VARIABLES$banking$vars, collapse = ", "), "\n")
  cat("  Level 3 - External:   ", paste(VARIABLES$external$vars, collapse = ", "), "\n")
  cat("  Level 2 - Endogenous: ", length(LEVEL2$endogenous$vars), " vars\n")
  cat("  Level 2 - Exogenous:  ", length(LEVEL2$exogenous$vars), " vars\n\n")
}


################################################################################
# 3. UTILITY FUNCTIONS
################################################################################

#' Calculate Year-over-Year growth rate
#' @param x Numeric vector (stock variable)
#' @return YoY growth rate in percentage
calculate_yoy_growth <- function(x) {
  (x / lag(x, 12) - 1) * 100
}

#' Apply sign adjustments to variables
#' @param data Data frame with variables
#' @param vars Vector of variable names
#' @param signs Vector of signs (+1 or -1)
#' @return Data frame with sign-adjusted variables
apply_signs <- function(data, vars, signs) {
  data_adj <- data
  for (i in seq_along(vars)) {
    if (vars[i] %in% names(data_adj)) {
      data_adj[[vars[i]]] <- data_adj[[vars[i]]] * signs[i]
    }
  }
  return(data_adj)
}

#' Rolling window standardization (Z-score)
#' @param data Data frame with numeric variables
#' @param window Rolling window size
#' @return Standardized data frame
rolling_standardize <- function(data, window = CONFIG$rolling_window) {
  data %>%
    mutate(across(where(is.numeric) & !matches("fecha"), ~ {
      rm <- zoo::rollapply(., width = window, FUN = mean, na.rm = TRUE,
                           fill = NA, align = "right", partial = TRUE)
      rsd <- zoo::rollapply(., width = window, FUN = sd, na.rm = TRUE,
                            fill = NA, align = "right", partial = TRUE)
      rsd <- pmax(rsd, 1e-10)
      (. - rm) / rsd
    }))
}


################################################################################
# 4. DATA LOADING AND PREPARATION
################################################################################

load_and_prepare_data <- function() {

  if (CONFIG$verbose) {
    cat("================================================================================\n")
    cat("LOADING AND PREPARING DATA\n")
    cat("================================================================================\n\n")
  }

  # Check file exists
  if (!file.exists(CONFIG$data_file)) {
    stop("Data file not found: ", CONFIG$data_file)
  }

  # Load raw data
  datos_raw <- read_excel(CONFIG$data_file, sheet = CONFIG$data_sheet)

  # Detect date column
  fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

  if (is.na(fecha_col)) {
    stop("No date column found in data")
  }

  # Prepare base data
  datos <- datos_raw %>%
    rename(fecha = !!sym(fecha_col)) %>%
    mutate(fecha = as.Date(fecha)) %>%
    arrange(fecha)

  # Truncate sample if CONFIG$sample_end is set (paper: Jan 1996 - Dec 2025)
  if (!is.null(CONFIG$sample_end)) {
    n_before <- nrow(datos)
    datos <- datos %>% filter(fecha <= CONFIG$sample_end)
    if (CONFIG$verbose) {
      cat(sprintf("Sample truncated to %s: %d -> %d observations\n",
                  format(CONFIG$sample_end, "%B %Y"), n_before, nrow(datos)))
    }
  }

  # ==========================================================================
  # CALCULATE CREDIT GROWTH FROM STOCK DATA
  # ==========================================================================
  # The database contains stock variables for credits.
  # We calculate YoY growth rate for use in the FCI.

  if ("Creditos_Sector_privado_totales" %in% names(datos)) {
    if (CONFIG$verbose) {
      cat("Calculating credit growth (YoY) from stock data...\n")
    }

    datos <- datos %>%
      arrange(fecha) %>%
      mutate(
        Crecimiento_creditos = calculate_yoy_growth(Creditos_Sector_privado_totales)
      )

    if (CONFIG$verbose) {
      credit_stats <- datos %>%
        filter(!is.na(Crecimiento_creditos)) %>%
        summarise(
          mean = mean(Crecimiento_creditos, na.rm = TRUE),
          sd = sd(Crecimiento_creditos, na.rm = TRUE),
          min = min(Crecimiento_creditos, na.rm = TRUE),
          max = max(Crecimiento_creditos, na.rm = TRUE)
        )

      cat(sprintf("  Credit growth (YoY): Mean: %.1f%%, SD: %.1f%%, Range: [%.1f%%, %.1f%%]\n",
                  credit_stats$mean, credit_stats$sd,
                  credit_stats$min, credit_stats$max))
    }
  }

  # Report data summary
  if (CONFIG$verbose) {
    cat("\nData loaded successfully:\n")
    cat("  Period:      ", format(min(datos$fecha)), "to", format(max(datos$fecha)), "\n")
    cat("  Observations:", nrow(datos), "\n\n")
  }

  # ==========================================================================
  # SAMPLE VERIFICATION (Phase 1.1: Critical Fix)
  # ==========================================================================
  # Export sample verification for transparency and report consistency
  sample_verification <- data.frame(
    execution_date = Sys.Date(),
    execution_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    sample_start = min(datos$fecha),
    sample_end = max(datos$fecha),
    n_observations = nrow(datos),
    n_variables = ncol(datos) - 1,  # Exclude fecha
    r_version = R.version.string
  )

  if (CONFIG$export_results) {
    write.csv(sample_verification,
              file.path(CONFIG$output_dir, "FCI_Sample_Verification.csv"),
              row.names = FALSE)

    if (CONFIG$verbose) {
      cat("================================================================================\n")
      cat("SAMPLE VERIFICATION EXPORTED\n")
      cat("================================================================================\n")
      cat("  Execution date:  ", format(Sys.Date(), "%Y-%m-%d"), "\n")
      cat("  Sample start:    ", format(min(datos$fecha), "%Y-%m"), "\n")
      cat("  Sample end:      ", format(max(datos$fecha), "%Y-%m"), "\n")
      cat("  N observations:  ", nrow(datos), "\n")
      cat("  Saved: FCI_Sample_Verification.csv\n\n")
    }
  }

  return(datos)
}


################################################################################
# 5. FCI CALCULATION FUNCTIONS
################################################################################

#' Calculate FCI using Z-Score method (equal-weighted average)
#' @param data Standardized data frame
#' @param vars Vector of variable names
#' @return Numeric vector with FCI values
calculate_zscore <- function(data, vars) {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) == 0) return(rep(NA, nrow(data)))

  fci <- rowMeans(data[, vars_available, drop = FALSE], na.rm = TRUE)
  return(fci)
}

#' Calculate FCI using PCA (first principal component)
#' @param data Standardized data frame
#' @param vars Vector of variable names
#' @return Numeric vector with FCI values
calculate_pca <- function(data, vars) {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) < 2) return(rep(NA, nrow(data)))

  datos_subset <- data[, vars_available, drop = FALSE]

  # Impute missing values with median for PCA
  datos_clean <- datos_subset
  for (col in names(datos_clean)) {
    median_val <- median(datos_clean[[col]], na.rm = TRUE)
    datos_clean[[col]][is.na(datos_clean[[col]])] <- median_val
  }

  tryCatch({
    pca_model <- PCA(datos_clean, ncp = 1, graph = FALSE)
    fci <- pca_model$ind$coord[, 1]

    # Sign correction: ensure positive correlation with average
    if (cor(fci, rowMeans(datos_clean), use = "complete.obs") < 0) {
      fci <- -fci
    }

    return(fci)
  }, error = function(e) {
    warning("PCA failed: ", e$message)
    return(rep(NA, nrow(data)))
  })
}

#' Calculate FCI using VAR model
#' @param data Standardized data frame
#' @param vars Vector of variable names
#' @return Numeric vector with FCI values
calculate_var <- function(data, vars) {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) < 2) return(rep(NA, nrow(data)))

  datos_subset <- data[, vars_available, drop = FALSE]

  # Impute missing values
  datos_clean <- datos_subset
  for (col in names(datos_clean)) {
    median_val <- median(datos_clean[[col]], na.rm = TRUE)
    datos_clean[[col]][is.na(datos_clean[[col]])] <- median_val
  }

  n_obs <- nrow(datos_clean)

  tryCatch({
    # Select optimal lag
    var_select <- VARselect(datos_clean,
                            lag.max = min(CONFIG$var_max_lags, floor(n_obs/10)),
                            type = "const")
    lag_opt <- max(1, min(var_select$selection["AIC(n)"], 6))

    # Estimate VAR
    var_model <- VAR(datos_clean, p = lag_opt, type = "const")

    # FCI = average of fitted values
    fci <- c(rep(NA, lag_opt), rowMeans(fitted(var_model)))

    return(fci)
  }, error = function(e) {
    warning("VAR failed: ", e$message)
    return(rep(NA, nrow(data)))
  })
}

#' Calculate FCI using Dynamic Factor Model
#' @param data Standardized data frame
#' @param vars Vector of variable names
#' @return Numeric vector with FCI values
calculate_dfm <- function(data, vars) {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) < 2) return(rep(NA, nrow(data)))

  datos_subset <- data[, vars_available, drop = FALSE]

  # Impute missing values
  datos_clean <- datos_subset
  for (col in names(datos_clean)) {
    median_val <- median(datos_clean[[col]], na.rm = TRUE)
    datos_clean[[col]][is.na(datos_clean[[col]])] <- median_val
  }

  n_vars <- ncol(datos_clean)

  tryCatch({
    Z_matrix <- matrix(as.list(paste0("z", 1:n_vars)), nrow = n_vars, ncol = 1)

    dfm_model <- MARSS(
      t(datos_clean),
      model = list(
        Z = Z_matrix,
        A = "zero",
        R = "diagonal and unequal",
        B = matrix(1),
        U = "zero",
        Q = matrix(1),
        x0 = matrix(0)
      ),
      silent = TRUE,
      control = list(maxit = CONFIG$dfm_max_iter)
    )

    fci <- as.numeric(dfm_model$states)

    # Sign correction
    if (cor(fci, rowMeans(datos_clean), use = "complete.obs") < 0) {
      fci <- -fci
    }

    return(fci)
  }, error = function(e) {
    warning("DFM failed: ", e$message)
    return(rep(NA, nrow(data)))
  })
}

#' Calculate FCI using all four methods
#' @param data Standardized data frame with fecha column
#' @param vars Vector of variable names
#' @param prefix Prefix for output column names
#' @return Data frame with FCI from all methods
calculate_fci_all_methods <- function(data, vars, prefix = "FCI") {

  vars_available <- intersect(vars, names(data))

  if (length(vars_available) < 2) {
    warning("Insufficient variables for ", prefix)
    return(NULL)
  }

  if (CONFIG$verbose) {
    cat("  Calculating ", prefix, " (", length(vars_available), " variables)...\n", sep = "")
  }

  results <- data.frame(fecha = data$fecha)

  # Z-Score
  results[[paste0(prefix, "_ZSCORE")]] <- calculate_zscore(data, vars)

  # PCA
  results[[paste0(prefix, "_PCA")]] <- calculate_pca(data, vars)

  # VAR
  results[[paste0(prefix, "_VAR")]] <- calculate_var(data, vars)

  # DFM
  results[[paste0(prefix, "_DFM")]] <- calculate_dfm(data, vars)

  # Normalize all methods
  methods <- c("ZSCORE", "PCA", "VAR", "DFM")
  for (m in methods) {
    col_name <- paste0(prefix, "_", m)
    if (col_name %in% names(results) && !all(is.na(results[[col_name]]))) {
      results[[paste0(col_name, "_norm")]] <- scale(results[[col_name]])[,1]
    }
  }

  # Average of normalized indices
  norm_cols <- paste0(prefix, "_", methods, "_norm")
  norm_cols <- norm_cols[norm_cols %in% names(results)]

  if (length(norm_cols) > 0) {
    results[[paste0(prefix, "_AVG")]] <- rowMeans(
      results[, norm_cols, drop = FALSE], na.rm = TRUE
    )
  }

  return(results)
}


################################################################################
# 6. ORTHOGONALIZATION
################################################################################

#' Orthogonalize domestic FCI with respect to external FCI
#' @param endo_fci Numeric vector of endogenous FCI
#' @param exo_fci Numeric vector of exogenous FCI
#' @param fecha Date vector
#' @return List with orthogonalized FCI and regression results
orthogonalize_fci <- function(endo_fci, exo_fci, fecha) {

  # Create data frame
  data <- data.frame(
    fecha = fecha,
    endo = endo_fci,
    exo = exo_fci
  ) %>% na.omit()

  if (nrow(data) < 30) {
    warning("Insufficient observations for orthogonalization")
    return(NULL)
  }

  # Run regression: endo = beta0 + beta1 * exo + epsilon
  model <- lm(endo ~ exo, data = data)

  # Orthogonalized = residuals (domestic net of external spillovers)
  ortho <- residuals(model)
  spillover <- fitted(model)

  if (CONFIG$verbose) {
    cat("\nOrthogonalization Results:\n")
    cat(sprintf("  FCI_Endo = %.3f + %.3f * FCI_Exo + epsilon\n",
                coef(model)[1], coef(model)[2]))
    cat(sprintf("  R-squared: %.3f (%.1f%% of domestic variation explained by external)\n",
                summary(model)$r.squared, summary(model)$r.squared * 100))
    cat(sprintf("  Spillover coefficient: %.3f (p = %.4f)\n",
                coef(model)[2], summary(model)$coefficients[2, 4]))
  }

  return(list(
    fecha = data$fecha,
    orthogonalized = ortho,
    spillover = spillover,
    model = model,
    r_squared = summary(model)$r.squared,
    beta = coef(model)[2]
  ))
}


################################################################################
# 7. VISUALIZATION FUNCTIONS
################################################################################

# Color palettes (colorblind-safe: avoids red-green combinations)
COLORS <- list(
  methods = c("ZSCORE" = "#0072B2", "PCA" = "#E69F00",
              "VAR" = "#009E73", "DFM" = "#CC79A7"),
  level2 = c("Endogenous" = "#D55E00", "Exogenous" = "#0072B2",
             "Orthogonalized" = "#CC79A7"),
  level3 = c("Rates" = "#D55E00", "Banking" = "#0072B2", "External" = "#009E73"),
  comprehensive = c("Comprehensive" = "#2C3E50")
)

#' Create time series comparison plot
create_comparison_plot <- function(data_long, title, subtitle, colors) {
  ggplot(data_long, aes(x = fecha, y = Value, color = Series)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    scale_color_manual(values = colors) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    theme_minimal(base_size = 12) +
    labs(title = title, subtitle = subtitle, x = NULL,
         y = "FCI (standard deviations)", color = NULL) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#' Create contribution area plot
create_contribution_plot <- function(data, title) {
  ggplot(data, aes(x = fecha, y = Contribution, fill = Channel)) +
    geom_area(alpha = 0.7, position = "stack") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    scale_fill_manual(values = COLORS$level3) +
    theme_minimal(base_size = 12) +
    labs(title = title, x = NULL, y = "Contribution to FCI", fill = "Channel") +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "bottom"
    )
}


################################################################################
# 8. MAIN EXECUTION
################################################################################

if (CONFIG$verbose) {
  cat("================================================================================\n")
  cat("STARTING FCI CALCULATION\n")
  cat("================================================================================\n\n")
}

# Create output directory
if (!dir.exists(CONFIG$output_dir)) {
  dir.create(CONFIG$output_dir, recursive = TRUE)
}

# Load and prepare data
datos <- load_and_prepare_data()

# Get all unique variables (12 variables)
all_vars_core <- unique(c(
  VARIABLES$rates$vars,
  VARIABLES$banking$vars,
  VARIABLES$external$vars
))

# For backward compatibility
all_vars <- all_vars_core

# =============================================================================
# PREPARE SIGN-ADJUSTED AND STANDARDIZED DATA
# =============================================================================

if (CONFIG$verbose) {
  cat("================================================================================\n")
  cat("PREPARING STANDARDIZED DATA\n")
  cat("================================================================================\n\n")
}

# Apply signs to all 12 variables
all_signs_core <- c(VARIABLES$rates$signs, VARIABLES$banking$signs, VARIABLES$external$signs)
datos_signed_core <- apply_signs(datos, all_vars_core, all_signs_core)

# Standardize data
datos_std <- datos_signed_core %>%
  dplyr::select(fecha, any_of(all_vars_core)) %>%
  rolling_standardize(CONFIG$rolling_window)

if (CONFIG$verbose) {
  cat("Data preparation complete.\n")
  cat("  Observations:", nrow(datos_std), "\n\n")
}

# =============================================================================
# LEVEL 1: COMPREHENSIVE FCI
# =============================================================================

if (CONFIG$verbose) {
  cat("================================================================================\n")
  cat("LEVEL 1: COMPREHENSIVE FCI\n")
  cat("================================================================================\n\n")
}

# Comprehensive FCI (12 variables)
fci_comprehensive <- calculate_fci_all_methods(datos_std, all_vars_core, "FCI_COMP")

# =============================================================================
# LEVEL 2: ENDOGENOUS VS EXOGENOUS
# =============================================================================

if (CONFIG$verbose) {
  cat("\n================================================================================\n")
  cat("LEVEL 2: ENDOGENOUS VS EXOGENOUS DECOMPOSITION\n")
  cat("================================================================================\n\n")
}

# Prepare endogenous data
datos_endo <- datos %>%
  dplyr::select(fecha, any_of(LEVEL2$endogenous$vars))
datos_endo <- apply_signs(datos_endo, LEVEL2$endogenous$vars, LEVEL2$endogenous$signs)
datos_endo <- rolling_standardize(datos_endo)

fci_endogenous <- calculate_fci_all_methods(datos_endo, LEVEL2$endogenous$vars, "FCI_ENDO")

# Prepare exogenous data
datos_exo <- datos %>%
  dplyr::select(fecha, any_of(LEVEL2$exogenous$vars))
datos_exo <- apply_signs(datos_exo, LEVEL2$exogenous$vars, LEVEL2$exogenous$signs)
datos_exo <- rolling_standardize(datos_exo)

fci_exogenous <- calculate_fci_all_methods(datos_exo, LEVEL2$exogenous$vars, "FCI_EXO")

# Orthogonalization
ortho_result <- NULL
if (!is.null(fci_endogenous) && !is.null(fci_exogenous)) {
  if ("FCI_ENDO_AVG" %in% names(fci_endogenous) && "FCI_EXO_AVG" %in% names(fci_exogenous)) {

    if (CONFIG$verbose) {
      cat("\nPerforming orthogonalization...\n")
    }

    ortho_result <- orthogonalize_fci(
      fci_endogenous$FCI_ENDO_AVG,
      fci_exogenous$FCI_EXO_AVG,
      fci_endogenous$fecha
    )
  }
}

# =============================================================================
# LEVEL 3: BY ECONOMIC CHANNEL
# =============================================================================

if (CONFIG$verbose) {
  cat("\n================================================================================\n")
  cat("LEVEL 3: BY ECONOMIC CHANNEL\n")
  cat("================================================================================\n\n")
}

# Rates channel
datos_rates <- datos %>%
  dplyr::select(fecha, any_of(VARIABLES$rates$vars))
datos_rates <- apply_signs(datos_rates, VARIABLES$rates$vars, VARIABLES$rates$signs)
datos_rates <- rolling_standardize(datos_rates)

fci_rates <- calculate_fci_all_methods(datos_rates, VARIABLES$rates$vars, "FCI_RATES")

# Banking channel
datos_banking <- datos %>%
  dplyr::select(fecha, any_of(VARIABLES$banking$vars))
datos_banking <- apply_signs(datos_banking, VARIABLES$banking$vars, VARIABLES$banking$signs)
datos_banking <- rolling_standardize(datos_banking)

fci_banking <- calculate_fci_all_methods(datos_banking, VARIABLES$banking$vars, "FCI_BANKING")

# External channel
datos_external <- datos %>%
  dplyr::select(fecha, any_of(VARIABLES$external$vars))
datos_external <- apply_signs(datos_external, VARIABLES$external$vars, VARIABLES$external$signs)
datos_external <- rolling_standardize(datos_external)

fci_external <- calculate_fci_all_methods(datos_external, VARIABLES$external$vars, "FCI_EXTERNAL")


# =============================================================================
# LEVEL 4: PURIFIED FCI (Phase 2.1: Orthogonalized to Monetary Policy)
# =============================================================================

if (CONFIG$verbose) {
  cat("\n================================================================================\n")
  cat("LEVEL 4: PURIFIED FCI (Orthogonalized to Monetary Policy)\n")
  cat("================================================================================\n\n")
}

# METHOD 1: FCI_exTPM - Simple exclusion of TPM
# This creates an FCI that excludes the policy rate entirely
vars_exTPM <- setdiff(all_vars_core, "TPM")
signs_exTPM <- all_signs_core[all_vars_core != "TPM"]

datos_exTPM <- datos %>%
  dplyr::select(fecha, any_of(vars_exTPM))
datos_exTPM <- apply_signs(datos_exTPM, vars_exTPM, signs_exTPM)
datos_exTPM <- rolling_standardize(datos_exTPM)

fci_exTPM <- calculate_fci_all_methods(datos_exTPM, vars_exTPM, "FCI_exTPM")

if (CONFIG$verbose) {
  cat("  FCI_exTPM: Excludes TPM (", length(vars_exTPM), " variables)\n", sep = "")
}

# METHOD 2: FCI_PURIFIED - Residuals after regressing on TPM + lags
# For each variable, regress on TPM(t), TPM(t-1), ..., TPM(t-k) and use residuals
# This preserves financial variation not driven by monetary policy

if (CONFIG$verbose) {
  cat("  Calculating FCI_PURIFIED (residuals after TPM regression)...\n")
}

purify_variable <- function(y, tpm, n_lags = 3) {
  # Create lag matrix
  df <- data.frame(y = y, tpm = tpm)
  for (i in 1:n_lags) {
    df[[paste0("tpm_lag", i)]] <- dplyr::lag(tpm, i)
  }

  df_complete <- na.omit(df)
  if (nrow(df_complete) < 50) return(y)  # Return original if insufficient data

  # Regress on TPM and lags
  formula_str <- paste("y ~ tpm +", paste(paste0("tpm_lag", 1:n_lags), collapse = " + "))
  model <- tryCatch({
    lm(as.formula(formula_str), data = df_complete)
  }, error = function(e) NULL)

  if (is.null(model)) return(y)

  # Get residuals
  residuals_full <- rep(NA, length(y))
  residuals_full[as.numeric(rownames(df_complete))] <- residuals(model)

  return(residuals_full)
}

# Create purified versions of all variables
datos_purified <- datos %>%
  dplyr::select(fecha, any_of(all_vars_core))

# Apply signs first
datos_purified_signed <- apply_signs(datos_purified, all_vars_core, all_signs_core)

# Purify each variable (except TPM itself)
tpm_series <- datos_purified_signed$TPM
vars_to_purify <- setdiff(all_vars_core, "TPM")

for (v in vars_to_purify) {
  if (v %in% names(datos_purified_signed) && !all(is.na(datos_purified_signed[[v]]))) {
    datos_purified_signed[[v]] <- purify_variable(datos_purified_signed[[v]], tpm_series)
  }
}

# Now standardize the purified data
datos_purified_std <- datos_purified_signed %>%
  dplyr::select(fecha, any_of(vars_to_purify)) %>%
  rolling_standardize(CONFIG$rolling_window)

# Calculate FCI_PURIFIED using all 4 methods
fci_purified <- calculate_fci_all_methods(datos_purified_std, vars_to_purify, "FCI_PURIFIED")

if (CONFIG$verbose) {
  cat("  FCI_PURIFIED: Residuals after TPM + 3 lags regression\n")
}

# Compare with baseline
if (!is.null(fci_purified) && "FCI_PURIFIED_AVG" %in% names(fci_purified)) {
  purified_comparison <- data.frame(
    fecha = fci_comprehensive$fecha,
    FCI_COMP_AVG = fci_comprehensive$FCI_COMP_AVG
  ) %>%
    left_join(fci_purified %>% dplyr::select(fecha, FCI_PURIFIED_AVG), by = "fecha") %>%
    left_join(fci_exTPM %>% dplyr::select(fecha, FCI_exTPM_AVG), by = "fecha") %>%
    na.omit()

  if (nrow(purified_comparison) > 30) {
    corr_purified <- cor(purified_comparison$FCI_COMP_AVG, purified_comparison$FCI_PURIFIED_AVG)
    corr_exTPM <- cor(purified_comparison$FCI_COMP_AVG, purified_comparison$FCI_exTPM_AVG)

    if (CONFIG$verbose) {
      cat("\nPurified FCI Comparison:\n")
      cat(sprintf("  Correlation FCI_COMP vs FCI_exTPM:    %.3f\n", corr_exTPM))
      cat(sprintf("  Correlation FCI_COMP vs FCI_PURIFIED: %.3f\n", corr_purified))

      if (corr_purified > 0.85) {
        cat("  >>> HIGH correlation: Monetary policy explains little FCI variation\n")
      } else if (corr_purified > 0.70) {
        cat("  >>> MODERATE correlation: Some FCI variation driven by policy\n")
      } else {
        cat("  >>> LOW correlation: Substantial FCI variation from policy\n")
      }
    }
  }
}


# =============================================================================
# LEVEL 5: ROBUSTNESS INDICES (Endogeneity Corrections)
# =============================================================================
# FCI_exCredit: Excludes Crecimiento_creditos (for credit-LHS regressions)
# FCI_ENDO_exCredit: Excludes Crecimiento_creditos from domestic component

if (CONFIG$verbose) {
  cat("\n================================================================================\n")
  cat("LEVEL 5: ROBUSTNESS INDICES (Endogeneity Corrections)\n")
  cat("================================================================================\n\n")
}

# --- FCI_exCredit: Comprehensive FCI excluding credit growth (11 vars) ---
vars_exCredit <- setdiff(all_vars_core, "Crecimiento_creditos")
signs_exCredit <- all_signs_core[all_vars_core != "Crecimiento_creditos"]

datos_exCredit <- datos %>%
  dplyr::select(fecha, any_of(vars_exCredit))
datos_exCredit <- apply_signs(datos_exCredit, vars_exCredit, signs_exCredit)
datos_exCredit <- rolling_standardize(datos_exCredit)

fci_exCredit <- calculate_fci_all_methods(datos_exCredit, vars_exCredit, "FCI_exCredit")

if (CONFIG$verbose) {
  cat("  FCI_exCredit: Excludes Crecimiento_creditos (", length(vars_exCredit), " variables)\n", sep = "")
}

# --- FCI_ENDO_exCredit: Endogenous FCI excluding credit growth (8 vars) ---
vars_ENDO_exCredit <- setdiff(LEVEL2$endogenous$vars, "Crecimiento_creditos")
signs_ENDO_exCredit <- LEVEL2$endogenous$signs[LEVEL2$endogenous$vars != "Crecimiento_creditos"]

datos_ENDO_exCredit <- datos %>%
  dplyr::select(fecha, any_of(vars_ENDO_exCredit))
datos_ENDO_exCredit <- apply_signs(datos_ENDO_exCredit, vars_ENDO_exCredit, signs_ENDO_exCredit)
datos_ENDO_exCredit <- rolling_standardize(datos_ENDO_exCredit)

fci_ENDO_exCredit <- calculate_fci_all_methods(datos_ENDO_exCredit, vars_ENDO_exCredit, "FCI_ENDO_exCredit")

if (CONFIG$verbose) {
  cat("  FCI_ENDO_exCredit: Excludes Crecimiento_creditos from ENDO (", length(vars_ENDO_exCredit), " variables)\n", sep = "")
}

# --- FCI_exNPL: Comprehensive FCI excluding Morosidad (11 vars) ---
vars_exNPL <- setdiff(all_vars_core, "Morosidad")
signs_exNPL <- all_signs_core[all_vars_core != "Morosidad"]

datos_exNPL <- datos %>%
  dplyr::select(fecha, any_of(vars_exNPL))
datos_exNPL <- apply_signs(datos_exNPL, vars_exNPL, signs_exNPL)
datos_exNPL <- rolling_standardize(datos_exNPL)

fci_exNPL <- calculate_fci_all_methods(datos_exNPL, vars_exNPL, "FCI_exNPL")

if (CONFIG$verbose) {
  cat("  FCI_exNPL: Excludes Morosidad (", length(vars_exNPL), " variables)\n\n", sep = "")
}

# --- FCI_ENDO_exCredit_exTCN: Endogenous FCI excluding credit growth AND TCN (7 vars) ---
vars_ENDO_exCredit_exTCN <- setdiff(vars_ENDO_exCredit, "TCN")
signs_ENDO_exCredit_exTCN <- signs_ENDO_exCredit[vars_ENDO_exCredit != "TCN"]

datos_ENDO_exCredit_exTCN <- datos %>%
  dplyr::select(fecha, any_of(vars_ENDO_exCredit_exTCN))
datos_ENDO_exCredit_exTCN <- apply_signs(datos_ENDO_exCredit_exTCN, vars_ENDO_exCredit_exTCN, signs_ENDO_exCredit_exTCN)
datos_ENDO_exCredit_exTCN <- rolling_standardize(datos_ENDO_exCredit_exTCN)

fci_ENDO_exCredit_exTCN <- calculate_fci_all_methods(datos_ENDO_exCredit_exTCN, vars_ENDO_exCredit_exTCN, "FCI_ENDO_exCredit_exTCN")

if (CONFIG$verbose) {
  cat("  FCI_ENDO_exCredit_exTCN: Excludes Crecimiento_creditos + TCN from ENDO (", length(vars_ENDO_exCredit_exTCN), " variables)\n\n", sep = "")
}


# =============================================================================
# CONSOLIDATE RESULTS
# =============================================================================

if (CONFIG$verbose) {
  cat("\n================================================================================\n")
  cat("CONSOLIDATING RESULTS\n")
  cat("================================================================================\n\n")
}

# Merge all results
all_results <- data.frame(fecha = datos$fecha)

if (!is.null(fci_comprehensive)) {
  all_results <- all_results %>% left_join(fci_comprehensive, by = "fecha")
}

if (!is.null(fci_endogenous)) {
  all_results <- all_results %>% left_join(fci_endogenous, by = "fecha")
}

if (!is.null(fci_exogenous)) {
  all_results <- all_results %>% left_join(fci_exogenous, by = "fecha")
}

if (!is.null(fci_rates)) {
  all_results <- all_results %>% left_join(fci_rates, by = "fecha")
}

if (!is.null(fci_banking)) {
  all_results <- all_results %>% left_join(fci_banking, by = "fecha")
}

if (!is.null(fci_external)) {
  all_results <- all_results %>% left_join(fci_external, by = "fecha")
}

# Add orthogonalized series
if (!is.null(ortho_result)) {
  ortho_df <- data.frame(
    fecha = ortho_result$fecha,
    FCI_ENDO_ORTHO = ortho_result$orthogonalized,
    FCI_SPILLOVER = ortho_result$spillover
  )
  all_results <- all_results %>% left_join(ortho_df, by = "fecha")
}

# Add purified FCI series (Level 4)
if (exists("fci_exTPM") && !is.null(fci_exTPM)) {
  all_results <- all_results %>% left_join(fci_exTPM, by = "fecha")
}

if (exists("fci_purified") && !is.null(fci_purified)) {
  all_results <- all_results %>% left_join(fci_purified, by = "fecha")
}

# Add robustness indices (Level 5)
if (exists("fci_exCredit") && !is.null(fci_exCredit)) {
  all_results <- all_results %>% left_join(fci_exCredit, by = "fecha")
}

if (exists("fci_ENDO_exCredit") && !is.null(fci_ENDO_exCredit)) {
  all_results <- all_results %>% left_join(fci_ENDO_exCredit, by = "fecha")
}

if (exists("fci_exNPL") && !is.null(fci_exNPL)) {
  all_results <- all_results %>% left_join(fci_exNPL, by = "fecha")
}

if (exists("fci_ENDO_exCredit_exTCN") && !is.null(fci_ENDO_exCredit_exTCN)) {
  all_results <- all_results %>% left_join(fci_ENDO_exCredit_exTCN, by = "fecha")
}

cat("Total indices calculated:", ncol(all_results) - 1, "\n\n")

# =============================================================================
# GENERATE VISUALIZATIONS
# =============================================================================

if (CONFIG$verbose) {
  cat("================================================================================\n")
  cat("GENERATING VISUALIZATIONS\n")
  cat("================================================================================\n\n")
}

plot_counter <- 1

# Plot 1: Comprehensive FCI - Method Comparison
if (!is.null(fci_comprehensive)) {
  plot_data <- fci_comprehensive %>%
    dplyr::select(fecha, ends_with("_norm")) %>%
    pivot_longer(-fecha, names_to = "Method", values_to = "Value") %>%
    mutate(Method = gsub("FCI_COMP_|_norm", "", Method))

  p1 <- create_comparison_plot(
    plot_data %>% rename(Series = Method),
    "Comprehensive FCI: Method Comparison",
    "All four methodologies applied to complete variable set",
    COLORS$methods
  ) +
    geom_vline(xintercept = as.Date("2011-05-01"), linetype = "dashed",
               color = "darkred", linewidth = 0.5, alpha = 0.7) +
    annotate("text", x = as.Date("2011-05-01"), y = Inf,
             label = "IT Adoption", hjust = -0.1, vjust = 1.5,
             size = 3, color = "darkred")

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Methods_Comparison.png", plot_counter)),
         p1, width = 12, height = 6, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_Methods_Comparison.png\n", plot_counter))
  plot_counter <- plot_counter + 1
}

# Plot 2: Level 2 - Endogenous vs Exogenous
if (!is.null(fci_endogenous) && !is.null(fci_exogenous)) {
  plot_data <- all_results %>%
    dplyr::select(fecha, FCI_ENDO_AVG, FCI_EXO_AVG, any_of("FCI_ENDO_ORTHO")) %>%
    pivot_longer(-fecha, names_to = "Series", values_to = "Value") %>%
    mutate(Series = case_when(
      Series == "FCI_ENDO_AVG" ~ "Endogenous",
      Series == "FCI_EXO_AVG" ~ "Exogenous",
      Series == "FCI_ENDO_ORTHO" ~ "Orthogonalized",
      TRUE ~ Series
    ))

  p2 <- create_comparison_plot(
    plot_data,
    "Level 2 Decomposition: Endogenous vs Exogenous FCI",
    "Endogenous = domestic (policy-influenced) | Exogenous = external (given)",
    COLORS$level2
  )

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Level2_Decomposition.png", plot_counter)),
         p2, width = 12, height = 6, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_Level2_Decomposition.png\n", plot_counter))
  plot_counter <- plot_counter + 1
}

# Plot 3: Level 3 - By Economic Channel
if (!is.null(fci_rates) && !is.null(fci_banking) && !is.null(fci_external)) {
  plot_data <- all_results %>%
    dplyr::select(fecha, FCI_RATES_AVG, FCI_BANKING_AVG, FCI_EXTERNAL_AVG) %>%
    pivot_longer(-fecha, names_to = "Series", values_to = "Value") %>%
    mutate(Series = case_when(
      Series == "FCI_RATES_AVG" ~ "Rates",
      Series == "FCI_BANKING_AVG" ~ "Banking",
      Series == "FCI_EXTERNAL_AVG" ~ "External",
      TRUE ~ Series
    ))

  p3 <- create_comparison_plot(
    plot_data,
    "Level 3 Decomposition: FCI by Economic Channel",
    "Rates = interest rate channel | Banking = credit channel | External = global conditions",
    COLORS$level3
  )

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Level3_Channels.png", plot_counter)),
         p3, width = 12, height = 6, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_Level3_Channels.png\n", plot_counter))
  plot_counter <- plot_counter + 1
}

# Plot 4: Channel Contributions
if (!is.null(fci_rates) && !is.null(fci_banking) && !is.null(fci_external)) {
  contrib_data <- all_results %>%
    dplyr::select(fecha, FCI_RATES_AVG, FCI_BANKING_AVG, FCI_EXTERNAL_AVG) %>%
    na.omit() %>%
    pivot_longer(-fecha, names_to = "Channel", values_to = "Contribution") %>%
    mutate(Channel = case_when(
      Channel == "FCI_RATES_AVG" ~ "Rates",
      Channel == "FCI_BANKING_AVG" ~ "Banking",
      Channel == "FCI_EXTERNAL_AVG" ~ "External",
      TRUE ~ Channel
    ))

  p4 <- create_contribution_plot(
    contrib_data,
    "Contribution to Financial Conditions by Channel"
  )

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Channel_Contributions.png", plot_counter)),
         p4, width = 12, height = 6, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_Channel_Contributions.png\n", plot_counter))
  plot_counter <- plot_counter + 1
}

# Plot 5: Spillover Analysis (if orthogonalization was performed)
if (!is.null(ortho_result)) {
  plot_data <- data.frame(
    fecha = ortho_result$fecha,
    Original = fci_endogenous$FCI_ENDO_AVG[fci_endogenous$fecha %in% ortho_result$fecha],
    Orthogonalized = ortho_result$orthogonalized,
    Spillover = ortho_result$spillover
  ) %>%
    pivot_longer(-fecha, names_to = "Series", values_to = "Value")

  p5 <- create_comparison_plot(
    plot_data,
    "Orthogonalization: Domestic FCI Before and After",
    sprintf("R-squared = %.3f | %.1f%% of domestic variation explained by external factors",
            ortho_result$r_squared, ortho_result$r_squared * 100),
    c("Original" = "#E41A1C", "Orthogonalized" = "#377EB8", "Spillover" = "#4DAF4A")
  )

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Spillover_Analysis.png", plot_counter)),
         p5, width = 12, height = 6, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_Spillover_Analysis.png\n", plot_counter))
  plot_counter <- plot_counter + 1
}

# =============================================================================
# ADDITIONAL PLOTS: VARIABLE CONTRIBUTIONS AND UNCERTAINTY RANGE
# =============================================================================

# Plot 6: Average FCI with Uncertainty Range (Min/Max across 4 methods)
if (!is.null(fci_comprehensive)) {
  # Get the 4 normalized method columns
  norm_cols <- c("FCI_COMP_ZSCORE_norm", "FCI_COMP_PCA_norm",
                 "FCI_COMP_VAR_norm", "FCI_COMP_DFM_norm")
  norm_cols <- norm_cols[norm_cols %in% names(fci_comprehensive)]

  if (length(norm_cols) >= 2) {
    uncertainty_data <- fci_comprehensive %>%
      dplyr::select(fecha, all_of(norm_cols), FCI_COMP_AVG) %>%
      rowwise() %>%
      mutate(
        FCI_min = min(c_across(all_of(norm_cols)), na.rm = TRUE),
        FCI_max = max(c_across(all_of(norm_cols)), na.rm = TRUE)
      ) %>%
      ungroup()

    p6 <- ggplot(uncertainty_data, aes(x = fecha)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_ribbon(aes(ymin = FCI_min, ymax = FCI_max),
                  fill = "#3498DB", alpha = 0.3) +
      geom_line(aes(y = FCI_COMP_AVG), color = "#2C3E50", linewidth = 1.2) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      theme_minimal(base_size = 12) +
      labs(
        title = "Comprehensive FCI: Average with Uncertainty Range",
        subtitle = "Solid line = Average of 4 methods | Shaded area = Min-Max range across methods",
        x = NULL, y = "FCI (standard deviations)"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Average_with_Uncertainty.png", plot_counter)),
           p6, width = 12, height = 6, dpi = 300)
    cat("Saved:", sprintf("%02d_FCI_Average_with_Uncertainty.png\n", plot_counter))
    plot_counter <- plot_counter + 1
  }
}

# Plot 7: ALL Variable Contributions to Comprehensive FCI (single chart)
# For Z-Score, contribution of each variable = standardized value / n_variables
{
  # Calculate contributions from the standardized data
  var_contributions <- datos_std %>%
    dplyr::select(fecha, any_of(all_vars)) %>%
    na.omit()

  n_vars <- length(intersect(all_vars, names(var_contributions))) - 1  # exclude fecha

  # Calculate contribution (each variable contributes its value / n_vars to the FCI)
  contrib_long <- var_contributions %>%
    pivot_longer(-fecha, names_to = "Variable", values_to = "Contribution") %>%
    mutate(Contribution = Contribution / n_vars)

  # Assign channel to each variable
  channel_map <- data.frame(
    Variable = c(VARIABLES$rates$vars, VARIABLES$banking$vars, VARIABLES$external$vars),
    Channel = c(rep("Rates", length(VARIABLES$rates$vars)),
                rep("Banking", length(VARIABLES$banking$vars)),
                rep("External", length(VARIABLES$external$vars)))
  )

  contrib_long <- contrib_long %>%
    left_join(channel_map, by = "Variable") %>%
    filter(!is.na(Channel))

  # Calculate the total FCI (sum of all contributions)
  total_fci <- contrib_long %>%
    group_by(fecha) %>%
    summarise(FCI_Total = sum(Contribution, na.rm = TRUE), .groups = "drop")

  # Create color palette for variables (grouped by channel)
  var_colors <- c(
    # Rates (reds/oranges)
    "TPM" = "#E41A1C", "Spread_activas_pasivas" = "#FC8D62", "Spread_mercado_TPM" = "#FDAE61",
    # Banking (blues)
    "Crecimiento_creditos" = "#377EB8", "Ratio_Cred_Depo" = "#4292C6",
    "Morosidad" = "#6BAED6", "Rentabilidad" = "#9ECAE1", "Liquidez" = "#C6DBEF",
    # External (greens)
    "TCN" = "#4DAF4A", "Commodities" = "#ADDD8E",
    "FFER" = "#D9F0A3", "VIX" = "#F7FCB4"
  )

  # Order variables by channel for better visualization
  var_order <- c(VARIABLES$rates$vars, VARIABLES$banking$vars, VARIABLES$external$vars)
  contrib_long <- contrib_long %>%
    mutate(Variable = factor(Variable, levels = var_order))

  # Plot 7a: All variables stacked with FCI line
  p7a <- ggplot(contrib_long, aes(x = fecha, y = Contribution, fill = Variable)) +
    geom_col(position = "stack", width = 28) +
    geom_line(data = total_fci, aes(x = fecha, y = FCI_Total, fill = NULL),
              color = "black", linewidth = 1.2) +
    geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = var_colors) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    theme_minimal(base_size = 12) +
    labs(
      title = "Variable Contributions to Comprehensive FCI",
      subtitle = "Stacked bars = Individual variable contributions | Black line = Total FCI (sum)",
      x = NULL, y = "FCI / Contribution", fill = "Variable"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Variable_Contributions_All.png", plot_counter)),
         p7a, width = 14, height = 7, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_Variable_Contributions_All.png\n", plot_counter))
  plot_counter <- plot_counter + 1

  # Plot 7b: Variable contributions by channel (faceted)
  # Calculate channel subtotals for the FCI line in each facet
  channel_totals <- contrib_long %>%
    group_by(fecha, Channel) %>%
    summarise(Channel_FCI = sum(Contribution, na.rm = TRUE), .groups = "drop")

  p7b <- ggplot(contrib_long, aes(x = fecha, y = Contribution, fill = Variable)) +
    geom_col(position = "stack", width = 28) +
    geom_line(data = channel_totals, aes(x = fecha, y = Channel_FCI, fill = NULL),
              color = "black", linewidth = 1) +
    geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = var_colors) +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    facet_wrap(~Channel, ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 11) +
    labs(
      title = "Variable Contributions by Economic Channel",
      subtitle = "Stacked bars = Variable contributions | Black line = Channel sub-index",
      x = NULL, y = "Contribution to FCI", fill = "Variable"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Variable_Contributions_byChannel.png", plot_counter)),
         p7b, width = 14, height = 10, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_Variable_Contributions_byChannel.png\n", plot_counter))
  plot_counter <- plot_counter + 1
}

# Plot 8: Channel Contributions to Integral FCI (sum of variables per channel)
{
  # Aggregate variable contributions by channel
  channel_contrib <- contrib_long %>%
    group_by(fecha, Channel) %>%
    summarise(Contribution = sum(Contribution, na.rm = TRUE), .groups = "drop")

  # Order channels
  channel_contrib <- channel_contrib %>%
    mutate(Channel = factor(Channel, levels = c("Rates", "Banking", "External")))

  p8 <- ggplot(channel_contrib, aes(x = fecha, y = Contribution, fill = Channel)) +
    geom_col(position = "stack", width = 28) +
    geom_line(data = total_fci, aes(x = fecha, y = FCI_Total, fill = NULL),
              color = "black", linewidth = 1.2) +
    geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = COLORS$level3) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    theme_minimal(base_size = 12) +
    labs(
      title = "Channel Contributions to Comprehensive FCI",
      subtitle = "Stacked bars = Channel contributions (sum of variables) | Black line = Total FCI",
      x = NULL, y = "FCI / Contribution", fill = "Channel"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Channel_Contributions_Integral.png", plot_counter)),
         p8, width = 12, height = 6, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_Channel_Contributions_Integral.png\n", plot_counter))
  plot_counter <- plot_counter + 1
}

# Plot 9: Variable Contributions Heatmap (last 24 months)
{
  recent_contrib <- contrib_long %>%
    filter(fecha >= max(fecha) - (365*2)) %>%
    group_by(Variable, Channel) %>%
    summarise(
      Mean_Contrib = mean(Contribution, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Channel, desc(abs(Mean_Contrib)))

  # Order variables by channel and contribution magnitude
  var_order <- recent_contrib %>%
    arrange(Channel, Mean_Contrib) %>%
    pull(Variable)

  recent_monthly <- contrib_long %>%
    filter(fecha >= max(fecha) - (365*2)) %>%
    mutate(Variable = factor(Variable, levels = var_order))

  p9 <- ggplot(recent_monthly, aes(x = fecha, y = Variable, fill = Contribution)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Contribution") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
    theme_minimal(base_size = 11) +
    labs(
      title = "Variable Contributions Heatmap (Last 24 Months)",
      subtitle = "Red = Tightening contribution | Blue = Loosening contribution",
      x = NULL, y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9)
    )

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_Variable_Heatmap.png", plot_counter)),
         p9, width = 14, height = 8, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_Variable_Heatmap.png\n", plot_counter))
  plot_counter <- plot_counter + 1
}

# Plot 10: Method Comparison with All Methods Visible
if (!is.null(fci_comprehensive)) {
  methods_comparison <- fci_comprehensive %>%
    dplyr::select(fecha, ends_with("_norm"), FCI_COMP_AVG) %>%
    pivot_longer(-fecha, names_to = "Method", values_to = "Value") %>%
    mutate(
      Method = gsub("FCI_COMP_|_norm", "", Method),
      Method = ifelse(Method == "AVG", "AVERAGE", Method),
      LineType = ifelse(Method == "AVERAGE", "Average", "Method"),
      LineSize = ifelse(Method == "AVERAGE", 1.3, 0.7)
    )

  method_colors <- c(COLORS$methods, "AVERAGE" = "#2C3E50")

  p10 <- ggplot(methods_comparison, aes(x = fecha, y = Value, color = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(aes(linewidth = LineType, alpha = LineType)) +
    scale_color_manual(values = method_colors) +
    scale_linewidth_manual(values = c("Average" = 1.3, "Method" = 0.6), guide = "none") +
    scale_alpha_manual(values = c("Average" = 1, "Method" = 0.7), guide = "none") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    theme_minimal(base_size = 12) +
    labs(
      title = "FCI: All Four Methods Plus Average",
      subtitle = "Z-Score, PCA, VAR, DFM | Bold line = Average",
      x = NULL, y = "FCI (standard deviations)", color = "Method"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(CONFIG$output_dir, sprintf("%02d_FCI_All_Methods.png", plot_counter)),
         p10, width = 12, height = 6, dpi = 300)
  cat("Saved:", sprintf("%02d_FCI_All_Methods.png\n", plot_counter))
  plot_counter <- plot_counter + 1
}



# =============================================================================
# EXPORT RESULTS
# =============================================================================

if (CONFIG$export_results) {
  if (CONFIG$verbose) {
    cat("\n================================================================================\n")
    cat("EXPORTING RESULTS\n")
    cat("================================================================================\n\n")
  }

  # Export all indices
  write.csv(all_results,
            file.path(CONFIG$output_dir, "FCI_Complete_Results.csv"),
            row.names = FALSE)
  cat("Saved: FCI_Complete_Results.csv\n")

  # Export summary statistics
  avg_cols <- grep("_AVG$", names(all_results), value = TRUE)

  if (length(avg_cols) > 0) {
    summary_stats <- data.frame(
      Index = character(),
      Mean = numeric(),
      SD = numeric(),
      Min = numeric(),
      Max = numeric(),
      Obs = integer(),
      stringsAsFactors = FALSE
    )

    for (col in avg_cols) {
      vals <- all_results[[col]]
      summary_stats <- rbind(summary_stats, data.frame(
        Index = col,
        Mean = round(mean(vals, na.rm = TRUE), 3),
        SD = round(sd(vals, na.rm = TRUE), 3),
        Min = round(min(vals, na.rm = TRUE), 2),
        Max = round(max(vals, na.rm = TRUE), 2),
        Obs = sum(!is.na(vals))
      ))
    }

    write.csv(summary_stats,
              file.path(CONFIG$output_dir, "FCI_Summary_Statistics.csv"),
              row.names = FALSE)
    cat("Saved: FCI_Summary_Statistics.csv\n")

    # Print summary
    if (CONFIG$verbose) {
      cat("\nSummary Statistics:\n")
      print(summary_stats, row.names = FALSE)
    }
  }

  # Export correlation matrix
  if (length(avg_cols) > 1) {
    cor_matrix <- cor(all_results[, avg_cols], use = "pairwise.complete.obs")
    write.csv(as.data.frame(cor_matrix),
              file.path(CONFIG$output_dir, "FCI_Correlation_Matrix.csv"))
    cat("Saved: FCI_Correlation_Matrix.csv\n")
  }

  # Export cross-method correlation matrix (COMP: Z-Score, PCA, VAR, DFM)
  method_cols <- grep("^FCI_COMP_(ZSCORE|PCA|VAR|DFM)_norm$",
                      names(all_results), value = TRUE)
  if (length(method_cols) == 4) {
    method_data <- na.omit(all_results[, method_cols])
    method_corr <- cor(method_data)
    write.csv(as.data.frame(method_corr),
              file.path(CONFIG$output_dir, "FCI_Method_Correlations.csv"))
    cat("Saved: FCI_Method_Correlations.csv\n")
    cat("Cross-method pairwise correlations:\n")
    print(round(method_corr, 3))
  }
}


# =============================================================================
# STORE RESULTS IN GLOBAL ENVIRONMENT
# =============================================================================

# Extract PCA loadings for comprehensive FCI (useful for weight reporting)
pca_loadings_comp <- tryCatch({
  datos_pca <- datos_std[, intersect(all_vars_core, names(datos_std)), drop = FALSE]
  for (col in names(datos_pca)) {
    median_val <- median(datos_pca[[col]], na.rm = TRUE)
    datos_pca[[col]][is.na(datos_pca[[col]])] <- median_val
  }
  pca_model <- PCA(datos_pca, ncp = 1, graph = FALSE)
  loadings <- pca_model$var$coord[, 1]
  if (cor(pca_model$ind$coord[, 1], rowMeans(datos_pca), use = "complete.obs") < 0) {
    loadings <- -loadings
  }
  loadings
}, error = function(e) NULL)

# Make results available to other scripts
resultado_fci <<- list(
  all_indices = all_results,
  comprehensive = fci_comprehensive,
  endogenous = fci_endogenous,
  exogenous = fci_exogenous,
  orthogonalization = ortho_result,
  rates = fci_rates,
  banking = fci_banking,
  external = fci_external,
  pca_loadings = pca_loadings_comp,
  config = CONFIG,
  variables = list(
    level2 = LEVEL2,
    level3 = VARIABLES,
    core_vars = all_vars_core
  )
)


# =============================================================================
# FINAL SUMMARY
# =============================================================================

if (CONFIG$verbose) {
  cat("\n================================================================================\n")
  cat("FCI ANALYSIS COMPLETE\n")
  cat("================================================================================\n\n")

  cat("Results Summary:\n")
  cat("  Total indices calculated: ", ncol(all_results) - 1, "\n")
  cat("  Visualizations generated: ", plot_counter - 1, "\n")
  cat("  Output directory:         ", CONFIG$output_dir, "\n\n")

  cat("Available FCI series:\n")
  cat("  Level 1: FCI_COMP_AVG (12 variables, full sample 1996-2025)\n")
  cat("  Level 2: FCI_ENDO_AVG, FCI_EXO_AVG, FCI_ENDO_ORTHO\n")
  cat("  Level 3: FCI_RATES_AVG, FCI_BANKING_AVG, FCI_EXTERNAL_AVG\n\n")

  cat("  ROBUSTNESS (endogeneity corrections):\n")
  cat("    FCI_exCredit_AVG: Excludes credit growth (11 vars) - for credit-LHS regressions\n")
  cat("    FCI_ENDO_exCredit_AVG: Excludes credit growth from ENDO (8 vars)\n")
  cat("    FCI_exNPL_AVG: Excludes Morosidad (11 vars) - for NPL-LHS regressions\n\n")

  cat("Access results via: resultado_fci$all_indices\n\n")
}


################################################################################
# END OF SCRIPT
################################################################################
