################################################################################
# FCI POST-IT SUBSAMPLE ANALYSIS
################################################################################
#
# Project:      Financial Conditions Index - Post-IT Regime Subsample
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Re-estimates core FCI specifications on the post-Inflation
#               Targeting subsample (May 2011 - present, ~175 obs). Addresses
#               the concern that full-sample FCI averages two distinct factor
#               structures (3 PCA loading sign reversals between pre-IT and
#               post-IT periods).
#
#   PART A: FCI Construction Diagnostics
#     - Within-period PCA loadings vs full-sample
#     - 4-method FCI (Z-Score, PCA, VAR, DFM) on post-IT data
#     - Cross-method correlations and variance decomposition
#
#   PART B: Granger Causality
#     - FCI_exCredit -> real credit, FCI_COMP -> IMAEP, FCI_EXO -> FCI_ENDO
#
#   PART C: Credit Local Projections
#     - FCI_exCredit -> credit by type; FCI decomposition
#
#   PART D: Sectoral Output
#     - FCI_COMP -> IMAEP, IMAEP_SANB, FBKf, Consumo
#
#   PART E: Asymmetric Effects
#     - Tightening vs easing on credit
#
#   PART F: Full-Sample vs Post-IT Comparison Dashboard
#     - Overlay LP IRFs, z-tests for equality, PCA loading comparison
#
# References:
#   - Jorda (2005) - Local Projections
#   - Brave & Butters (2011) - Chicago Fed NFCI
#
################################################################################


################################################################################
# 1. SETUP AND CONFIGURATION
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(FactoMineR)
  library(vars)
  library(MARSS)
  library(ggplot2)
  library(gridExtra)
  library(sandwich)
  library(lmtest)
})

POSTIT_CONFIG <- list(
  data_file = "../data/FCI_data_1.xlsx",
  macro_sheet = "Datos_macro",
  IT_START = as.Date("2011-05-01"),
  rolling_window = 60,
  var_max_lags = 6,
  dfm_max_iter = 1000,
  max_horizon = 18,
  horizons_report = c(3, 6, 12, 18),
  n_lags = 2,
  confidence_level = 0.90,
  output_dir = "../output",
  verbose = TRUE
)

set.seed(20260210)

cat("\n################################################################################\n")
cat("FCI POST-IT SUBSAMPLE ANALYSIS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Post-IT start:", format(POSTIT_CONFIG$IT_START, "%Y-%m"), "\n")
cat("Max LP horizon:", POSTIT_CONFIG$max_horizon, "\n\n")


################################################################################
# 2. DATA PREPARATION
################################################################################

cat("================================================================================\n")
cat("DATA PREPARATION\n")
cat("================================================================================\n\n")

# Load FCI from main script (with CONFIG save/restore)
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  CONFIG_SAVE <- if (exists("CONFIG")) CONFIG else NULL
  source("01_FCI_Complete.R")
  if (!is.null(CONFIG_SAVE)) CONFIG <- CONFIG_SAVE
  rm(CONFIG_SAVE)
}

# Extract full-sample FCI for later comparison
fullsample_fci <- resultado_fci$all_indices %>%
  dplyr::select(fecha,
                FCI_COMP_full = FCI_COMP_AVG,
                FCI_ENDO_full = FCI_ENDO_AVG,
                FCI_EXO_full = FCI_EXO_AVG)

# Full-sample PCA loadings
fullsample_pca_loadings <- resultado_fci$pca_loadings

# ---- Variable definitions (from script 01) ----
VARIABLES <- list(
  rates = list(
    vars = c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM"),
    signs = c(+1, +1, +1)
  ),
  banking = list(
    vars = c("Crecimiento_creditos", "Ratio_Cred_Depo", "Morosidad",
             "Rentabilidad", "Liquidez"),
    signs = c(-1, -1, +1, -1, -1)
  ),
  external = list(
    vars = c("TCN", "Commodities", "FFER", "VIX"),
    signs = c(+1, -1, +1, +1)
  )
)

LEVEL2 <- list(
  endogenous = list(
    vars = c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM",
             "Crecimiento_creditos", "Ratio_Cred_Depo", "Morosidad",
             "Rentabilidad", "Liquidez", "TCN"),
    signs = c(+1, +1, +1, -1, -1, +1, -1, -1, +1)
  ),
  exogenous = list(
    vars = c("FFER", "VIX", "Commodities"),
    signs = c(+1, +1, -1)
  )
)

all_vars_core <- unique(c(VARIABLES$rates$vars, VARIABLES$banking$vars,
                           VARIABLES$external$vars))
all_signs_core <- c(VARIABLES$rates$signs, VARIABLES$banking$signs,
                     VARIABLES$external$signs)

# ---- Load raw data, filter to post-IT ----
datos_raw <- read_excel(POSTIT_CONFIG$data_file)
fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

datos <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Calculate credit growth BEFORE filtering (need 12 lags)
if ("Creditos_Sector_privado_totales" %in% names(datos)) {
  datos <- datos %>%
    mutate(Crecimiento_creditos = (Creditos_Sector_privado_totales /
                                     lag(Creditos_Sector_privado_totales, 12) - 1) * 100)
}

# Filter to post-IT
datos <- datos %>% filter(fecha >= POSTIT_CONFIG$IT_START)

cat("Post-IT raw data:\n")
cat("  Period:", format(min(datos$fecha)), "to", format(max(datos$fecha)), "\n")
cat("  Observations:", nrow(datos), "\n\n")

# ---- Load macro data ----
macro_raw <- tryCatch({
  read_excel(POSTIT_CONFIG$data_file, sheet = POSTIT_CONFIG$macro_sheet)
}, error = function(e) {
  cat("Error loading macro sheet:", e$message, "\n")
  NULL
})

macro_vars <- c("IMAEP", "IPC", "IMAEP_SANB", "FBKf", "Consumo")
fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

macro_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Calculate YoY growth for macro variables (before filtering)
for (var in macro_vars) {
  if (var %in% names(macro_data)) {
    macro_data <- macro_data %>%
      mutate(!!paste0(var, "_yoy") := (!!sym(var) / lag(!!sym(var), 12) - 1) * 100)
  }
}

# Get IPC and deflated credit
ipc_credit_data <- macro_data %>%
  mutate(
    IPC_yoy = (IPC / lag(IPC, 12) - 1) * 100,
    Cred_Real_yoy = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100
  ) %>%
  dplyr::select(fecha, IPC_yoy, Cred_Real_yoy)

# Credit data by type
credit_data <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  left_join(
    macro_raw %>%
      rename(fecha = !!sym(fecha_col_macro)) %>%
      mutate(fecha = as.Date(fecha)) %>%
      dplyr::select(fecha, IPC, Creditos_deflactados),
    by = "fecha"
  ) %>%
  mutate(
    Cred_Total = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100,
    Cred_USD = (Creditos_Sector_privado_USD_equivalente /
                  lag(Creditos_Sector_privado_USD_equivalente, 12) - 1) * 100,
    Cred_Real_MN = ((Creditos_Sector_privado_MN / IPC) /
                      lag(Creditos_Sector_privado_MN / IPC, 12) - 1) * 100
  ) %>%
  dplyr::select(fecha, Cred_Total, Cred_USD, Cred_Real_MN)

# Filter macro and credit to post-IT
macro_data <- macro_data %>% filter(fecha >= POSTIT_CONFIG$IT_START)
credit_data <- credit_data %>% filter(fecha >= POSTIT_CONFIG$IT_START)
ipc_credit_data <- ipc_credit_data %>% filter(fecha >= POSTIT_CONFIG$IT_START)

cat("Post-IT macro observations:", nrow(macro_data), "\n")
cat("Post-IT credit observations:", nrow(credit_data), "\n\n")


################################################################################
# 3. HELPER FUNCTIONS
################################################################################

# -- FCI calculation functions (from script 01) --

apply_signs <- function(data, vars, signs) {
  data_adj <- data
  for (i in seq_along(vars)) {
    if (vars[i] %in% names(data_adj)) {
      data_adj[[vars[i]]] <- data_adj[[vars[i]]] * signs[i]
    }
  }
  return(data_adj)
}

rolling_standardize <- function(data, window = POSTIT_CONFIG$rolling_window) {
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

calculate_zscore <- function(data, vars) {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) == 0) return(rep(NA, nrow(data)))
  fci <- rowMeans(data[, vars_available, drop = FALSE], na.rm = TRUE)
  return(fci)
}

calculate_pca <- function(data, vars) {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) < 2) return(rep(NA, nrow(data)))

  datos_subset <- data[, vars_available, drop = FALSE]
  datos_clean <- datos_subset
  for (col in names(datos_clean)) {
    median_val <- median(datos_clean[[col]], na.rm = TRUE)
    datos_clean[[col]][is.na(datos_clean[[col]])] <- median_val
  }

  tryCatch({
    pca_model <- PCA(datos_clean, ncp = 1, graph = FALSE)
    fci <- pca_model$ind$coord[, 1]
    if (cor(fci, rowMeans(datos_clean), use = "complete.obs") < 0) {
      fci <- -fci
    }
    return(fci)
  }, error = function(e) {
    warning("PCA failed: ", e$message)
    return(rep(NA, nrow(data)))
  })
}

calculate_var <- function(data, vars) {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) < 2) return(rep(NA, nrow(data)))

  datos_subset <- data[, vars_available, drop = FALSE]
  datos_clean <- datos_subset
  for (col in names(datos_clean)) {
    median_val <- median(datos_clean[[col]], na.rm = TRUE)
    datos_clean[[col]][is.na(datos_clean[[col]])] <- median_val
  }

  n_obs <- nrow(datos_clean)

  tryCatch({
    var_select <- VARselect(datos_clean,
                            lag.max = min(POSTIT_CONFIG$var_max_lags, floor(n_obs / 10)),
                            type = "const")
    lag_opt <- max(1, min(var_select$selection["AIC(n)"], POSTIT_CONFIG$var_max_lags))

    var_model <- VAR(datos_clean, p = lag_opt, type = "const")
    fci <- c(rep(NA, lag_opt), rowMeans(fitted(var_model)))
    return(fci)
  }, error = function(e) {
    warning("VAR failed: ", e$message)
    return(rep(NA, nrow(data)))
  })
}

calculate_dfm <- function(data, vars) {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) < 2) return(rep(NA, nrow(data)))

  datos_subset <- data[, vars_available, drop = FALSE]
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
      control = list(maxit = POSTIT_CONFIG$dfm_max_iter)
    )

    fci <- as.numeric(dfm_model$states)
    if (cor(fci, rowMeans(datos_clean), use = "complete.obs") < 0) {
      fci <- -fci
    }
    return(fci)
  }, error = function(e) {
    warning("DFM failed: ", e$message)
    return(rep(NA, nrow(data)))
  })
}

calculate_fci_all_methods <- function(data, vars, prefix = "FCI") {
  vars_available <- intersect(vars, names(data))

  if (length(vars_available) < 2) {
    warning("Insufficient variables for ", prefix)
    return(NULL)
  }

  cat("  Calculating ", prefix, " (", length(vars_available), " variables)...\n", sep = "")

  results <- data.frame(fecha = data$fecha)

  results[[paste0(prefix, "_ZSCORE")]] <- calculate_zscore(data, vars)
  results[[paste0(prefix, "_PCA")]] <- calculate_pca(data, vars)
  results[[paste0(prefix, "_VAR")]] <- calculate_var(data, vars)
  results[[paste0(prefix, "_DFM")]] <- calculate_dfm(data, vars)

  methods <- c("ZSCORE", "PCA", "VAR", "DFM")
  for (m in methods) {
    col_name <- paste0(prefix, "_", m)
    if (col_name %in% names(results) && !all(is.na(results[[col_name]]))) {
      results[[paste0(col_name, "_norm")]] <- scale(results[[col_name]])[, 1]
    }
  }

  norm_cols <- paste0(prefix, "_", methods, "_norm")
  norm_cols <- norm_cols[norm_cols %in% names(results)]

  if (length(norm_cols) > 0) {
    results[[paste0(prefix, "_AVG")]] <- rowMeans(
      results[, norm_cols, drop = FALSE], na.rm = TRUE
    )
  }

  return(results)
}

# -- Check method convergence --
check_convergence <- function(results, prefix) {
  methods <- c("ZSCORE", "PCA", "VAR", "DFM")
  status <- sapply(methods, function(m) {
    col <- paste0(prefix, "_", m)
    col %in% names(results) && !all(is.na(results[[col]]))
  })
  names(status) <- methods
  converged <- names(status[status])
  failed <- names(status[!status])
  cat(sprintf("    CONVERGED: %s", paste(converged, collapse = ", ")))
  if (length(failed) > 0) cat(sprintf(" | FAILED: %s", paste(failed, collapse = ", ")))
  cat(sprintf(" [%d/4 methods]\n", length(converged)))
  return(list(converged = converged, failed = failed, n_converged = length(converged)))
}

# -- LP functions (from script 05) --

run_lp_standard <- function(data, y_var, fci_var, max_h, n_lags = 2, control_vars = NULL) {
  results <- data.frame()
  z_crit <- qnorm(1 - (1 - POSTIT_CONFIG$confidence_level) / 2)

  for (h in 1:max_h) {
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

    control_cols <- c()
    if (!is.null(control_vars)) {
      for (cv in control_vars) {
        if (cv %in% names(data)) {
          data_h <- data_h %>% mutate(!!cv := !!sym(cv))
          control_cols <- c(control_cols, cv)
          for (j in 1:n_lags) {
            cv_lag_name <- paste0(cv, "_lag", j)
            data_h <- data_h %>% mutate(!!cv_lag_name := lag(!!sym(cv), j))
            control_cols <- c(control_cols, cv_lag_name)
          }
        }
      }
    }

    formula_str <- paste("y_fwd ~", fci_var, "+", paste(c(lag_vars, control_cols), collapse = " + "))

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
      p_value = coef_test[idx, 4],
      ci_lower = coef_test[idx, 1] - z_crit * coef_test[idx, 2],
      ci_upper = coef_test[idx, 1] + z_crit * coef_test[idx, 2],
      n_obs = nrow(reg_data)
    ))
  }

  return(results)
}

run_lp_asymmetric <- function(data, y_var, fci_var, max_h, n_lags = 2, control_vars = NULL) {
  results <- data.frame()
  z_crit <- qnorm(1 - (1 - POSTIT_CONFIG$confidence_level) / 2)

  fci_pos <- paste0(fci_var, "_pos")
  fci_neg <- paste0(fci_var, "_neg")

  for (h in 1:max_h) {
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

    control_cols <- c()
    if (!is.null(control_vars)) {
      for (cv in control_vars) {
        if (cv %in% names(data)) {
          data_h <- data_h %>% mutate(!!cv := !!sym(cv))
          control_cols <- c(control_cols, cv)
          for (j in 1:n_lags) {
            cv_lag_name <- paste0(cv, "_lag", j)
            data_h <- data_h %>% mutate(!!cv_lag_name := lag(!!sym(cv), j))
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


################################################################################
# 4. PART A: FCI CONSTRUCTION DIAGNOSTICS
################################################################################

cat("================================================================================\n")
cat("PART A: FCI CONSTRUCTION DIAGNOSTICS (Post-IT Subsample)\n")
cat("================================================================================\n\n")

# ---- A1: Within-period PCA ----
cat("A1: Within-period PCA loadings...\n")

# Prepare post-IT data for PCA
datos_pca_postit <- datos %>%
  dplyr::select(fecha, any_of(all_vars_core))

datos_pca_postit <- apply_signs(datos_pca_postit, all_vars_core, all_signs_core)

# Scale post-IT data for PCA
pca_data <- datos_pca_postit %>% dplyr::select(-fecha)
pca_data_clean <- pca_data
for (col in names(pca_data_clean)) {
  median_val <- median(pca_data_clean[[col]], na.rm = TRUE)
  pca_data_clean[[col]][is.na(pca_data_clean[[col]])] <- median_val
}
pca_data_scaled <- as.data.frame(scale(pca_data_clean))

postit_pca <- PCA(pca_data_scaled, ncp = 5, graph = FALSE)
postit_loadings_raw <- postit_pca$var$coord[, 1]

# Sign correction
if (cor(postit_pca$ind$coord[, 1], rowMeans(pca_data_scaled), use = "complete.obs") < 0) {
  postit_loadings_raw <- -postit_loadings_raw
}

postit_var_explained <- postit_pca$eig[1, 2]  # % variance by PC1

cat(sprintf("  Post-IT PC1 variance explained: %.1f%%\n", postit_var_explained))

# Compare with full-sample loadings
loading_comparison <- data.frame(
  variable = names(postit_loadings_raw),
  postit_loading = as.numeric(postit_loadings_raw),
  fullsample_loading = if (!is.null(fullsample_pca_loadings)) {
    as.numeric(fullsample_pca_loadings[names(postit_loadings_raw)])
  } else {
    rep(NA, length(postit_loadings_raw))
  }
)
loading_comparison$sign_consistent <- sign(loading_comparison$postit_loading) ==
  sign(loading_comparison$fullsample_loading)

cat("  Sign reversals vs full-sample: ",
    sum(!loading_comparison$sign_consistent, na.rm = TRUE), "/",
    nrow(loading_comparison), " variables\n\n")

# ---- A2: 4-method FCI on post-IT data ----
cat("A2: Constructing post-IT FCI using 4 methods...\n\n")

# Prepare standardized data for post-IT
datos_signed <- apply_signs(datos, all_vars_core, all_signs_core)
datos_std_postit <- datos_signed %>%
  dplyr::select(fecha, any_of(all_vars_core)) %>%
  rolling_standardize(POSTIT_CONFIG$rolling_window)

methods_status <- list()

# COMP (12 vars)
fci_comp <- calculate_fci_all_methods(datos_std_postit, all_vars_core, "FCI_COMP")
methods_status$COMP <- check_convergence(fci_comp, "FCI_COMP")

# ENDO (9 vars)
datos_endo <- apply_signs(datos, LEVEL2$endogenous$vars, LEVEL2$endogenous$signs)
datos_endo_std <- datos_endo %>%
  dplyr::select(fecha, any_of(LEVEL2$endogenous$vars)) %>%
  rolling_standardize(POSTIT_CONFIG$rolling_window)
fci_endo <- calculate_fci_all_methods(datos_endo_std, LEVEL2$endogenous$vars, "FCI_ENDO")
methods_status$ENDO <- check_convergence(fci_endo, "FCI_ENDO")

# EXO (3 vars)
datos_exo <- apply_signs(datos, LEVEL2$exogenous$vars, LEVEL2$exogenous$signs)
datos_exo_std <- datos_exo %>%
  dplyr::select(fecha, any_of(LEVEL2$exogenous$vars)) %>%
  rolling_standardize(POSTIT_CONFIG$rolling_window)
fci_exo <- calculate_fci_all_methods(datos_exo_std, LEVEL2$exogenous$vars, "FCI_EXO")
methods_status$EXO <- check_convergence(fci_exo, "FCI_EXO")

# exCredit (11 vars)
vars_exCredit <- setdiff(all_vars_core, "Crecimiento_creditos")
signs_exCredit <- all_signs_core[all_vars_core != "Crecimiento_creditos"]
datos_exCred <- apply_signs(datos, vars_exCredit, signs_exCredit)
datos_exCred_std <- datos_exCred %>%
  dplyr::select(fecha, any_of(vars_exCredit)) %>%
  rolling_standardize(POSTIT_CONFIG$rolling_window)
fci_exCredit <- calculate_fci_all_methods(datos_exCred_std, vars_exCredit, "FCI_exCredit")
methods_status$exCredit <- check_convergence(fci_exCredit, "FCI_exCredit")

# ENDO_exCredit (8 vars)
vars_ENDO_exCredit <- setdiff(LEVEL2$endogenous$vars, "Crecimiento_creditos")
signs_ENDO_exCredit <- LEVEL2$endogenous$signs[LEVEL2$endogenous$vars != "Crecimiento_creditos"]
datos_endoExCred <- apply_signs(datos, vars_ENDO_exCredit, signs_ENDO_exCredit)
datos_endoExCred_std <- datos_endoExCred %>%
  dplyr::select(fecha, any_of(vars_ENDO_exCredit)) %>%
  rolling_standardize(POSTIT_CONFIG$rolling_window)
fci_endo_exCredit <- calculate_fci_all_methods(datos_endoExCred_std, vars_ENDO_exCredit,
                                                "FCI_ENDO_exCredit")
methods_status$ENDO_exCredit <- check_convergence(fci_endo_exCredit, "FCI_ENDO_exCredit")

# ---- Merge all post-IT FCI ----
postit_fci <- fci_comp %>%
  dplyr::select(fecha, FCI_COMP = FCI_COMP_AVG)

if (!is.null(fci_endo) && "FCI_ENDO_AVG" %in% names(fci_endo)) {
  postit_fci <- postit_fci %>%
    left_join(fci_endo %>% dplyr::select(fecha, FCI_ENDO = FCI_ENDO_AVG), by = "fecha")
}
if (!is.null(fci_exo) && "FCI_EXO_AVG" %in% names(fci_exo)) {
  postit_fci <- postit_fci %>%
    left_join(fci_exo %>% dplyr::select(fecha, FCI_EXO = FCI_EXO_AVG), by = "fecha")
}
if (!is.null(fci_exCredit) && "FCI_exCredit_AVG" %in% names(fci_exCredit)) {
  postit_fci <- postit_fci %>%
    left_join(fci_exCredit %>% dplyr::select(fecha, FCI_exCredit = FCI_exCredit_AVG), by = "fecha")
}
if (!is.null(fci_endo_exCredit) && "FCI_ENDO_exCredit_AVG" %in% names(fci_endo_exCredit)) {
  postit_fci <- postit_fci %>%
    left_join(fci_endo_exCredit %>% dplyr::select(fecha, FCI_ENDO_exCredit = FCI_ENDO_exCredit_AVG),
              by = "fecha")
}

# ---- A3: Cross-method correlations ----
cat("\nA3: Cross-method correlations (post-IT)...\n")

method_names <- c("ZSCORE", "PCA", "VAR", "DFM")
comp_methods <- fci_comp %>%
  dplyr::select(any_of(paste0("FCI_COMP_", method_names)))
comp_methods <- comp_methods[, colSums(!is.na(comp_methods)) > 10, drop = FALSE]

if (ncol(comp_methods) >= 2) {
  cross_method_cor <- cor(comp_methods, use = "pairwise.complete.obs")
  cat("  Cross-method correlation range: [",
      sprintf("%.3f", min(cross_method_cor[lower.tri(cross_method_cor)])), ", ",
      sprintf("%.3f", max(cross_method_cor[lower.tri(cross_method_cor)])), "]\n")
} else {
  cross_method_cor <- NULL
}

# ---- A4: Variance decomposition ----
cat("\nA4: Variance decomposition (post-IT)...\n")

var_decomp <- tryCatch({
  df_vd <- postit_fci %>% dplyr::select(FCI_ENDO, FCI_EXO) %>% na.omit()
  if (nrow(df_vd) > 30) {
    model_vd <- lm(FCI_ENDO ~ FCI_EXO, data = df_vd)
    r2 <- summary(model_vd)$r.squared
    cat(sprintf("  External share (R2): %.1f%% (full-sample: 9.1%%)\n", r2 * 100))
    r2
  } else {
    NA
  }
}, error = function(e) {
  cat("  Variance decomposition failed:", e$message, "\n")
  NA
})

# ---- A5: Post-IT vs full-sample FCI correlation ----
cat("\nA5: Post-IT FCI vs full-sample FCI correlation...\n")

fci_compare <- postit_fci %>%
  inner_join(fullsample_fci, by = "fecha") %>%
  na.omit()

if (nrow(fci_compare) > 10) {
  corr_comp <- cor(fci_compare$FCI_COMP, fci_compare$FCI_COMP_full)
  cat(sprintf("  Post-IT vs full-sample FCI_COMP: r = %.3f (n=%d)\n",
              corr_comp, nrow(fci_compare)))
} else {
  corr_comp <- NA
}

# ---- Plot 170: FCI Construction Diagnostics ----
cat("\nGenerating plots...\n")

# Panel A: PCA loading comparison
p170a <- ggplot(loading_comparison, aes(x = reorder(variable, postit_loading))) +
  geom_col(aes(y = postit_loading, fill = "Post-IT"), alpha = 0.8, width = 0.4,
           position = position_nudge(x = -0.2)) +
  geom_col(aes(y = fullsample_loading, fill = "Full Sample"), alpha = 0.8, width = 0.4,
           position = position_nudge(x = 0.2)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Post-IT" = "#E41A1C", "Full Sample" = "#377EB8")) +
  coord_flip() +
  labs(title = "PCA Loadings", subtitle = sprintf("PC1: %.1f%% var (post-IT)", postit_var_explained),
       x = NULL, y = "Loading", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Panel B: Post-IT FCI time series
fci_ts_long <- postit_fci %>%
  dplyr::select(fecha, FCI_COMP, FCI_ENDO, FCI_EXO) %>%
  pivot_longer(-fecha, names_to = "index", values_to = "value") %>%
  na.omit()

p170b <- ggplot(fci_ts_long, aes(x = fecha, y = value, color = index)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("FCI_COMP" = "#377EB8", "FCI_ENDO" = "#4DAF4A",
                                 "FCI_EXO" = "#FF7F00"),
                     labels = c("Comprehensive", "Endogenous", "Exogenous")) +
  labs(title = "Post-IT FCI Indices",
       subtitle = sprintf("Average of %d/4 methods", methods_status$COMP$n_converged),
       x = NULL, y = "FCI (std. units)", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Panel C: Cross-method correlation heatmap
if (!is.null(cross_method_cor)) {
  cor_df <- expand.grid(Var1 = rownames(cross_method_cor),
                        Var2 = colnames(cross_method_cor), stringsAsFactors = FALSE)
  cor_df$value <- as.vector(cross_method_cor)
  # Clean names for display
  cor_df$Var1 <- gsub("FCI_COMP_", "", cor_df$Var1)
  cor_df$Var2 <- gsub("FCI_COMP_", "", cor_df$Var2)

  p170c <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0.5, limits = c(0, 1)) +
    labs(title = "Cross-Method Correlation",
         subtitle = sprintf("%d converged methods", methods_status$COMP$n_converged),
         x = NULL, y = NULL, fill = "r") +
    theme_minimal()
} else {
  p170c <- ggplot() + theme_void() + labs(title = "Cross-Method Correlation (insufficient data)")
}

# Panel D: Post-IT vs full-sample scatter
if (nrow(fci_compare) > 10) {
  p170d <- ggplot(fci_compare, aes(x = FCI_COMP_full, y = FCI_COMP)) +
    geom_point(alpha = 0.4, size = 1.5, color = "#377EB8") +
    geom_smooth(method = "lm", se = TRUE, color = "#E41A1C", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Post-IT vs Full-Sample FCI",
         subtitle = sprintf("r = %.3f, n = %d", corr_comp, nrow(fci_compare)),
         x = "Full-Sample FCI_COMP", y = "Post-IT FCI_COMP") +
    theme_minimal()
} else {
  p170d <- ggplot() + theme_void() + labs(title = "Post-IT vs Full-Sample (insufficient data)")
}

p170 <- grid.arrange(
  p170a, p170b, p170c, p170d,
  ncol = 2,
  top = grid::textGrob("Post-IT FCI Construction Diagnostics",
                        gp = grid::gpar(fontsize = 16, fontface = "bold"))
)

ggsave(file.path(POSTIT_CONFIG$output_dir, "170_PostIT_FCI_Construction.png"), p170,
       width = 14, height = 12, dpi = 300)
cat("Saved: 170_PostIT_FCI_Construction.png\n")

# Export diagnostics CSV
diag_export <- loading_comparison
diag_export$postit_var_explained_pct <- postit_var_explained
diag_export$external_share_r2 <- ifelse(is.na(var_decomp), NA, var_decomp)
diag_export$fci_fullsample_corr <- ifelse(is.na(corr_comp), NA, corr_comp)
diag_export$n_methods_COMP <- methods_status$COMP$n_converged
diag_export$methods_converged_COMP <- paste(methods_status$COMP$converged, collapse = ",")
diag_export$n_obs_postit <- nrow(datos)

write.csv(diag_export, file.path(POSTIT_CONFIG$output_dir, "PostIT_FCI_Diagnostics.csv"),
          row.names = FALSE)
cat("Saved: PostIT_FCI_Diagnostics.csv\n\n")


################################################################################
# 5. PART B: GRANGER CAUSALITY
################################################################################

cat("================================================================================\n")
cat("PART B: GRANGER CAUSALITY (Post-IT)\n")
cat("================================================================================\n\n")

# Merge FCI with macro/credit data for Granger tests
granger_data <- postit_fci %>%
  inner_join(credit_data, by = "fecha") %>%
  inner_join(ipc_credit_data, by = "fecha") %>%
  left_join(macro_data %>% dplyr::select(fecha, IMAEP_yoy), by = "fecha") %>%
  arrange(fecha) %>%
  na.omit()

cat("Granger data: n =", nrow(granger_data), "\n\n")

granger_results <- data.frame()
granger_tests <- list(
  list(y = "Cred_Total", x = "FCI_exCredit", label = "FCI_exCredit -> Real Credit"),
  list(y = "IMAEP_yoy", x = "FCI_COMP", label = "FCI_COMP -> IMAEP"),
  list(y = "FCI_ENDO", x = "FCI_EXO", label = "FCI_EXO -> FCI_ENDO")
)

for (test in granger_tests) {
  for (lag in c(3, 6, 12)) {
    if (lag >= nrow(granger_data) / 3) next

    tryCatch({
      gt <- grangertest(
        as.formula(paste(test$y, "~", test$x)),
        order = lag,
        data = granger_data
      )

      granger_results <- rbind(granger_results, data.frame(
        direction = test$label,
        lags = lag,
        F_stat = gt$F[2],
        p_value = gt$`Pr(>F)`[2],
        n_obs = nrow(granger_data)
      ))

      stars <- ifelse(gt$`Pr(>F)`[2] < 0.01, "***",
                       ifelse(gt$`Pr(>F)`[2] < 0.05, "**",
                              ifelse(gt$`Pr(>F)`[2] < 0.10, "*", "")))
      cat(sprintf("  %s (lag=%d): F=%.2f, p=%.3f %s\n",
                  test$label, lag, gt$F[2], gt$`Pr(>F)`[2], stars))
    }, error = function(e) {
      cat(sprintf("  %s (lag=%d): FAILED - %s\n", test$label, lag, e$message))
    })
  }
}

write.csv(granger_results, file.path(POSTIT_CONFIG$output_dir, "PostIT_Granger.csv"),
          row.names = FALSE)
cat("\nSaved: PostIT_Granger.csv\n\n")


################################################################################
# 6. PART C: CREDIT LOCAL PROJECTIONS
################################################################################

cat("================================================================================\n")
cat("PART C: CREDIT LOCAL PROJECTIONS (Post-IT)\n")
cat("================================================================================\n\n")

# Build analysis dataset for LPs
lp_data <- postit_fci %>%
  inner_join(credit_data, by = "fecha") %>%
  inner_join(ipc_credit_data, by = "fecha") %>%
  left_join(macro_data %>% dplyr::select(fecha, IMAEP_yoy), by = "fecha") %>%
  arrange(fecha) %>%
  mutate(
    FCI_exCredit_pos = pmax(FCI_exCredit, 0),
    FCI_exCredit_neg = abs(pmin(FCI_exCredit, 0)),
    FCI_ENDO_exCredit_pos = pmax(FCI_ENDO_exCredit, 0),
    FCI_ENDO_exCredit_neg = abs(pmin(FCI_ENDO_exCredit, 0))
  ) %>%
  na.omit()

cat("LP data: n =", nrow(lp_data), "\n\n")

# C1: Credit by type (using FCI_exCredit for endogeneity correction)
credit_lp_results <- list()
credit_types <- c("Cred_Total", "Cred_Real_MN", "Cred_USD")
credit_labels <- c("Cred_Total" = "Real Total Credit",
                    "Cred_Real_MN" = "Real MN Credit",
                    "Cred_USD" = "USD Credit")

for (cred in credit_types) {
  cat(sprintf("  FCI_exCredit -> %s ... ", credit_labels[cred]))
  credit_lp_results[[cred]] <- run_lp_standard(
    lp_data, cred, "FCI_exCredit", POSTIT_CONFIG$max_horizon, POSTIT_CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  ) %>% mutate(credit_type = cred, fci_type = "FCI_exCredit")
  cat("done\n")
}

all_credit_lp <- bind_rows(credit_lp_results)

# C1b: Credit by type using FCI_ENDO_exCredit (domestic component, for currency decomposition)
cat("\n  Currency decomposition with FCI_ENDO_exCredit:\n")
credit_endo_lp_results <- list()
for (cred in credit_types) {
  cat(sprintf("  FCI_ENDO_exCredit -> %s ... ", credit_labels[cred]))
  credit_endo_lp_results[[cred]] <- run_lp_standard(
    lp_data, cred, "FCI_ENDO_exCredit", POSTIT_CONFIG$max_horizon, POSTIT_CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  ) %>% mutate(credit_type = cred, fci_type = "FCI_ENDO_exCredit")
  cat("done\n")
}
all_credit_lp <- bind_rows(all_credit_lp, bind_rows(credit_endo_lp_results))

# C2: FCI decomposition (exCredit, ENDO_exCredit, EXO) -> Total credit
fci_decomp_results <- list()
fci_decomp_types <- c("FCI_exCredit", "FCI_ENDO_exCredit", "FCI_EXO")
fci_decomp_labels <- c("FCI_exCredit" = "Comprehensive (exCredit)",
                        "FCI_ENDO_exCredit" = "Endogenous (exCredit)",
                        "FCI_EXO" = "Exogenous")

for (fci_type in fci_decomp_types) {
  cat(sprintf("  %s -> Cred_Total ... ", fci_decomp_labels[fci_type]))
  fci_decomp_results[[fci_type]] <- run_lp_standard(
    lp_data, "Cred_Total", fci_type, POSTIT_CONFIG$max_horizon, POSTIT_CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  ) %>% mutate(fci_type = fci_type, credit_type = "Cred_Total")
  cat("done\n")
}

all_decomp_lp <- bind_rows(fci_decomp_results)

# Combine for export
all_credit_export <- bind_rows(all_credit_lp, all_decomp_lp)
write.csv(all_credit_export, file.path(POSTIT_CONFIG$output_dir, "PostIT_LP_Credit.csv"),
          row.names = FALSE)
cat("Saved: PostIT_LP_Credit.csv\n")

# Print summary at key horizons
cat("\nCredit LP Summary (Post-IT):\n")
cat(sprintf("  %3s  %12s  %8s  %8s  %6s\n", "h", "Credit Type", "Coef", "SE", ""))
for (cred in credit_types) {
  for (h in POSTIT_CONFIG$horizons_report) {
    r <- all_credit_lp %>% filter(credit_type == cred, horizon == h)
    if (nrow(r) > 0) {
      stars <- ifelse(r$p_value < 0.01, "***", ifelse(r$p_value < 0.05, "**",
                                                       ifelse(r$p_value < 0.10, "*", "")))
      cat(sprintf("  %3d  %12s  %8.2f  %8.2f  %s\n", h, cred, r$coef, r$se, stars))
    }
  }
}

# ---- Plot 171: Credit LP by type ----
p171 <- ggplot(all_credit_lp, aes(x = horizon, y = coef)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "#E41A1C") +
  geom_line(color = "#E41A1C", linewidth = 0.8) +
  geom_point(color = "#E41A1C", size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~credit_type, scales = "free_y", ncol = 3,
             labeller = labeller(credit_type = credit_labels)) +
  labs(title = "Post-IT Credit Local Projections",
       subtitle = "Response to 1 SD FCI_exCredit tightening | 90% CI (Newey-West)",
       x = "Horizon (months)", y = "Effect (pp)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold"))

ggsave(file.path(POSTIT_CONFIG$output_dir, "171_PostIT_LP_Credit.png"), p171,
       width = 14, height = 5, dpi = 300)
cat("\nSaved: 171_PostIT_LP_Credit.png\n")

# ---- Plot 172: Credit decomposition ----
p172 <- ggplot(all_decomp_lp, aes(x = horizon, y = coef, color = fci_type)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = fci_type), alpha = 0.15) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("FCI_exCredit" = "#377EB8", "FCI_ENDO_exCredit" = "#4DAF4A",
                                 "FCI_EXO" = "#FF7F00"),
                     labels = fci_decomp_labels) +
  scale_fill_manual(values = c("FCI_exCredit" = "#377EB8", "FCI_ENDO_exCredit" = "#4DAF4A",
                                "FCI_EXO" = "#FF7F00"),
                    labels = fci_decomp_labels) +
  labs(title = "Post-IT Credit Channel Decomposition",
       subtitle = "FCI components -> Real Total Credit | 90% CI",
       x = "Horizon (months)", y = "Effect (pp)", color = NULL, fill = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom")

ggsave(file.path(POSTIT_CONFIG$output_dir, "172_PostIT_LP_Credit_Decomposition.png"), p172,
       width = 10, height = 6, dpi = 300)
cat("Saved: 172_PostIT_LP_Credit_Decomposition.png\n\n")


################################################################################
# 7. PART D: SECTORAL OUTPUT
################################################################################

cat("================================================================================\n")
cat("PART D: SECTORAL OUTPUT LOCAL PROJECTIONS (Post-IT)\n")
cat("================================================================================\n\n")

# Add macro variables to LP data
output_vars_yoy <- c("IMAEP_yoy", "IMAEP_SANB_yoy", "FBKf_yoy", "Consumo_yoy")
output_labels <- c("IMAEP_yoy" = "Aggregate GDP",
                    "IMAEP_SANB_yoy" = "Non-Agro GDP",
                    "FBKf_yoy" = "Investment (FBKf)",
                    "Consumo_yoy" = "Consumption")

lp_output_data <- postit_fci %>%
  inner_join(ipc_credit_data, by = "fecha") %>%
  left_join(macro_data %>% dplyr::select(fecha, any_of(output_vars_yoy)), by = "fecha") %>%
  arrange(fecha) %>%
  na.omit()

cat("Output LP data: n =", nrow(lp_output_data), "\n\n")

output_lp_results <- list()
available_output_vars <- intersect(output_vars_yoy, names(lp_output_data))

for (var in available_output_vars) {
  cat(sprintf("  FCI_COMP -> %s ... ", output_labels[var]))
  output_lp_results[[var]] <- run_lp_standard(
    lp_output_data, var, "FCI_COMP", POSTIT_CONFIG$max_horizon, POSTIT_CONFIG$n_lags,
    control_vars = c("IPC_yoy", "Cred_Real_yoy")
  ) %>% mutate(variable = var)
  cat("done\n")
}

all_output_lp <- bind_rows(output_lp_results)

write.csv(all_output_lp, file.path(POSTIT_CONFIG$output_dir, "PostIT_LP_Output.csv"),
          row.names = FALSE)
cat("Saved: PostIT_LP_Output.csv\n")

# Print summary
cat("\nOutput LP Summary (Post-IT):\n")
for (var in available_output_vars) {
  r12 <- all_output_lp %>% filter(variable == var, horizon == 12)
  if (nrow(r12) > 0) {
    stars <- ifelse(r12$p_value < 0.01, "***", ifelse(r12$p_value < 0.05, "**",
                                                       ifelse(r12$p_value < 0.10, "*", "")))
    cat(sprintf("  %s (h=12): %.2f (p=%.3f) %s\n", output_labels[var], r12$coef, r12$p_value, stars))
  }
}

# ---- Plot 173: Output LP ----
p173 <- ggplot(all_output_lp, aes(x = horizon, y = coef)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "#377EB8") +
  geom_line(color = "#377EB8", linewidth = 0.8) +
  geom_point(color = "#377EB8", size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~variable, scales = "free_y", ncol = 2,
             labeller = labeller(variable = output_labels)) +
  labs(title = "Post-IT Output Local Projections",
       subtitle = "Response to 1 SD FCI_COMP tightening | 90% CI (Newey-West)",
       x = "Horizon (months)", y = "Effect (pp)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold"))

ggsave(file.path(POSTIT_CONFIG$output_dir, "173_PostIT_LP_Output.png"), p173,
       width = 12, height = 8, dpi = 300)
cat("\nSaved: 173_PostIT_LP_Output.png\n\n")


################################################################################
# 8. PART E: ASYMMETRIC EFFECTS
################################################################################

cat("================================================================================\n")
cat("PART E: ASYMMETRIC EFFECTS (Post-IT)\n")
cat("================================================================================\n\n")

asym_results <- list()

for (cred in credit_types) {
  cat(sprintf("  Asymmetric FCI_exCredit -> %s ... ", credit_labels[cred]))
  asym_results[[cred]] <- run_lp_asymmetric(
    lp_data, cred, "FCI_exCredit", POSTIT_CONFIG$max_horizon, POSTIT_CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  ) %>% mutate(credit_type = cred)
  cat("done\n")
}

all_asym <- bind_rows(asym_results)

write.csv(all_asym, file.path(POSTIT_CONFIG$output_dir, "PostIT_LP_Asymmetric.csv"),
          row.names = FALSE)
cat("Saved: PostIT_LP_Asymmetric.csv\n")

# ---- Plot 174: Asymmetric effects ----
asym_long <- all_asym %>%
  dplyr::select(horizon, credit_type,
                coef_tight, ci_tight_lo, ci_tight_hi,
                coef_ease, ci_ease_lo, ci_ease_hi) %>%
  pivot_longer(
    cols = c(coef_tight, coef_ease),
    names_to = "direction",
    values_to = "coef"
  ) %>%
  mutate(
    ci_lower = ifelse(direction == "coef_tight", ci_tight_lo, ci_ease_lo),
    ci_upper = ifelse(direction == "coef_tight", ci_tight_hi, ci_ease_hi),
    direction_label = ifelse(direction == "coef_tight", "Tightening", "Easing")
  )

p174 <- ggplot(asym_long, aes(x = horizon, y = coef, color = direction_label)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = direction_label), alpha = 0.15) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~credit_type, scales = "free_y", ncol = 3,
             labeller = labeller(credit_type = credit_labels)) +
  scale_color_manual(values = c("Tightening" = "#E41A1C", "Easing" = "#377EB8")) +
  scale_fill_manual(values = c("Tightening" = "#E41A1C", "Easing" = "#377EB8")) +
  labs(title = "Post-IT Asymmetric Credit Effects",
       subtitle = "Tightening vs Easing | 90% CI (Newey-West)",
       x = "Horizon (months)", y = "Effect (pp)", color = NULL, fill = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(POSTIT_CONFIG$output_dir, "174_PostIT_LP_Asymmetric.png"), p174,
       width = 14, height = 5, dpi = 300)
cat("Saved: 174_PostIT_LP_Asymmetric.png\n\n")


################################################################################
# 9. PART F: FULL-SAMPLE vs POST-IT COMPARISON DASHBOARD
################################################################################

cat("================================================================================\n")
cat("PART F: FULL-SAMPLE vs POST-IT COMPARISON\n")
cat("================================================================================\n\n")

# ---- F1: Load full-sample LP results ----
cat("Loading full-sample results for comparison...\n")

fullsample_credit <- tryCatch(
  read.csv(file.path(POSTIT_CONFIG$output_dir, "LP_Credit_Standard.csv")),
  error = function(e) { cat("  LP_Credit_Standard.csv not found\n"); NULL }
)
fullsample_macro <- tryCatch(
  read.csv(file.path(POSTIT_CONFIG$output_dir, "LP_Macro_Standard.csv")),
  error = function(e) { cat("  LP_Macro_Standard.csv not found\n"); NULL }
)

# ---- F2: Credit overlay ----
if (!is.null(fullsample_credit)) {
  cat("Creating credit LP overlay...\n")

  # Full-sample credit uses fci_type column - get FCI_COMP (which maps to FCI_exCredit)
  full_credit_comp <- fullsample_credit %>%
    filter(fci_type == "FCI_COMP") %>%
    mutate(sample = "Full Sample",
           credit_type_std = paste0("Cred_", credit_type))

  postit_credit_comp <- all_credit_lp %>%
    filter(fci_type == "FCI_exCredit") %>%
    mutate(sample = "Post-IT",
           credit_type_std = credit_type)

  # Ensure matching columns
  cols_common <- c("horizon", "coef", "se", "p_value", "ci_lower", "ci_upper",
                    "credit_type_std", "sample")
  full_sub <- full_credit_comp %>% dplyr::select(any_of(cols_common))
  postit_sub <- postit_credit_comp %>% dplyr::select(any_of(cols_common))

  overlay_credit <- bind_rows(full_sub, postit_sub)

  p175 <- ggplot(overlay_credit, aes(x = horizon, y = coef, color = sample)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = sample), alpha = 0.12) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~credit_type_std, scales = "free_y", ncol = 3,
               labeller = labeller(credit_type_std = c("Cred_Total" = "Real Total Credit",
                                                        "Cred_Real_MN" = "Real MN Credit",
                                                        "Cred_USD" = "USD Credit"))) +
    scale_color_manual(values = c("Full Sample" = "#377EB8", "Post-IT" = "#E41A1C")) +
    scale_fill_manual(values = c("Full Sample" = "#377EB8", "Post-IT" = "#E41A1C")) +
    labs(title = "Credit LP: Full Sample vs Post-IT",
         subtitle = "FCI_exCredit -> Credit Growth | 90% CI (Newey-West)",
         x = "Horizon (months)", y = "Effect (pp)", color = NULL, fill = NULL) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(POSTIT_CONFIG$output_dir, "175_PostIT_LP_vs_FullSample.png"), p175,
         width = 14, height = 5, dpi = 300)
  cat("Saved: 175_PostIT_LP_vs_FullSample.png\n")
}

# ---- F3: Output overlay ----
if (!is.null(fullsample_macro)) {
  cat("Creating output LP overlay...\n")

  full_macro_sub <- fullsample_macro %>%
    filter(variable %in% c("IMAEP", "Creditos")) %>%
    mutate(sample = "Full Sample",
           variable = ifelse(variable == "IMAEP", "IMAEP_yoy", "Cred_Total_from_macro"))

  postit_macro_sub <- all_output_lp %>%
    filter(variable %in% c("IMAEP_yoy", "IMAEP_SANB_yoy")) %>%
    mutate(sample = "Post-IT")

  # Only overlay IMAEP
  full_imaep <- fullsample_macro %>%
    filter(variable == "IMAEP") %>%
    mutate(sample = "Full Sample", variable = "IMAEP_yoy")

  postit_imaep <- all_output_lp %>%
    filter(variable == "IMAEP_yoy") %>%
    mutate(sample = "Post-IT")

  overlay_output <- bind_rows(
    full_imaep %>% dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, variable, sample),
    postit_imaep %>% dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, variable, sample)
  )

  # Also show post-IT IMAEP_SANB for comparison
  postit_sanb <- all_output_lp %>%
    filter(variable == "IMAEP_SANB_yoy") %>%
    mutate(sample = "Post-IT (SANB)")

  if (nrow(postit_sanb) > 0) {
    overlay_output <- bind_rows(
      overlay_output,
      postit_sanb %>% dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, variable, sample) %>%
        mutate(variable = "IMAEP_yoy")  # Same facet
    )
  }

  p176 <- ggplot(overlay_output, aes(x = horizon, y = coef, color = sample)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = sample), alpha = 0.12) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Full Sample" = "#377EB8", "Post-IT" = "#E41A1C",
                                   "Post-IT (SANB)" = "#4DAF4A")) +
    scale_fill_manual(values = c("Full Sample" = "#377EB8", "Post-IT" = "#E41A1C",
                                  "Post-IT (SANB)" = "#4DAF4A")) +
    labs(title = "Output LP: Full Sample vs Post-IT",
         subtitle = "FCI_COMP -> Output Growth | 90% CI (Newey-West)",
         x = "Horizon (months)", y = "Effect (pp)", color = NULL, fill = NULL) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")

  ggsave(file.path(POSTIT_CONFIG$output_dir, "176_PostIT_Output_vs_FullSample.png"), p176,
         width = 10, height = 6, dpi = 300)
  cat("Saved: 176_PostIT_Output_vs_FullSample.png\n")
}

# ---- F4: PCA loading 3-way bar chart ----
cat("Creating PCA loading comparison...\n")

# Load pre-IT loadings from regime PCA comparison
regime_pca <- tryCatch(
  read.csv(file.path(POSTIT_CONFIG$output_dir, "FCI_Regime_PCA_Comparison.csv")),
  error = function(e) { cat("  FCI_Regime_PCA_Comparison.csv not found\n"); NULL }
)

if (!is.null(regime_pca) && "pre_IT" %in% names(regime_pca)) {
  loading_3way <- data.frame(
    variable = loading_comparison$variable,
    `Pre-IT` = regime_pca$pre_IT[match(loading_comparison$variable, regime_pca$variable)],
    `Post-IT` = loading_comparison$postit_loading,
    `Full Sample` = loading_comparison$fullsample_loading,
    check.names = FALSE
  )
} else {
  loading_3way <- data.frame(
    variable = loading_comparison$variable,
    `Post-IT` = loading_comparison$postit_loading,
    `Full Sample` = loading_comparison$fullsample_loading,
    check.names = FALSE
  )
}

loading_3way_long <- loading_3way %>%
  pivot_longer(-variable, names_to = "period", values_to = "loading") %>%
  filter(!is.na(loading))

# Order variables by post-IT loading
var_order <- loading_comparison %>% arrange(postit_loading) %>% pull(variable)
loading_3way_long$variable <- factor(loading_3way_long$variable, levels = var_order)

# Order periods
period_levels <- c("Pre-IT", "Post-IT", "Full Sample")
loading_3way_long$period <- factor(loading_3way_long$period,
                                    levels = intersect(period_levels, unique(loading_3way_long$period)))

p177 <- ggplot(loading_3way_long, aes(x = variable, y = loading, fill = period)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Pre-IT" = "#FF7F00", "Post-IT" = "#E41A1C",
                                "Full Sample" = "#377EB8")) +
  coord_flip() +
  labs(title = "PCA Loading Comparison Across Periods",
       subtitle = "First Principal Component Loadings (sign-corrected)",
       x = NULL, y = "Loading", fill = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom")

ggsave(file.path(POSTIT_CONFIG$output_dir, "177_PostIT_PCA_Loading_Comparison.png"), p177,
       width = 10, height = 8, dpi = 300)
cat("Saved: 177_PostIT_PCA_Loading_Comparison.png\n")

# ---- F5: Summary comparison table + z-test ----
cat("Creating comparison summary...\n")

comparison_rows <- list()

# Compare credit LP at key horizons
if (!is.null(fullsample_credit)) {
  for (h in POSTIT_CONFIG$horizons_report) {
    for (cred in c("Total", "Real_MN", "USD")) {
      full_r <- fullsample_credit %>%
        filter(fci_type == "FCI_COMP", credit_type == cred, horizon == h)
      postit_r <- all_credit_lp %>%
        filter(fci_type == "FCI_exCredit", credit_type == paste0("Cred_", cred), horizon == h)

      if (nrow(full_r) > 0 && nrow(postit_r) > 0) {
        z_stat <- abs(full_r$coef - postit_r$coef) / sqrt(full_r$se^2 + postit_r$se^2)
        z_p <- 2 * (1 - pnorm(z_stat))

        comparison_rows[[length(comparison_rows) + 1]] <- data.frame(
          variable = paste0("Credit_", cred),
          horizon = h,
          full_coef = full_r$coef,
          full_se = full_r$se,
          full_p = full_r$p_value,
          postit_coef = postit_r$coef,
          postit_se = postit_r$se,
          postit_p = postit_r$p_value,
          z_stat_equality = z_stat,
          z_p_equality = z_p,
          ci_overlap = !(full_r$ci_upper < postit_r$ci_lower | postit_r$ci_upper < full_r$ci_lower)
        )
      }
    }
  }
}

# Compare output LP at key horizons
if (!is.null(fullsample_macro)) {
  for (h in POSTIT_CONFIG$horizons_report) {
    full_r <- fullsample_macro %>% filter(variable == "IMAEP", horizon == h)
    postit_r <- all_output_lp %>% filter(variable == "IMAEP_yoy", horizon == h)

    if (nrow(full_r) > 0 && nrow(postit_r) > 0) {
      z_stat <- abs(full_r$coef - postit_r$coef) / sqrt(full_r$se^2 + postit_r$se^2)
      z_p <- 2 * (1 - pnorm(z_stat))

      comparison_rows[[length(comparison_rows) + 1]] <- data.frame(
        variable = "IMAEP",
        horizon = h,
        full_coef = full_r$coef,
        full_se = full_r$se,
        full_p = full_r$p_value,
        postit_coef = postit_r$coef,
        postit_se = postit_r$se,
        postit_p = postit_r$p_value,
        z_stat_equality = z_stat,
        z_p_equality = z_p,
        ci_overlap = !(full_r$ci_upper < postit_r$ci_lower | postit_r$ci_upper < full_r$ci_lower)
      )
    }
  }
}

comparison_summary <- bind_rows(comparison_rows)

if (nrow(comparison_summary) > 0) {
  write.csv(comparison_summary, file.path(POSTIT_CONFIG$output_dir, "PostIT_Comparison_Summary.csv"),
            row.names = FALSE)
  cat("Saved: PostIT_Comparison_Summary.csv\n")

  # Print comparison
  cat("\nFull-Sample vs Post-IT Comparison (z-test for equality):\n")
  cat(sprintf("  %15s  h  %8s  %8s  %8s  %8s  %5s  %8s\n",
              "Variable", "Full", "PostIT", "z-stat", "z-p", "Overlap", "Status"))
  cat(paste(rep("-", 90), collapse = ""), "\n")

  for (i in 1:nrow(comparison_summary)) {
    r <- comparison_summary[i, ]
    status <- ifelse(r$z_p_equality < 0.05, "DIFFER", "EQUAL")
    cat(sprintf("  %15s  %2d  %8.2f  %8.2f  %8.2f  %8.3f  %5s  %8s\n",
                r$variable, r$horizon, r$full_coef, r$postit_coef,
                r$z_stat_equality, r$z_p_equality,
                ifelse(r$ci_overlap, "Yes", "No"), status))
  }
}

# ---- F6: Summary Dashboard (4-panel) ----
cat("\nCreating summary dashboard...\n")

# Panel A: Credit LP overlay (simplified)
p178a <- if (!is.null(fullsample_credit)) {
  overlay_total <- overlay_credit %>% filter(credit_type_std == "Cred_Total")
  ggplot(overlay_total, aes(x = horizon, y = coef, color = sample)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = sample), alpha = 0.12) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Full Sample" = "#377EB8", "Post-IT" = "#E41A1C")) +
    scale_fill_manual(values = c("Full Sample" = "#377EB8", "Post-IT" = "#E41A1C")) +
    labs(title = "Credit LP: Total", x = NULL, y = "Effect (pp)", color = NULL, fill = NULL) +
    theme_minimal() + theme(legend.position = "bottom")
} else {
  ggplot() + theme_void() + labs(title = "Credit LP (N/A)")
}

# Panel B: PCA loadings
p178b <- p170a + labs(title = "PCA Loadings Comparison")

# Panel C: Post-IT FCI time series
p178c <- p170b + labs(title = "Post-IT FCI Indices")

# Panel D: Granger p-values
if (nrow(granger_results) > 0) {
  granger_plot_data <- granger_results %>%
    mutate(
      sig = ifelse(p_value < 0.01, "***",
                   ifelse(p_value < 0.05, "**",
                          ifelse(p_value < 0.10, "*", "n.s."))),
      label = paste0("F=", sprintf("%.1f", F_stat), " ", sig)
    )

  p178d <- ggplot(granger_plot_data, aes(x = factor(lags), y = direction)) +
    geom_tile(aes(fill = -log10(p_value + 0.001)), color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient(low = "white", high = "#E41A1C",
                        name = "-log10(p)") +
    labs(title = "Granger Causality (Post-IT)", x = "Lags", y = NULL) +
    theme_minimal()
} else {
  p178d <- ggplot() + theme_void() + labs(title = "Granger Causality (N/A)")
}

p178 <- grid.arrange(
  p178a, p178b, p178c, p178d,
  ncol = 2,
  top = grid::textGrob("Post-IT Subsample Analysis Summary",
                        gp = grid::gpar(fontsize = 16, fontface = "bold"))
)

ggsave(file.path(POSTIT_CONFIG$output_dir, "178_PostIT_Summary_Dashboard.png"), p178,
       width = 14, height = 12, dpi = 300)
cat("Saved: 178_PostIT_Summary_Dashboard.png\n")


################################################################################
# 10. FINAL SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("POST-IT SUBSAMPLE ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Sample: ", format(POSTIT_CONFIG$IT_START, "%Y-%m"), " to ",
    format(max(datos$fecha), "%Y-%m"), "\n")
cat("Observations (raw):", nrow(datos), "\n")
cat("Observations (LP data):", nrow(lp_data), "\n\n")

cat("Method convergence:\n")
for (idx_name in names(methods_status)) {
  ms <- methods_status[[idx_name]]
  cat(sprintf("  %s: %d/4 (%s)\n", idx_name, ms$n_converged,
              paste(ms$converged, collapse = ", ")))
}

cat("\nKey results (Post-IT):\n")
for (cred in credit_types) {
  r12 <- all_credit_lp %>% filter(credit_type == cred, horizon == 12)
  if (nrow(r12) > 0) {
    stars <- ifelse(r12$p_value < 0.01, "***", ifelse(r12$p_value < 0.05, "**",
                                                       ifelse(r12$p_value < 0.10, "*", "")))
    cat(sprintf("  FCI_exCredit -> %s (h=12): %.2f pp (p=%.3f) %s\n",
                credit_labels[cred], r12$coef, r12$p_value, stars))
  }
}

cat("\nOutput files:\n")
cat("  PNG: 170-178 (9 charts)\n")
cat("  CSV: PostIT_FCI_Diagnostics.csv\n")
cat("       PostIT_Granger.csv\n")
cat("       PostIT_LP_Credit.csv\n")
cat("       PostIT_LP_Output.csv\n")
cat("       PostIT_LP_Asymmetric.csv\n")
cat("       PostIT_Comparison_Summary.csv\n")

cat("\n[DONE] Post-IT Subsample Analysis\n\n")
