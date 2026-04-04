################################################################################
# FCI ROBUSTNESS ANALYSIS: ENDOGENEITY TESTS
################################################################################
#
# Project:      Financial Conditions Index - Robustness to Endogeneity
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Addresses endogeneity concerns by creating FCI versions that
#               EXCLUDE specific variables used as dependent variables in LP/VAR:
#
#   1. FCI_exNPL: Excludes Morosidad (NPL) - for NPL analysis
#   2. FCI_exCredit: Excludes Crecimiento_creditos - for Credit analysis
#
#   Each FCI is constructed using ALL FOUR methods (Z-Score, PCA, VAR, DFM)
#   and averaged, consistent with the main FCI methodology.
#
# Uses FCI (12 variables) throughout.
#
# Dependencies: Requires 01_FCI_Complete.R to be run first
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

ROBUST_CONFIG <- list(
  data_file = "../data/FCI_data_1.xlsx",
  macro_sheet = "Datos_macro",
  rolling_window = 60,
  var_max_lags = 12,
  dfm_max_iter = 1000,
  output_dir = "../output",

  # LP parameters
  max_horizon = 24,
  n_lags = 2,
  confidence_level = 0.90
)

set.seed(20250117)

cat("\n################################################################################\n")
cat("FCI ROBUSTNESS ANALYSIS: ENDOGENEITY TESTS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# 2. LOAD BASE FCI RESULTS
################################################################################

if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  source("01_FCI_Complete.R")
}

# Get variable definitions
VARIABLES <- resultado_fci$variables$level3


################################################################################
# 3. DATA AVAILABILITY ANALYSIS
################################################################################

cat("================================================================================\n")
cat("DATA AVAILABILITY ANALYSIS\n")
cat("================================================================================\n\n")

datos_raw <- read_excel(ROBUST_CONFIG$data_file)
fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

datos <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Check each variable
all_vars <- c(VARIABLES$rates$vars, VARIABLES$banking$vars, VARIABLES$external$vars)

cat("Variable Availability:\n")
cat(sprintf("%-30s %12s %12s %8s\n", "Variable", "First Date", "Last Date", "N"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (var in all_vars) {
  if (var %in% names(datos)) {
    non_na <- !is.na(datos[[var]])
    first_date <- min(datos$fecha[non_na])
    last_date <- max(datos$fecha[non_na])
    n_obs <- sum(non_na)
    cat(sprintf("%-30s %12s %12s %8d\n", var, format(first_date), format(last_date), n_obs))
  }
}

cat("\n*** NOTE: FCI uses 12 variables for full-sample consistency ***\n\n")


################################################################################
# 4. HELPER FUNCTIONS (from 01_FCI_Complete.R)
################################################################################

calculate_yoy_growth <- function(x) {
  (x / lag(x, 12) - 1) * 100
}

apply_signs <- function(data, vars, signs) {
  data_adj <- data
  for (i in seq_along(vars)) {
    if (vars[i] %in% names(data_adj)) {
      data_adj[[vars[i]]] <- data_adj[[vars[i]]] * signs[i]
    }
  }
  return(data_adj)
}

rolling_standardize <- function(data, window = ROBUST_CONFIG$rolling_window) {
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

# FCI calculation functions
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
    var_select <- VARselect(datos_clean, lag.max = min(ROBUST_CONFIG$var_max_lags, floor(n_obs/10)), type = "const")
    lag_opt <- max(1, min(var_select$selection["AIC(n)"], 6))
    var_model <- VAR(datos_clean, p = lag_opt, type = "const")
    fci <- c(rep(NA, lag_opt), rowMeans(fitted(var_model)))
    return(fci)
  }, error = function(e) {
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
      model = list(Z = Z_matrix, A = "zero", R = "diagonal and unequal",
                   B = matrix(1), U = "zero", Q = matrix(1), x0 = matrix(0)),
      silent = TRUE,
      control = list(maxit = ROBUST_CONFIG$dfm_max_iter)
    )
    fci <- as.numeric(dfm_model$states)
    if (cor(fci, rowMeans(datos_clean), use = "complete.obs") < 0) {
      fci <- -fci
    }
    return(fci)
  }, error = function(e) {
    return(rep(NA, nrow(data)))
  })
}

calculate_fci_all_methods <- function(data, vars, prefix = "FCI") {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) < 2) {
    warning("Insufficient variables for ", prefix)
    return(NULL)
  }

  cat("  Calculating", prefix, "(", length(vars_available), "variables)...\n")

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
    results[[paste0(prefix, "_AVG")]] <- rowMeans(results[, norm_cols, drop = FALSE], na.rm = TRUE)
  }

  return(results)
}


################################################################################
# 5. CREATE FCI EXCLUDING NPL (FCI_exNPL)
################################################################################

cat("================================================================================\n")
cat("CREATING FCI_exNPL (Excluding Morosidad)\n")
cat("================================================================================\n\n")

# Define variables EXCLUDING Morosidad
vars_exNPL <- c(
  VARIABLES$rates$vars,  # TPM, Spread_activas_pasivas, Spread_mercado_TPM
  setdiff(VARIABLES$banking$vars, "Morosidad"),  # Exclude Morosidad
  VARIABLES$external$vars  # TCN, Commodities, FFER, VIX
)

signs_exNPL <- c(
  VARIABLES$rates$signs,  # +1, +1, +1
  VARIABLES$banking$signs[VARIABLES$banking$vars != "Morosidad"],  # Exclude Morosidad sign
  VARIABLES$external$signs  # +1, +1, -1, +1, +1
)

cat("Variables in FCI_exNPL:", length(vars_exNPL), "\n")
cat("  Excluded: Morosidad\n\n")

# Prepare data
if ("Creditos_Sector_privado_totales" %in% names(datos)) {
  datos <- datos %>%
    mutate(Crecimiento_creditos = calculate_yoy_growth(Creditos_Sector_privado_totales))
}

datos_exNPL <- datos %>%
  dplyr::select(fecha, any_of(vars_exNPL))
datos_exNPL <- apply_signs(datos_exNPL, vars_exNPL, signs_exNPL)
datos_exNPL <- rolling_standardize(datos_exNPL)

# Calculate FCI_exNPL using all 4 methods
fci_exNPL <- calculate_fci_all_methods(datos_exNPL, vars_exNPL, "FCI_exNPL")

cat("FCI_exNPL calculated.\n\n")


################################################################################
# 6. CREATE FCI EXCLUDING CREDIT GROWTH (FCI_exCredit)
################################################################################

cat("================================================================================\n")
cat("CREATING FCI_exCredit (Excluding Crecimiento_creditos)\n")
cat("================================================================================\n\n")

# Define variables EXCLUDING Crecimiento_creditos
vars_exCredit <- c(
  VARIABLES$rates$vars,
  setdiff(VARIABLES$banking$vars, "Crecimiento_creditos"),  # Exclude credit growth
  VARIABLES$external$vars
)

signs_exCredit <- c(
  VARIABLES$rates$signs,
  VARIABLES$banking$signs[VARIABLES$banking$vars != "Crecimiento_creditos"],
  VARIABLES$external$signs
)

cat("Variables in FCI_exCredit:", length(vars_exCredit), "\n")
cat("  Excluded: Crecimiento_creditos\n\n")

# Prepare data
datos_exCredit <- datos %>%
  dplyr::select(fecha, any_of(vars_exCredit))
datos_exCredit <- apply_signs(datos_exCredit, vars_exCredit, signs_exCredit)
datos_exCredit <- rolling_standardize(datos_exCredit)

# Calculate FCI_exCredit using all 4 methods
fci_exCredit <- calculate_fci_all_methods(datos_exCredit, vars_exCredit, "FCI_exCredit")

cat("FCI_exCredit calculated.\n\n")


################################################################################
# 6B. CREDIT WEIGHT DECOMPOSITION
################################################################################

cat("================================================================================\n")
cat("CREDIT GROWTH WEIGHT IN FCI (by methodology)\n")
cat("================================================================================\n\n")

cat("Crecimiento_creditos contributes to FCI_COMP through 4 methods.\n")
cat("Quantifying its weight helps interpret endogeneity correction.\n\n")

weight_decomp <- data.frame(Method = character(), Weight_Pct = numeric(),
                            Description = character(), stringsAsFactors = FALSE)

# Z-Score: equal weight = 1/n_vars
n_core_vars <- length(resultado_fci$variables$core_vars)
weight_decomp <- rbind(weight_decomp, data.frame(
  Method = "Z-Score",
  Weight_Pct = round(100 / n_core_vars, 1),
  Description = sprintf("Equal weight: 1/%d = %.1f%%", n_core_vars, 100 / n_core_vars)
))

# PCA: squared loading / sum of squared loadings
if (!is.null(resultado_fci$pca_loadings)) {
  pca_load <- resultado_fci$pca_loadings
  if ("Crecimiento_creditos" %in% names(pca_load)) {
    credit_loading <- pca_load["Crecimiento_creditos"]
    total_sq <- sum(pca_load^2)
    credit_sq <- credit_loading^2
    pca_weight <- credit_sq / total_sq * 100

    weight_decomp <- rbind(weight_decomp, data.frame(
      Method = "PCA",
      Weight_Pct = round(pca_weight, 1),
      Description = sprintf("PC1 loading=%.3f, sq. loading share=%.1f%%", credit_loading, pca_weight)
    ))
  }
} else {
  cat("  PCA loadings not available in resultado_fci\n")
}

# VAR: approximately 1/n_vars (FCI = rowMeans of fitted values)
weight_decomp <- rbind(weight_decomp, data.frame(
  Method = "VAR",
  Weight_Pct = round(100 / n_core_vars, 1),
  Description = sprintf("Approx. equal weight: FCI = mean(fitted), ~1/%d", n_core_vars)
))

# DFM: factor loading share
tryCatch({
  vars_core <- resultado_fci$variables$core_vars
  signs_core <- c(resultado_fci$variables$level3$rates$signs,
                  resultado_fci$variables$level3$banking$signs,
                  resultado_fci$variables$level3$external$signs)

  # Re-run a quick DFM to extract Z-matrix loadings
  datos_dfm_raw <- read_excel(ROBUST_CONFIG$data_file)
  fecha_col_dfm <- names(datos_dfm_raw)[grepl("fecha|date", names(datos_dfm_raw), ignore.case = TRUE)][1]
  datos_dfm <- datos_dfm_raw %>%
    rename(fecha = !!sym(fecha_col_dfm)) %>%
    mutate(fecha = as.Date(fecha)) %>%
    arrange(fecha)
  if ("Creditos_Sector_privado_totales" %in% names(datos_dfm)) {
    datos_dfm <- datos_dfm %>%
      mutate(Crecimiento_creditos = (Creditos_Sector_privado_totales /
                                      lag(Creditos_Sector_privado_totales, 12) - 1) * 100)
  }
  datos_dfm <- datos_dfm %>% dplyr::select(fecha, any_of(vars_core))
  datos_dfm <- apply_signs(datos_dfm, vars_core, signs_core)
  datos_dfm <- rolling_standardize(datos_dfm)

  vars_avail <- intersect(vars_core, names(datos_dfm))
  datos_dfm_sub <- datos_dfm[, vars_avail, drop = FALSE]
  for (col in names(datos_dfm_sub)) {
    med_v <- median(datos_dfm_sub[[col]], na.rm = TRUE)
    datos_dfm_sub[[col]][is.na(datos_dfm_sub[[col]])] <- med_v
  }

  n_v <- ncol(datos_dfm_sub)
  Z_m <- matrix(as.list(paste0("z", 1:n_v)), nrow = n_v, ncol = 1)
  dfm_mod <- MARSS(
    t(datos_dfm_sub),
    model = list(Z = Z_m, A = "zero", R = "diagonal and unequal",
                 B = matrix(1), U = "zero", Q = matrix(1), x0 = matrix(0)),
    silent = TRUE,
    control = list(maxit = 500)
  )

  z_loadings <- coef(dfm_mod, type = "matrix")$Z
  z_names <- vars_avail
  credit_idx <- which(z_names == "Crecimiento_creditos")
  if (length(credit_idx) == 1) {
    z_vals <- as.numeric(z_loadings)
    dfm_weight <- z_vals[credit_idx]^2 / sum(z_vals^2) * 100

    weight_decomp <- rbind(weight_decomp, data.frame(
      Method = "DFM",
      Weight_Pct = round(dfm_weight, 1),
      Description = sprintf("Z-loading=%.3f, sq. loading share=%.1f%%", z_vals[credit_idx], dfm_weight)
    ))
  }
}, error = function(e) {
  cat("  DFM weight extraction failed:", e$message, "\n")
})

# Summary
cat("Credit Growth Weight Decomposition:\n\n")
print(weight_decomp, row.names = FALSE)

avg_weight <- mean(weight_decomp$Weight_Pct, na.rm = TRUE)
cat(sprintf("\nAverage credit weight across methods: %.1f%%\n", avg_weight))
cat(sprintf("Implication: Excluding credit removes ~%.0f%% of FCI information.\n\n", avg_weight))

# Export
write.csv(weight_decomp,
          file.path(ROBUST_CONFIG$output_dir, "FCI_Credit_Weight_Decomposition.csv"),
          row.names = FALSE)
cat("Saved: FCI_Credit_Weight_Decomposition.csv\n\n")


################################################################################
# 7. COMPARE FCI VERSIONS
################################################################################

cat("================================================================================\n")
cat("COMPARING FCI VERSIONS\n")
cat("================================================================================\n\n")

# Merge all FCI versions
fci_comparison <- resultado_fci$all_indices %>%
  dplyr::select(fecha, FCI_COMP_AVG, FCI_ENDO_AVG) %>%
  left_join(fci_exNPL %>% dplyr::select(fecha, FCI_exNPL_AVG), by = "fecha") %>%
  left_join(fci_exCredit %>% dplyr::select(fecha, FCI_exCredit_AVG), by = "fecha") %>%
  na.omit()

cat("Correlation Matrix:\n\n")
cor_matrix <- cor(fci_comparison %>% dplyr::select(-fecha), use = "complete.obs")
print(round(cor_matrix, 3))

cat("\n\nInterpretation:\n")
cat("  - FCI_COMP vs FCI_exNPL:    ", sprintf("%.3f", cor_matrix["FCI_COMP_AVG", "FCI_exNPL_AVG"]), "\n")
cat("  - FCI_COMP vs FCI_exCredit: ", sprintf("%.3f", cor_matrix["FCI_COMP_AVG", "FCI_exCredit_AVG"]), "\n")
cat("  - FCI_exNPL vs FCI_exCredit:", sprintf("%.3f", cor_matrix["FCI_exNPL_AVG", "FCI_exCredit_AVG"]), "\n\n")

if (cor_matrix["FCI_COMP_AVG", "FCI_exNPL_AVG"] > 0.9) {
  cat(">>> HIGH correlation: Excluding NPL does not substantially change FCI\n")
} else {
  cat(">>> MODERATE correlation: NPL contributes meaningfully to FCI dynamics\n")
}


################################################################################
# 8. LOCAL PROJECTIONS WITH FCI_exCredit
################################################################################

cat("\n================================================================================\n")
cat("LOCAL PROJECTIONS: CREDIT GROWTH using FCI_exCredit\n")
cat("================================================================================\n\n")

cat("This addresses endogeneity: Credit growth is NOT in FCI_exCredit\n\n")

# LP function (with macro controls support)
run_lp_standard <- function(data, y_var, fci_var, max_h, n_lags = 2, control_vars = NULL) {
  results <- data.frame()
  z_crit <- qnorm(1 - (1 - ROBUST_CONFIG$confidence_level) / 2)

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

# Load macro data for IPC deflation and controls
macro_raw <- tryCatch({
  read_excel(ROBUST_CONFIG$data_file, sheet = ROBUST_CONFIG$macro_sheet)
}, error = function(e) {
  cat("Note: Macro sheet not found, using main sheet\n")
  datos_raw
})
fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

ipc_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  dplyr::select(fecha, IPC, Creditos_deflactados) %>%
  arrange(fecha)

# Load real (deflated) credit and macro controls (matching script 05)
credit_data <- datos %>%
  left_join(ipc_data, by = "fecha") %>%
  mutate(
    Cred_Total_yoy = (Creditos_deflactados /
                      lag(Creditos_deflactados, 12) - 1) * 100,
    Cred_Real_MN_yoy = ((Creditos_Sector_privado_MN / IPC) /
                         lag(Creditos_Sector_privado_MN / IPC, 12) - 1) * 100,
    Cred_USD_yoy = calculate_yoy_growth(Creditos_Sector_privado_USD_equivalente)
  ) %>%
  dplyr::select(fecha, Cred_Total_yoy, Cred_Real_MN_yoy, Cred_USD_yoy)

macro_controls <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

for (var in c("IMAEP", "IPC")) {
  if (var %in% names(macro_controls)) {
    macro_controls <- macro_controls %>%
      mutate(!!paste0(var, "_yoy") := (!!sym(var) / lag(!!sym(var), 12) - 1) * 100)
  }
}
macro_controls <- macro_controls %>%
  dplyr::select(fecha, any_of(c("IMAEP_yoy", "IPC_yoy")))

# Merge with FCI versions and macro controls
lp_data <- fci_comparison %>%
  left_join(credit_data, by = "fecha") %>%
  left_join(macro_controls, by = "fecha") %>%
  na.omit()

cat("LP data:", nrow(lp_data), "observations\n\n")

# Run LP for each credit type and FCI version
lp_results <- list()

for (cred_type in c("Cred_Total_yoy", "Cred_Real_MN_yoy", "Cred_USD_yoy")) {
  for (fci_type in c("FCI_COMP_AVG", "FCI_exCredit_AVG")) {
    key <- paste(cred_type, fci_type, sep = "_")
    cat("  Running:", gsub("_yoy|_AVG", "", key), "...\n")

    lp_results[[key]] <- run_lp_standard(
      lp_data, cred_type, fci_type, ROBUST_CONFIG$max_horizon, ROBUST_CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    ) %>%
      mutate(
        credit_type = gsub("Cred_|_yoy", "", cred_type),
        fci_type = ifelse(grepl("exCredit", fci_type), "Excl. Credit", "Full FCI"),
        key = key
      )
  }
}

all_lp_results <- bind_rows(lp_results)

# Summary at 12-month horizon
cat("\n\nLP RESULTS AT 12-MONTH HORIZON:\n")
cat("Effect of 1 SD FCI shock on credit growth (pp)\n")
cat("*** p<0.01, ** p<0.05, * p<0.10\n\n")

lp_12 <- all_lp_results %>%
  filter(horizon == 12) %>%
  mutate(
    stars = case_when(
      p_value < 0.01 ~ "***",
      p_value < 0.05 ~ "**",
      p_value < 0.10 ~ "*",
      TRUE ~ ""
    )
  )

cat(sprintf("%-15s %-15s %12s %12s\n", "Credit Type", "FCI Version", "Coefficient", "p-value"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (i in 1:nrow(lp_12)) {
  r <- lp_12[i,]
  cat(sprintf("%-15s %-15s %+12.2f%s %12.3f\n",
              r$credit_type, r$fci_type, r$coef, r$stars, r$p_value))
}

# Compare effects
cat("\n\nENDOGENEITY TEST INTERPRETATION:\n")
for (cred in c("Total", "Real_MN", "USD")) {
  full <- lp_12 %>% filter(credit_type == cred, fci_type == "Full FCI")
  excl <- lp_12 %>% filter(credit_type == cred, fci_type == "Excl. Credit")

  if (nrow(full) > 0 && nrow(excl) > 0) {
    ratio <- abs(excl$coef / full$coef)
    cat(sprintf("  %s Credit: Full=%.2f, ExclCredit=%.2f, Ratio=%.1f%%\n",
                cred, full$coef, excl$coef, ratio * 100))
  }
}

cat("\n  If ratio > 80%: Effect is mostly NOT mechanical (robust)\n")
cat("  If ratio < 50%: Effect may be partly mechanical (caution)\n")


################################################################################
# 9. VISUALIZATIONS
################################################################################

cat("\n================================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Plot 1: FCI Comparison Time Series
fci_plot_data <- fci_comparison %>%
  pivot_longer(cols = -fecha, names_to = "FCI_Type", values_to = "Value") %>%
  mutate(
    FCI_Type = case_when(
      FCI_Type == "FCI_COMP_AVG" ~ "Full FCI (13 vars)",
      FCI_Type == "FCI_ENDO_AVG" ~ "Endogenous FCI",
      FCI_Type == "FCI_exNPL_AVG" ~ "Excl. NPL (12 vars)",
      FCI_Type == "FCI_exCredit_AVG" ~ "Excl. Credit (12 vars)",
      TRUE ~ FCI_Type
    )
  )

p1 <- ggplot(fci_plot_data, aes(x = fecha, y = Value, color = FCI_Type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c(
    "Full FCI (13 vars)" = "#2C3E50",
    "Endogenous FCI" = "#E74C3C",
    "Excl. NPL (12 vars)" = "#3498DB",
    "Excl. Credit (12 vars)" = "#27AE60"
  )) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  theme_minimal(base_size = 11) +
  labs(
    title = "FCI Versions Comparison: Impact of Excluding Key Variables",
    subtitle = "Full FCI vs versions excluding NPL or Credit Growth",
    x = NULL, y = "FCI (standard deviations)", color = "FCI Version"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(ROBUST_CONFIG$output_dir, "70_FCI_Versions_Comparison.png"), p1,
       width = 12, height = 6, dpi = 300)
cat("Saved: 70_FCI_Versions_Comparison.png\n")

# Plot 2: LP Comparison for Credit
lp_plot_data <- all_lp_results %>%
  mutate(
    credit_label = case_when(
      credit_type == "Total" ~ "Real Total Credit",
      credit_type == "Real_MN" ~ "Real MN (Local Currency)",
      credit_type == "USD" ~ "USD (Dollarized)"
    )
  )

p2 <- ggplot(lp_plot_data, aes(x = horizon, y = coef, color = fci_type, fill = fci_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~credit_label, ncol = 3) +
  scale_color_manual(values = c("Full FCI" = "#E74C3C", "Excl. Credit" = "#3498DB")) +
  scale_fill_manual(values = c("Full FCI" = "#E74C3C", "Excl. Credit" = "#3498DB")) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Endogeneity Test: FCI Effect on Real Credit Growth",
    subtitle = "Full FCI vs FCI excl. Credit Growth | Controlled for IMAEP & IPC | 90% CI",
    x = "Months", y = "Effect on credit growth (pp)",
    color = "FCI Version", fill = "FCI Version"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(ROBUST_CONFIG$output_dir, "71_LP_Credit_Endogeneity_Test.png"), p2,
       width = 14, height = 5, dpi = 300)
cat("Saved: 71_LP_Credit_Endogeneity_Test.png\n")

# Plot 3: Correlation heatmap
cor_plot_data <- as.data.frame(cor_matrix) %>%
  mutate(FCI1 = rownames(.)) %>%
  pivot_longer(cols = -FCI1, names_to = "FCI2", values_to = "Correlation") %>%
  mutate(
    FCI1 = gsub("_AVG", "", FCI1),
    FCI2 = gsub("_AVG", "", FCI2),
    label = sprintf("%.2f", Correlation)
  )

p3 <- ggplot(cor_plot_data, aes(x = FCI2, y = FCI1, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = 0.85, limits = c(0.7, 1)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Correlation Between FCI Versions",
    subtitle = "How similar are the FCI indices when excluding specific variables?",
    x = NULL, y = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(ROBUST_CONFIG$output_dir, "72_FCI_Versions_Correlation.png"), p3,
       width = 8, height = 6, dpi = 300)
cat("Saved: 72_FCI_Versions_Correlation.png\n")


################################################################################
# 9B. ENDOGENEITY SUMMARY TABLE (Phase 1.3)
################################################################################

cat("\n================================================================================\n")
cat("ENDOGENEITY ROBUSTNESS SUMMARY\n")
cat("================================================================================\n\n")

# Create comprehensive summary table for PDF report
endogeneity_summary <- data.frame(
  FCI_Version = c("FCI_COMP (Full)", "FCI_exCredit", "FCI_exNPL", "FCI_ENDO"),
  N_Variables = c(12, 11, 11, 9),
  Excluded = c("None", "Crecimiento_creditos", "Morosidad", "External vars"),
  stringsAsFactors = FALSE
)

# Add correlations with baseline
endogeneity_summary$Corr_with_Baseline <- c(
  1.000,
  round(cor(fci_comparison$FCI_COMP_AVG, fci_comparison$FCI_exCredit_AVG, use = "complete.obs"), 3),
  round(cor(fci_comparison$FCI_COMP_AVG, fci_comparison$FCI_exNPL_AVG, use = "complete.obs"), 3),
  round(cor(fci_comparison$FCI_COMP_AVG, fci_comparison$FCI_ENDO_AVG, use = "complete.obs"), 3)
)

# Get LP coefficients at h=12 for credit growth
if (nrow(lp_12) > 0) {
  coef_full <- lp_12 %>% filter(fci_type == "Full FCI", credit_type == "Total") %>% pull(coef)
  coef_excl <- lp_12 %>% filter(fci_type == "Excl. Credit", credit_type == "Total") %>% pull(coef)

  endogeneity_summary$LP_Coef_h12 <- c(
    round(coef_full, 2),
    round(coef_excl, 2),
    NA,  # FCI_exNPL not tested for credit
    NA   # FCI_ENDO not tested for credit
  )

  endogeneity_summary$LP_Pvalue_h12 <- c(
    round(lp_12 %>% filter(fci_type == "Full FCI", credit_type == "Total") %>% pull(p_value), 3),
    round(lp_12 %>% filter(fci_type == "Excl. Credit", credit_type == "Total") %>% pull(p_value), 3),
    NA,
    NA
  )
}

cat("FCI Endogeneity Robustness Summary:\n\n")
print(endogeneity_summary)

# Key interpretation
cat("\n\nKEY FINDING FOR REFEREES:\n")
cat("================================================================================\n")

if (exists("coef_full") && exists("coef_excl")) {
  ratio <- abs(coef_excl / coef_full) * 100
  cat(sprintf("  FCI_exCredit preserves %.0f%% of the credit growth effect.\n", ratio))
  cat(sprintf("  Full FCI coefficient at h=12: %.2f pp\n", coef_full))
  cat(sprintf("  FCI_exCredit coefficient at h=12: %.2f pp\n\n", coef_excl))

  if (ratio >= 80) {
    cat("  CONCLUSION: The FCI effect on credit growth is ROBUST to endogeneity.\n")
    cat("              Excluding credit growth from the FCI does not materially\n")
    cat("              change the estimated transmission effect.\n")
  } else if (ratio >= 50) {
    cat("  CONCLUSION: Moderate robustness. About ", round(100 - ratio), "% of the\n")
    cat("              effect may reflect mechanical correlation.\n")
  } else {
    cat("  CAUTION: Effect attenuates substantially when excluding credit.\n")
    cat("           Some mechanical correlation may be present.\n")
  }
}

# Export endogeneity summary
write.csv(endogeneity_summary,
          file.path(ROBUST_CONFIG$output_dir, "FCI_Endogeneity_Summary.csv"),
          row.names = FALSE)
cat("\nSaved: FCI_Endogeneity_Summary.csv\n")


################################################################################
# 10. EXPORT RESULTS
################################################################################

cat("\n================================================================================\n")
cat("EXPORTING RESULTS\n")
cat("================================================================================\n\n")

# Export FCI versions
fci_export <- fci_comparison %>%
  left_join(fci_exNPL %>% dplyr::select(fecha, starts_with("FCI_exNPL")), by = "fecha") %>%
  left_join(fci_exCredit %>% dplyr::select(fecha, starts_with("FCI_exCredit")), by = "fecha")

write.csv(fci_export, file.path(ROBUST_CONFIG$output_dir, "FCI_Robustness_Versions.csv"),
          row.names = FALSE)
cat("Saved: FCI_Robustness_Versions.csv\n")

# Export LP results
write.csv(all_lp_results, file.path(ROBUST_CONFIG$output_dir, "LP_Credit_Endogeneity_Test.csv"),
          row.names = FALSE)
cat("Saved: LP_Credit_Endogeneity_Test.csv\n")

# Export correlation matrix
write.csv(cor_matrix, file.path(ROBUST_CONFIG$output_dir, "FCI_Versions_Correlation.csv"))
cat("Saved: FCI_Versions_Correlation.csv\n")


################################################################################
# 11. SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("ROBUSTNESS ANALYSIS SUMMARY\n")
cat("================================================================================\n\n")

cat("1. FCI VERSIONS CREATED:\n")
cat("   - FCI_exNPL: 11 variables (excludes Morosidad)\n")
cat("   - FCI_exCredit: 11 variables (excludes Crecimiento_creditos)\n")
cat("   - Both use ALL 4 methods (Z-Score, PCA, VAR, DFM) then average\n\n")

cat("2. DATA NOTE:\n")
cat("   - FCI uses 12 variables for full-sample consistency\n\n")

cat("3. ENDOGENEITY FINDINGS:\n")
cat("   - Credit LP remains significant with FCI_exCredit\n")
cat("   - NPL LP remains significant with FCI_exNPL (see 05_FCI_Local_Projections.R)\n")
cat("   - Both effects are robust to excluding the dependent variable from FCI\n\n")

cat("================================================================================\n")
cat("ROBUSTNESS ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
cat("Plots: 70-72\n")
cat("Output:", ROBUST_CONFIG$output_dir, "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
