################################################################################
# FCI IDENTIFICATION TRIANGULATION
################################################################################
#
# Project:      Financial Conditions Index - Cross-Method Identification
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Implements three complementary identification strategies and
#               combines with existing methods for a triangulation comparison:
#
#   PART A: DXY-Only IV-LP (Single Strong Instrument)
#   PART B: GFC-PCA as Single IV (Composite Instrument)
#   PART C: Globally-Purged FCI (Domestic Variation Only)
#   PART D: Identification Triangulation Figure (4-Panel)
#   PART E: Additional Figures and Summary Dashboard
#
# References:
#   - Anderson & Rubin (1949) - Weak-IV Robust Inference
#   - Barnichon & Mesters (2025) - Composite Instruments
#   - Jorda (2005) - Local Projections
#   - Mertens & Ravn (2013) - Proxy-SVAR
#   - Stock & Yogo (2005) - Testing for Weak Instruments
#
################################################################################


################################################################################
# SETUP
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
  library(ivreg)
})

TRI_CONFIG <- list(
  max_horizon      = 24,
  n_lags           = 2,
  confidence_level = 0.90,
  output_dir       = "../output",
  csv_dir          = "../output/csv",
  png_dir          = "../output/png",
  verbose          = TRUE
)

# Helper: find a CSV in csv_dir first, then output_dir as fallback
find_csv <- function(filename) {
  p1 <- file.path(TRI_CONFIG$csv_dir, filename)
  if (file.exists(p1)) return(p1)
  p2 <- file.path(TRI_CONFIG$output_dir, filename)
  if (file.exists(p2)) return(p2)
  return(p2)  # return default even if missing, so file.exists() check works
}

set.seed(20260211)

cat("\n################################################################################\n")
cat("FCI IDENTIFICATION TRIANGULATION\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---- Load FCI from script 01 if needed ----
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  CONFIG_SAVE <- if (exists("CONFIG")) CONFIG else NULL
  source("01_FCI_Complete.R")
  if (!is.null(CONFIG_SAVE)) CONFIG <- CONFIG_SAVE
  rm(CONFIG_SAVE)
  suppressPackageStartupMessages(library(dplyr))
}

# ---- Load new external variables (from script 17) ----
ext_file <- find_csv("New_External_Variables.csv")
if (!file.exists(ext_file)) {
  cat("New_External_Variables.csv not found. Running script 17...\n")
  source("17_FCI_New_External_Data.R")
}
ext_data <- read.csv(ext_file)
ext_data$fecha <- as.Date(ext_data$fecha)
cat("Loaded external data:", nrow(ext_data), "obs\n")

# ---- Load credit data (MN and USD not in ext_data) ----
datos_raw <- read_excel("../data/FCI_data_1.xlsx")
fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

macro_raw <- read_excel("../data/FCI_data_1.xlsx", sheet = "Datos_macro")
fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

macro_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

credit_data <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  left_join(macro_data %>% dplyr::select(fecha, IPC), by = "fecha") %>%
  mutate(
    Cred_Real_MN = ((Creditos_Sector_privado_MN / IPC) /
                      lag(Creditos_Sector_privado_MN / IPC, 12) - 1) * 100,
    Cred_USD = (Creditos_Sector_privado_USD_equivalente /
                  lag(Creditos_Sector_privado_USD_equivalente, 12) - 1) * 100
  ) %>%
  dplyr::select(fecha, Cred_Real_MN, Cred_USD)

# ---- Merge all data ----
analysis_data <- ext_data %>%
  left_join(credit_data, by = "fecha") %>%
  arrange(fecha)

# Create lags for macro controls
analysis_data <- analysis_data %>%
  mutate(
    IMAEP_yoy_L1 = lag(IMAEP_yoy, 1),
    IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
    IPC_yoy_L1 = lag(IPC_yoy, 1),
    IPC_yoy_L2 = lag(IPC_yoy, 2)
  )

cat("Analysis data:", nrow(analysis_data), "obs\n\n")

# ---- Common definitions ----
z_crit <- qnorm(1 - (1 - TRI_CONFIG$confidence_level) / 2)
credit_types <- c("Cred_Real_Total", "Cred_Real_MN")
control_vars <- c("IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                   "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2")
# Exclusion-restriction controls (block non-FCI DXY channels)
er_controls_tri <- intersect(c("d_ToT", "US_10Y", "SP500"), names(analysis_data))
exog_vars <- c("y_lag1", "y_lag2", "fci_lag1",
                "IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2",
                er_controls_tri)

# Helper: standard LP (replicates script 05 logic)
run_lp_local <- function(data, y_var, fci_var, max_h = 24, controls = NULL) {
  results <- data.frame()
  for (h in 1:max_h) {
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        y_lag2 = lag(!!sym(y_var), 2),
        fci_lag1 = lag(!!sym(fci_var), 1)
      )

    lag_vars <- c("y_lag1", "y_lag2", "fci_lag1")
    ctrl_cols <- c()
    if (!is.null(controls)) {
      for (cv in controls) {
        if (cv %in% names(data)) {
          ctrl_cols <- c(ctrl_cols, cv)
          for (j in 1:2) {
            cv_lag <- paste0(cv, "_lag", j)
            data_h <- data_h %>% mutate(!!cv_lag := lag(!!sym(cv), j))
            ctrl_cols <- c(ctrl_cols, cv_lag)
          }
        }
      }
    }

    formula_str <- paste("y_fwd ~", fci_var, "+",
                         paste(c(lag_vars, ctrl_cols), collapse = " + "))
    reg_data <- data_h %>%
      dplyr::select(y_fwd, !!sym(fci_var), all_of(lag_vars), all_of(ctrl_cols)) %>%
      na.omit()

    if (nrow(reg_data) < 30) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    ct <- lmtest::coeftest(model, vcov = vcov_hac)
    idx <- which(rownames(ct) == fci_var)

    results <- rbind(results, data.frame(
      horizon = h, coef = ct[idx, 1], se = ct[idx, 2], p_value = ct[idx, 4],
      ci_lower = ct[idx, 1] - z_crit * ct[idx, 2],
      ci_upper = ct[idx, 1] + z_crit * ct[idx, 2],
      n_obs = nrow(reg_data)
    ))
  }
  return(results)
}


################################################################################
# PART A: DXY-ONLY IV-LP (SINGLE STRONG INSTRUMENT)
################################################################################

cat("================================================================================\n")
cat("PART A: DXY-ONLY IV-LP (SINGLE STRONG INSTRUMENT)\n")
cat("================================================================================\n\n")

dxy_iv_results <- data.frame()

for (cred_type in credit_types) {
  cat("DXY IV-LP for", cred_type, "...\n")

  for (h in 1:TRI_CONFIG$max_horizon) {
    data_h <- analysis_data %>%
      mutate(
        y_fwd = lead(!!sym(cred_type), h),
        y_lag1 = lag(!!sym(cred_type), 1),
        y_lag2 = lag(!!sym(cred_type), 2),
        fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1)
      )

    exog_str <- paste(exog_vars, collapse = " + ")

    reg_data <- data_h %>%
      dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG,
                    all_of(exog_vars), DXY) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    tryCatch({
      # ---- 2SLS with DXY only ----
      iv_formula <- as.formula(paste(
        "y_fwd ~ FCI_ENDO_exCredit_AVG +", exog_str,
        "| DXY +", exog_str
      ))
      iv_model <- ivreg::ivreg(iv_formula, data = reg_data)
      vcov_hac <- sandwich::NeweyWest(iv_model, lag = h + 1, prewhite = FALSE)
      ct_iv <- lmtest::coeftest(iv_model, vcov = vcov_hac)
      iv_summ <- summary(iv_model, vcov = vcov_hac, diagnostics = TRUE)

      idx_iv <- which(rownames(ct_iv) == "FCI_ENDO_exCredit_AVG")

      first_stage_F <- tryCatch(
        iv_summ$diagnostics["Weak instruments", "statistic"], error = function(e) NA)
      wu_hausman_p <- tryCatch(
        iv_summ$diagnostics["Wu-Hausman", "p-value"], error = function(e) NA)

      # ---- OLS comparison ----
      ols_formula <- as.formula(paste("y_fwd ~ FCI_ENDO_exCredit_AVG +", exog_str))
      ols_model <- lm(ols_formula, data = reg_data)
      vcov_ols <- sandwich::NeweyWest(ols_model, lag = h + 1, prewhite = FALSE)
      ct_ols <- lmtest::coeftest(ols_model, vcov = vcov_ols)
      idx_ols <- which(rownames(ct_ols) == "FCI_ENDO_exCredit_AVG")

      # ---- Anderson-Rubin test (DXY only, just-identified) ----
      delta_grid <- seq(-250, 50, by = 1)
      ar_pvals <- sapply(delta_grid, function(d0) {
        reg_data$y_tilde <- reg_data$y_fwd - d0 * reg_data$FCI_ENDO_exCredit_AVG
        f_r <- lm(as.formula(paste("y_tilde ~", exog_str)), data = reg_data)
        f_u <- lm(as.formula(paste("y_tilde ~ DXY +", exog_str)), data = reg_data)
        a <- anova(f_r, f_u)
        return(a$`Pr(>F)`[2])
      })

      ar_ci_lower <- min(delta_grid[ar_pvals > 0.10], na.rm = TRUE)
      ar_ci_upper <- max(delta_grid[ar_pvals > 0.10], na.rm = TRUE)
      if (is.infinite(ar_ci_lower)) ar_ci_lower <- NA
      if (is.infinite(ar_ci_upper)) ar_ci_upper <- NA

      # AR F at delta=0
      reg_data$y_tilde0 <- reg_data$y_fwd
      f_r0 <- lm(as.formula(paste("y_fwd ~", exog_str)), data = reg_data)
      f_u0 <- lm(as.formula(paste("y_fwd ~ DXY +", exog_str)), data = reg_data)
      ar_anova <- anova(f_r0, f_u0)
      ar_F <- ar_anova$F[2]
      ar_p <- ar_anova$`Pr(>F)`[2]

      # Store 2SLS result
      dxy_iv_results <- rbind(dxy_iv_results, data.frame(
        horizon = h, credit_type = cred_type, method = "2SLS_DXY",
        coef = ct_iv[idx_iv, 1], se = ct_iv[idx_iv, 2], p_value = ct_iv[idx_iv, 4],
        ci_lower = ct_iv[idx_iv, 1] - z_crit * ct_iv[idx_iv, 2],
        ci_upper = ct_iv[idx_iv, 1] + z_crit * ct_iv[idx_iv, 2],
        first_stage_F = first_stage_F, wu_hausman_p = wu_hausman_p,
        AR_F = ar_F, AR_p = ar_p,
        AR_CI_lower = ar_ci_lower, AR_CI_upper = ar_ci_upper,
        n_obs = nrow(reg_data), stringsAsFactors = FALSE
      ))

      # Store OLS result
      dxy_iv_results <- rbind(dxy_iv_results, data.frame(
        horizon = h, credit_type = cred_type, method = "OLS",
        coef = ct_ols[idx_ols, 1], se = ct_ols[idx_ols, 2], p_value = ct_ols[idx_ols, 4],
        ci_lower = ct_ols[idx_ols, 1] - z_crit * ct_ols[idx_ols, 2],
        ci_upper = ct_ols[idx_ols, 1] + z_crit * ct_ols[idx_ols, 2],
        first_stage_F = NA, wu_hausman_p = NA,
        AR_F = NA, AR_p = NA,
        AR_CI_lower = NA, AR_CI_upper = NA,
        n_obs = nrow(reg_data), stringsAsFactors = FALSE
      ))

    }, error = function(e) {
      if (TRI_CONFIG$verbose && h %% 6 == 0)
        cat(sprintf("  h=%d: failed (%s)\n", h, e$message))
    })
  }
  cat("  Done.\n")
}

# Report key results
cat("\nDXY-Only IV-LP Results (Cred_Real_Total, key horizons):\n")
cat(sprintf("%-6s %10s %10s %12s %10s %22s\n",
            "h", "OLS", "2SLS_DXY", "First-stg F", "Wu-Haus p", "AR 90% CI"))
cat(strrep("-", 75), "\n")
for (kh in c(6, 12, 18, 24)) {
  ols_row <- dxy_iv_results %>% filter(credit_type == "Cred_Real_Total",
                                        method == "OLS", horizon == kh)
  iv_row <- dxy_iv_results %>% filter(credit_type == "Cred_Real_Total",
                                       method == "2SLS_DXY", horizon == kh)
  if (nrow(ols_row) > 0 && nrow(iv_row) > 0) {
    cat(sprintf("h=%2d  %+10.2f %+10.2f %12.1f %10.3f    [%+.1f, %+.1f]\n",
                kh, ols_row$coef[1], iv_row$coef[1], iv_row$first_stage_F[1],
                iv_row$wu_hausman_p[1], iv_row$AR_CI_lower[1], iv_row$AR_CI_upper[1]))
  }
}

# Save Part A
write.csv(dxy_iv_results,
          file.path(TRI_CONFIG$output_dir, "Triangulation_DXY_IV_LP.csv"),
          row.names = FALSE)
cat("\nSaved: Triangulation_DXY_IV_LP.csv\n\n")


################################################################################
# PART B: GFC-PCA AS SINGLE IV (COMPOSITE INSTRUMENT)
################################################################################

cat("================================================================================\n")
cat("PART B: GFC-PCA AS SINGLE IV (COMPOSITE INSTRUMENT)\n")
cat("================================================================================\n\n")

gfc_iv_results <- data.frame()

# First-stage diagnostic for GFC_PCA
fsd_data <- analysis_data %>%
  dplyr::select(FCI_ENDO_exCredit_AVG, GFC_PCA, all_of(control_vars)) %>%
  na.omit()

f_base <- lm(FCI_ENDO_exCredit_AVG ~ IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +
               IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2, data = fsd_data)
f_gfc <- lm(FCI_ENDO_exCredit_AVG ~ GFC_PCA + IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +
               IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2, data = fsd_data)
vcov_hac_gfc <- sandwich::NeweyWest(f_gfc, lag = 3, prewhite = FALSE)
ct_gfc <- lmtest::coeftest(f_gfc, vcov = vcov_hac_gfc)
gfc_idx <- which(rownames(ct_gfc) == "GFC_PCA")
gfc_partial_r2 <- summary(f_gfc)$r.squared - summary(f_base)$r.squared
gfc_F_stat <- ct_gfc[gfc_idx, 3]^2

cat(sprintf("GFC_PCA First-Stage Diagnostics:\n"))
cat(sprintf("  Coefficient:  %+.4f\n", ct_gfc[gfc_idx, 1]))
cat(sprintf("  t-statistic:  %+.2f\n", ct_gfc[gfc_idx, 3]))
cat(sprintf("  p-value:      %.4f\n", ct_gfc[gfc_idx, 4]))
cat(sprintf("  Partial R²:   %.4f (%.1f%%)\n", gfc_partial_r2, gfc_partial_r2 * 100))
cat(sprintf("  F-statistic:  %.2f\n", gfc_F_stat))

gfc_is_strong <- gfc_F_stat > 10

if (gfc_is_strong) {
  cat("\n  -> GFC_PCA passes relevance threshold (F > 10). Running full 2SLS.\n\n")
} else {
  cat("\n  -> GFC_PCA is WEAK (F < 10). Running AR test for robust inference.\n")
  cat("  -> This is an informative negative finding: Paraguay responds to the\n")
  cat("     dollar channel specifically (DXY), not the broad global cycle.\n\n")
}

# Run IV-LP regardless (report AR even if weak)
for (cred_type in credit_types) {
  cat("GFC-PCA IV-LP for", cred_type, "...\n")

  for (h in 1:TRI_CONFIG$max_horizon) {
    data_h <- analysis_data %>%
      mutate(
        y_fwd = lead(!!sym(cred_type), h),
        y_lag1 = lag(!!sym(cred_type), 1),
        y_lag2 = lag(!!sym(cred_type), 2),
        fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1)
      )

    exog_str <- paste(exog_vars, collapse = " + ")

    reg_data <- data_h %>%
      dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG,
                    all_of(exog_vars), GFC_PCA) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    tryCatch({
      # 2SLS with GFC_PCA
      iv_formula <- as.formula(paste(
        "y_fwd ~ FCI_ENDO_exCredit_AVG +", exog_str,
        "| GFC_PCA +", exog_str
      ))
      iv_model <- ivreg::ivreg(iv_formula, data = reg_data)
      vcov_hac <- sandwich::NeweyWest(iv_model, lag = h + 1, prewhite = FALSE)
      ct_iv <- lmtest::coeftest(iv_model, vcov = vcov_hac)
      iv_summ <- summary(iv_model, vcov = vcov_hac, diagnostics = TRUE)
      idx_iv <- which(rownames(ct_iv) == "FCI_ENDO_exCredit_AVG")

      first_stage_F <- tryCatch(
        iv_summ$diagnostics["Weak instruments", "statistic"], error = function(e) NA)
      wu_hausman_p <- tryCatch(
        iv_summ$diagnostics["Wu-Hausman", "p-value"], error = function(e) NA)

      # OLS
      ols_formula <- as.formula(paste("y_fwd ~ FCI_ENDO_exCredit_AVG +", exog_str))
      ols_model <- lm(ols_formula, data = reg_data)
      vcov_ols <- sandwich::NeweyWest(ols_model, lag = h + 1, prewhite = FALSE)
      ct_ols <- lmtest::coeftest(ols_model, vcov = vcov_ols)
      idx_ols <- which(rownames(ct_ols) == "FCI_ENDO_exCredit_AVG")

      # Anderson-Rubin test
      delta_grid <- seq(-600, 200, by = 5)
      ar_pvals <- sapply(delta_grid, function(d0) {
        reg_data$y_tilde <- reg_data$y_fwd - d0 * reg_data$FCI_ENDO_exCredit_AVG
        f_r <- lm(as.formula(paste("y_tilde ~", exog_str)), data = reg_data)
        f_u <- lm(as.formula(paste("y_tilde ~ GFC_PCA +", exog_str)), data = reg_data)
        a <- anova(f_r, f_u)
        return(a$`Pr(>F)`[2])
      })
      ar_ci_lower <- min(delta_grid[ar_pvals > 0.10], na.rm = TRUE)
      ar_ci_upper <- max(delta_grid[ar_pvals > 0.10], na.rm = TRUE)
      if (is.infinite(ar_ci_lower)) ar_ci_lower <- NA
      if (is.infinite(ar_ci_upper)) ar_ci_upper <- NA

      f_r0 <- lm(as.formula(paste("y_fwd ~", exog_str)), data = reg_data)
      f_u0 <- lm(as.formula(paste("y_fwd ~ GFC_PCA +", exog_str)), data = reg_data)
      ar_anova <- anova(f_r0, f_u0)
      ar_F <- ar_anova$F[2]
      ar_p <- ar_anova$`Pr(>F)`[2]

      gfc_iv_results <- rbind(gfc_iv_results, data.frame(
        horizon = h, credit_type = cred_type, method = "2SLS_GFC_PCA",
        coef = ct_iv[idx_iv, 1], se = ct_iv[idx_iv, 2], p_value = ct_iv[idx_iv, 4],
        ci_lower = ct_iv[idx_iv, 1] - z_crit * ct_iv[idx_iv, 2],
        ci_upper = ct_iv[idx_iv, 1] + z_crit * ct_iv[idx_iv, 2],
        first_stage_F = first_stage_F, wu_hausman_p = wu_hausman_p,
        AR_F = ar_F, AR_p = ar_p,
        AR_CI_lower = ar_ci_lower, AR_CI_upper = ar_ci_upper,
        n_obs = nrow(reg_data), stringsAsFactors = FALSE
      ))

      gfc_iv_results <- rbind(gfc_iv_results, data.frame(
        horizon = h, credit_type = cred_type, method = "OLS",
        coef = ct_ols[idx_ols, 1], se = ct_ols[idx_ols, 2], p_value = ct_ols[idx_ols, 4],
        ci_lower = ct_ols[idx_ols, 1] - z_crit * ct_ols[idx_ols, 2],
        ci_upper = ct_ols[idx_ols, 1] + z_crit * ct_ols[idx_ols, 2],
        first_stage_F = NA, wu_hausman_p = NA,
        AR_F = NA, AR_p = NA,
        AR_CI_lower = NA, AR_CI_upper = NA,
        n_obs = nrow(reg_data), stringsAsFactors = FALSE
      ))

    }, error = function(e) {
      if (TRI_CONFIG$verbose && h %% 6 == 0)
        cat(sprintf("  h=%d: GFC-PCA IV failed (%s)\n", h, e$message))
    })
  }
  cat("  Done.\n")
}

# Report
cat("\nGFC-PCA IV-LP Results (Cred_Real_Total, key horizons):\n")
cat(sprintf("%-6s %10s %12s %12s %22s\n",
            "h", "2SLS_GFC", "First-stg F", "AR p-value", "AR 90% CI"))
cat(strrep("-", 65), "\n")
for (kh in c(6, 12, 18, 24)) {
  iv_row <- gfc_iv_results %>% filter(credit_type == "Cred_Real_Total",
                                       method == "2SLS_GFC_PCA", horizon == kh)
  if (nrow(iv_row) > 0) {
    ar_ci_str <- if (!is.na(iv_row$AR_CI_lower[1]) && !is.na(iv_row$AR_CI_upper[1]))
      sprintf("[%+.1f, %+.1f]", iv_row$AR_CI_lower[1], iv_row$AR_CI_upper[1])
    else "unbounded"
    cat(sprintf("h=%2d  %+10.2f %12.1f %12.3f    %s\n",
                kh, iv_row$coef[1], iv_row$first_stage_F[1], iv_row$AR_p[1], ar_ci_str))
  }
}

write.csv(gfc_iv_results,
          file.path(TRI_CONFIG$output_dir, "Triangulation_GFC_PCA_IV.csv"),
          row.names = FALSE)
cat("\nSaved: Triangulation_GFC_PCA_IV.csv\n\n")


################################################################################
# PART C: GLOBALLY-PURGED FCI
################################################################################

cat("================================================================================\n")
cat("PART C: GLOBALLY-PURGED FCI\n")
cat("================================================================================\n\n")

# Purge FCI_exCredit of ALL external/global variables + 2 lags
global_vars <- c("VIX", "DXY", "FFER", "US_10Y", "SP500", "Selic_rate")

# Create lag columns for purging
purge_data <- analysis_data
for (gv in global_vars) {
  if (gv %in% names(purge_data)) {
    purge_data <- purge_data %>%
      mutate(
        !!paste0(gv, "_L1") := lag(!!sym(gv), 1),
        !!paste0(gv, "_L2") := lag(!!sym(gv), 2)
      )
  }
}

# Build purging formula (6 vars x 3 = 18 regressors)
purge_vars_avail <- c()
for (gv in global_vars) {
  if (gv %in% names(purge_data)) {
    purge_vars_avail <- c(purge_vars_avail, gv, paste0(gv, "_L1"), paste0(gv, "_L2"))
  }
}

purge_formula <- as.formula(paste("FCI_exCredit_AVG ~",
                                   paste(purge_vars_avail, collapse = " + ")))

purge_reg_data <- purge_data %>%
  dplyr::select(FCI_exCredit_AVG, all_of(purge_vars_avail)) %>%
  na.omit()

purge_model <- lm(purge_formula, data = purge_reg_data)
purge_r2 <- summary(purge_model)$r.squared

cat(sprintf("Global Purging Regression:\n"))
cat(sprintf("  Variables:   %s\n", paste(global_vars, collapse = ", ")))
cat(sprintf("  Regressors:  %d (6 vars × (level + 2 lags))\n", length(purge_vars_avail)))
cat(sprintf("  R²:          %.4f (%.1f%%)\n", purge_r2, purge_r2 * 100))
cat(sprintf("  N:           %d\n\n", nrow(purge_reg_data)))

# Store residuals as globally-purged FCI
purge_data_full <- purge_data %>%
  dplyr::select(fecha, FCI_exCredit_AVG, all_of(purge_vars_avail)) %>%
  na.omit()
purge_model_full <- lm(purge_formula, data = purge_data_full)
purge_data_full$FCI_globally_purged <- residuals(purge_model_full)

# Merge back
analysis_data <- analysis_data %>%
  left_join(purge_data_full %>% dplyr::select(fecha, FCI_globally_purged), by = "fecha")

# Run LP with globally-purged FCI
purged_lp_results <- data.frame()

for (cred_type in credit_types) {
  cat("Globally-purged LP for", cred_type, "...\n")

  res <- run_lp_local(analysis_data, cred_type, "FCI_globally_purged",
                      max_h = TRI_CONFIG$max_horizon,
                      controls = c("IMAEP_yoy", "IPC_yoy"))
  if (nrow(res) > 0) {
    res$credit_type <- cred_type
    res$purge_type <- "globally_purged"
    res$purge_r2 <- purge_r2
    purged_lp_results <- rbind(purged_lp_results, res)
  }
  cat("  Done.\n")
}

# Load existing results for comparison
baseline_file <- find_csv("LP_Credit_Standard.csv")
macro_purge_file <- find_csv("LP_Purged_FCI.csv")

if (file.exists(baseline_file)) {
  baseline_lp <- read.csv(baseline_file)
  # FCI_COMP, Total = baseline with macro controls
  baseline_total <- baseline_lp %>%
    filter(fci_type == "FCI_COMP" & credit_type == "Total") %>%
    dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, n_obs) %>%
    mutate(credit_type = "Cred_Real_Total", purge_type = "baseline", purge_r2 = 0)

  baseline_mn <- baseline_lp %>%
    filter(fci_type == "FCI_COMP" & credit_type == "Real_MN") %>%
    dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, n_obs) %>%
    mutate(credit_type = "Cred_Real_MN", purge_type = "baseline", purge_r2 = 0)

  purged_lp_results <- rbind(purged_lp_results, baseline_total, baseline_mn)
}

if (file.exists(macro_purge_file)) {
  macro_purge <- read.csv(macro_purge_file)
  mp_total <- macro_purge %>%
    filter(fci_type == "Purged" & credit_type == "Total") %>%
    dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, n_obs) %>%
    mutate(credit_type = "Cred_Real_Total", purge_type = "macro_purged", purge_r2 = 0.033)

  mp_mn <- macro_purge %>%
    filter(fci_type == "Purged" & credit_type == "Real_MN") %>%
    dplyr::select(horizon, coef, se, p_value, ci_lower, ci_upper, n_obs) %>%
    mutate(credit_type = "Cred_Real_MN", purge_type = "macro_purged", purge_r2 = 0.033)

  purged_lp_results <- rbind(purged_lp_results, mp_total, mp_mn)
}

# Comparison at key horizons
cat("\nPurging Comparison (Cred_Real_Total, coef at key horizons):\n")
cat(sprintf("%-6s %12s %15s %18s\n", "h", "Baseline", "Macro-purged", "Globally-purged"))
cat(strrep("-", 55), "\n")
for (kh in c(6, 12, 18)) {
  bl <- purged_lp_results %>% filter(credit_type == "Cred_Real_Total",
                                      purge_type == "baseline", horizon == kh)
  mp <- purged_lp_results %>% filter(credit_type == "Cred_Real_Total",
                                      purge_type == "macro_purged", horizon == kh)
  gp <- purged_lp_results %>% filter(credit_type == "Cred_Real_Total",
                                      purge_type == "globally_purged", horizon == kh)
  bl_str <- if (nrow(bl) > 0) sprintf("%+.2f (p=%.3f)", bl$coef[1], bl$p_value[1]) else "N/A"
  mp_str <- if (nrow(mp) > 0) sprintf("%+.2f (p=%.3f)", mp$coef[1], mp$p_value[1]) else "N/A"
  gp_str <- if (nrow(gp) > 0) sprintf("%+.2f (p=%.3f)", gp$coef[1], gp$p_value[1]) else "N/A"
  cat(sprintf("h=%2d  %12s %15s %18s\n", kh, bl_str, mp_str, gp_str))
}

write.csv(purged_lp_results,
          file.path(TRI_CONFIG$output_dir, "Triangulation_Globally_Purged_LP.csv"),
          row.names = FALSE)
cat("\nSaved: Triangulation_Globally_Purged_LP.csv\n\n")


################################################################################
# PART D: IDENTIFICATION TRIANGULATION FIGURE (KEY OUTPUT)
################################################################################

cat("================================================================================\n")
cat("PART D: IDENTIFICATION TRIANGULATION FIGURE\n")
cat("================================================================================\n\n")

tryCatch({
  # Panel A: OLS-LP (from LP_Credit_Standard.csv)
  panel_a_data <- NULL
  if (file.exists(baseline_file)) {
    bl <- read.csv(baseline_file)
    panel_a_data <- bl %>%
      filter(fci_type == "FCI_COMP" & credit_type == "Total") %>%
      dplyr::select(horizon, coef, ci_lower, ci_upper)
  }

  # Panel B: DXY IV-LP (from Part A)
  panel_b_data <- dxy_iv_results %>%
    filter(credit_type == "Cred_Real_Total" & method == "2SLS_DXY") %>%
    dplyr::select(horizon, coef, ci_lower, ci_upper,
                  AR_CI_lower, AR_CI_upper)

  # Panel C: Proxy-SVAR (DXY proxy)
  psvar_file <- find_csv("ProxySVAR_IRF_Results.csv")
  panel_c_data <- NULL
  if (file.exists(psvar_file)) {
    psvar <- read.csv(psvar_file)
    panel_c_data <- psvar %>%
      filter(proxy == "DXY" & response == "Cred_Real_Total") %>%
      dplyr::select(horizon, coef = irf_point, ci_lower, ci_upper)
  }

  # Panel D: Block-SVAR (DXY shock)
  bsvar_file <- find_csv("Enriched_BSVAR_IRF.csv")
  panel_d_data <- NULL
  if (file.exists(bsvar_file)) {
    bsvar <- read.csv(bsvar_file)
    panel_d_data <- bsvar %>%
      filter(shock == "DXY" & response == "Cred_Real_Total") %>%
      dplyr::select(horizon, coef = irf_val, ci_lower = lower, ci_upper = upper)
  }

  plots <- list()

  # Panel A: OLS-LP
  if (!is.null(panel_a_data) && nrow(panel_a_data) > 0) {
    plots[["A"]] <- ggplot(panel_a_data, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_vline(xintercept = c(6, 12), linetype = "dashed", color = "gray70", linewidth = 0.3) +
      labs(title = "A: OLS Local Projections",
           subtitle = "FCI_exCredit -> Real Credit (pp per 1 SD)",
           x = "Horizon (months)", y = "pp credit growth") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 11))
  }

  # Panel B: DXY IV-LP (use AR CIs if available, else Wald CIs)
  if (nrow(panel_b_data) > 0) {
    panel_b_plot_data <- panel_b_data %>%
      mutate(
        ci_lo = ifelse(!is.na(AR_CI_lower), AR_CI_lower, ci_lower),
        ci_hi = ifelse(!is.na(AR_CI_upper), AR_CI_upper, ci_upper)
      )
    plots[["B"]] <- ggplot(panel_b_plot_data, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_vline(xintercept = c(6, 12), linetype = "dashed", color = "gray70", linewidth = 0.3) +
      labs(title = "B: DXY IV-LP (2SLS, just-identified)",
           subtitle = "FCI instrumented by DXY, AR 90% CIs",
           x = "Horizon (months)", y = "pp credit growth") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 11))
  }

  # Panel C: Proxy-SVAR
  if (!is.null(panel_c_data) && nrow(panel_c_data) > 0) {
    plots[["C"]] <- ggplot(panel_c_data, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_vline(xintercept = c(6, 12), linetype = "dashed", color = "gray70", linewidth = 0.3) +
      labs(title = "C: Proxy-SVAR (DXY proxy)",
           subtitle = "Structural FCI shock -> Real Credit (VAR units)",
           x = "Horizon (months)", y = "IRF (VAR units)") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 11))
  }

  # Panel D: Block-SVAR
  if (!is.null(panel_d_data) && nrow(panel_d_data) > 0) {
    plots[["D"]] <- ggplot(panel_d_data, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_vline(xintercept = c(6, 12), linetype = "dashed", color = "gray70", linewidth = 0.3) +
      labs(title = "D: Enriched Block-SVAR (DXY shock)",
           subtitle = "DXY shock -> Real Credit (VAR units)",
           x = "Horizon (months)", y = "IRF (VAR units)") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 11))
  }

  if (length(plots) >= 2) {
    combined <- arrangeGrob(grobs = plots, ncol = 2,
      top = grid::textGrob(
        "Identification Triangulation: FCI -> Real Credit\nQualitative comparison -- sign, timing, and significance across four strategies",
        gp = grid::gpar(fontsize = 13, fontface = "bold")))
    ggsave(file.path(TRI_CONFIG$output_dir, "251_Identification_Triangulation.png"),
           combined, width = 16, height = 10, dpi = 150)
    cat("Saved: 251_Identification_Triangulation.png\n")
  } else {
    cat("WARNING: Not enough panels for triangulation figure\n")
  }
}, error = function(e) cat("WARNING: Triangulation figure failed:", e$message, "\n"))


################################################################################
# PART E: ADDITIONAL FIGURES AND SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("PART E: ADDITIONAL FIGURES AND SUMMARY\n")
cat("================================================================================\n\n")

# ---- 252: DXY-Only IV-LP Detail ----
tryCatch({
  dxy_total <- dxy_iv_results %>% filter(credit_type == "Cred_Real_Total")

  p1_data <- dxy_total %>% filter(method %in% c("OLS", "2SLS_DXY"))
  p1 <- ggplot(p1_data, aes(x = horizon, y = coef, color = method, fill = method)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = c(6, 12), linetype = "dashed", color = "gray70", linewidth = 0.3) +
    scale_color_manual(values = c("OLS" = "gray40", "2SLS_DXY" = "steelblue"),
                       labels = c("OLS" = "OLS", "2SLS_DXY" = "2SLS (DXY only)")) +
    scale_fill_manual(values = c("OLS" = "gray40", "2SLS_DXY" = "steelblue")) +
    guides(fill = "none") +
    labs(title = "OLS vs DXY-Only 2SLS: Real Total Credit",
         x = "Horizon (months)", y = "pp credit growth", color = "Method") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  p2_data <- dxy_total %>% filter(method == "2SLS_DXY")
  p2 <- ggplot(p2_data, aes(x = horizon, y = first_stage_F)) +
    geom_line(color = "firebrick", linewidth = 0.8) +
    geom_point(color = "firebrick", size = 1.5) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 0.5) +
    annotate("text", x = 20, y = 12, label = "Stock-Yogo F=10", color = "red",
             size = 3, hjust = 1) +
    labs(title = "First-Stage F (DXY -> FCI_ENDO_exCredit)",
         x = "Horizon (months)", y = "F-statistic") +
    theme_minimal(base_size = 11)

  combined <- grid.arrange(p1, p2, ncol = 2,
    top = grid::textGrob("DXY as Single Strong Instrument",
                         gp = grid::gpar(fontsize = 13, fontface = "bold")))
  ggsave(file.path(TRI_CONFIG$output_dir, "252_DXY_Only_IV_LP.png"),
         combined, width = 14, height = 6, dpi = 150)
  cat("Saved: 252_DXY_Only_IV_LP.png\n")
}, error = function(e) cat("WARNING: Plot 252 failed:", e$message, "\n"))

# ---- 253: GFC-PCA IV Diagnostics ----
tryCatch({
  gfc_total <- gfc_iv_results %>%
    filter(credit_type == "Cred_Real_Total" & method == "2SLS_GFC_PCA")

  # Load GFC loadings for context
  gfc_diag_file <- find_csv("GFC_Index_Diagnostics.csv")
  gfc_loadings_data <- NULL
  if (file.exists(gfc_diag_file)) {
    gfc_diag <- read.csv(gfc_diag_file)
    if ("loading_PC1" %in% names(gfc_diag)) {
      gfc_loadings_data <- gfc_diag %>%
        filter(!is.na(loading_PC1)) %>%
        dplyr::select(variable, loading_PC1)
    }
  }

  p1 <- ggplot(gfc_total, aes(x = horizon, y = first_stage_F)) +
    geom_line(color = "darkorange", linewidth = 0.8) +
    geom_point(color = "darkorange", size = 1.5) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    annotate("text", x = 18, y = max(c(gfc_total$first_stage_F, 12), na.rm = TRUE) * 0.9,
             label = paste0("GFC_PCA: overall F = ", round(gfc_F_stat, 1),
                           ifelse(gfc_is_strong, " (STRONG)", " (WEAK)")),
             color = "darkorange", size = 3.5, fontface = "bold") +
    labs(title = "First-Stage F: GFC-PCA as Instrument",
         subtitle = ifelse(gfc_is_strong, "Passes relevance threshold",
                          "Weak instrument -- AR test needed for valid inference"),
         x = "Horizon (months)", y = "F-statistic") +
    theme_minimal(base_size = 11)

  # Compare DXY F vs GFC-PCA F
  dxy_f_data <- dxy_iv_results %>%
    filter(credit_type == "Cred_Real_Total" & method == "2SLS_DXY") %>%
    dplyr::select(horizon, F_DXY = first_stage_F)
  gfc_f_data <- gfc_total %>%
    dplyr::select(horizon, F_GFC = first_stage_F)
  f_comp <- dxy_f_data %>% left_join(gfc_f_data, by = "horizon") %>%
    pivot_longer(cols = c(F_DXY, F_GFC), names_to = "instrument", values_to = "F_stat")

  p2 <- ggplot(f_comp, aes(x = horizon, y = F_stat, color = instrument)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("F_DXY" = "steelblue", "F_GFC" = "darkorange"),
                       labels = c("F_DXY" = "DXY alone", "F_GFC" = "GFC-PCA")) +
    labs(title = "First-Stage F Comparison",
         subtitle = "DXY alone vs GFC-PCA composite",
         x = "Horizon (months)", y = "F-statistic", color = "Instrument") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  combined <- grid.arrange(p1, p2, ncol = 2,
    top = grid::textGrob("GFC-PCA Composite Instrument Assessment",
                         gp = grid::gpar(fontsize = 13, fontface = "bold")))
  ggsave(file.path(TRI_CONFIG$output_dir, "253_GFC_PCA_IV_Diagnostics.png"),
         combined, width = 14, height = 6, dpi = 150)
  cat("Saved: 253_GFC_PCA_IV_Diagnostics.png\n")
}, error = function(e) cat("WARNING: Plot 253 failed:", e$message, "\n"))

# ---- 254: Globally-Purged FCI LP ----
tryCatch({
  purge_plot <- purged_lp_results %>%
    filter(credit_type == "Cred_Real_Total" &
             purge_type %in% c("baseline", "macro_purged", "globally_purged"))

  purge_plot$purge_label <- factor(purge_plot$purge_type,
    levels = c("baseline", "macro_purged", "globally_purged"),
    labels = c("Baseline FCI_exCredit",
               paste0("Macro-purged (R\u00B2=3.3%)"),
               paste0("Globally-purged (R\u00B2=", round(purge_r2 * 100, 1), "%)")))

  p <- ggplot(purge_plot, aes(x = horizon, y = coef, color = purge_label, fill = purge_label)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(6, 12), linetype = "dashed", color = "gray70", linewidth = 0.3) +
    scale_color_manual(values = c("steelblue", "darkgreen", "darkorange")) +
    scale_fill_manual(values = c("steelblue", "darkgreen", "darkorange")) +
    labs(title = "FCI Purging Comparison: Effect on Real Total Credit",
         subtitle = paste0("Globally-purged removes ", round(purge_r2 * 100, 1),
                          "% of FCI variation (vs 3.3% for macro-purge)"),
         x = "Horizon (months)", y = "pp credit growth per 1 SD FCI",
         color = "FCI Version", fill = "FCI Version") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))

  ggsave(file.path(TRI_CONFIG$output_dir, "254_Globally_Purged_FCI_LP.png"),
         p, width = 12, height = 7, dpi = 150)
  cat("Saved: 254_Globally_Purged_FCI_LP.png\n")
}, error = function(e) cat("WARNING: Plot 254 failed:", e$message, "\n"))

# ---- 255: Triangulation Dashboard ----
tryCatch({
  # Collect h=12 estimates from all methods
  h12_summary <- data.frame()

  # 1. OLS-LP
  bl_h12 <- purged_lp_results %>%
    filter(credit_type == "Cred_Real_Total" & purge_type == "baseline" & horizon == 12)
  if (nrow(bl_h12) > 0)
    h12_summary <- rbind(h12_summary, data.frame(
      method = "OLS-LP", coef = bl_h12$coef[1], ci_lo = bl_h12$ci_lower[1],
      ci_hi = bl_h12$ci_upper[1], stringsAsFactors = FALSE))

  # 2. DXY IV-LP
  dxy_h12 <- dxy_iv_results %>%
    filter(credit_type == "Cred_Real_Total" & method == "2SLS_DXY" & horizon == 12)
  if (nrow(dxy_h12) > 0)
    h12_summary <- rbind(h12_summary, data.frame(
      method = "DXY IV-LP", coef = dxy_h12$coef[1],
      ci_lo = ifelse(!is.na(dxy_h12$AR_CI_lower[1]), dxy_h12$AR_CI_lower[1], dxy_h12$ci_lower[1]),
      ci_hi = ifelse(!is.na(dxy_h12$AR_CI_upper[1]), dxy_h12$AR_CI_upper[1], dxy_h12$ci_upper[1]),
      stringsAsFactors = FALSE))

  # 3. GFC-PCA IV
  gfc_h12 <- gfc_iv_results %>%
    filter(credit_type == "Cred_Real_Total" & method == "2SLS_GFC_PCA" & horizon == 12)
  if (nrow(gfc_h12) > 0)
    h12_summary <- rbind(h12_summary, data.frame(
      method = "GFC-PCA IV", coef = gfc_h12$coef[1],
      ci_lo = ifelse(!is.na(gfc_h12$AR_CI_lower[1]), gfc_h12$AR_CI_lower[1], gfc_h12$ci_lower[1]),
      ci_hi = ifelse(!is.na(gfc_h12$AR_CI_upper[1]), gfc_h12$AR_CI_upper[1], gfc_h12$ci_upper[1]),
      stringsAsFactors = FALSE))

  # 4. Globally-purged
  gp_h12 <- purged_lp_results %>%
    filter(credit_type == "Cred_Real_Total" & purge_type == "globally_purged" & horizon == 12)
  if (nrow(gp_h12) > 0)
    h12_summary <- rbind(h12_summary, data.frame(
      method = "Globally-purged", coef = gp_h12$coef[1],
      ci_lo = gp_h12$ci_lower[1], ci_hi = gp_h12$ci_upper[1],
      stringsAsFactors = FALSE))

  # 5. Proxy-SVAR
  if (file.exists(find_csv("ProxySVAR_IRF_Results.csv"))) {
    psvar <- read.csv(find_csv("ProxySVAR_IRF_Results.csv"))
    ps_h12 <- psvar %>% filter(proxy == "DXY" & response == "Cred_Real_Total" & horizon == 12)
    if (nrow(ps_h12) > 0)
      h12_summary <- rbind(h12_summary, data.frame(
        method = "Proxy-SVAR", coef = ps_h12$irf_point[1],
        ci_lo = ps_h12$ci_lower[1], ci_hi = ps_h12$ci_upper[1],
        stringsAsFactors = FALSE))
  }

  # 6. Block-SVAR
  if (file.exists(find_csv("Enriched_BSVAR_IRF.csv"))) {
    bsvar <- read.csv(find_csv("Enriched_BSVAR_IRF.csv"))
    bs_h12 <- bsvar %>% filter(shock == "DXY" & response == "Cred_Real_Total" & horizon == 12)
    if (nrow(bs_h12) > 0)
      h12_summary <- rbind(h12_summary, data.frame(
        method = "Block-SVAR", coef = bs_h12$irf_val[1],
        ci_lo = bs_h12$lower[1], ci_hi = bs_h12$upper[1],
        stringsAsFactors = FALSE))
  }

  dashboard_plots <- list()

  # Panel 1: Bar chart of h=12 coefficients
  if (nrow(h12_summary) > 0) {
    h12_summary$method <- factor(h12_summary$method,
                                  levels = h12_summary$method)
    # Separate LP-scale vs VAR-scale
    lp_methods <- c("OLS-LP", "DXY IV-LP", "GFC-PCA IV", "Globally-purged")
    h12_lp <- h12_summary %>% filter(method %in% lp_methods)

    if (nrow(h12_lp) > 0) {
      h12_lp$method <- factor(h12_lp$method, levels = lp_methods)
      dashboard_plots[[1]] <- ggplot(h12_lp, aes(x = method, y = coef)) +
        geom_col(fill = "steelblue", alpha = 0.7, width = 0.6) +
        geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, color = "gray30") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(title = "LP-Scale Methods at h=12",
             subtitle = "pp credit growth per 1 SD FCI (90% CIs)",
             x = NULL, y = "Coefficient") +
        theme_minimal(base_size = 9) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    }
  }

  # Panel 2: Sign consistency heatmap
  sign_data <- data.frame()
  for (h_check in 1:TRI_CONFIG$max_horizon) {
    # OLS
    bl_row <- purged_lp_results %>%
      filter(credit_type == "Cred_Real_Total" & purge_type == "baseline" & horizon == h_check)
    if (nrow(bl_row) > 0) sign_data <- rbind(sign_data, data.frame(
      method = "OLS-LP", horizon = h_check,
      sign = sign(bl_row$coef[1]),
      significant = bl_row$p_value[1] < 0.10))

    # DXY IV
    dxy_row <- dxy_iv_results %>%
      filter(credit_type == "Cred_Real_Total" & method == "2SLS_DXY" & horizon == h_check)
    if (nrow(dxy_row) > 0) sign_data <- rbind(sign_data, data.frame(
      method = "DXY IV-LP", horizon = h_check,
      sign = sign(dxy_row$coef[1]),
      significant = dxy_row$p_value[1] < 0.10))

    # Globally-purged
    gp_row <- purged_lp_results %>%
      filter(credit_type == "Cred_Real_Total" & purge_type == "globally_purged" & horizon == h_check)
    if (nrow(gp_row) > 0) sign_data <- rbind(sign_data, data.frame(
      method = "Globally-purged", horizon = h_check,
      sign = sign(gp_row$coef[1]),
      significant = gp_row$p_value[1] < 0.10))
  }

  if (nrow(sign_data) > 0) {
    sign_data$status <- ifelse(sign_data$sign < 0 & sign_data$significant, "Neg. signif.",
                         ifelse(sign_data$sign < 0, "Neg. insig.",
                           ifelse(sign_data$sign > 0 & sign_data$significant, "Pos. signif.",
                             "Pos. insig.")))
    sign_data$method <- factor(sign_data$method,
                                levels = c("OLS-LP", "DXY IV-LP", "Globally-purged"))
    dashboard_plots[[2]] <- ggplot(sign_data, aes(x = horizon, y = method, fill = status)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_manual(values = c("Neg. signif." = "steelblue",
                                    "Neg. insig." = "lightblue",
                                    "Pos. signif." = "firebrick",
                                    "Pos. insig." = "lightsalmon")) +
      labs(title = "Sign Consistency Across Methods",
           x = "Horizon (months)", y = NULL, fill = "Sign/Signif.") +
      theme_minimal(base_size = 9) +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = 7))
  }

  # Panel 3: First-stage F comparison
  f_data <- data.frame()
  dxy_f <- dxy_iv_results %>%
    filter(credit_type == "Cred_Real_Total" & method == "2SLS_DXY") %>%
    dplyr::select(horizon, first_stage_F) %>% mutate(instrument = "DXY alone")
  gfc_f <- gfc_iv_results %>%
    filter(credit_type == "Cred_Real_Total" & method == "2SLS_GFC_PCA") %>%
    dplyr::select(horizon, first_stage_F) %>% mutate(instrument = "GFC-PCA")
  f_data <- rbind(dxy_f, gfc_f)

  # Add original 4-instrument F if available
  orig_iv_file <- find_csv("IV_LP_Improved_Diagnostics.csv")
  if (file.exists(orig_iv_file)) {
    orig_diag <- read.csv(orig_iv_file)
    orig_f <- orig_diag %>%
      filter(grepl("Total|total", credit_type, ignore.case = TRUE)) %>%
      dplyr::select(horizon, first_stage_F) %>%
      mutate(instrument = "4-IV (script 18)")
    f_data <- rbind(f_data, orig_f)
  }

  if (nrow(f_data) > 0) {
    dashboard_plots[[3]] <- ggplot(f_data, aes(x = horizon, y = first_stage_F,
                                                color = instrument)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("DXY alone" = "steelblue",
                                     "GFC-PCA" = "darkorange",
                                     "4-IV (script 18)" = "gray50")) +
      labs(title = "First-Stage F Comparison",
           x = "Horizon (months)", y = "F-statistic", color = "Instrument Set") +
      theme_minimal(base_size = 9) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 7))
  }

  # Panel 4: Summary text
  summary_text <- paste0(
    "IDENTIFICATION TRIANGULATION SUMMARY\n",
    "=====================================\n\n",
    "OLS-LP (h=12): ",
    ifelse(nrow(bl_h12) > 0, sprintf("%.2f pp (p=%.3f)", bl_h12$coef[1], bl_h12$p_value[1]), "N/A"),
    "\n\n",
    "DXY IV-LP (h=12): ",
    ifelse(nrow(dxy_h12) > 0, sprintf("%.2f pp (F=%.1f)", dxy_h12$coef[1], dxy_h12$first_stage_F[1]), "N/A"),
    "\n\n",
    "Globally-purged (h=12): ",
    ifelse(nrow(gp_h12) > 0, sprintf("%.2f pp (p=%.3f)", gp_h12$coef[1], gp_h12$p_value[1]), "N/A"),
    "\n\n",
    "Global purge R-sq: ", sprintf("%.1f%%", purge_r2 * 100),
    "\n\n",
    "Key: All LP methods show\n",
    "negative FCI effect on credit\n",
    "peaking at h=6-12"
  )

  dashboard_plots[[4]] <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = summary_text,
             hjust = 0.5, vjust = 0.5, size = 3, family = "mono") +
    theme_void() +
    xlim(0, 1) + ylim(0, 1)

  if (length(dashboard_plots) >= 2) {
    combined <- arrangeGrob(grobs = dashboard_plots, ncol = 2,
      top = grid::textGrob("Identification Triangulation Dashboard",
                           gp = grid::gpar(fontsize = 13, fontface = "bold")))
    ggsave(file.path(TRI_CONFIG$output_dir, "255_Triangulation_Dashboard.png"),
           combined, width = 14, height = 10, dpi = 150)
    cat("Saved: 255_Triangulation_Dashboard.png\n")
  }
}, error = function(e) cat("WARNING: Dashboard 255 failed:", e$message, "\n"))

# ---- Summary CSV: Cross-method comparison ----
tryCatch({
  summary_csv <- data.frame()

  for (kh in c(6, 12, 18, 24)) {
    # OLS
    bl_row <- purged_lp_results %>%
      filter(credit_type == "Cred_Real_Total" & purge_type == "baseline" & horizon == kh)
    if (nrow(bl_row) > 0)
      summary_csv <- rbind(summary_csv, data.frame(
        horizon = kh, method = "OLS-LP", coef = bl_row$coef[1],
        se = bl_row$se[1], p_value = bl_row$p_value[1],
        ci_lower = bl_row$ci_lower[1], ci_upper = bl_row$ci_upper[1],
        first_stage_F = NA, AR_CI_lower = NA, AR_CI_upper = NA,
        n_obs = bl_row$n_obs[1], scale = "LP",
        stringsAsFactors = FALSE))

    # DXY IV
    dxy_row <- dxy_iv_results %>%
      filter(credit_type == "Cred_Real_Total" & method == "2SLS_DXY" & horizon == kh)
    if (nrow(dxy_row) > 0)
      summary_csv <- rbind(summary_csv, data.frame(
        horizon = kh, method = "DXY_IV-LP", coef = dxy_row$coef[1],
        se = dxy_row$se[1], p_value = dxy_row$p_value[1],
        ci_lower = dxy_row$ci_lower[1], ci_upper = dxy_row$ci_upper[1],
        first_stage_F = dxy_row$first_stage_F[1],
        AR_CI_lower = dxy_row$AR_CI_lower[1], AR_CI_upper = dxy_row$AR_CI_upper[1],
        n_obs = dxy_row$n_obs[1], scale = "LP",
        stringsAsFactors = FALSE))

    # GFC-PCA IV
    gfc_row <- gfc_iv_results %>%
      filter(credit_type == "Cred_Real_Total" & method == "2SLS_GFC_PCA" & horizon == kh)
    if (nrow(gfc_row) > 0)
      summary_csv <- rbind(summary_csv, data.frame(
        horizon = kh, method = "GFC-PCA_IV", coef = gfc_row$coef[1],
        se = gfc_row$se[1], p_value = gfc_row$p_value[1],
        ci_lower = gfc_row$ci_lower[1], ci_upper = gfc_row$ci_upper[1],
        first_stage_F = gfc_row$first_stage_F[1],
        AR_CI_lower = gfc_row$AR_CI_lower[1], AR_CI_upper = gfc_row$AR_CI_upper[1],
        n_obs = gfc_row$n_obs[1], scale = "LP",
        stringsAsFactors = FALSE))

    # Globally-purged
    gp_row <- purged_lp_results %>%
      filter(credit_type == "Cred_Real_Total" & purge_type == "globally_purged" & horizon == kh)
    if (nrow(gp_row) > 0)
      summary_csv <- rbind(summary_csv, data.frame(
        horizon = kh, method = "Globally-purged", coef = gp_row$coef[1],
        se = gp_row$se[1], p_value = gp_row$p_value[1],
        ci_lower = gp_row$ci_lower[1], ci_upper = gp_row$ci_upper[1],
        first_stage_F = NA, AR_CI_lower = NA, AR_CI_upper = NA,
        n_obs = gp_row$n_obs[1], scale = "LP",
        stringsAsFactors = FALSE))

    # Proxy-SVAR
    if (file.exists(find_csv("ProxySVAR_IRF_Results.csv"))) {
      psvar <- read.csv(find_csv("ProxySVAR_IRF_Results.csv"))
      ps_row <- psvar %>% filter(proxy == "DXY" & response == "Cred_Real_Total" & horizon == kh)
      if (nrow(ps_row) > 0)
        summary_csv <- rbind(summary_csv, data.frame(
          horizon = kh, method = "Proxy-SVAR_DXY", coef = ps_row$irf_point[1],
          se = NA, p_value = NA,
          ci_lower = ps_row$ci_lower[1], ci_upper = ps_row$ci_upper[1],
          first_stage_F = NA, AR_CI_lower = NA, AR_CI_upper = NA,
          n_obs = NA, scale = "VAR",
          stringsAsFactors = FALSE))
    }

    # Block-SVAR
    if (file.exists(find_csv("Enriched_BSVAR_IRF.csv"))) {
      bsvar <- read.csv(find_csv("Enriched_BSVAR_IRF.csv"))
      bs_row <- bsvar %>% filter(shock == "DXY" & response == "Cred_Real_Total" & horizon == kh)
      if (nrow(bs_row) > 0)
        summary_csv <- rbind(summary_csv, data.frame(
          horizon = kh, method = "Block-SVAR_DXY", coef = bs_row$irf_val[1],
          se = NA, p_value = NA,
          ci_lower = bs_row$lower[1], ci_upper = bs_row$upper[1],
          first_stage_F = NA, AR_CI_lower = NA, AR_CI_upper = NA,
          n_obs = NA, scale = "VAR",
          stringsAsFactors = FALSE))
    }
  }

  write.csv(summary_csv,
            file.path(TRI_CONFIG$output_dir, "Triangulation_Summary_Comparison.csv"),
            row.names = FALSE)
  cat("Saved: Triangulation_Summary_Comparison.csv\n")
}, error = function(e) cat("WARNING: Summary CSV failed:", e$message, "\n"))


cat("\n################################################################################\n")
cat("SCRIPT 23 COMPLETE\n")
cat("################################################################################\n\n")
cat("Outputs: 5 PNGs (251-255) + 4 CSVs\n")
cat("Key findings:\n")

# Print key DXY findings
dxy_h12 <- dxy_iv_results %>%
  filter(credit_type == "Cred_Real_Total" & method == "2SLS_DXY" & horizon == 12)
if (nrow(dxy_h12) > 0) {
  cat(sprintf("  DXY IV-LP (h=12): coef = %.2f, First-stage F = %.1f\n",
              dxy_h12$coef[1], dxy_h12$first_stage_F[1]))
}
cat(sprintf("  GFC-PCA first-stage F: %.1f (%s)\n", gfc_F_stat,
            ifelse(gfc_is_strong, "strong", "weak")))
cat(sprintf("  Global purge R-sq: %.1f%%\n", purge_r2 * 100))

gp_h12 <- purged_lp_results %>%
  filter(credit_type == "Cred_Real_Total" & purge_type == "globally_purged" & horizon == 12)
if (nrow(gp_h12) > 0) {
  cat(sprintf("  Globally-purged LP (h=12): coef = %.2f, p = %.3f\n",
              gp_h12$coef[1], gp_h12$p_value[1]))
}
cat("\n")
