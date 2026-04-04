################################################################################
# FCI IMPROVED INSTRUMENTAL VARIABLES LOCAL PROJECTIONS
################################################################################
#
# Project:      Financial Conditions Index - Expanded IV Analysis
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Uses expanded instrument set (DXY, US10Y, SP500, Selic, ToT)
#               to address weak-IV problem from original IV-LP (first-stage
#               F=0.9-1.6 with VIX/FFER/Commodities). Implements weak-IV
#               robust inference (Anderson-Rubin test).
#
#   SECTION A: First-Stage Diagnostic Battery
#   SECTION B: Instrument Selection
#   SECTION C: Formal IV-LP with Improved Instruments
#   SECTION D: Weak-IV Robust Inference (Anderson-Rubin)
#   SECTION E: Comparison with Original IV-LP
#
# References:
#   - Stock & Yogo (2005) - Testing for Weak Instruments
#   - Anderson & Rubin (1949) - Weak-IV Robust Inference
#   - Lee et al. (2022) - tF Procedure
#   - Jorda (2005) - Local Projections
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

IV_CONFIG <- list(
  max_horizon     = 24,
  n_lags          = 2,
  confidence_level = 0.90,
  output_dir      = "../output",
  verbose         = TRUE
)

set.seed(20260211)

cat("\n################################################################################\n")
cat("FCI IMPROVED INSTRUMENTAL VARIABLES LOCAL PROJECTIONS\n")
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
ext_file <- file.path(IV_CONFIG$output_dir, "New_External_Variables.csv")
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

# ---- Merge all data (ext_data already has Cred_Real_Total) ----
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


################################################################################
# SECTION A: FIRST-STAGE DIAGNOSTIC BATTERY
################################################################################

cat("================================================================================\n")
cat("SECTION A: FIRST-STAGE DIAGNOSTIC BATTERY\n")
cat("================================================================================\n\n")

# Full instrument set
all_instruments <- c("VIX", "FFER", "DXY", "US_10Y", "SP500", "Selic_rate", "d_ToT")
available_instruments <- intersect(all_instruments, names(analysis_data))
cat("Available instruments:", paste(available_instruments, collapse = ", "), "\n\n")

# Controls (same as LP)
control_vars <- c("IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                   "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2")

# Dependent variable: FCI_ENDO_exCredit_AVG
first_stage_data <- analysis_data %>%
  dplyr::select(FCI_ENDO_exCredit_AVG, all_of(available_instruments), all_of(control_vars)) %>%
  na.omit()

cat("First-stage sample:", nrow(first_stage_data), "obs\n\n")

# ---- Individual instrument diagnostics ----
instrument_diag <- data.frame()

for (inst in available_instruments) {
  # Individual regression
  f_indiv <- as.formula(paste("FCI_ENDO_exCredit_AVG ~", inst, "+",
                               paste(control_vars, collapse = " + ")))
  m_indiv <- lm(f_indiv, data = first_stage_data)
  vcov_hac <- sandwich::NeweyWest(m_indiv, lag = 3, prewhite = FALSE)
  ct <- lmtest::coeftest(m_indiv, vcov = vcov_hac)
  idx <- which(rownames(ct) == inst)

  # Partial R² (incremental from adding this instrument to controls-only model)
  f_base <- as.formula(paste("FCI_ENDO_exCredit_AVG ~", paste(control_vars, collapse = " + ")))
  m_base <- lm(f_base, data = first_stage_data)
  partial_r2 <- summary(m_indiv)$r.squared - summary(m_base)$r.squared

  # Individual F-statistic
  f_stat <- ct[idx, 3]^2  # t² ≈ F for single restriction

  instrument_diag <- rbind(instrument_diag, data.frame(
    Instrument = inst,
    Coefficient = ct[idx, 1],
    Std_Error = ct[idx, 2],
    t_stat = ct[idx, 3],
    p_value = ct[idx, 4],
    Partial_R2 = partial_r2,
    F_stat = f_stat,
    stringsAsFactors = FALSE
  ))
}

cat("Individual Instrument Diagnostics:\n")
cat(sprintf("%-15s %10s %10s %10s %10s\n", "Instrument", "t-stat", "p-value", "Partial R²", "F-stat"))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(instrument_diag))) {
  cat(sprintf("%-15s %+10.2f %10.3f %10.4f %10.2f\n",
              instrument_diag$Instrument[i],
              instrument_diag$t_stat[i],
              instrument_diag$p_value[i],
              instrument_diag$Partial_R2[i],
              instrument_diag$F_stat[i]))
}

# ---- Joint F-statistics ----
# Full set
f_full <- as.formula(paste("FCI_ENDO_exCredit_AVG ~",
                            paste(c(available_instruments, control_vars), collapse = " + ")))
m_full <- lm(f_full, data = first_stage_data)
f_full_restricted <- lm(as.formula(paste("FCI_ENDO_exCredit_AVG ~",
                                          paste(control_vars, collapse = " + "))),
                         data = first_stage_data)
joint_F_full <- anova(f_full_restricted, m_full)$F[2]

cat(sprintf("\nJoint F-statistic (all %d instruments): %.2f\n",
            length(available_instruments), joint_F_full))

# Theory-preferred subset: Selic, DXY, SP500
theory_subset <- intersect(c("Selic_rate", "DXY", "SP500"), available_instruments)
if (length(theory_subset) > 0) {
  f_theory <- as.formula(paste("FCI_ENDO_exCredit_AVG ~",
                                paste(c(theory_subset, control_vars), collapse = " + ")))
  m_theory <- lm(f_theory, data = first_stage_data)
  joint_F_theory <- anova(f_full_restricted, m_theory)$F[2]
  cat(sprintf("Joint F-statistic (theory subset: %s): %.2f\n",
              paste(theory_subset, collapse = ", "), joint_F_theory))
}
cat("\n")


################################################################################
# SECTION B: INSTRUMENT SELECTION
################################################################################

cat("================================================================================\n")
cat("SECTION B: INSTRUMENT SELECTION\n")
cat("================================================================================\n\n")

# Selection criteria: |t| > 2.0 OR partial R² > 0.01 (relaxed since instruments are for SVAR-style)
selected_instruments <- instrument_diag %>%
  filter(abs(t_stat) > 2.0 | Partial_R2 > 0.01) %>%
  pull(Instrument)

cat("Instruments passing selection criteria (|t|>2 or R²>0.01):\n")
if (length(selected_instruments) > 0) {
  cat("  ", paste(selected_instruments, collapse = ", "), "\n\n")
} else {
  cat("  None meet strict criteria. Using theory-preferred subset.\n")
  selected_instruments <- theory_subset
  cat("  Theory-preferred:", paste(selected_instruments, collapse = ", "), "\n\n")
}

# If still empty, use all available
if (length(selected_instruments) == 0) {
  selected_instruments <- available_instruments
  cat("  Fallback: using all instruments\n\n")
}

# Exclusion restriction assessment
cat("Exclusion Restriction Assessment:\n")
er_assessment <- data.frame(
  Instrument = c("Selic_rate", "DXY", "SP500", "US_10Y", "VIX", "FFER", "d_ToT"),
  Validity = c("Strong: regional spillover, no direct PY credit channel",
               "Likely valid: affects PY through FCI, not directly",
               "Likely valid: wealth effect implausible for Paraguay",
               "Possibly valid: may affect long-term investment",
               "Possibly valid: risk appetite channel",
               "Possibly valid: direct dollar credit cost channel",
               "Problematic: affects agricultural income -> credit demand"),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(er_assessment))) {
  if (er_assessment$Instrument[i] %in% available_instruments) {
    cat(sprintf("  %-12s: %s\n", er_assessment$Instrument[i], er_assessment$Validity[i]))
  }
}
cat("\n")


################################################################################
# SECTION C: FORMAL IV-LP WITH IMPROVED INSTRUMENTS
################################################################################

cat("================================================================================\n")
cat("SECTION C: IV-LP WITH IMPROVED INSTRUMENTS\n")
cat("================================================================================\n\n")

z_crit <- qnorm(1 - (1 - IV_CONFIG$confidence_level) / 2)
credit_types <- c("Cred_Real_Total", "Cred_Real_MN", "Cred_USD")

iv_results_all <- data.frame()
iv_diagnostics_all <- data.frame()

for (cred_type in credit_types) {
  cat("Running IV-LP for", cred_type, "...\n")

  for (h in 1:IV_CONFIG$max_horizon) {
    # Prepare data
    data_h <- analysis_data %>%
      mutate(
        y_fwd = lead(!!sym(cred_type), h),
        y_lag1 = lag(!!sym(cred_type), 1),
        y_lag2 = lag(!!sym(cred_type), 2),
        fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1)
      )

    # Exclusion-restriction controls: block non-FCI channels through which DXY affects credit
    er_controls <- c("d_ToT", "US_10Y", "SP500")
    er_controls_avail <- intersect(er_controls, names(analysis_data))

    exog_vars <- c("y_lag1", "y_lag2", "fci_lag1",
                    "IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                    "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2",
                    er_controls_avail)

    # Excluded instruments: only those NOT used as controls
    excluded_instruments <- setdiff(selected_instruments, er_controls_avail)
    if (length(excluded_instruments) == 0) excluded_instruments <- selected_instruments[1]

    reg_data <- data_h %>%
      dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG,
                    all_of(exog_vars),
                    all_of(excluded_instruments)) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    exog_str <- paste(exog_vars, collapse = " + ")
    inst_str <- paste(excluded_instruments, collapse = " + ")

    # ---- 2SLS IV regression ----
    iv_formula <- as.formula(paste(
      "y_fwd ~ FCI_ENDO_exCredit_AVG +", exog_str,
      "|", inst_str, "+", exog_str
    ))

    tryCatch({
      iv_model <- ivreg::ivreg(iv_formula, data = reg_data)
      vcov_hac <- sandwich::NeweyWest(iv_model, lag = h + 1, prewhite = FALSE)
      ct_iv <- lmtest::coeftest(iv_model, vcov = vcov_hac)
      iv_summ <- summary(iv_model, vcov = vcov_hac, diagnostics = TRUE)

      idx <- which(rownames(ct_iv) == "FCI_ENDO_exCredit_AVG")

      # ---- OLS for comparison ----
      ols_formula <- as.formula(paste("y_fwd ~ FCI_ENDO_exCredit_AVG +", exog_str))
      ols_model <- lm(ols_formula, data = reg_data)
      vcov_ols <- sandwich::NeweyWest(ols_model, lag = h + 1, prewhite = FALSE)
      ct_ols <- lmtest::coeftest(ols_model, vcov = vcov_ols)
      idx_ols <- which(rownames(ct_ols) == "FCI_ENDO_exCredit_AVG")

      iv_results_all <- rbind(iv_results_all, data.frame(
        credit_type = cred_type,
        horizon = h,
        method = "2SLS_improved",
        coef = ct_iv[idx, 1],
        se = ct_iv[idx, 2],
        p_value = ct_iv[idx, 4],
        ci_lower = ct_iv[idx, 1] - z_crit * ct_iv[idx, 2],
        ci_upper = ct_iv[idx, 1] + z_crit * ct_iv[idx, 2],
        n_obs = nrow(reg_data),
        stringsAsFactors = FALSE
      ))

      iv_results_all <- rbind(iv_results_all, data.frame(
        credit_type = cred_type,
        horizon = h,
        method = "OLS",
        coef = ct_ols[idx_ols, 1],
        se = ct_ols[idx_ols, 2],
        p_value = ct_ols[idx_ols, 4],
        ci_lower = ct_ols[idx_ols, 1] - z_crit * ct_ols[idx_ols, 2],
        ci_upper = ct_ols[idx_ols, 1] + z_crit * ct_ols[idx_ols, 2],
        n_obs = nrow(reg_data),
        stringsAsFactors = FALSE
      ))

      # Diagnostics
      first_stage_F <- tryCatch(iv_summ$diagnostics["Weak instruments", "statistic"], error = function(e) NA)
      wu_hausman_p <- tryCatch(iv_summ$diagnostics["Wu-Hausman", "p-value"], error = function(e) NA)
      sargan_p <- tryCatch(iv_summ$diagnostics["Sargan", "p-value"], error = function(e) NA)

      iv_diagnostics_all <- rbind(iv_diagnostics_all, data.frame(
        credit_type = cred_type,
        horizon = h,
        first_stage_F = first_stage_F,
        wu_hausman_p = wu_hausman_p,
        sargan_p = sargan_p,
        n_instruments = length(excluded_instruments),
        instruments = paste(excluded_instruments, collapse = ";"),
        er_controls = paste(er_controls_avail, collapse = ";"),
        n_obs = nrow(reg_data),
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      if (IV_CONFIG$verbose) cat(sprintf("  h=%d: IV failed (%s)\n", h, e$message))
    })
  }
  cat("  Done.\n")
}
cat("\n")


################################################################################
# SECTION D: WEAK-IV ROBUST INFERENCE (ANDERSON-RUBIN)
################################################################################

cat("================================================================================\n")
cat("SECTION D: ANDERSON-RUBIN WEAK-IV ROBUST INFERENCE\n")
cat("================================================================================\n\n")

ar_results <- data.frame()

for (cred_type in credit_types) {
  cat("Anderson-Rubin test for", cred_type, "...\n")

  for (h in 1:IV_CONFIG$max_horizon) {
    data_h <- analysis_data %>%
      mutate(
        y_fwd = lead(!!sym(cred_type), h),
        y_lag1 = lag(!!sym(cred_type), 1),
        y_lag2 = lag(!!sym(cred_type), 2),
        fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1)
      )

    # Same exclusion-restriction controls as in Section C
    er_controls_ar <- intersect(c("d_ToT", "US_10Y", "SP500"), names(analysis_data))
    exog_vars_ar <- c("y_lag1", "y_lag2", "fci_lag1",
                    "IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                    "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2",
                    er_controls_ar)
    excluded_instruments_ar <- setdiff(selected_instruments, er_controls_ar)
    if (length(excluded_instruments_ar) == 0) excluded_instruments_ar <- selected_instruments[1]

    reg_data <- data_h %>%
      dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG,
                    all_of(exog_vars_ar),
                    all_of(excluded_instruments_ar)) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    tryCatch({
      # Anderson-Rubin test: regress y on instruments + controls, test instruments jointly
      # Under H0: delta=0, the instruments should not predict y_fwd (after controls)
      exog_str <- paste(exog_vars_ar, collapse = " + ")

      # Restricted model (controls only)
      f_restricted <- as.formula(paste("y_fwd ~", exog_str))
      m_restricted <- lm(f_restricted, data = reg_data)

      # Unrestricted model (controls + instruments)
      f_unrestricted <- as.formula(paste("y_fwd ~",
                                          paste(excluded_instruments_ar, collapse = " + "), "+",
                                          exog_str))
      m_unrestricted <- lm(f_unrestricted, data = reg_data)

      # F-test for joint significance of instruments
      ar_anova <- anova(m_restricted, m_unrestricted)
      ar_F <- ar_anova$F[2]
      ar_p <- ar_anova$`Pr(>F)`[2]

      # Anderson-Rubin confidence interval (grid search)
      # Invert the AR test by testing H0: delta = delta_0 for grid of delta_0
      delta_grid <- seq(-40, 20, by = 0.5)
      ar_reject <- sapply(delta_grid, function(d0) {
        # Under H0: delta=d0, construct y_tilde = y_fwd - d0 * FCI
        reg_data$y_tilde <- reg_data$y_fwd - d0 * reg_data$FCI_ENDO_exCredit_AVG

        f_r <- lm(as.formula(paste("y_tilde ~", exog_str)), data = reg_data)
        f_u <- lm(as.formula(paste("y_tilde ~",
                                    paste(excluded_instruments_ar, collapse = " + "), "+",
                                    exog_str)), data = reg_data)
        a <- anova(f_r, f_u)
        return(a$`Pr(>F)`[2])
      })

      # CI = values where we fail to reject
      ar_ci_lower <- min(delta_grid[ar_reject > 0.10], na.rm = TRUE)
      ar_ci_upper <- max(delta_grid[ar_reject > 0.10], na.rm = TRUE)
      if (is.infinite(ar_ci_lower)) ar_ci_lower <- NA
      if (is.infinite(ar_ci_upper)) ar_ci_upper <- NA

      ar_results <- rbind(ar_results, data.frame(
        credit_type = cred_type,
        horizon = h,
        AR_F = ar_F,
        AR_p_value = ar_p,
        AR_CI_lower = ar_ci_lower,
        AR_CI_upper = ar_ci_upper,
        n_obs = nrow(reg_data),
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      if (IV_CONFIG$verbose && h %% 6 == 0)
        cat(sprintf("  h=%d: AR test failed (%s)\n", h, e$message))
    })
  }
  cat("  Done.\n")
}

# Report key results
cat("\nAnderson-Rubin Test Results (Cred_Real_Total, key horizons):\n")
ar_total <- ar_results %>% filter(credit_type == "Cred_Real_Total")
key_h <- c(6, 12, 18, 24)
for (kh in key_h) {
  row <- ar_total %>% filter(horizon == kh)
  if (nrow(row) > 0) {
    cat(sprintf("  h=%2d: AR F=%.2f, p=%.3f, 90%% CI=[%.1f, %.1f]\n",
                kh, row$AR_F[1], row$AR_p_value[1],
                row$AR_CI_lower[1], row$AR_CI_upper[1]))
  }
}
cat("\n")


################################################################################
# SECTION E: COMPARISON WITH ORIGINAL IV-LP
################################################################################

cat("================================================================================\n")
cat("SECTION E: COMPARISON WITH ORIGINAL IV-LP\n")
cat("================================================================================\n\n")

# Load original IV results if available
orig_iv_file <- file.path(IV_CONFIG$output_dir, "LP_IV_Diagnostics.csv")
if (file.exists(orig_iv_file)) {
  orig_diag <- read.csv(orig_iv_file)
  cat("Original IV-LP diagnostics loaded.\n")

  # Compare first-stage F at key horizons
  cat("\nFirst-Stage F Comparison (Cred_Real_Total):\n")
  cat(sprintf("%-8s %15s %15s\n", "Horizon", "Original (3 IV)", "Improved"))
  cat(strrep("-", 45), "\n")
  new_diag <- iv_diagnostics_all %>% filter(credit_type == "Cred_Real_Total")
  for (kh in key_h) {
    orig_row <- orig_diag %>%
      filter(horizon == kh & grepl("Total|total", credit_type, ignore.case = TRUE))
    new_row <- new_diag %>% filter(horizon == kh)
    orig_F <- if (nrow(orig_row) > 0) sprintf("%.2f", orig_row$first_stage_F[1]) else "N/A"
    new_F <- if (nrow(new_row) > 0) sprintf("%.2f", new_row$first_stage_F[1]) else "N/A"
    cat(sprintf("h=%2d    %15s %15s\n", kh, orig_F, new_F))
  }
} else {
  cat("Original IV-LP diagnostics not found (LP_IV_Diagnostics.csv).\n")
}
cat("\n")


################################################################################
# SAVE OUTPUTS
################################################################################

cat("================================================================================\n")
cat("SAVING OUTPUTS\n")
cat("================================================================================\n\n")

# ---- CSVs ----
write.csv(instrument_diag, file.path(IV_CONFIG$output_dir, "IV_First_Stage_Battery.csv"),
          row.names = FALSE)
cat("Saved: IV_First_Stage_Battery.csv\n")

write.csv(iv_results_all, file.path(IV_CONFIG$output_dir, "IV_LP_Improved_Results.csv"),
          row.names = FALSE)
cat("Saved: IV_LP_Improved_Results.csv\n")

write.csv(iv_diagnostics_all, file.path(IV_CONFIG$output_dir, "IV_LP_Improved_Diagnostics.csv"),
          row.names = FALSE)
cat("Saved: IV_LP_Improved_Diagnostics.csv\n")

write.csv(ar_results, file.path(IV_CONFIG$output_dir, "IV_Weak_Robust_Inference.csv"),
          row.names = FALSE)
cat("Saved: IV_Weak_Robust_Inference.csv\n")


################################################################################
# VISUALIZATIONS
################################################################################

cat("\n--- Generating Visualizations ---\n\n")

# ---- 196: First-Stage Diagnostics ----
tryCatch({
  p1 <- ggplot(instrument_diag, aes(x = reorder(Instrument, abs(t_stat)),
                                     y = abs(t_stat), fill = abs(t_stat) > 2)) +
    geom_col(show.legend = FALSE) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = "First-Stage |t-statistics| by Instrument",
         subtitle = "Dashed line: |t|=2 threshold",
         x = NULL, y = "|t-statistic|") +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
    theme_minimal(base_size = 11)

  p2 <- ggplot(instrument_diag, aes(x = reorder(Instrument, Partial_R2),
                                     y = Partial_R2 * 100, fill = Partial_R2 > 0.01)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    labs(title = "Partial R² (%) by Instrument",
         subtitle = "Incremental variance explained beyond controls",
         x = NULL, y = "Partial R² (%)") +
    scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "gray60")) +
    theme_minimal(base_size = 11)

  combined <- grid.arrange(p1, p2, ncol = 2)
  ggsave(file.path(IV_CONFIG$output_dir, "196_First_Stage_Diagnostics.png"),
         combined, width = 14, height = 6, dpi = 150)
  cat("Saved: 196_First_Stage_Diagnostics.png\n")
}, error = function(e) cat("  WARNING: Plot 196 failed:", e$message, "\n"))

# ---- 197: IV-LP IRFs for real credit ----
tryCatch({
  iv_2sls <- iv_results_all %>%
    filter(method == "2SLS_improved")

  if (nrow(iv_2sls) > 0) {
    iv_2sls$credit_label <- factor(iv_2sls$credit_type,
                                    levels = c("Cred_Real_Total", "Cred_Real_MN", "Cred_USD"),
                                    labels = c("Real Total", "Real MN", "USD"))

    p <- ggplot(iv_2sls, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~credit_label, scales = "free_y") +
      labs(title = "IV-LP: FCI → Credit (Improved Instruments)",
           subtitle = paste("Instruments:", paste(selected_instruments, collapse = ", ")),
           x = "Horizon (months)", y = "pp response to 1 SD FCI tightening") +
      theme_minimal(base_size = 11)

    ggsave(file.path(IV_CONFIG$output_dir, "197_IV_LP_Improved_Credit.png"),
           p, width = 14, height = 5, dpi = 150)
    cat("Saved: 197_IV_LP_Improved_Credit.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 197 failed:", e$message, "\n"))

# ---- 198: OLS vs 2SLS Comparison ----
tryCatch({
  comp_data <- iv_results_all %>%
    filter(credit_type == "Cred_Real_Total")

  if (nrow(comp_data) > 0) {
    p <- ggplot(comp_data, aes(x = horizon, y = coef, color = method, fill = method)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "OLS vs 2SLS: FCI → Real Credit (Total)",
           x = "Horizon (months)", y = "pp response") +
      scale_color_manual(values = c("OLS" = "gray40", "2SLS_improved" = "steelblue"),
                         labels = c("OLS" = "OLS (baseline)", "2SLS_improved" = "2SLS (improved IV)")) +
      scale_fill_manual(values = c("OLS" = "gray40", "2SLS_improved" = "steelblue")) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")

    ggsave(file.path(IV_CONFIG$output_dir, "198_IV_LP_OLS_vs_2SLS_Comparison.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 198_IV_LP_OLS_vs_2SLS_Comparison.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 198 failed:", e$message, "\n"))

# ---- 199: Anderson-Rubin Test ----
tryCatch({
  ar_total <- ar_results %>% filter(credit_type == "Cred_Real_Total")
  if (nrow(ar_total) > 0) {
    p1 <- ggplot(ar_total, aes(x = horizon, y = AR_p_value)) +
      geom_line(color = "firebrick", linewidth = 0.8) +
      geom_point(color = "firebrick") +
      geom_hline(yintercept = 0.10, linetype = "dashed", color = "gray50") +
      labs(title = "Anderson-Rubin Test p-values",
           subtitle = "H0: FCI has no effect on credit (weak-IV robust)",
           x = "Horizon (months)", y = "p-value") +
      theme_minimal(base_size = 11)

    p2 <- ggplot(ar_total, aes(x = horizon)) +
      geom_ribbon(aes(ymin = AR_CI_lower, ymax = AR_CI_upper), alpha = 0.2, fill = "firebrick") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Anderson-Rubin 90% Confidence Intervals",
           subtitle = "Weak-IV robust CI for FCI effect on real credit",
           x = "Horizon (months)", y = "Effect (pp)") +
      theme_minimal(base_size = 11)

    combined <- grid.arrange(p1, p2, ncol = 2)
    ggsave(file.path(IV_CONFIG$output_dir, "199_Anderson_Rubin_Test.png"),
           combined, width = 14, height = 6, dpi = 150)
    cat("Saved: 199_Anderson_Rubin_Test.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 199 failed:", e$message, "\n"))

# ---- 200: Dashboard ----
tryCatch({
  plots <- list()

  # First-stage F comparison
  diag_total <- iv_diagnostics_all %>% filter(credit_type == "Cred_Real_Total")
  if (nrow(diag_total) > 0) {
    plots[[1]] <- ggplot(diag_total, aes(x = horizon, y = first_stage_F)) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_point(color = "steelblue", size = 1.5) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
      labs(title = "First-Stage F (Improved)", x = "Horizon", y = "F-statistic") +
      theme_minimal(base_size = 9)
  }

  # Wu-Hausman p-values
  if (nrow(diag_total) > 0 && any(!is.na(diag_total$wu_hausman_p))) {
    plots[[2]] <- ggplot(diag_total, aes(x = horizon, y = wu_hausman_p)) +
      geom_line(color = "darkgreen", linewidth = 0.8) +
      geom_point(color = "darkgreen", size = 1.5) +
      geom_hline(yintercept = 0.10, linetype = "dashed", color = "gray50") +
      labs(title = "Wu-Hausman p-values", x = "Horizon", y = "p-value") +
      theme_minimal(base_size = 9)
  }

  # 2SLS coefficients
  iv_total_2sls <- iv_results_all %>%
    filter(credit_type == "Cred_Real_Total" & method == "2SLS_improved")
  if (nrow(iv_total_2sls) > 0) {
    plots[[3]] <- ggplot(iv_total_2sls, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "2SLS Coefficient Path", x = "Horizon", y = "pp") +
      theme_minimal(base_size = 9)
  }

  # AR test
  if (nrow(ar_total) > 0) {
    plots[[4]] <- ggplot(ar_total, aes(x = horizon, y = AR_p_value)) +
      geom_line(color = "firebrick", linewidth = 0.8) +
      geom_hline(yintercept = 0.10, linetype = "dashed") +
      labs(title = "AR Test p-values", x = "Horizon", y = "p-value") +
      theme_minimal(base_size = 9)
  }

  if (length(plots) > 0) {
    combined <- do.call(grid.arrange, c(plots, ncol = 2))
    ggsave(file.path(IV_CONFIG$output_dir, "200_IV_Improvement_Dashboard.png"),
           combined, width = 12, height = 8, dpi = 150)
    cat("Saved: 200_IV_Improvement_Dashboard.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 200 failed:", e$message, "\n"))


################################################################################
# SECTION F: DXY MAGNITUDE & exTCN ROBUSTNESS
################################################################################

cat("================================================================================\n")
cat("SECTION F: DXY MAGNITUDE COMPUTATION & exTCN ROBUSTNESS\n")
cat("================================================================================\n\n")

# ---- Compute SD(DXY) for correct DXY shock magnitude ----
if ("DXY" %in% names(analysis_data)) {
  dxy_sd <- sd(analysis_data$DXY, na.rm = TRUE)
  dxy_mean <- mean(analysis_data$DXY, na.rm = TRUE)
  cat(sprintf("DXY in analysis_data: Mean = %.2f, SD = %.2f\n", dxy_mean, dxy_sd))

  # First-stage coefficient from Table 3: 0.042 per DXY index point
  first_stage_beta <- instrument_diag$Coefficient[instrument_diag$Instrument == "DXY"]
  if (length(first_stage_beta) > 0) {
    fci_per_sd_dxy <- first_stage_beta * dxy_sd
    cat(sprintf("First-stage beta (DXY): %.4f per index point\n", first_stage_beta))
    cat(sprintf("1 SD DXY (%.1f points) -> %.3f SD FCI tightening\n", dxy_sd, fci_per_sd_dxy))

    # OLS credit effect at h=12
    ols_h12 <- iv_results_all %>%
      filter(credit_type == "Cred_Real_Total", method == "OLS", horizon == 12)
    if (nrow(ols_h12) > 0) {
      credit_via_dxy <- fci_per_sd_dxy * ols_h12$coef[1]
      cat(sprintf("1 SD DXY -> credit effect: %.3f * %.2f = %.1f pp at h=12\n",
                  fci_per_sd_dxy, ols_h12$coef[1], credit_via_dxy))
    }
  }
}

# ---- exTCN first-stage robustness ----
cat("\n--- FCI_ENDO_exCredit_exTCN first-stage check ---\n")
if ("FCI_ENDO_exCredit_exTCN_AVG" %in% names(analysis_data)) {
  exTCN_fs <- analysis_data %>%
    dplyr::select(FCI_ENDO_exCredit_exTCN_AVG, DXY, all_of(control_vars)) %>%
    na.omit()

  if (nrow(exTCN_fs) > 30) {
    f_base_tcn <- lm(FCI_ENDO_exCredit_exTCN_AVG ~ .,
                      data = exTCN_fs %>% dplyr::select(-DXY))
    f_full_tcn <- lm(FCI_ENDO_exCredit_exTCN_AVG ~ .,
                      data = exTCN_fs)
    partial_r2_tcn <- summary(f_full_tcn)$r.squared - summary(f_base_tcn)$r.squared
    f_stat_tcn <- anova(f_base_tcn, f_full_tcn)$F[2]
    cat(sprintf("  DXY -> FCI_ENDO_exCredit_exTCN: partial R² = %.4f, F = %.2f\n",
                partial_r2_tcn, f_stat_tcn))
    cat(sprintf("  Compare with DXY -> FCI_ENDO_exCredit: partial R² = %.4f, F = %.2f\n",
                instrument_diag$Partial_R2[instrument_diag$Instrument == "DXY"],
                instrument_diag$F_stat[instrument_diag$Instrument == "DXY"]))
    if (f_stat_tcn < instrument_diag$F_stat[instrument_diag$Instrument == "DXY"] * 0.5) {
      cat("  >>> TCN is a major channel for DXY -> FCI (F dropped by >50%)\n")
    } else {
      cat("  >>> Non-TCN domestic variables also respond to DXY (F preserved)\n")
    }
  }
} else {
  cat("  FCI_ENDO_exCredit_exTCN_AVG not available in analysis_data.\n")
  cat("  Ensure script 01 was run with the exTCN variant.\n")
}

cat("\n################################################################################\n")
cat("SCRIPT 18 COMPLETE\n")
cat("################################################################################\n\n")
cat("Outputs: 5 PNGs (196-200) + 4 CSVs + DXY magnitude diagnostics\n")
cat("Key finding: First-stage F with enriched controls (see diagnostics)\n\n")
