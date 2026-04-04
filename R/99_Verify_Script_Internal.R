################################################################################
# VERIFICATION SCRIPT: Generate CSV Evidence for Script-Internal Computations
################################################################################
#
# Project:      Financial Conditions Index for Paraguay
# Purpose:      Generate verifiable CSV output for numbers that previously
#               existed only as script-internal computations. This ensures
#               every number in the paper and appendix has a traceable CSV.
#
# Items verified:
#   1. Post-IT unconditional first-stage F (paper line ~1056, App. K.7)
#   2. Post-IT conditional first-stage F (paper line ~1059, App. K.8)
#   3. Chow structural break test (paper lines 1081-1083)
#   4. Sample end-date sensitivity — June 2025 (Appendix N.5)
#
# NOTE: This script reads pre-computed FCI series from CSV files.
#       It does NOT re-estimate any FCI or re-run the full pipeline.
#       Existing CSV files are NOT modified.
#
# Output: output/csv/Verification_Script_Internal.csv
#
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(sandwich)
  library(lmtest)
  library(ivreg)
})

set.seed(20260328)

cat("\n################################################################################\n")
cat("VERIFICATION SCRIPT: Script-Internal Computations\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# DATA LOADING — from pre-computed pipeline outputs + raw data
# ==============================================================================

csv_dir <- "../output/csv"

# New_External_Variables.csv already contains FCI series, DXY, ToT,
# macro controls, and credit data — use it directly as the analysis dataset
ext_data <- read.csv(file.path(csv_dir, "New_External_Variables.csv"))
ext_data$fecha <- as.Date(ext_data$fecha)
cat("Loaded New_External_Variables.csv:", nrow(ext_data), "obs\n")

analysis_data <- ext_data %>% arrange(fecha)

# Create macro control lags
analysis_data <- analysis_data %>%
  mutate(
    IMAEP_yoy_L1 = lag(IMAEP_yoy, 1),
    IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
    IPC_yoy_L1 = lag(IPC_yoy, 1),
    IPC_yoy_L2 = lag(IPC_yoy, 2),
    regime_IT = as.integer(fecha >= as.Date("2011-05-01"))
  )

cat("Merged analysis data:", nrow(analysis_data), "obs\n\n")

# Control variables used throughout
control_vars <- c("IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                   "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2")
control_str <- paste(control_vars, collapse = " + ")

# Container for all results
all_results <- data.frame()

# ==============================================================================
# ITEM 1: POST-IT UNCONDITIONAL FIRST-STAGE F
# Paper: F = 0.14, coef = -0.003, t = -0.38, p = 0.707, R² = 0.19%
# Appendix K.7: same values
# ==============================================================================

cat("=" |> rep(80) |> paste(collapse=""), "\n")
cat("ITEM 1: Post-IT Unconditional First-Stage\n")
cat("=" |> rep(80) |> paste(collapse=""), "\n\n")

post_it <- analysis_data %>% filter(fecha >= as.Date("2011-05-01"))

fs_data <- post_it %>%
  dplyr::select(FCI_ENDO_exCredit_AVG, DXY, all_of(control_vars)) %>%
  na.omit()

cat("Post-IT first-stage sample:", nrow(fs_data), "obs\n")

# Full model with DXY
f_full <- as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY +", control_str))
m_full <- lm(f_full, data = fs_data)
vcov_hac <- sandwich::NeweyWest(m_full, lag = 3, prewhite = FALSE)
ct <- lmtest::coeftest(m_full, vcov = vcov_hac)
idx_dxy <- which(rownames(ct) == "DXY")

# Controls-only model for partial R²
f_base <- as.formula(paste("FCI_ENDO_exCredit_AVG ~", control_str))
m_base <- lm(f_base, data = fs_data)
partial_r2 <- summary(m_full)$r.squared - summary(m_base)$r.squared

# F-stat = t²
t_stat <- ct[idx_dxy, 3]
f_stat <- t_stat^2

cat(sprintf("  DXY coefficient:  %+.6f\n", ct[idx_dxy, 1]))
cat(sprintf("  t-statistic:      %+.4f\n", t_stat))
cat(sprintf("  p-value:          %.6f\n", ct[idx_dxy, 4]))
cat(sprintf("  Partial R²:       %.6f (%.2f%%)\n", partial_r2, partial_r2 * 100))
cat(sprintf("  F-stat (t²):      %.4f\n", f_stat))

all_results <- rbind(all_results, data.frame(
  item = "PostIT_unconditional_FS",
  horizon = NA_integer_,
  statistic = c("DXY_coef", "t_stat", "p_value", "partial_R2_pct", "F_stat"),
  value = c(ct[idx_dxy, 1], t_stat, ct[idx_dxy, 4], partial_r2 * 100, f_stat),
  n_obs = nrow(fs_data),
  paper_value = c("-0.003", "-0.38", "0.707", "0.19", "0.14"),
  stringsAsFactors = FALSE
))

cat("\n")

# ==============================================================================
# ITEM 2: POST-IT CONDITIONAL FIRST-STAGE F
# Paper: F = 6.5-6.7, Appendix K.8: h=6 F=6.47, h=12 F=6.74, h=18 F=6.58
# ==============================================================================

cat("=" |> rep(80) |> paste(collapse=""), "\n")
cat("ITEM 2: Post-IT Conditional First-Stage\n")
cat("=" |> rep(80) |> paste(collapse=""), "\n\n")

er_controls <- intersect(c("d_ToT", "US_10Y", "SP500"), names(analysis_data))

for (h in c(6, 12, 18)) {
  data_h <- post_it %>%
    mutate(
      y_fwd = lead(Cred_Real_Total, h),
      y_lag1 = lag(Cred_Real_Total, 1),
      y_lag2 = lag(Cred_Real_Total, 2),
      fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1)
    )

  exog_vars <- c("y_lag1", "y_lag2", "fci_lag1", control_vars, er_controls)
  exog_str <- paste(exog_vars, collapse = " + ")

  reg_data <- data_h %>%
    dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG, all_of(exog_vars), DXY) %>%
    na.omit()

  # IV regression for conditional first-stage F
  iv_formula <- as.formula(paste(
    "y_fwd ~ FCI_ENDO_exCredit_AVG +", exog_str, "| DXY +", exog_str
  ))
  iv_model <- ivreg::ivreg(iv_formula, data = reg_data)
  vcov_hac_iv <- sandwich::NeweyWest(iv_model, lag = h + 1, prewhite = FALSE)
  iv_summ <- summary(iv_model, vcov = vcov_hac_iv, diagnostics = TRUE)
  cond_F <- iv_summ$diagnostics["Weak instruments", "statistic"]

  # Manual first-stage for t-stat and partial R²
  fs_full <- lm(as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY +", exog_str)), data = reg_data)
  fs_base <- lm(as.formula(paste("FCI_ENDO_exCredit_AVG ~", exog_str)), data = reg_data)
  vcov_fs <- sandwich::NeweyWest(fs_full, lag = h + 1, prewhite = FALSE)
  ct_fs <- lmtest::coeftest(fs_full, vcov = vcov_fs)
  idx_d <- which(rownames(ct_fs) == "DXY")
  pr2 <- summary(fs_full)$r.squared - summary(fs_base)$r.squared

  cat(sprintf("  h=%d: N=%d, Cond F=%.4f, DXY t=%+.4f, p=%.6f, partial R²=%.4f%%\n",
              h, nrow(reg_data), cond_F, ct_fs[idx_d, 3], ct_fs[idx_d, 4], pr2 * 100))

  all_results <- rbind(all_results, data.frame(
    item = "PostIT_conditional_FS",
    horizon = h,
    statistic = c("cond_F", "DXY_t", "DXY_p", "partial_R2_pct"),
    value = c(cond_F, ct_fs[idx_d, 3], ct_fs[idx_d, 4], pr2 * 100),
    n_obs = nrow(reg_data),
    paper_value = c(
      ifelse(h == 6, "6.47", ifelse(h == 12, "6.74", "6.58")),
      ifelse(h == 6, "2.64", ifelse(h == 12, "2.78", "2.86")),
      ifelse(h == 6, "0.009", ifelse(h == 12, "0.006", "0.005")),
      ifelse(h == 6, "0.4", ifelse(h == 12, "0.4", "0.5"))
    ),
    stringsAsFactors = FALSE
  ))
}

cat("\n")

# ==============================================================================
# ITEM 3: CHOW TEST FOR DXY→FCI STRUCTURAL BREAK AT MAY 2011
# Paper: Chow F = 82.19, HAC Wald F = 29.38, interaction = -0.064, t = -6.49
# ==============================================================================

cat("=" |> rep(80) |> paste(collapse=""), "\n")
cat("ITEM 3: Chow Test for Structural Break\n")
cat("=" |> rep(80) |> paste(collapse=""), "\n\n")

chow_data <- analysis_data %>%
  dplyr::select(fecha, FCI_ENDO_exCredit_AVG, DXY, regime_IT,
                all_of(control_vars)) %>%
  na.omit()

chow_data$DXY_x_IT <- chow_data$DXY * chow_data$regime_IT

cat("Chow test sample:", nrow(chow_data), "obs\n")
cat("  Pre-IT:", sum(chow_data$regime_IT == 0), "obs\n")
cat("  Post-IT:", sum(chow_data$regime_IT == 1), "obs\n\n")

# Restricted model (no break)
m_restricted <- lm(as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY +", control_str)),
                    data = chow_data)

# Unrestricted model (with break)
m_unrestricted <- lm(
  as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY + regime_IT + DXY_x_IT +", control_str)),
  data = chow_data
)

# Standard Chow F-test
anova_result <- anova(m_restricted, m_unrestricted)
chow_F <- anova_result$F[2]
chow_p <- anova_result$`Pr(>F)`[2]

# HAC-robust coefficients
vcov_unr <- sandwich::NeweyWest(m_unrestricted, lag = 3, prewhite = FALSE)
ct_unr <- lmtest::coeftest(m_unrestricted, vcov = vcov_unr)

# HAC-robust Wald test
suppressPackageStartupMessages(library(car))
wald_result <- linearHypothesis(m_unrestricted, c("regime_IT = 0", "DXY_x_IT = 0"),
                                 vcov. = vcov_unr)
wald_F <- wald_result$F[2]
wald_p <- wald_result$`Pr(>F)`[2]

# Extract interaction coefficient
idx_interact <- which(rownames(ct_unr) == "DXY_x_IT")
interact_coef <- ct_unr[idx_interact, 1]
interact_t <- ct_unr[idx_interact, 3]
interact_p <- ct_unr[idx_interact, 4]

cat(sprintf("  Chow F-stat:      %.4f (p = %.2e)\n", chow_F, chow_p))
cat(sprintf("  HAC Wald F-stat:  %.4f (p = %.2e)\n", wald_F, wald_p))
cat(sprintf("  DXY×IT coef:      %+.6f\n", interact_coef))
cat(sprintf("  DXY×IT t-stat:    %+.4f\n", interact_t))
cat(sprintf("  DXY×IT p-value:   %.2e\n", interact_p))

all_results <- rbind(all_results, data.frame(
  item = "Chow_test",
  horizon = NA_integer_,
  statistic = c("Chow_F", "Chow_p", "HAC_Wald_F", "HAC_Wald_p",
                "DXY_x_IT_coef", "DXY_x_IT_t", "DXY_x_IT_p"),
  value = c(chow_F, chow_p, wald_F, wald_p,
            interact_coef, interact_t, interact_p),
  n_obs = nrow(chow_data),
  paper_value = c("82.19", "<0.001", "29.38", "<0.001",
                  "-0.064", "-6.49", "<0.001"),
  stringsAsFactors = FALSE
))

cat("\n")

# ==============================================================================
# ITEM 4: SAMPLE END-DATE SENSITIVITY (JUNE 2025 vs DECEMBER 2025)
# Appendix Table N.5
# ==============================================================================

cat("=" |> rep(80) |> paste(collapse=""), "\n")
cat("ITEM 4: Sample End-Date Sensitivity (June 2025 Truncation)\n")
cat("=" |> rep(80) |> paste(collapse=""), "\n\n")

# Truncate to June 2025
trunc_data <- analysis_data %>% filter(fecha <= as.Date("2025-06-01"))
cat("Truncated sample:", nrow(trunc_data), "obs (to", as.character(max(trunc_data$fecha)), ")\n\n")

# ---- 4a: Full-sample LP — FCI_exCredit → Real Total Credit, h=12 ----
cat("--- 4a: Full-sample LP (FCI_exCredit → Credit, h=12) ---\n")

h <- 12
data_h <- trunc_data %>%
  mutate(
    y_fwd = lead(Cred_Real_Total, h),
    y_lag1 = lag(Cred_Real_Total, 1),
    y_lag2 = lag(Cred_Real_Total, 2),
    fci_lag1 = lag(FCI_exCredit_AVG, 1)
  )

lp_vars <- c("y_fwd", "FCI_exCredit_AVG", "y_lag1", "y_lag2", "fci_lag1", control_vars)
reg_data_4a <- data_h %>% dplyr::select(all_of(lp_vars)) %>% na.omit()

lp_formula <- y_fwd ~ FCI_exCredit_AVG + y_lag1 + y_lag2 + fci_lag1 +
  IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +
  IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2

m_4a <- lm(lp_formula, data = reg_data_4a)
vcov_4a <- sandwich::NeweyWest(m_4a, lag = h + 1, prewhite = FALSE)
ct_4a <- lmtest::coeftest(m_4a, vcov = vcov_4a)
idx_fci <- which(rownames(ct_4a) == "FCI_exCredit_AVG")

cat(sprintf("  N=%d, coef=%.4f, p=%.6f\n", nrow(reg_data_4a), ct_4a[idx_fci, 1], ct_4a[idx_fci, 4]))

all_results <- rbind(all_results, data.frame(
  item = "SampleSensitivity_Jun2025",
  horizon = 12L,
  statistic = c("credit_LP_coef", "credit_LP_p"),
  value = c(ct_4a[idx_fci, 1], ct_4a[idx_fci, 4]),
  n_obs = nrow(reg_data_4a),
  paper_value = c("-8.05", "0.005"),
  stringsAsFactors = FALSE
))

# ---- 4b: Post-IT LP — FCI_ENDO_exCredit → Real Total Credit, h=12 ----
cat("--- 4b: Post-IT LP (FCI_ENDO_exCredit → Credit, h=12) ---\n")

trunc_postit <- trunc_data %>% filter(fecha >= as.Date("2011-05-01"))

data_h_b <- trunc_postit %>%
  mutate(
    y_fwd = lead(Cred_Real_Total, h),
    y_lag1 = lag(Cred_Real_Total, 1),
    y_lag2 = lag(Cred_Real_Total, 2),
    fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1)
  )

lp_vars_b <- c("y_fwd", "FCI_ENDO_exCredit_AVG", "y_lag1", "y_lag2", "fci_lag1", control_vars)
reg_data_4b <- data_h_b %>% dplyr::select(all_of(lp_vars_b)) %>% na.omit()

lp_formula_b <- y_fwd ~ FCI_ENDO_exCredit_AVG + y_lag1 + y_lag2 + fci_lag1 +
  IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +
  IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2

m_4b <- lm(lp_formula_b, data = reg_data_4b)
vcov_4b <- sandwich::NeweyWest(m_4b, lag = h + 1, prewhite = FALSE)
ct_4b <- lmtest::coeftest(m_4b, vcov = vcov_4b)
idx_fci_b <- which(rownames(ct_4b) == "FCI_ENDO_exCredit_AVG")

cat(sprintf("  N=%d, coef=%.4f, p=%.6f\n", nrow(reg_data_4b), ct_4b[idx_fci_b, 1], ct_4b[idx_fci_b, 4]))

all_results <- rbind(all_results, data.frame(
  item = "SampleSensitivity_Jun2025",
  horizon = 12L,
  statistic = c("postIT_LP_coef", "postIT_LP_p"),
  value = c(ct_4b[idx_fci_b, 1], ct_4b[idx_fci_b, 4]),
  n_obs = nrow(reg_data_4b),
  paper_value = c("-4.78", "0.008"),
  stringsAsFactors = FALSE
))

# ---- 4c: DXY First-Stage F (truncated full sample) ----
cat("--- 4c: DXY First-Stage F (truncated full sample) ---\n")

fs_data_trunc <- trunc_data %>%
  dplyr::select(FCI_ENDO_exCredit_AVG, DXY, all_of(control_vars)) %>%
  na.omit()

m_full_trunc <- lm(as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY +", control_str)),
                    data = fs_data_trunc)
vcov_fs_trunc <- sandwich::NeweyWest(m_full_trunc, lag = 3, prewhite = FALSE)
ct_fs_trunc <- lmtest::coeftest(m_full_trunc, vcov = vcov_fs_trunc)
idx_d_trunc <- which(rownames(ct_fs_trunc) == "DXY")
t_trunc <- ct_fs_trunc[idx_d_trunc, 3]
f_trunc <- t_trunc^2

cat(sprintf("  N=%d, F=%.4f, p=%.2e\n", nrow(fs_data_trunc), f_trunc, ct_fs_trunc[idx_d_trunc, 4]))

all_results <- rbind(all_results, data.frame(
  item = "SampleSensitivity_Jun2025",
  horizon = NA_integer_,
  statistic = c("DXY_FS_F", "DXY_FS_p"),
  value = c(f_trunc, ct_fs_trunc[idx_d_trunc, 4]),
  n_obs = nrow(fs_data_trunc),
  paper_value = c("35.52", "<0.001"),
  stringsAsFactors = FALSE
))

# ---- 4d: Anderson-Rubin F at h=12 (truncated full sample) ----
cat("--- 4d: Anderson-Rubin F at h=12 (truncated full sample) ---\n")

er_controls <- intersect(c("d_ToT", "US_10Y", "SP500"), names(trunc_data))
exog_vars_ar <- c("y_lag1", "y_lag2", "fci_lag1", control_vars, er_controls)
exog_str_ar <- paste(exog_vars_ar, collapse = " + ")

data_h_d <- trunc_data %>%
  mutate(
    y_fwd = lead(Cred_Real_Total, h),
    y_lag1 = lag(Cred_Real_Total, 1),
    y_lag2 = lag(Cred_Real_Total, 2),
    fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1)
  )

reg_data_4d <- data_h_d %>%
  dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG, all_of(exog_vars_ar), DXY) %>%
  na.omit()

# AR test at delta=0
f_r0 <- lm(as.formula(paste("y_fwd ~", exog_str_ar)), data = reg_data_4d)
f_u0 <- lm(as.formula(paste("y_fwd ~ DXY +", exog_str_ar)), data = reg_data_4d)
ar_anova <- anova(f_r0, f_u0)
ar_F <- ar_anova$F[2]
ar_p <- ar_anova$`Pr(>F)`[2]

cat(sprintf("  N=%d, AR F=%.4f, p=%.2e\n", nrow(reg_data_4d), ar_F, ar_p))

all_results <- rbind(all_results, data.frame(
  item = "SampleSensitivity_Jun2025",
  horizon = 12L,
  statistic = c("AR_F", "AR_p"),
  value = c(ar_F, ar_p),
  n_obs = nrow(reg_data_4d),
  paper_value = c("141.03", "<0.001"),
  stringsAsFactors = FALSE
))

# ---- 4e: FCI×ToT interaction at h=12 (truncated full sample) ----
cat("--- 4e: FCI×ToT interaction at h=12 (truncated full sample) ---\n")

data_h_e <- trunc_data %>%
  mutate(
    y_fwd = lead(Cred_Real_Total, h),
    y_lag1 = lag(Cred_Real_Total, 1),
    y_lag2 = lag(Cred_Real_Total, 2),
    fci_lag1 = lag(FCI_exCredit_AVG, 1),
    FCI_x_ToT = FCI_exCredit_AVG * d_ToT
  )

reg_data_4e <- data_h_e %>%
  dplyr::select(y_fwd, FCI_exCredit_AVG, d_ToT, FCI_x_ToT,
                y_lag1, y_lag2, fci_lag1, all_of(control_vars)) %>%
  na.omit()

f_interact <- y_fwd ~ FCI_exCredit_AVG + d_ToT + FCI_x_ToT +
  y_lag1 + y_lag2 + fci_lag1 +
  IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +
  IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2

m_4e <- lm(f_interact, data = reg_data_4e)
vcov_4e <- sandwich::NeweyWest(m_4e, lag = h + 1, prewhite = FALSE)
ct_4e <- lmtest::coeftest(m_4e, vcov = vcov_4e)
idx_inter <- which(rownames(ct_4e) == "FCI_x_ToT")

cat(sprintf("  N=%d, FCI×ToT coef=%.6f, p=%.6f\n",
            nrow(reg_data_4e), ct_4e[idx_inter, 1], ct_4e[idx_inter, 4]))

all_results <- rbind(all_results, data.frame(
  item = "SampleSensitivity_Jun2025",
  horizon = 12L,
  statistic = c("FCI_x_ToT_coef", "FCI_x_ToT_p"),
  value = c(ct_4e[idx_inter, 1], ct_4e[idx_inter, 4]),
  n_obs = nrow(reg_data_4e),
  paper_value = c("-0.229", "0.003"),
  stringsAsFactors = FALSE
))

cat("\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

write.csv(all_results, file.path(csv_dir, "Verification_Script_Internal.csv"),
          row.names = FALSE)

cat("=" |> rep(80) |> paste(collapse=""), "\n")
cat("SAVED: Verification_Script_Internal.csv\n")
cat("=" |> rep(80) |> paste(collapse=""), "\n\n")

# ==============================================================================
# VERIFICATION SUMMARY
# ==============================================================================

cat("VERIFICATION SUMMARY:\n\n")

for (itm in unique(all_results$item)) {
  cat(sprintf("  %s:\n", itm))
  subset <- all_results[all_results$item == itm, ]
  for (i in seq_len(nrow(subset))) {
    r <- subset[i, ]
    h_str <- if (is.na(r$horizon)) "" else sprintf(" (h=%d)", r$horizon)
    match_status <- ""
    pv <- as.character(r$paper_value)
    if (grepl("<", pv)) {
      match_status <- if (r$value < 0.001) "MATCH" else "MISMATCH"
    } else {
      pv_num <- as.numeric(pv)
      # Check rounding match
      n_dec <- nchar(sub(".*\\.", "", pv))
      rounded <- round(r$value, n_dec)
      match_status <- if (abs(rounded - pv_num) < 10^(-n_dec)) "MATCH" else "MISMATCH"
    }
    cat(sprintf("    %-20s%s: computed=%.6f, paper=%s → %s\n",
                r$statistic, h_str, r$value, pv, match_status))
  }
  cat("\n")
}

cat("Done.\n")
