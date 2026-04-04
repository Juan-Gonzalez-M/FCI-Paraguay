################################################################################
# EXPERT FEEDBACK COMPUTATIONS
################################################################################
#
# Project:      Financial Conditions Index - Expert Reviewer Response
# Institution:  Banco Central del Paraguay
#
# Description:  Three targeted computations for expert feedback response:
#   PART A: IV-LP with FCI_ENDO_exCredit_exTCN + DXY (TCN conduit test)
#   PART B: Chow test for DXY→FCI structural break at May 2011
#   PART C: ToT interaction with FCI_ENDO_exCredit_exTCN (domestic-only)
#
# Output:
#   - output/csv/Expert_IV_LP_exTCN.csv
#   - output/csv/Expert_Chow_Test.csv
#   - output/csv/Expert_ToT_Interaction_Domestic.csv
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
  library(sandwich)
  library(lmtest)
  library(ivreg)
})

EXP_CONFIG <- list(
  max_horizon = 24,
  confidence_level = 0.90,
  csv_dir    = "../output/csv",
  verbose    = TRUE
)

find_csv <- function(filename) {
  p1 <- file.path(EXP_CONFIG$csv_dir, filename)
  if (file.exists(p1)) return(p1)
  p2 <- file.path("../output", filename)
  if (file.exists(p2)) return(p2)
  return(p2)
}

set.seed(20260317)
z_crit <- qnorm(1 - (1 - EXP_CONFIG$confidence_level) / 2)

cat("\n################################################################################\n")
cat("EXPERT FEEDBACK COMPUTATIONS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---- Load data (same pattern as scripts 23/25) ----
ext_file <- find_csv("New_External_Variables.csv")
if (!file.exists(ext_file)) {
  stop("New_External_Variables.csv not found. Run script 17 first.")
}
ext_data <- read.csv(ext_file)
ext_data$fecha <- as.Date(ext_data$fecha)
cat("Loaded external data:", nrow(ext_data), "obs\n")

# Load FCI_ENDO_exCredit_exTCN_AVG from FCI_Complete_Results.csv
fci_file <- find_csv("FCI_Complete_Results.csv")
fci_data <- read.csv(fci_file)
fci_data$fecha <- as.Date(fci_data$fecha)
cat("Loaded FCI data:", nrow(fci_data), "obs,", ncol(fci_data), "columns\n")

# Load credit data (MN and USD)
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

# Merge all data
analysis_data <- ext_data %>%
  left_join(fci_data %>% dplyr::select(fecha, FCI_ENDO_exCredit_exTCN_AVG),
            by = "fecha") %>%
  left_join(credit_data, by = "fecha") %>%
  arrange(fecha)

# Create macro lags
analysis_data <- analysis_data %>%
  mutate(
    IMAEP_yoy_L1 = lag(IMAEP_yoy, 1),
    IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
    IPC_yoy_L1 = lag(IPC_yoy, 1),
    IPC_yoy_L2 = lag(IPC_yoy, 2)
  )

# Load regime indicators
regime_file <- find_csv("FCI_Regime_Indicators.csv")
regime_data <- read.csv(regime_file)
regime_data$fecha <- as.Date(regime_data$fecha)
analysis_data <- analysis_data %>%
  left_join(regime_data %>% dplyr::select(fecha, regime_IT), by = "fecha")

cat("Full analysis data:", nrow(analysis_data), "obs\n")
cat("FCI_ENDO_exCredit_exTCN_AVG available:", sum(!is.na(analysis_data$FCI_ENDO_exCredit_exTCN_AVG)), "obs\n\n")

# Exclusion-restriction controls (same as script 23)
er_controls <- intersect(c("d_ToT", "US_10Y", "SP500"), names(analysis_data))
exog_vars <- c("y_lag1", "y_lag2", "fci_lag1",
                "IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2",
                er_controls)


################################################################################
# PART A: IV-LP WITH FCI_ENDO_exCredit_exTCN + DXY
################################################################################

cat("================================================================================\n")
cat("PART A: IV-LP WITH FCI_ENDO_exCredit_exTCN (7-var domestic, excl TCN) + DXY\n")
cat("================================================================================\n\n")

fci_var_a <- "FCI_ENDO_exCredit_exTCN_AVG"
iv_results_a <- data.frame()

for (h in 1:EXP_CONFIG$max_horizon) {
  data_h <- analysis_data %>%
    mutate(
      y_fwd = lead(Cred_Real_Total, h),
      y_lag1 = lag(Cred_Real_Total, 1),
      y_lag2 = lag(Cred_Real_Total, 2),
      fci_lag1 = lag(!!sym(fci_var_a), 1)
    )

  exog_str <- paste(exog_vars, collapse = " + ")

  reg_data <- data_h %>%
    dplyr::select(y_fwd, !!sym(fci_var_a),
                  all_of(exog_vars), DXY) %>%
    na.omit()

  if (nrow(reg_data) < 50) next

  tryCatch({
    # ---- 2SLS with DXY ----
    iv_formula <- as.formula(paste(
      "y_fwd ~", fci_var_a, "+", exog_str,
      "| DXY +", exog_str
    ))
    iv_model <- ivreg::ivreg(iv_formula, data = reg_data)
    vcov_hac <- sandwich::NeweyWest(iv_model, lag = h + 1, prewhite = FALSE)
    ct_iv <- lmtest::coeftest(iv_model, vcov = vcov_hac)
    iv_summ <- summary(iv_model, vcov = vcov_hac, diagnostics = TRUE)

    idx_iv <- which(rownames(ct_iv) == fci_var_a)

    first_stage_F <- tryCatch(
      iv_summ$diagnostics["Weak instruments", "statistic"], error = function(e) NA)
    wu_hausman_p <- tryCatch(
      iv_summ$diagnostics["Wu-Hausman", "p-value"], error = function(e) NA)

    # ---- OLS comparison ----
    ols_formula <- as.formula(paste("y_fwd ~", fci_var_a, "+", exog_str))
    ols_model <- lm(ols_formula, data = reg_data)
    vcov_ols <- sandwich::NeweyWest(ols_model, lag = h + 1, prewhite = FALSE)
    ct_ols <- lmtest::coeftest(ols_model, vcov = vcov_ols)
    idx_ols <- which(rownames(ct_ols) == fci_var_a)

    # ---- Anderson-Rubin test ----
    delta_grid <- seq(-300, 50, by = 1)
    ar_pvals <- sapply(delta_grid, function(d0) {
      reg_data$y_tilde <- reg_data$y_fwd - d0 * reg_data[[fci_var_a]]
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
    f_r0 <- lm(as.formula(paste("y_fwd ~", exog_str)), data = reg_data)
    f_u0 <- lm(as.formula(paste("y_fwd ~ DXY +", exog_str)), data = reg_data)
    ar_anova <- anova(f_r0, f_u0)
    ar_F <- ar_anova$F[2]
    ar_p <- ar_anova$`Pr(>F)`[2]

    # Store 2SLS
    iv_results_a <- rbind(iv_results_a, data.frame(
      horizon = h, method = "2SLS_DXY",
      coef = ct_iv[idx_iv, 1], se = ct_iv[idx_iv, 2], p_value = ct_iv[idx_iv, 4],
      ci_lower = ct_iv[idx_iv, 1] - z_crit * ct_iv[idx_iv, 2],
      ci_upper = ct_iv[idx_iv, 1] + z_crit * ct_iv[idx_iv, 2],
      first_stage_F = first_stage_F, wu_hausman_p = wu_hausman_p,
      AR_F = ar_F, AR_p = ar_p,
      AR_CI_lower = ar_ci_lower, AR_CI_upper = ar_ci_upper,
      n_obs = nrow(reg_data), stringsAsFactors = FALSE
    ))

    # Store OLS
    iv_results_a <- rbind(iv_results_a, data.frame(
      horizon = h, method = "OLS",
      coef = ct_ols[idx_ols, 1], se = ct_ols[idx_ols, 2], p_value = ct_ols[idx_ols, 4],
      ci_lower = ct_ols[idx_ols, 1] - z_crit * ct_ols[idx_ols, 2],
      ci_upper = ct_ols[idx_ols, 1] + z_crit * ct_ols[idx_ols, 2],
      first_stage_F = NA, wu_hausman_p = NA,
      AR_F = NA, AR_p = NA,
      AR_CI_lower = NA, AR_CI_upper = NA,
      n_obs = nrow(reg_data), stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    if (h %% 6 == 0) cat(sprintf("  h=%d: failed (%s)\n", h, e$message))
  })
}

# Report key results
cat("\nIV-LP Results: FCI_ENDO_exCredit_exTCN (7-var domestic, no TCN) + DXY\n")
cat(sprintf("%-6s %10s %10s %12s %10s %22s\n",
            "h", "OLS", "2SLS_DXY", "First-stg F", "Wu-Haus p", "AR 90% CI"))
cat(strrep("-", 75), "\n")
for (kh in c(6, 12, 18, 24)) {
  ols_row <- iv_results_a %>% filter(method == "OLS", horizon == kh)
  iv_row <- iv_results_a %>% filter(method == "2SLS_DXY", horizon == kh)
  if (nrow(ols_row) > 0 && nrow(iv_row) > 0) {
    cat(sprintf("h=%2d  %+10.2f %+10.2f %12.1f %10.3f    [%+.1f, %+.1f]\n",
                kh, ols_row$coef[1], iv_row$coef[1], iv_row$first_stage_F[1],
                iv_row$wu_hausman_p[1], iv_row$AR_CI_lower[1], iv_row$AR_CI_upper[1]))
  }
}

write.csv(iv_results_a, file.path(EXP_CONFIG$csv_dir, "Expert_IV_LP_exTCN.csv"),
          row.names = FALSE)
cat("\nSaved: Expert_IV_LP_exTCN.csv\n\n")


################################################################################
# PART B: CHOW TEST FOR DXY→FCI BREAK AT MAY 2011
################################################################################

cat("================================================================================\n")
cat("PART B: CHOW TEST FOR DXY→FCI STRUCTURAL BREAK AT MAY 2011\n")
cat("================================================================================\n\n")

# Use FCI_ENDO_exCredit_AVG as dependent variable (matches paper's first-stage)
chow_data <- analysis_data %>%
  dplyr::select(fecha, FCI_ENDO_exCredit_AVG, DXY, regime_IT,
                IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
  na.omit()

cat("Chow test sample:", nrow(chow_data), "obs\n")
cat("  Pre-IT:", sum(chow_data$regime_IT == 0), "obs\n")
cat("  Post-IT:", sum(chow_data$regime_IT == 1), "obs\n\n")

control_str <- "IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 + IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2"

# Restricted model: no structural break
f_restricted <- as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY +", control_str))
m_restricted <- lm(f_restricted, data = chow_data)

# Unrestricted model: with break
chow_data$DXY_x_IT <- chow_data$DXY * chow_data$regime_IT
f_unrestricted <- as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY + regime_IT + DXY_x_IT +", control_str))
m_unrestricted <- lm(f_unrestricted, data = chow_data)

# F-test
anova_result <- anova(m_restricted, m_unrestricted)
chow_F <- anova_result$F[2]
chow_p <- anova_result$`Pr(>F)`[2]

# Extract interaction coefficient with HAC standard errors
vcov_unr <- sandwich::NeweyWest(m_unrestricted, lag = 3, prewhite = FALSE)
ct_unr <- lmtest::coeftest(m_unrestricted, vcov = vcov_unr)

# Report
cat("CHOW TEST RESULTS:\n")
cat(sprintf("  F-statistic:      %.3f\n", chow_F))
cat(sprintf("  p-value:          %.4f\n", chow_p))
cat(sprintf("  Result:           %s at 5%% level\n",
            ifelse(chow_p < 0.05, "REJECTS stability", "FAILS TO REJECT stability")))

cat("\nUnrestricted model coefficients (HAC SEs):\n")
for (vn in c("DXY", "regime_IT", "DXY_x_IT")) {
  idx <- which(rownames(ct_unr) == vn)
  if (length(idx) > 0) {
    cat(sprintf("  %-12s: coef = %+.4f, SE = %.4f, t = %+.2f, p = %.4f\n",
                vn, ct_unr[idx, 1], ct_unr[idx, 2], ct_unr[idx, 3], ct_unr[idx, 4]))
  }
}

# Also run the Chow test with HAC-based Wald test
# F-test for joint significance of regime_IT + DXY_x_IT
library(car)
tryCatch({
  wald_result <- linearHypothesis(m_unrestricted, c("regime_IT = 0", "DXY_x_IT = 0"),
                                   vcov. = vcov_unr)
  wald_F <- wald_result$F[2]
  wald_p <- wald_result$`Pr(>F)`[2]
  cat(sprintf("\nHAC-robust Wald test (joint regime_IT + DXY×IT = 0):\n"))
  cat(sprintf("  F-statistic:      %.3f\n", wald_F))
  cat(sprintf("  p-value:          %.4f\n", wald_p))
}, error = function(e) {
  wald_F <- NA; wald_p <- NA
  cat(sprintf("\nWald test failed: %s\n", e$message))
})

# Save
chow_results <- data.frame(
  test = c("Chow_F_test", "HAC_Wald_test"),
  F_stat = c(chow_F, ifelse(exists("wald_F"), wald_F, NA)),
  p_value = c(chow_p, ifelse(exists("wald_p"), wald_p, NA)),
  DXY_coef_preIT = ct_unr[which(rownames(ct_unr) == "DXY"), 1],
  DXY_coef_interaction = ct_unr[which(rownames(ct_unr) == "DXY_x_IT"), 1],
  DXY_interaction_p = ct_unr[which(rownames(ct_unr) == "DXY_x_IT"), 4],
  n_obs = nrow(chow_data),
  n_preIT = sum(chow_data$regime_IT == 0),
  n_postIT = sum(chow_data$regime_IT == 1),
  stringsAsFactors = FALSE
)

write.csv(chow_results, file.path(EXP_CONFIG$csv_dir, "Expert_Chow_Test.csv"),
          row.names = FALSE)
cat("\nSaved: Expert_Chow_Test.csv\n\n")


################################################################################
# PART C: ToT INTERACTION WITH FCI_ENDO_exCredit_exTCN (DOMESTIC ONLY)
################################################################################

cat("================================================================================\n")
cat("PART C: ToT INTERACTION WITH DOMESTIC-ONLY FCI (excl commodities, TCN, credit)\n")
cat("================================================================================\n\n")

fci_var_c <- "FCI_ENDO_exCredit_exTCN_AVG"
interaction_results <- data.frame()

for (h in 1:EXP_CONFIG$max_horizon) {
  data_h <- analysis_data %>%
    mutate(
      y_fwd = lead(Cred_Real_Total, h),
      y_lag1 = lag(Cred_Real_Total, 1),
      y_lag2 = lag(Cred_Real_Total, 2),
      fci_lag1 = lag(!!sym(fci_var_c), 1),
      FCI_x_ToT = !!sym(fci_var_c) * d_ToT
    )

  reg_data <- data_h %>%
    dplyr::select(y_fwd, !!sym(fci_var_c), d_ToT, FCI_x_ToT,
                  y_lag1, y_lag2, fci_lag1,
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    na.omit()

  if (nrow(reg_data) < 50) next

  f <- as.formula(paste("y_fwd ~", fci_var_c, "+ d_ToT + FCI_x_ToT +",
                         "y_lag1 + y_lag2 + fci_lag1 +",
                         "IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +",
                         "IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2"))

  model <- lm(f, data = reg_data)
  vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
  ct <- lmtest::coeftest(model, vcov = vcov_hac)

  for (var_name in c(fci_var_c, "d_ToT", "FCI_x_ToT")) {
    idx <- which(rownames(ct) == var_name)
    if (length(idx) > 0) {
      interaction_results <- rbind(interaction_results, data.frame(
        horizon = h,
        variable = var_name,
        coef = ct[idx, 1],
        se = ct[idx, 2],
        p_value = ct[idx, 4],
        ci_lower = ct[idx, 1] - z_crit * ct[idx, 2],
        ci_upper = ct[idx, 1] + z_crit * ct[idx, 2],
        n_obs = nrow(reg_data),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Report key results
cat("ToT Interaction LP with DOMESTIC-ONLY FCI (7 vars, no commodities/TCN/credit):\n\n")
cat(sprintf("%-30s %8s %10s %10s\n", "Variable", "Horizon", "Coef", "p-value"))
cat(strrep("-", 62), "\n")
for (kh in c(6, 12, 18)) {
  for (vn in c(fci_var_c, "d_ToT", "FCI_x_ToT")) {
    row <- interaction_results %>% filter(horizon == kh & variable == vn)
    if (nrow(row) > 0) {
      label <- ifelse(vn == fci_var_c, "FCI_ENDO_exCred_exTCN",
                      ifelse(vn == "FCI_x_ToT", "FCI × ToT", vn))
      cat(sprintf("%-30s %8d %+10.3f %10.3f\n", label, kh, row$coef[1], row$p_value[1]))
    }
  }
  cat("\n")
}

write.csv(interaction_results, file.path(EXP_CONFIG$csv_dir, "Expert_ToT_Interaction_Domestic.csv"),
          row.names = FALSE)
cat("Saved: Expert_ToT_Interaction_Domestic.csv\n\n")


################################################################################
# SUMMARY
################################################################################

cat("================================================================================\n")
cat("EXPERT FEEDBACK COMPUTATIONS COMPLETE\n")
cat("================================================================================\n\n")

cat("Outputs:\n")
cat("  1. Expert_IV_LP_exTCN.csv           — IV-LP with TCN-excluded FCI\n")
cat("  2. Expert_Chow_Test.csv             — Structural break test at IT adoption\n")
cat("  3. Expert_ToT_Interaction_Domestic.csv — Commodity interaction with domestic FCI\n\n")

# Quick summary for paper integration
cat("KEY RESULTS FOR PAPER:\n\n")

# Part A summary
iv12 <- iv_results_a %>% filter(method == "2SLS_DXY", horizon == 12)
ols12 <- iv_results_a %>% filter(method == "OLS", horizon == 12)
if (nrow(iv12) > 0) {
  cat(sprintf("Part A: First-stage F (h=12) = %.2f %s\n",
              iv12$first_stage_F,
              ifelse(iv12$first_stage_F > 10, "(STRONG: DXY predicts FCI beyond TCN)",
                     ifelse(iv12$first_stage_F > 3, "(MODERATE: some power beyond TCN)",
                            "(WEAK: DXY operates mainly through TCN)"))))
  if (nrow(ols12) > 0) {
    cat(sprintf("  OLS: %.2f pp (p=%.3f)  |  2SLS: %.2f pp\n",
                ols12$coef, ols12$p_value, iv12$coef))
  }
}

# Part B summary
cat(sprintf("\nPart B: Chow test F = %.3f, p = %.4f → %s\n",
            chow_F, chow_p,
            ifelse(chow_p < 0.05, "STRUCTURAL BREAK CONFIRMED",
                   "No formal break detected")))

# Part C summary
int12 <- interaction_results %>% filter(horizon == 12, variable == "FCI_x_ToT")
if (nrow(int12) > 0) {
  cat(sprintf("\nPart C: FCI×ToT interaction (h=12) = %.3f (p=%.3f) → %s\n",
              int12$coef, int12$p_value,
              ifelse(int12$p_value < 0.10, "INTERACTION SURVIVES with domestic-only FCI",
                     "Interaction weakens (partly reflects commodity composition)")))
}

cat("\n")
