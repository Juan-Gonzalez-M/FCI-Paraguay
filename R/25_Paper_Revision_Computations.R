################################################################################
# POST-IT FIRST-STAGE DIAGNOSTICS & ADDITIONAL ROBUSTNESS CHECKS
################################################################################
#
# Project:      Financial Conditions Index - Paper Revision Computations
# Institution:  Banco Central del Paraguay
#
# Description:  Computes two items needed for the paper revision:
#   PART A: Post-IT first-stage F for DXY → FCI_ENDO_exCredit
#   PART B: Additional robustness checks (1-lag, ex-COVID, ex-GFC)
#
# Output:
#   - output/csv/PostIT_First_Stage_Diagnostics.csv
#   - output/csv/Robustness_Additional_Checks.csv
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

REV_CONFIG <- list(
  output_dir = "../output",
  csv_dir    = "../output/csv",
  verbose    = TRUE
)

find_csv <- function(filename) {
  p1 <- file.path(REV_CONFIG$csv_dir, filename)
  if (file.exists(p1)) return(p1)
  p2 <- file.path(REV_CONFIG$output_dir, filename)
  if (file.exists(p2)) return(p2)
  return(p2)
}

set.seed(20260316)

cat("\n################################################################################\n")
cat("POST-IT FIRST-STAGE DIAGNOSTICS & ROBUSTNESS CHECKS\n")
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

# ---- Load credit data ----
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

cat("Full analysis data:", nrow(analysis_data), "obs\n\n")


################################################################################
# PART A: POST-IT FIRST-STAGE DIAGNOSTICS
################################################################################

cat("================================================================================\n")
cat("PART A: POST-IT FIRST-STAGE DIAGNOSTICS\n")
cat("================================================================================\n\n")

# Post-IT: May 2011+
post_it_data <- analysis_data %>% filter(fecha >= as.Date("2011-05-01"))
cat("Post-IT sample:", nrow(post_it_data), "obs (from",
    as.character(min(post_it_data$fecha)), "to",
    as.character(max(post_it_data$fecha)), ")\n\n")

control_vars <- c("IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                   "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2")

# ---- A1: Unconditional first-stage (DXY only, with macro controls) ----
cat("--- A1: Unconditional First-Stage (Post-IT) ---\n\n")

fs_data <- post_it_data %>%
  dplyr::select(FCI_ENDO_exCredit_AVG, DXY, all_of(control_vars)) %>%
  na.omit()

cat("First-stage sample (post-IT):", nrow(fs_data), "obs\n")

# Full model with DXY
f_full <- as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY +",
                            paste(control_vars, collapse = " + ")))
m_full <- lm(f_full, data = fs_data)
vcov_hac <- sandwich::NeweyWest(m_full, lag = 3, prewhite = FALSE)
ct <- lmtest::coeftest(m_full, vcov = vcov_hac)
idx_dxy <- which(rownames(ct) == "DXY")

# Controls-only model for partial R²
f_base <- as.formula(paste("FCI_ENDO_exCredit_AVG ~",
                            paste(control_vars, collapse = " + ")))
m_base <- lm(f_base, data = fs_data)
partial_r2 <- summary(m_full)$r.squared - summary(m_base)$r.squared

# F-stat (from ANOVA)
f_stat_anova <- anova(m_base, m_full)$F[2]

# Also compute t²
t_stat <- ct[idx_dxy, 3]
f_stat_t2 <- t_stat^2

cat(sprintf("  DXY coefficient:  %+.4f\n", ct[idx_dxy, 1]))
cat(sprintf("  t-statistic:      %+.2f\n", t_stat))
cat(sprintf("  p-value:          %.4f\n", ct[idx_dxy, 4]))
cat(sprintf("  Partial R²:       %.4f (%.1f%%)\n", partial_r2, partial_r2 * 100))
cat(sprintf("  F-stat (t²):      %.2f\n", f_stat_t2))
cat(sprintf("  F-stat (ANOVA):   %.2f\n", f_stat_anova))
cat("\n")

# ---- A2: Conditional first-stage at h=6,12,18 (in full LP spec) ----
cat("--- A2: Conditional First-Stage at Key Horizons (Post-IT) ---\n\n")

er_controls <- intersect(c("d_ToT", "US_10Y", "SP500"), names(analysis_data))
horizons_check <- c(6, 12, 18)

conditional_results <- data.frame()

for (h in horizons_check) {
  data_h <- post_it_data %>%
    mutate(
      y_fwd = lead(Cred_Real_Total, h),
      y_lag1 = lag(Cred_Real_Total, 1),
      y_lag2 = lag(Cred_Real_Total, 2),
      fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1)
    )

  exog_vars <- c("y_lag1", "y_lag2", "fci_lag1",
                  "IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                  "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2",
                  er_controls)
  exog_str <- paste(exog_vars, collapse = " + ")

  reg_data <- data_h %>%
    dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG,
                  all_of(exog_vars), DXY) %>%
    na.omit()

  # IV regression to get conditional first-stage F
  iv_formula <- as.formula(paste(
    "y_fwd ~ FCI_ENDO_exCredit_AVG +", exog_str,
    "| DXY +", exog_str
  ))

  tryCatch({
    iv_model <- ivreg::ivreg(iv_formula, data = reg_data)
    vcov_hac_iv <- sandwich::NeweyWest(iv_model, lag = h + 1, prewhite = FALSE)
    iv_summ <- summary(iv_model, vcov = vcov_hac_iv, diagnostics = TRUE)

    cond_F <- tryCatch(
      iv_summ$diagnostics["Weak instruments", "statistic"], error = function(e) NA)

    # Also run manual first-stage for t-stat and partial R²
    fs_formula_full <- as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY +", exog_str))
    fs_formula_base <- as.formula(paste("FCI_ENDO_exCredit_AVG ~", exog_str))
    m_fs_full <- lm(fs_formula_full, data = reg_data)
    m_fs_base <- lm(fs_formula_base, data = reg_data)
    vcov_fs <- sandwich::NeweyWest(m_fs_full, lag = h + 1, prewhite = FALSE)
    ct_fs <- lmtest::coeftest(m_fs_full, vcov = vcov_fs)
    idx_d <- which(rownames(ct_fs) == "DXY")
    pr2 <- summary(m_fs_full)$r.squared - summary(m_fs_base)$r.squared

    conditional_results <- rbind(conditional_results, data.frame(
      horizon = h,
      n_obs = nrow(reg_data),
      cond_F = cond_F,
      DXY_coef = ct_fs[idx_d, 1],
      DXY_t = ct_fs[idx_d, 3],
      DXY_p = ct_fs[idx_d, 4],
      partial_R2 = pr2,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  h=%d: N=%d, Cond F=%.2f, DXY t=%.2f, partial R²=%.3f\n",
                h, nrow(reg_data), cond_F, ct_fs[idx_d, 3], pr2))

  }, error = function(e) {
    cat(sprintf("  h=%d: ERROR: %s\n", h, e$message))
  })
}

cat("\n")

# ---- Save Part A results ----
post_it_diag <- data.frame(
  type = c("unconditional", rep("conditional", nrow(conditional_results))),
  horizon = c(NA, conditional_results$horizon),
  sample = "post-IT (May 2011+)",
  n_obs = c(nrow(fs_data), conditional_results$n_obs),
  DXY_coef = c(ct[idx_dxy, 1], conditional_results$DXY_coef),
  DXY_t = c(t_stat, conditional_results$DXY_t),
  DXY_p = c(ct[idx_dxy, 4], conditional_results$DXY_p),
  partial_R2 = c(partial_r2, conditional_results$partial_R2),
  F_stat = c(f_stat_t2, conditional_results$cond_F),
  stringsAsFactors = FALSE
)

write.csv(post_it_diag, file.path(REV_CONFIG$csv_dir, "PostIT_First_Stage_Diagnostics.csv"),
          row.names = FALSE)
cat("Saved: PostIT_First_Stage_Diagnostics.csv\n\n")


################################################################################
# PART B: ADDITIONAL ROBUSTNESS CHECKS
################################################################################

cat("================================================================================\n")
cat("PART B: ADDITIONAL ROBUSTNESS CHECKS\n")
cat("================================================================================\n\n")

# Common LP function
run_robustness_lp <- function(data, fci_var, y_var, n_lags, max_h = 18,
                               label = "baseline") {
  z_crit <- qnorm(0.95)
  results <- data.frame()

  for (h in 1:max_h) {
    data_h <- data %>%
      mutate(y_fwd = lead(!!sym(y_var), h))

    # Create lags dynamically
    for (j in 1:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", j) := lag(!!sym(y_var), j))
    }
    data_h <- data_h %>%
      mutate(fci_lag1 = lag(!!sym(fci_var), 1))

    lag_vars <- c(paste0("y_lag", 1:n_lags), "fci_lag1")
    ctrl_vars <- c("IMAEP_yoy", "IPC_yoy")
    # Add macro control lags
    for (cv in ctrl_vars) {
      for (j in 1:n_lags) {
        cv_lag <- paste0(cv, "_L", j)
        if (cv_lag %in% names(data_h)) {
          lag_vars <- c(lag_vars, cv_lag)
        } else {
          data_h <- data_h %>% mutate(!!cv_lag := lag(!!sym(cv), j))
          lag_vars <- c(lag_vars, cv_lag)
        }
      }
    }

    all_rhs <- c(fci_var, lag_vars)
    formula_str <- paste("y_fwd ~", paste(all_rhs, collapse = " + "))

    reg_data <- data_h %>%
      dplyr::select(y_fwd, all_of(all_rhs)) %>%
      na.omit()

    if (nrow(reg_data) < 30) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    ct <- lmtest::coeftest(model, vcov = vcov_hac)
    idx <- which(rownames(ct) == fci_var)

    if (length(idx) == 0) next

    results <- rbind(results, data.frame(
      specification = label,
      horizon = h,
      coef = ct[idx, 1],
      se = ct[idx, 2],
      p_value = ct[idx, 4],
      ci_lower = ct[idx, 1] - z_crit * ct[idx, 2],
      ci_upper = ct[idx, 1] + z_crit * ct[idx, 2],
      n_obs = nrow(reg_data),
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}

robustness_all <- data.frame()

# Use FCI_exCredit_AVG (comprehensive 11-var, excl credit) as baseline
fci_var_use <- "FCI_exCredit_AVG"
if (!fci_var_use %in% names(analysis_data)) {
  fci_var_use <- "FCI_ENDO_exCredit_AVG"
  cat("  Fallback: using FCI_ENDO_exCredit_AVG\n")
}
cat("  FCI variable:", fci_var_use, "\n")

# ---- B1: Baseline (2 lags, full sample) ----
cat("B1: Baseline (2 lags, full sample)...\n")
r_baseline <- run_robustness_lp(analysis_data, fci_var_use,
                                 "Cred_Real_Total", n_lags = 2, label = "Baseline (2 lags)")
robustness_all <- rbind(robustness_all, r_baseline)
if (nrow(r_baseline) > 0) {
  h12 <- r_baseline %>% filter(horizon == 12)
  if (nrow(h12) > 0) cat(sprintf("  h=12: %.2f pp, p=%.3f, N=%d\n", h12$coef, h12$p_value, h12$n_obs))
}

# ---- B2: Alternative lag (1 lag instead of 2) ----
cat("B2: Alternative lag (1 lag)...\n")
r_1lag <- run_robustness_lp(analysis_data, fci_var_use,
                             "Cred_Real_Total", n_lags = 1, label = "1 lag (instead of 2)")
robustness_all <- rbind(robustness_all, r_1lag)
if (nrow(r_1lag) > 0) {
  h12 <- r_1lag %>% filter(horizon == 12)
  if (nrow(h12) > 0) cat(sprintf("  h=12: %.2f pp, p=%.3f, N=%d\n", h12$coef, h12$p_value, h12$n_obs))
}

# ---- B3: Exclude COVID (Mar 2020 – Dec 2021) ----
cat("B3: Exclude COVID (Mar 2020 - Dec 2021)...\n")
data_excovid <- analysis_data %>%
  filter(!(fecha >= as.Date("2020-03-01") & fecha <= as.Date("2021-12-31")))
cat("  Obs after excluding COVID:", nrow(data_excovid), "\n")

r_excovid <- run_robustness_lp(data_excovid, fci_var_use,
                                "Cred_Real_Total", n_lags = 2, label = "Exclude COVID")
robustness_all <- rbind(robustness_all, r_excovid)
if (nrow(r_excovid) > 0) {
  h12 <- r_excovid %>% filter(horizon == 12)
  if (nrow(h12) > 0) cat(sprintf("  h=12: %.2f pp, p=%.3f, N=%d\n", h12$coef, h12$p_value, h12$n_obs))
}

# ---- B4: Exclude GFC (Sep 2008 – Jun 2009) ----
cat("B4: Exclude GFC (Sep 2008 - Jun 2009)...\n")
data_exgfc <- analysis_data %>%
  filter(!(fecha >= as.Date("2008-09-01") & fecha <= as.Date("2009-06-30")))
cat("  Obs after excluding GFC:", nrow(data_exgfc), "\n")

r_exgfc <- run_robustness_lp(data_exgfc, fci_var_use,
                              "Cred_Real_Total", n_lags = 2, label = "Exclude GFC")
robustness_all <- rbind(robustness_all, r_exgfc)
if (nrow(r_exgfc) > 0) {
  h12 <- r_exgfc %>% filter(horizon == 12)
  if (nrow(h12) > 0) cat(sprintf("  h=12: %.2f pp, p=%.3f, N=%d\n", h12$coef, h12$p_value, h12$n_obs))
}

cat("\n")

# ---- Summary at h=12 ----
cat("ROBUSTNESS SUMMARY AT h=12:\n")
cat(sprintf("%-30s %10s %10s %6s\n", "Specification", "Coeff (pp)", "p-value", "N"))
cat(strrep("-", 60), "\n")
for (spec in unique(robustness_all$specification)) {
  h12 <- robustness_all %>% filter(specification == spec, horizon == 12)
  if (nrow(h12) > 0) {
    cat(sprintf("%-30s %+10.2f %10.3f %6d\n",
                spec, h12$coef, h12$p_value, h12$n_obs))
  }
}
cat("\n")

# ---- Save Part B results ----
write.csv(robustness_all, file.path(REV_CONFIG$csv_dir, "Robustness_Additional_Checks.csv"),
          row.names = FALSE)
cat("Saved: Robustness_Additional_Checks.csv\n\n")

cat("################################################################################\n")
cat("REVISION COMPUTATIONS COMPLETE\n")
cat("################################################################################\n\n")
