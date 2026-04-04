################################################################################
# FCI REGIONAL SPILLOVER ANALYSIS
################################################################################
#
# Project:      Financial Conditions Index - Brazil Selic Regional Spillovers
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Tests whether Brazil's monetary policy (Selic) creates financial
#               contagion in Paraguay beyond common global shocks. Distinguishes
#               regional contagion from common exposure to global factors.
#
#   SECTION A: Direct Spillover Test
#   SECTION B: Contagion vs Common Shocks
#   SECTION C: Granger Causality
#   SECTION D: Asymmetric Regional Spillovers
#   SECTION E: Regional Spillover LP
#   SECTION F: Bivariate & Controlled VAR
#
# References:
#   - Forbes & Rigobon (2002) - No Contagion, Only Interdependence
#   - Kaminsky, Reinhart & Vegh (2003) - The Unholy Trinity of Financial Contagion
#   - Canova (2005) - The Transmission of US Shocks to Latin America
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
  library(vars)
  library(ggplot2)
  library(gridExtra)
  library(sandwich)
  library(lmtest)
})

REG_CONFIG <- list(
  max_horizon     = 24,
  max_var_lags    = 6,
  n_lags          = 2,
  confidence_level = 0.90,
  n_bootstrap     = 500,
  output_dir      = "../output",
  verbose         = TRUE
)

set.seed(20260211)

cat("\n################################################################################\n")
cat("FCI REGIONAL SPILLOVER ANALYSIS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---- Load FCI ----
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  CONFIG_SAVE <- if (exists("CONFIG")) CONFIG else NULL
  source("01_FCI_Complete.R")
  if (!is.null(CONFIG_SAVE)) CONFIG <- CONFIG_SAVE
  rm(CONFIG_SAVE)
  suppressPackageStartupMessages(library(dplyr))
}

# ---- Load external data ----
ext_file <- file.path(REG_CONFIG$output_dir, "New_External_Variables.csv")
if (!file.exists(ext_file)) {
  cat("New_External_Variables.csv not found. Running script 17...\n")
  source("17_FCI_New_External_Data.R")
}
ext_data <- read.csv(ext_file)
ext_data$fecha <- as.Date(ext_data$fecha)

# ---- Load credit data ----
macro_raw <- read_excel("../data/FCI_data_1.xlsx", sheet = "Datos_macro")
fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]
macro_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  mutate(
    IMAEP_yoy = (IMAEP / lag(IMAEP, 12) - 1) * 100,
    IPC_yoy   = (IPC / lag(IPC, 12) - 1) * 100,
    Cred_Real_Total = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100
  )

analysis_data <- ext_data %>%
  left_join(macro_data %>% dplyr::select(fecha, any_of(c("IMAEP_yoy", "IPC_yoy", "Cred_Real_Total"))),
            by = "fecha", suffix = c("", ".y")) %>%
  arrange(fecha)

# Clean duplicates
for (col in names(analysis_data)) {
  if (grepl("\\.y$", col)) {
    base_col <- sub("\\.y$", "", col)
    if (base_col %in% names(analysis_data)) {
      na_mask <- is.na(analysis_data[[base_col]])
      analysis_data[[base_col]][na_mask] <- analysis_data[[col]][na_mask]
    }
    analysis_data[[col]] <- NULL
  }
}

# Create lags
analysis_data <- analysis_data %>%
  mutate(
    IMAEP_yoy_L1 = lag(IMAEP_yoy, 1),
    IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
    IPC_yoy_L1 = lag(IPC_yoy, 1),
    IPC_yoy_L2 = lag(IPC_yoy, 2),
    d_Selic = Selic_rate - lag(Selic_rate, 1),
    d_Selic_pos = pmax(d_Selic, 0),
    d_Selic_neg = pmin(d_Selic, 0)
  )

z_crit <- qnorm(1 - (1 - REG_CONFIG$confidence_level) / 2)
cat("Analysis data:", nrow(analysis_data), "obs\n\n")


################################################################################
# SECTION A: DIRECT SPILLOVER TEST
################################################################################

cat("================================================================================\n")
cat("SECTION A: DIRECT SPILLOVER TEST\n")
cat("================================================================================\n\n")

control_str <- "IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 + IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2"

# Simple: Selic + VIX + FFER
spillover_results <- data.frame()

if ("Selic_rate" %in% names(analysis_data) && "VIX" %in% names(analysis_data)) {
  f1 <- as.formula(paste("FCI_ENDO_AVG ~ Selic_rate + VIX + FFER +", control_str))
  reg1 <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, Selic_rate, VIX, FFER,
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    na.omit()

  if (nrow(reg1) >= 30) {
    m1 <- lm(f1, data = reg1)
    vcov1 <- sandwich::NeweyWest(m1, lag = 3, prewhite = FALSE)
    ct1 <- lmtest::coeftest(m1, vcov = vcov1)

    cat("Direct Spillover Test (FCI ~ Selic + VIX + FFER + controls):\n")
    for (v in c("Selic_rate", "VIX", "FFER")) {
      idx <- which(rownames(ct1) == v)
      cat(sprintf("  %-12s: coef = %+.4f, t = %+.2f, p = %.3f\n",
                  v, ct1[idx, 1], ct1[idx, 3], ct1[idx, 4]))
      spillover_results <- rbind(spillover_results, data.frame(
        spec = "Direct (3 global)", variable = v,
        coef = ct1[idx, 1], se = ct1[idx, 2], t_stat = ct1[idx, 3], p_value = ct1[idx, 4],
        n_obs = nrow(reg1), stringsAsFactors = FALSE
      ))
    }
    cat(sprintf("  R² = %.4f\n\n", summary(m1)$r.squared))
  }
}


################################################################################
# SECTION B: CONTAGION vs COMMON SHOCKS
################################################################################

cat("================================================================================\n")
cat("SECTION B: CONTAGION vs COMMON SHOCKS\n")
cat("================================================================================\n\n")

# Full specification with ALL global variables
global_vars_seq <- c("VIX", "FFER", "DXY", "SP500", "US_10Y")
available_global <- intersect(global_vars_seq, names(analysis_data))

# Sequential addition to track Selic coefficient stability
selic_stability <- data.frame()

for (k in 0:length(available_global)) {
  if (k == 0) {
    glob_str <- ""
    spec_label <- "Selic only"
    f <- as.formula(paste("FCI_ENDO_AVG ~ Selic_rate +", control_str))
  } else {
    glob_subset <- available_global[1:k]
    glob_str <- paste(glob_subset, collapse = " + ")
    spec_label <- paste("+ ", paste(glob_subset, collapse = ", "))
    f <- as.formula(paste("FCI_ENDO_AVG ~ Selic_rate +", glob_str, "+", control_str))
  }

  reg_vars <- c("FCI_ENDO_AVG", "Selic_rate",
                 if (k > 0) available_global[1:k] else NULL,
                 "IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                 "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2")
  reg_vars <- intersect(reg_vars, names(analysis_data))

  reg_k <- analysis_data %>%
    dplyr::select(all_of(reg_vars)) %>%
    na.omit()

  if (nrow(reg_k) < 30) next

  m_k <- lm(f, data = reg_k)
  vcov_k <- sandwich::NeweyWest(m_k, lag = 3, prewhite = FALSE)
  ct_k <- lmtest::coeftest(m_k, vcov = vcov_k)
  idx <- which(rownames(ct_k) == "Selic_rate")

  selic_stability <- rbind(selic_stability, data.frame(
    step = k,
    spec = spec_label,
    selic_coef = ct_k[idx, 1],
    selic_se = ct_k[idx, 2],
    selic_t = ct_k[idx, 3],
    selic_p = ct_k[idx, 4],
    R2 = summary(m_k)$r.squared,
    n_obs = nrow(reg_k),
    stringsAsFactors = FALSE
  ))
}

cat("Selic Coefficient Stability (Sequential Global Controls):\n")
cat(sprintf("%-35s %10s %10s %10s\n", "Specification", "Selic β", "t-stat", "p-value"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(selic_stability))) {
  cat(sprintf("%-35s %+10.4f %10.2f %10.3f\n",
              selic_stability$spec[i],
              selic_stability$selic_coef[i],
              selic_stability$selic_t[i],
              selic_stability$selic_p[i]))
}

# Key test: is Selic significant in full specification?
full_spec <- selic_stability[nrow(selic_stability), ]
if (full_spec$selic_p < 0.10) {
  cat("\n→ CONTAGION: Selic retains significance after ALL global controls (p=",
      round(full_spec$selic_p, 3), ")\n\n")
} else {
  cat("\n→ COMMON SHOCKS: Selic loses significance with full global controls (p=",
      round(full_spec$selic_p, 3), ")\n\n")
}


################################################################################
# SECTION C: GRANGER CAUSALITY
################################################################################

cat("================================================================================\n")
cat("SECTION C: GRANGER CAUSALITY\n")
cat("================================================================================\n\n")

granger_results <- data.frame()

# ---- C1: Bivariate Granger ----
bivar_data <- analysis_data %>%
  dplyr::select(FCI_ENDO_AVG, Selic_rate) %>%
  na.omit()

for (n_lag in c(3, 6, 12)) {
  if (nrow(bivar_data) <= n_lag + 10) next
  tryCatch({
    var_g <- VAR(bivar_data, p = n_lag, type = "const")

    # Selic → FCI
    g1 <- causality(var_g, cause = "Selic_rate")
    granger_results <- rbind(granger_results, data.frame(
      test = "Bivariate", direction = "Selic → FCI",
      lag = n_lag, F_stat = g1$Granger$statistic, p_value = g1$Granger$p.value,
      n_obs = nrow(bivar_data) - n_lag, stringsAsFactors = FALSE
    ))

    # FCI → Selic (should be insignificant)
    g2 <- causality(var_g, cause = "FCI_ENDO_AVG")
    granger_results <- rbind(granger_results, data.frame(
      test = "Bivariate", direction = "FCI → Selic",
      lag = n_lag, F_stat = g2$Granger$statistic, p_value = g2$Granger$p.value,
      n_obs = nrow(bivar_data) - n_lag, stringsAsFactors = FALSE
    ))
  }, error = function(e) NULL)
}

# ---- C2: Conditional Granger (controlling for global) ----
cond_vars <- intersect(c("FCI_ENDO_AVG", "Selic_rate", "VIX", "FFER", "DXY"),
                        names(analysis_data))
cond_data <- analysis_data %>%
  dplyr::select(all_of(cond_vars)) %>%
  na.omit()

for (n_lag in c(3, 6, 12)) {
  if (nrow(cond_data) <= n_lag * length(cond_vars) + 10) next
  tryCatch({
    var_cond <- VAR(cond_data, p = n_lag, type = "const")
    g_cond <- causality(var_cond, cause = "Selic_rate")
    granger_results <- rbind(granger_results, data.frame(
      test = "Conditional", direction = "Selic → FCI | Global",
      lag = n_lag, F_stat = g_cond$Granger$statistic, p_value = g_cond$Granger$p.value,
      n_obs = nrow(cond_data) - n_lag, stringsAsFactors = FALSE
    ))
  }, error = function(e) NULL)
}

cat("Granger Causality Results:\n")
cat(sprintf("%-15s %-25s %5s %10s %10s\n", "Test", "Direction", "Lag", "F-stat", "p-value"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(granger_results))) {
  stars <- ifelse(granger_results$p_value[i] < 0.01, "***",
                  ifelse(granger_results$p_value[i] < 0.05, "**",
                         ifelse(granger_results$p_value[i] < 0.10, "*", "")))
  cat(sprintf("%-15s %-25s %5d %10.2f %10.3f %s\n",
              granger_results$test[i], granger_results$direction[i],
              granger_results$lag[i], granger_results$F_stat[i],
              granger_results$p_value[i], stars))
}
cat("\n")


################################################################################
# SECTION D: ASYMMETRIC REGIONAL SPILLOVERS
################################################################################

cat("================================================================================\n")
cat("SECTION D: ASYMMETRIC REGIONAL SPILLOVERS\n")
cat("================================================================================\n\n")

asym_results <- data.frame()

if ("d_Selic_pos" %in% names(analysis_data) && "d_Selic_neg" %in% names(analysis_data)) {
  # Include global controls
  glob_vars <- intersect(c("VIX", "FFER", "DXY"), names(analysis_data))
  glob_str <- paste(glob_vars, collapse = " + ")

  f_asym <- as.formula(paste("FCI_ENDO_AVG ~ d_Selic_pos + I(abs(d_Selic_neg)) +",
                              glob_str, "+", control_str))

  reg_asym <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, d_Selic_pos, d_Selic_neg,
                  all_of(glob_vars),
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    mutate(abs_d_Selic_neg = abs(d_Selic_neg)) %>%
    na.omit()

  if (nrow(reg_asym) >= 30) {
    # Redo with explicit abs column
    f_asym2 <- as.formula(paste("FCI_ENDO_AVG ~ d_Selic_pos + abs_d_Selic_neg +",
                                 glob_str, "+", control_str))
    m_asym <- lm(f_asym2, data = reg_asym)
    vcov_asym <- sandwich::NeweyWest(m_asym, lag = 3, prewhite = FALSE)
    ct_asym <- lmtest::coeftest(m_asym, vcov = vcov_asym)

    idx_pos <- which(rownames(ct_asym) == "d_Selic_pos")
    idx_neg <- which(rownames(ct_asym) == "abs_d_Selic_neg")

    cat("Asymmetric Selic Spillovers:\n")
    cat(sprintf("  Tightening (ΔSelic > 0): β⁺ = %+.4f, t = %.2f, p = %.3f\n",
                ct_asym[idx_pos, 1], ct_asym[idx_pos, 3], ct_asym[idx_pos, 4]))
    cat(sprintf("  Easing (ΔSelic < 0):     β⁻ = %+.4f, t = %.2f, p = %.3f\n",
                ct_asym[idx_neg, 1], ct_asym[idx_neg, 3], ct_asym[idx_neg, 4]))

    # Wald test: H0: β⁺ = β⁻
    R <- matrix(0, 1, length(coef(m_asym)))
    R[1, idx_pos] <- 1
    R[1, idx_neg] <- -1
    wald_stat <- (R %*% coef(m_asym))^2 / (R %*% vcov_asym %*% t(R))
    wald_p <- 1 - pchisq(wald_stat, 1)
    cat(sprintf("  Wald test (β⁺ = β⁻): χ² = %.2f, p = %.3f\n\n",
                wald_stat[1, 1], wald_p[1, 1]))

    asym_results <- data.frame(
      beta_pos = ct_asym[idx_pos, 1],
      beta_pos_p = ct_asym[idx_pos, 4],
      beta_neg = ct_asym[idx_neg, 1],
      beta_neg_p = ct_asym[idx_neg, 4],
      wald_stat = wald_stat[1, 1],
      wald_p = wald_p[1, 1]
    )
  }
}


################################################################################
# SECTION E: REGIONAL SPILLOVER LP
################################################################################

cat("================================================================================\n")
cat("SECTION E: REGIONAL SPILLOVER LP\n")
cat("================================================================================\n\n")

spillover_lp <- data.frame()

# LP: ΔSelic → FCI_ENDO and → Credit
targets <- c("FCI_ENDO_AVG", "Cred_Real_Total")
glob_vars <- intersect(c("VIX", "FFER", "DXY"), names(analysis_data))

for (target in targets) {
  cat("LP: ΔSelic →", target, "\n")

  for (h in 1:REG_CONFIG$max_horizon) {
    data_h <- analysis_data %>%
      mutate(
        y_fwd = lead(!!sym(target), h),
        y_lag1 = lag(!!sym(target), 1),
        y_lag2 = lag(!!sym(target), 2),
        selic_lag1 = lag(d_Selic, 1)
      )

    reg_vars <- c("y_fwd", "d_Selic", "y_lag1", "y_lag2", "selic_lag1",
                   glob_vars,
                   "IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                   "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2")
    available_reg <- intersect(reg_vars, names(data_h))

    reg_data <- data_h %>%
      dplyr::select(all_of(available_reg)) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    f <- as.formula(paste("y_fwd ~ d_Selic + y_lag1 + y_lag2 + selic_lag1 +",
                           paste(glob_vars, collapse = " + "), "+",
                           "IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +",
                           "IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2"))

    model <- lm(f, data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    ct <- lmtest::coeftest(model, vcov = vcov_hac)
    idx <- which(rownames(ct) == "d_Selic")

    if (length(idx) > 0) {
      spillover_lp <- rbind(spillover_lp, data.frame(
        target = target,
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
  }
  cat("  Done.\n")
}

# Report key results
cat("\nSpillover LP Results at key horizons:\n")
for (target in targets) {
  cat(sprintf("\n  %s:\n", target))
  for (kh in c(6, 12, 18, 24)) {
    row <- spillover_lp %>% filter(target == !!target & horizon == kh)
    if (nrow(row) > 0) {
      cat(sprintf("    h=%2d: β = %+.3f, p = %.3f\n", kh, row$coef[1], row$p_value[1]))
    }
  }
}
cat("\n")


################################################################################
# SECTION F: BIVARIATE & CONTROLLED VAR
################################################################################

cat("================================================================================\n")
cat("SECTION F: VAR ANALYSIS\n")
cat("================================================================================\n\n")

var_irf_results <- data.frame()

# ---- F1: Bivariate VAR ----
cat("Bivariate VAR: [Selic, FCI_ENDO]\n")
bivar_vars <- c("Selic_rate", "FCI_ENDO_AVG")
bivar_data <- analysis_data %>%
  dplyr::select(fecha, all_of(bivar_vars)) %>%
  na.omit()

bivar_ts <- bivar_data[, bivar_vars]
bivar_select <- VARselect(bivar_ts, lag.max = REG_CONFIG$max_var_lags, type = "const")
bivar_lag <- max(1, min(bivar_select$selection["AIC(n)"], REG_CONFIG$max_var_lags))
bivar_model <- VAR(bivar_ts, p = bivar_lag, type = "const")

tryCatch({
  irf_bivar <- irf(bivar_model, impulse = "Selic_rate", response = "FCI_ENDO_AVG",
                     n.ahead = REG_CONFIG$max_horizon, boot = TRUE,
                     ci = REG_CONFIG$confidence_level, runs = REG_CONFIG$n_bootstrap)

  bivar_irf_df <- data.frame(
    specification = "Bivariate",
    horizon = 0:REG_CONFIG$max_horizon,
    irf_val = as.vector(irf_bivar$irf$Selic_rate),
    lower = as.vector(irf_bivar$Lower$Selic_rate),
    upper = as.vector(irf_bivar$Upper$Selic_rate),
    stringsAsFactors = FALSE
  )
  var_irf_results <- rbind(var_irf_results, bivar_irf_df)

  peak_idx <- which.max(abs(bivar_irf_df$irf_val))
  cat(sprintf("  Peak: %.4f at h=%d\n", bivar_irf_df$irf_val[peak_idx], bivar_irf_df$horizon[peak_idx]))

  # FEVD
  fevd_bivar <- fevd(bivar_model, n.ahead = REG_CONFIG$max_horizon)
  fci_fevd <- as.data.frame(fevd_bivar$FCI_ENDO_AVG)
  cat(sprintf("  FCI variance from Selic: h=6: %.1f%%, h=12: %.1f%%, h=24: %.1f%%\n",
              fci_fevd$Selic_rate[6] * 100, fci_fevd$Selic_rate[12] * 100,
              fci_fevd$Selic_rate[min(24, nrow(fci_fevd))] * 100))
}, error = function(e) cat("  IRF failed:", e$message, "\n"))

# ---- F2: Controlled VAR ----
cat("\nControlled VAR: [Selic, VIX, FFER, FCI_ENDO]\n")
ctrl_vars <- intersect(c("Selic_rate", "VIX", "FFER", "FCI_ENDO_AVG"), names(analysis_data))
ctrl_data <- analysis_data %>%
  dplyr::select(fecha, all_of(ctrl_vars)) %>%
  na.omit()

ctrl_ts <- ctrl_data[, ctrl_vars]

tryCatch({
  ctrl_select <- VARselect(ctrl_ts, lag.max = REG_CONFIG$max_var_lags, type = "const")
  ctrl_lag <- max(1, min(ctrl_select$selection["AIC(n)"], REG_CONFIG$max_var_lags))
  ctrl_model <- VAR(ctrl_ts, p = ctrl_lag, type = "const")

  irf_ctrl <- irf(ctrl_model, impulse = "Selic_rate", response = "FCI_ENDO_AVG",
                    n.ahead = REG_CONFIG$max_horizon, boot = TRUE,
                    ci = REG_CONFIG$confidence_level, runs = REG_CONFIG$n_bootstrap)

  ctrl_irf_df <- data.frame(
    specification = "Controlled",
    horizon = 0:REG_CONFIG$max_horizon,
    irf_val = as.vector(irf_ctrl$irf$Selic_rate),
    lower = as.vector(irf_ctrl$Lower$Selic_rate),
    upper = as.vector(irf_ctrl$Upper$Selic_rate),
    stringsAsFactors = FALSE
  )
  var_irf_results <- rbind(var_irf_results, ctrl_irf_df)

  peak_idx_c <- which.max(abs(ctrl_irf_df$irf_val))
  cat(sprintf("  Peak: %.4f at h=%d\n", ctrl_irf_df$irf_val[peak_idx_c], ctrl_irf_df$horizon[peak_idx_c]))

  # FEVD
  fevd_ctrl <- fevd(ctrl_model, n.ahead = REG_CONFIG$max_horizon)
  fci_fevd_c <- as.data.frame(fevd_ctrl$FCI_ENDO_AVG)
  if ("Selic_rate" %in% names(fci_fevd_c)) {
    cat(sprintf("  FCI variance from Selic (controlled): h=6: %.1f%%, h=12: %.1f%%\n",
                fci_fevd_c$Selic_rate[6] * 100, fci_fevd_c$Selic_rate[12] * 100))
  }
}, error = function(e) cat("  Controlled VAR failed:", e$message, "\n"))
cat("\n")


################################################################################
# SAVE OUTPUTS
################################################################################

cat("================================================================================\n")
cat("SAVING OUTPUTS\n")
cat("================================================================================\n\n")

# Combine spillover regression results
all_regression_results <- rbind(
  spillover_results,
  if (nrow(selic_stability) > 0) {
    data.frame(
      spec = paste("Sequential:", selic_stability$spec),
      variable = "Selic_rate",
      coef = selic_stability$selic_coef,
      se = selic_stability$selic_se,
      t_stat = selic_stability$selic_t,
      p_value = selic_stability$selic_p,
      n_obs = selic_stability$n_obs,
      stringsAsFactors = FALSE
    )
  } else data.frame()
)

write.csv(all_regression_results, file.path(REG_CONFIG$output_dir, "Regional_Spillover_Regressions.csv"),
          row.names = FALSE)
cat("Saved: Regional_Spillover_Regressions.csv\n")

write.csv(granger_results, file.path(REG_CONFIG$output_dir, "Regional_Granger_Results.csv"),
          row.names = FALSE)
cat("Saved: Regional_Granger_Results.csv\n")

write.csv(spillover_lp, file.path(REG_CONFIG$output_dir, "Regional_Spillover_LP.csv"),
          row.names = FALSE)
cat("Saved: Regional_Spillover_LP.csv\n")

write.csv(var_irf_results, file.path(REG_CONFIG$output_dir, "Regional_Spillover_VAR_IRF.csv"),
          row.names = FALSE)
cat("Saved: Regional_Spillover_VAR_IRF.csv\n")


################################################################################
# VISUALIZATIONS
################################################################################

cat("\n--- Generating Visualizations ---\n\n")

# ---- 241: Direct Spillover with sequential controls ----
tryCatch({
  if (nrow(selic_stability) > 0) {
    selic_stability$spec_short <- factor(
      seq_len(nrow(selic_stability)),
      labels = c("Selic only", paste0("+", available_global[1:min(nrow(selic_stability)-1, length(available_global))]))
    )

    p <- ggplot(selic_stability, aes(x = spec_short, y = selic_coef)) +
      geom_col(aes(fill = selic_p < 0.10), show.legend = FALSE) +
      geom_errorbar(aes(ymin = selic_coef - z_crit * selic_se,
                         ymax = selic_coef + z_crit * selic_se), width = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Selic Coefficient with Sequential Global Controls",
           subtitle = "Testing regional contagion vs common global exposure",
           x = "Specification", y = "Selic coefficient on FCI") +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    ggsave(file.path(REG_CONFIG$output_dir, "241_Regional_Direct_Spillover.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 241_Regional_Direct_Spillover.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 241 failed:", e$message, "\n"))

# ---- 242: Contagion vs Common Shocks ----
tryCatch({
  if (nrow(selic_stability) > 0) {
    p <- ggplot(selic_stability, aes(x = step, y = selic_coef)) +
      geom_ribbon(aes(ymin = selic_coef - z_crit * selic_se,
                       ymax = selic_coef + z_crit * selic_se), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_point(aes(color = selic_p < 0.10), size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Selic Coefficient Stability",
           subtitle = "Steps = sequential addition of global controls",
           x = "Number of Global Controls Added", y = "Selic β") +
      scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60"),
                         name = "p < 0.10") +
      theme_minimal(base_size = 11)

    ggsave(file.path(REG_CONFIG$output_dir, "242_Contagion_vs_CommonShocks.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 242_Contagion_vs_CommonShocks.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 242 failed:", e$message, "\n"))

# ---- 243: Regional Granger Causality ----
tryCatch({
  if (nrow(granger_results) > 0) {
    granger_results$label <- paste(granger_results$test, granger_results$direction, sep = ": ")

    p <- ggplot(granger_results, aes(x = factor(lag), y = p_value, fill = direction)) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 0.10, linetype = "dashed", color = "red") +
      facet_wrap(~test) +
      labs(title = "Regional Granger Causality",
           subtitle = "Selic → FCI vs FCI → Selic | Dashed = 10% significance",
           x = "Lag Order", y = "p-value") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")

    ggsave(file.path(REG_CONFIG$output_dir, "243_Regional_Granger_Causality.png"),
           p, width = 12, height = 6, dpi = 150)
    cat("Saved: 243_Regional_Granger_Causality.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 243 failed:", e$message, "\n"))

# ---- 244: Asymmetric Selic Spillover ----
tryCatch({
  if (nrow(asym_results) > 0) {
    asym_plot <- data.frame(
      Direction = c("Tightening (ΔSelic > 0)", "Easing (ΔSelic < 0)"),
      Coef = c(asym_results$beta_pos, asym_results$beta_neg),
      p_value = c(asym_results$beta_pos_p, asym_results$beta_neg_p)
    )

    p <- ggplot(asym_plot, aes(x = Direction, y = Coef, fill = p_value < 0.10)) +
      geom_col(show.legend = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Asymmetric Selic Spillovers to Paraguay FCI",
           subtitle = sprintf("Wald test (β⁺ = β⁻): p = %.3f", asym_results$wald_p),
           x = NULL, y = "Coefficient") +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
      theme_minimal(base_size = 11)

    ggsave(file.path(REG_CONFIG$output_dir, "244_Asymmetric_Selic_Spillover.png"),
           p, width = 8, height = 6, dpi = 150)
    cat("Saved: 244_Asymmetric_Selic_Spillover.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 244 failed:", e$message, "\n"))

# ---- 245: Regional Spillover LP → FCI ----
tryCatch({
  lp_fci <- spillover_lp %>% filter(target == "FCI_ENDO_AVG")
  if (nrow(lp_fci) > 0) {
    p <- ggplot(lp_fci, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "darkgreen") +
      geom_line(color = "darkgreen", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Regional Spillover LP: ΔSelic → FCI_ENDO",
           subtitle = "With global controls (VIX, FFER, DXY) | 90% CI",
           x = "Horizon (months)", y = "FCI response to 1pp ΔSelic") +
      theme_minimal(base_size = 11)

    ggsave(file.path(REG_CONFIG$output_dir, "245_Regional_Spillover_LP_FCI.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 245_Regional_Spillover_LP_FCI.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 245 failed:", e$message, "\n"))

# ---- 246: Regional Spillover LP → Credit ----
tryCatch({
  lp_credit <- spillover_lp %>% filter(target == "Cred_Real_Total")
  if (nrow(lp_credit) > 0) {
    p <- ggplot(lp_credit, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "firebrick") +
      geom_line(color = "firebrick", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Regional Spillover LP: ΔSelic → Real Credit",
           subtitle = "Indirect channel (via FCI) | With global controls | 90% CI",
           x = "Horizon (months)", y = "Credit response to 1pp ΔSelic") +
      theme_minimal(base_size = 11)

    ggsave(file.path(REG_CONFIG$output_dir, "246_Regional_Spillover_LP_Credit.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 246_Regional_Spillover_LP_Credit.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 246 failed:", e$message, "\n"))

# ---- 247: VAR IRF Comparison ----
tryCatch({
  if (nrow(var_irf_results) > 0) {
    p <- ggplot(var_irf_results, aes(x = horizon, y = irf_val,
                                      color = specification, fill = specification)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = NA) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "VAR IRF: Selic → FCI_ENDO",
           subtitle = "Bivariate vs Controlled (VIX, FFER partialed out)",
           x = "Horizon (months)", y = "FCI response") +
      scale_color_manual(values = c("Bivariate" = "steelblue", "Controlled" = "firebrick")) +
      scale_fill_manual(values = c("Bivariate" = "steelblue", "Controlled" = "firebrick")) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")

    ggsave(file.path(REG_CONFIG$output_dir, "247_Regional_VAR_IRF.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 247_Regional_VAR_IRF.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 247 failed:", e$message, "\n"))

# ---- 248: Dashboard ----
tryCatch({
  plots <- list()

  # Selic stability
  if (nrow(selic_stability) > 0) {
    plots[[1]] <- ggplot(selic_stability, aes(x = step, y = selic_coef)) +
      geom_line(color = "steelblue") +
      geom_point(aes(color = selic_p < 0.10), size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Selic β Stability", x = "# Controls", y = "β") +
      scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
      theme_minimal(base_size = 8) + theme(legend.position = "none")
  }

  # Granger
  granger_selic_fci <- granger_results %>% filter(grepl("Selic.*FCI", direction))
  if (nrow(granger_selic_fci) > 0) {
    plots[[2]] <- ggplot(granger_selic_fci, aes(x = factor(lag), y = p_value, fill = test)) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 0.10, linetype = "dashed", color = "red") +
      labs(title = "Granger: Selic→FCI", x = "Lag", y = "p-value") +
      theme_minimal(base_size = 8) + theme(legend.position = "bottom")
  }

  # LP FCI
  lp_fci <- spillover_lp %>% filter(target == "FCI_ENDO_AVG")
  if (nrow(lp_fci) > 0) {
    plots[[3]] <- ggplot(lp_fci, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "darkgreen") +
      geom_line(color = "darkgreen") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "LP: ΔSelic→FCI", x = "Horizon", y = "β") +
      theme_minimal(base_size = 8)
  }

  # VAR IRF
  if (nrow(var_irf_results) > 0) {
    plots[[4]] <- ggplot(var_irf_results, aes(x = horizon, y = irf_val, color = specification)) +
      geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "VAR: Selic→FCI", x = "Horizon", y = "IRF") +
      scale_color_manual(values = c("Bivariate" = "steelblue", "Controlled" = "firebrick")) +
      theme_minimal(base_size = 8) + theme(legend.position = "bottom")
  }

  if (length(plots) > 0) {
    combined <- do.call(grid.arrange, c(plots, ncol = 2))
    ggsave(file.path(REG_CONFIG$output_dir, "248_Regional_Spillovers_Dashboard.png"),
           combined, width = 12, height = 8, dpi = 150)
    cat("Saved: 248_Regional_Spillovers_Dashboard.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 248 failed:", e$message, "\n"))


cat("\n################################################################################\n")
cat("SCRIPT 22 COMPLETE\n")
cat("################################################################################\n\n")
cat("Outputs: 8 PNGs (241-248) + 4 CSVs\n\n")
