################################################################################
# FCI COMMODITY PRICE PUZZLE INVESTIGATION
################################################################################
#
# Project:      Financial Conditions Index - Commodity Puzzle Diagnosis
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Investigates why aggregate commodity prices don't predict FCI
#               despite Paraguay being a major commodity exporter. Tests
#               disaggregated effects, Terms of Trade, offsetting effects
#               (agricultural exports vs oil imports), lagged/nonlinear effects,
#               and commodity-credit interactions.
#
#   SECTION A: Disaggregated Commodity → FCI Regressions
#   SECTION B: Terms of Trade vs Individual Prices
#   SECTION C: Offsetting Effects Test
#   SECTION D: Lagged & Nonlinear Effects
#   SECTION E: Commodity-Credit Interaction LP
#   SECTION F: Granger Causality
#
# References:
#   - Chen, Rogoff & Rossi (2010) - Can Exchange Rates Forecast Commodity Prices?
#   - Fernandez, Schmitt-Grohe & Uribe (2017) - World Shocks, World Prices, and
#     Business Cycles
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
  library(vars)
})

COMM_CONFIG <- list(
  max_horizon     = 24,
  n_lags          = 2,
  confidence_level = 0.90,
  output_dir      = "../output",
  verbose         = TRUE
)

set.seed(20260211)

cat("\n################################################################################\n")
cat("FCI COMMODITY PRICE PUZZLE INVESTIGATION\n")
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
ext_file <- file.path(COMM_CONFIG$output_dir, "New_External_Variables.csv")
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

# Create macro lags
analysis_data <- analysis_data %>%
  mutate(
    IMAEP_yoy_L1 = lag(IMAEP_yoy, 1),
    IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
    IPC_yoy_L1 = lag(IPC_yoy, 1),
    IPC_yoy_L2 = lag(IPC_yoy, 2)
  )

commodity_vars <- c("Soybean", "Soybean_oil", "Soybean_flour", "Cotton",
                    "Corn", "Wheat", "Meat", "Sugar", "Oil_Brent")
commodity_yoy <- paste0("d_", commodity_vars)

cat("Analysis data:", nrow(analysis_data), "obs\n\n")


################################################################################
# SECTION A: DISAGGREGATED COMMODITY → FCI REGRESSIONS
################################################################################

cat("================================================================================\n")
cat("SECTION A: DISAGGREGATED COMMODITY → FCI REGRESSIONS\n")
cat("================================================================================\n\n")

control_str <- "IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 + IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2"
z_crit <- qnorm(1 - (1 - COMM_CONFIG$confidence_level) / 2)

# ---- Individual commodity regressions ----
indiv_results <- data.frame()

for (cv in commodity_yoy) {
  if (!(cv %in% names(analysis_data))) next

  f <- as.formula(paste("FCI_ENDO_AVG ~", cv, "+", control_str))
  reg_data <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, !!sym(cv),
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    na.omit()

  if (nrow(reg_data) < 30) next

  model <- lm(f, data = reg_data)
  vcov_hac <- sandwich::NeweyWest(model, lag = 3, prewhite = FALSE)
  ct <- lmtest::coeftest(model, vcov = vcov_hac)
  idx <- which(rownames(ct) == cv)

  indiv_results <- rbind(indiv_results, data.frame(
    commodity = cv,
    coef = ct[idx, 1],
    se = ct[idx, 2],
    t_stat = ct[idx, 3],
    p_value = ct[idx, 4],
    ci_lower = ct[idx, 1] - z_crit * ct[idx, 2],
    ci_upper = ct[idx, 1] + z_crit * ct[idx, 2],
    n_obs = nrow(reg_data),
    stringsAsFactors = FALSE
  ))
}

cat("Individual Commodity Regressions (FCI_ENDO ~ commodity + controls):\n")
cat(sprintf("%-20s %10s %10s %10s\n", "Commodity", "Coef", "t-stat", "p-value"))
cat(strrep("-", 55), "\n")
for (i in seq_len(nrow(indiv_results))) {
  stars <- ifelse(indiv_results$p_value[i] < 0.01, "***",
                  ifelse(indiv_results$p_value[i] < 0.05, "**",
                         ifelse(indiv_results$p_value[i] < 0.10, "*", "")))
  cat(sprintf("%-20s %+10.4f %10.2f %10.3f %s\n",
              indiv_results$commodity[i],
              indiv_results$coef[i],
              indiv_results$t_stat[i],
              indiv_results$p_value[i],
              stars))
}

# ---- Multivariate: all commodities jointly ----
available_comm_yoy <- intersect(commodity_yoy, names(analysis_data))
f_joint <- as.formula(paste("FCI_ENDO_AVG ~",
                             paste(available_comm_yoy, collapse = " + "), "+",
                             control_str))
reg_joint <- analysis_data %>%
  dplyr::select(FCI_ENDO_AVG, all_of(available_comm_yoy),
                IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
  na.omit()

if (nrow(reg_joint) >= 30) {
  model_joint <- lm(f_joint, data = reg_joint)
  vcov_joint <- sandwich::NeweyWest(model_joint, lag = 3, prewhite = FALSE)
  ct_joint <- lmtest::coeftest(model_joint, vcov = vcov_joint)

  # F-test for joint significance of all commodities
  f_restricted <- as.formula(paste("FCI_ENDO_AVG ~", control_str))
  model_restricted <- lm(f_restricted, data = reg_joint)
  joint_F <- anova(model_restricted, model_joint)
  cat(sprintf("\nJoint F-test (all commodities): F=%.2f, p=%.4f\n",
              joint_F$F[2], joint_F$`Pr(>F)`[2]))
  cat(sprintf("Joint model R²: %.4f (vs controls-only R²: %.4f)\n",
              summary(model_joint)$r.squared, summary(model_restricted)$r.squared))
}
cat("\n")


################################################################################
# SECTION A.2: SOYBEAN ORTHOGONALIZED ON DXY
################################################################################
# Paper claim: residualizing soybean on DXY yields weaker coefficient,
# confirming marginal information beyond DXY is negligible.

cat("================================================================================\n")
cat("SECTION A.2: SOYBEAN ORTHOGONALIZED ON DXY\n")
cat("================================================================================\n\n")

if ("DXY" %in% names(analysis_data) && "d_Soybean" %in% names(analysis_data)) {
  orth_data <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, d_Soybean, DXY,
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    na.omit()

  cat(sprintf("Orthogonalization sample: %d observations\n", nrow(orth_data)))

  # Step 1: Residualize soybean on DXY
  soy_on_dxy <- lm(d_Soybean ~ DXY, data = orth_data)
  orth_data$d_Soybean_orth <- residuals(soy_on_dxy)
  cat(sprintf("  Soybean ~ DXY: R² = %.4f\n", summary(soy_on_dxy)$r.squared))

  # Step 2: Regress FCI on orthogonalized soybean (with macro controls)
  f_orth <- as.formula(paste("FCI_ENDO_AVG ~ d_Soybean_orth +", control_str))
  model_orth <- lm(f_orth, data = orth_data)
  vcov_orth <- sandwich::NeweyWest(model_orth, lag = 3, prewhite = FALSE)
  ct_orth <- lmtest::coeftest(model_orth, vcov = vcov_orth)

  soy_orth_coef <- ct_orth["d_Soybean_orth", "Estimate"]
  soy_orth_se   <- ct_orth["d_Soybean_orth", "Std. Error"]
  soy_orth_t    <- ct_orth["d_Soybean_orth", "t value"]
  soy_orth_p    <- ct_orth["d_Soybean_orth", "Pr(>|t|)"]

  cat(sprintf("  Soybean (orth. on DXY) → FCI: coef = %.5f, t = %.2f, p = %.3f\n",
              soy_orth_coef, soy_orth_t, soy_orth_p))

  # Compare with direct soybean coefficient
  f_direct <- as.formula(paste("FCI_ENDO_AVG ~ d_Soybean +", control_str))
  model_direct <- lm(f_direct, data = orth_data)
  vcov_direct <- sandwich::NeweyWest(model_direct, lag = 3, prewhite = FALSE)
  ct_direct <- lmtest::coeftest(model_direct, vcov = vcov_direct)
  cat(sprintf("  Soybean (direct) → FCI:       coef = %.5f, t = %.2f, p = %.3f\n",
              ct_direct["d_Soybean", "Estimate"],
              ct_direct["d_Soybean", "t value"],
              ct_direct["d_Soybean", "Pr(>|t|)"]))

  # Append to individual results for export
  orth_row <- data.frame(
    commodity = "d_Soybean_orth_DXY",
    coef = soy_orth_coef,
    se = soy_orth_se,
    t_stat = soy_orth_t,
    p_value = soy_orth_p,
    ci_lower = soy_orth_coef - z_crit * soy_orth_se,
    ci_upper = soy_orth_coef + z_crit * soy_orth_se,
    n_obs = nrow(orth_data),
    stringsAsFactors = FALSE
  )
  indiv_results <- rbind(indiv_results, orth_row)

  cat("\n")
} else {
  cat("DXY not available in analysis_data — skipping orthogonalization.\n")
  cat("Ensure 17_FCI_New_External_Data.R has been sourced first.\n\n")
}


################################################################################
# SECTION B: TERMS OF TRADE vs INDIVIDUAL PRICES
################################################################################

cat("================================================================================\n")
cat("SECTION B: TERMS OF TRADE vs INDIVIDUAL PRICES\n")
cat("================================================================================\n\n")

model_comparison <- data.frame()

# Model 1: Aggregate Commodities
if ("Commodities" %in% names(analysis_data)) {
  f1 <- as.formula(paste("FCI_ENDO_AVG ~ Commodities +", control_str))
  reg1 <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, Commodities,
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    na.omit()
  m1 <- lm(f1, data = reg1)
  vcov1 <- sandwich::NeweyWest(m1, lag = 3, prewhite = FALSE)
  ct1 <- lmtest::coeftest(m1, vcov = vcov1)
  idx1 <- which(rownames(ct1) == "Commodities")

  model_comparison <- rbind(model_comparison, data.frame(
    Model = "M1: Aggregate Commodities",
    Coef = ct1[idx1, 1], SE = ct1[idx1, 2], p_value = ct1[idx1, 4],
    R2 = summary(m1)$r.squared, AIC = AIC(m1), n_obs = nrow(reg1),
    stringsAsFactors = FALSE
  ))
}

# Model 2: Paraguay ToT
if ("d_ToT" %in% names(analysis_data)) {
  f2 <- as.formula(paste("FCI_ENDO_AVG ~ d_ToT +", control_str))
  reg2 <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, d_ToT,
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    na.omit()
  m2 <- lm(f2, data = reg2)
  vcov2 <- sandwich::NeweyWest(m2, lag = 3, prewhite = FALSE)
  ct2 <- lmtest::coeftest(m2, vcov = vcov2)
  idx2 <- which(rownames(ct2) == "d_ToT")

  model_comparison <- rbind(model_comparison, data.frame(
    Model = "M2: Paraguay ToT",
    Coef = ct2[idx2, 1], SE = ct2[idx2, 2], p_value = ct2[idx2, 4],
    R2 = summary(m2)$r.squared, AIC = AIC(m2), n_obs = nrow(reg2),
    stringsAsFactors = FALSE
  ))
}

# Model 3: All disaggregated
if (nrow(reg_joint) >= 30) {
  model_comparison <- rbind(model_comparison, data.frame(
    Model = "M3: All Disaggregated",
    Coef = NA, SE = NA, p_value = joint_F$`Pr(>F)`[2],
    R2 = summary(model_joint)$r.squared, AIC = AIC(model_joint), n_obs = nrow(reg_joint),
    stringsAsFactors = FALSE
  ))
}

cat("Model Comparison:\n")
print(model_comparison, row.names = FALSE)
cat("\n")


################################################################################
# SECTION C: OFFSETTING EFFECTS TEST
################################################################################

cat("================================================================================\n")
cat("SECTION C: OFFSETTING EFFECTS TEST\n")
cat("================================================================================\n\n")

# Construct export commodity index (equal-weighted agricultural)
ag_commodities <- paste0("d_", c("Soybean", "Soybean_oil", "Soybean_flour",
                                  "Cotton", "Corn", "Wheat", "Meat", "Sugar"))
available_ag <- intersect(ag_commodities, names(analysis_data))

if (length(available_ag) > 0 && "d_Oil_Brent" %in% names(analysis_data)) {
  analysis_data$d_Export_index <- rowMeans(analysis_data[, available_ag, drop = FALSE], na.rm = TRUE)
  analysis_data$d_Export_index[rowSums(!is.na(analysis_data[, available_ag])) < 3] <- NA

  # Offsetting test: FCI ~ Export_index + Oil_Brent
  f_offset <- as.formula(paste("FCI_ENDO_AVG ~ d_Export_index + d_Oil_Brent +", control_str))
  reg_offset <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, d_Export_index, d_Oil_Brent,
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    na.omit()

  if (nrow(reg_offset) >= 30) {
    m_offset <- lm(f_offset, data = reg_offset)
    vcov_offset <- sandwich::NeweyWest(m_offset, lag = 3, prewhite = FALSE)
    ct_offset <- lmtest::coeftest(m_offset, vcov = vcov_offset)

    cat("Offsetting Effects Test:\n")
    cat(sprintf("  Export index: coef = %+.4f, t = %.2f, p = %.3f\n",
                ct_offset["d_Export_index", 1], ct_offset["d_Export_index", 3],
                ct_offset["d_Export_index", 4]))
    cat(sprintf("  Oil (Brent):  coef = %+.4f, t = %.2f, p = %.3f\n",
                ct_offset["d_Oil_Brent", 1], ct_offset["d_Oil_Brent", 3],
                ct_offset["d_Oil_Brent", 4]))

    # Interpretation
    if (ct_offset["d_Export_index", 1] < 0 && ct_offset["d_Oil_Brent", 1] > 0) {
      cat("  → OFFSETTING pattern confirmed: Exports loosen, Oil tightens\n")
    } else {
      cat("  → No clear offsetting pattern\n")
    }

    # Net effect test
    analysis_data$d_Net_commodity <- analysis_data$d_Export_index - analysis_data$d_Oil_Brent
    f_net <- as.formula(paste("FCI_ENDO_AVG ~ d_Net_commodity +", control_str))
    reg_net <- analysis_data %>%
      dplyr::select(FCI_ENDO_AVG, d_Net_commodity,
                    IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                    IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
      na.omit()
    m_net <- lm(f_net, data = reg_net)
    vcov_net <- sandwich::NeweyWest(m_net, lag = 3, prewhite = FALSE)
    ct_net <- lmtest::coeftest(m_net, vcov = vcov_net)
    cat(sprintf("  Net effect (Export - Oil): coef = %+.4f, t = %.2f, p = %.3f\n\n",
                ct_net["d_Net_commodity", 1], ct_net["d_Net_commodity", 3],
                ct_net["d_Net_commodity", 4]))
  }
} else {
  cat("Insufficient data for offsetting test\n\n")
}


################################################################################
# SECTION D: LAGGED & NONLINEAR EFFECTS
################################################################################

cat("================================================================================\n")
cat("SECTION D: LAGGED & NONLINEAR EFFECTS\n")
cat("================================================================================\n\n")

# ---- D1: Cumulative effects at different lag structures ----
lag_results <- data.frame()

if ("d_ToT" %in% names(analysis_data)) {
  for (L in c(0, 1, 3, 6, 12)) {
    # Create cumulative lags
    lag_terms <- character(0)
    for (k in 0:L) {
      lag_name <- if (k == 0) "d_ToT" else paste0("d_ToT_L", k)
      if (k > 0) {
        analysis_data[[lag_name]] <- lag(analysis_data$d_ToT, k)
      }
      lag_terms <- c(lag_terms, lag_name)
    }

    available_lags <- intersect(lag_terms, names(analysis_data))
    f_lag <- as.formula(paste("FCI_ENDO_AVG ~",
                               paste(available_lags, collapse = " + "), "+",
                               control_str))
    reg_lag <- analysis_data %>%
      dplyr::select(FCI_ENDO_AVG, all_of(available_lags),
                    IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                    IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
      na.omit()

    if (nrow(reg_lag) >= 30) {
      m_lag <- lm(f_lag, data = reg_lag)
      vcov_lag <- sandwich::NeweyWest(m_lag, lag = 3, prewhite = FALSE)
      ct_lag <- lmtest::coeftest(m_lag, vcov = vcov_lag)

      # Cumulative effect
      tot_coefs <- ct_lag[available_lags, 1]
      cum_effect <- sum(tot_coefs)

      # Joint F-test
      m_restricted_lag <- lm(as.formula(paste("FCI_ENDO_AVG ~", control_str)), data = reg_lag)
      f_test <- anova(m_restricted_lag, m_lag)

      lag_results <- rbind(lag_results, data.frame(
        max_lag = L,
        cumulative_effect = cum_effect,
        F_stat = f_test$F[2],
        F_p_value = f_test$`Pr(>F)`[2],
        n_obs = nrow(reg_lag),
        stringsAsFactors = FALSE
      ))
    }
  }

  cat("Cumulative ToT Effects at Different Lag Structures:\n")
  cat(sprintf("%-10s %15s %10s %10s\n", "Max Lag", "Cum. Effect", "F-stat", "p-value"))
  cat(strrep("-", 50), "\n")
  for (i in seq_len(nrow(lag_results))) {
    cat(sprintf("L=%2d      %+15.4f %10.2f %10.3f\n",
                lag_results$max_lag[i], lag_results$cumulative_effect[i],
                lag_results$F_stat[i], lag_results$F_p_value[i]))
  }
  cat("\n")
}

# ---- D2: Nonlinear effects ----
nonlinear_results <- data.frame()

if ("d_ToT" %in% names(analysis_data)) {
  # Threshold test
  tot_sd <- sd(analysis_data$d_ToT, na.rm = TRUE)
  analysis_data$d_ToT_large <- analysis_data$d_ToT * (abs(analysis_data$d_ToT) > tot_sd)
  analysis_data$d_ToT_small <- analysis_data$d_ToT * (abs(analysis_data$d_ToT) <= tot_sd)

  f_thresh <- as.formula(paste("FCI_ENDO_AVG ~ d_ToT_large + d_ToT_small +", control_str))
  reg_thresh <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, d_ToT_large, d_ToT_small,
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    na.omit()

  if (nrow(reg_thresh) >= 30) {
    m_thresh <- lm(f_thresh, data = reg_thresh)
    vcov_thresh <- sandwich::NeweyWest(m_thresh, lag = 3, prewhite = FALSE)
    ct_thresh <- lmtest::coeftest(m_thresh, vcov = vcov_thresh)

    cat("Threshold Test (|ToT change| > 1 SD):\n")
    cat(sprintf("  Large moves: coef = %+.4f, t = %.2f, p = %.3f\n",
                ct_thresh["d_ToT_large", 1], ct_thresh["d_ToT_large", 3],
                ct_thresh["d_ToT_large", 4]))
    cat(sprintf("  Small moves: coef = %+.4f, t = %.2f, p = %.3f\n",
                ct_thresh["d_ToT_small", 1], ct_thresh["d_ToT_small", 3],
                ct_thresh["d_ToT_small", 4]))
  }

  # Quadratic test
  analysis_data$d_ToT_sq <- analysis_data$d_ToT^2
  f_quad <- as.formula(paste("FCI_ENDO_AVG ~ d_ToT + d_ToT_sq +", control_str))
  reg_quad <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, d_ToT, d_ToT_sq,
                  IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                  IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
    na.omit()

  if (nrow(reg_quad) >= 30) {
    m_quad <- lm(f_quad, data = reg_quad)
    vcov_quad <- sandwich::NeweyWest(m_quad, lag = 3, prewhite = FALSE)
    ct_quad <- lmtest::coeftest(m_quad, vcov = vcov_quad)

    cat(sprintf("\nQuadratic Test:\n  Linear: coef = %+.4f, p = %.3f\n  Squared: coef = %+.6f, p = %.3f\n\n",
                ct_quad["d_ToT", 1], ct_quad["d_ToT", 4],
                ct_quad["d_ToT_sq", 1], ct_quad["d_ToT_sq", 4]))
  }
}


################################################################################
# SECTION E: COMMODITY-CREDIT INTERACTION LP
################################################################################

cat("================================================================================\n")
cat("SECTION E: COMMODITY-CREDIT INTERACTION LP\n")
cat("================================================================================\n\n")

interaction_lp <- data.frame()
marginal_cov <- data.frame()

if ("d_ToT" %in% names(analysis_data) && "FCI_exCredit_AVG" %in% names(analysis_data) &&
    "Cred_Real_Total" %in% names(analysis_data)) {

  cat("Running interaction LP: Credit ~ FCI + ToT + FCI×ToT + controls...\n")

  for (h in 1:COMM_CONFIG$max_horizon) {
    data_h <- analysis_data %>%
      mutate(
        y_fwd = lead(Cred_Real_Total, h),
        y_lag1 = lag(Cred_Real_Total, 1),
        y_lag2 = lag(Cred_Real_Total, 2),
        fci_lag1 = lag(FCI_exCredit_AVG, 1),
        FCI_x_ToT = FCI_exCredit_AVG * d_ToT
      )

    reg_data <- data_h %>%
      dplyr::select(y_fwd, FCI_exCredit_AVG, d_ToT, FCI_x_ToT,
                    y_lag1, y_lag2, fci_lag1,
                    IMAEP_yoy, IMAEP_yoy_L1, IMAEP_yoy_L2,
                    IPC_yoy, IPC_yoy_L1, IPC_yoy_L2) %>%
      na.omit()

    if (nrow(reg_data) < 50) next

    f <- y_fwd ~ FCI_exCredit_AVG + d_ToT + FCI_x_ToT +
      y_lag1 + y_lag2 + fci_lag1 +
      IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +
      IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2

    model <- lm(f, data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    ct <- lmtest::coeftest(model, vcov = vcov_hac)

    # Store covariance info for marginal effects
    if ("FCI_exCredit_AVG" %in% names(coef(model)) && "FCI_x_ToT" %in% names(coef(model))) {
      marginal_cov <- rbind(marginal_cov, data.frame(
        horizon = h,
        beta1 = coef(model)["FCI_exCredit_AVG"],
        beta3 = coef(model)["FCI_x_ToT"],
        var_beta1 = vcov_hac["FCI_exCredit_AVG", "FCI_exCredit_AVG"],
        var_beta3 = vcov_hac["FCI_x_ToT", "FCI_x_ToT"],
        cov_beta1_beta3 = vcov_hac["FCI_exCredit_AVG", "FCI_x_ToT"],
        stringsAsFactors = FALSE
      ))
    }

    for (var_name in c("FCI_exCredit_AVG", "d_ToT", "FCI_x_ToT")) {
      idx <- which(rownames(ct) == var_name)
      if (length(idx) > 0) {
        interaction_lp <- rbind(interaction_lp, data.frame(
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

  cat("Interaction LP Results at key horizons:\n")
  cat(sprintf("%-20s %8s %10s %10s\n", "Variable", "Horizon", "Coef", "p-value"))
  cat(strrep("-", 52), "\n")
  for (kh in c(6, 12, 18, 24)) {
    for (vn in c("FCI_exCredit_AVG", "d_ToT", "FCI_x_ToT")) {
      row <- interaction_lp %>% filter(horizon == kh & variable == vn)
      if (nrow(row) > 0) {
        cat(sprintf("%-20s %8d %+10.3f %10.3f\n", vn, kh, row$coef[1], row$p_value[1]))
      }
    }
    cat("\n")
  }
}


################################################################################
# SECTION E2: MARGINAL EFFECTS OF FCI BY ToT REGIME
################################################################################

cat("================================================================================\n")
cat("SECTION E2: MARGINAL EFFECTS OF FCI BY ToT REGIME\n")
cat("================================================================================\n\n")

if (nrow(marginal_cov) > 0 && "d_ToT" %in% names(analysis_data)) {

  # ToT standard deviation for grid
  tot_sd <- sd(analysis_data$d_ToT, na.rm = TRUE)
  tot_grid <- seq(-2 * tot_sd, 2 * tot_sd, length.out = 41)
  tot_grid_sd <- tot_grid / tot_sd

  # Compute marginal effects for h = 6, 12, 18
  marginal_effects <- data.frame()

  for (hh in c(6, 12, 18)) {
    row_h <- marginal_cov[marginal_cov$horizon == hh, ]
    if (nrow(row_h) == 0) next

    b1 <- row_h$beta1
    b3 <- row_h$beta3
    v1 <- row_h$var_beta1
    v3 <- row_h$var_beta3
    c13 <- row_h$cov_beta1_beta3

    for (i in seq_along(tot_grid)) {
      tot_val <- tot_grid[i]
      me <- b1 + b3 * tot_val
      se_me <- sqrt(v1 + tot_val^2 * v3 + 2 * tot_val * c13)

      marginal_effects <- rbind(marginal_effects, data.frame(
        horizon = hh,
        tot_value = tot_val,
        tot_sd_units = tot_grid_sd[i],
        marginal_effect = me,
        se = se_me,
        ci_lower = me - z_crit * se_me,
        ci_upper = me + z_crit * se_me,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Print key values
  cat("Marginal Effects at Key ToT Values (h=12):\n")
  me_h12 <- marginal_effects[marginal_effects$horizon == 12, ]
  for (sd_val in c(-2, -1, 0, 1, 2)) {
    closest <- me_h12[which.min(abs(me_h12$tot_sd_units - sd_val)), ]
    cat(sprintf("  ToT = %+d SD: ME = %+.2f pp [%.2f, %.2f]\n",
                sd_val, closest$marginal_effect, closest$ci_lower, closest$ci_upper))
  }

  # Save CSV
  write.csv(marginal_effects, file.path(COMM_CONFIG$output_dir, "Marginal_Effects_FCI_ToT.csv"),
            row.names = FALSE)
  cat("\nSaved: Marginal_Effects_FCI_ToT.csv\n")

  # ---- 234: Marginal Effects Plot ----
  tryCatch({
    marginal_effects$horizon_label <- factor(
      paste0("h = ", marginal_effects$horizon),
      levels = c("h = 6", "h = 12", "h = 18")
    )

    p234 <- ggplot(marginal_effects, aes(x = tot_sd_units, y = marginal_effect)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
      geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray70") +
      facet_wrap(~horizon_label, ncol = 3) +
      scale_x_continuous(
        name = "Terms-of-Trade Growth (SD units)",
        sec.axis = sec_axis(~ . * tot_sd, name = "ToT Growth (pp)")
      ) +
      labs(
        title = "Marginal Effect of FCI Tightening on Credit by ToT Regime",
        subtitle = expression(paste("ME(FCI) = ", beta[1], " + ", beta[3],
                                     " \u00D7 ToT  |  Delta-method 90% CI from HAC vcov")),
        y = "Marginal FCI Effect on Credit (pp)"
      ) +
      theme_minimal(base_size = 11) +
      theme(strip.text = element_text(face = "bold"))

    ggsave(file.path(COMM_CONFIG$output_dir, "234_Marginal_Effects_FCI_Credit.png"),
           p234, width = 14, height = 5, dpi = 300)
    cat("Saved: 234_Marginal_Effects_FCI_Credit.png\n")
  }, error = function(e) cat("  WARNING: Plot 234 failed:", e$message, "\n"))

} else {
  cat("Marginal effects skipped: insufficient covariance data\n\n")
}


################################################################################
# SECTION F: GRANGER CAUSALITY
################################################################################

cat("================================================================================\n")
cat("SECTION F: GRANGER CAUSALITY\n")
cat("================================================================================\n\n")

granger_results <- data.frame()

# Bivariate Granger: each commodity → FCI_ENDO
for (cv in c(commodity_yoy, "d_ToT")) {
  if (!(cv %in% names(analysis_data))) next

  granger_data <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, !!sym(cv)) %>%
    na.omit()

  if (nrow(granger_data) < 50) next

  for (n_lag in c(3, 6, 12)) {
    tryCatch({
      var_granger <- VAR(granger_data, p = n_lag, type = "const")
      g_test <- causality(var_granger, cause = cv)

      granger_results <- rbind(granger_results, data.frame(
        cause = cv,
        effect = "FCI_ENDO_AVG",
        lag = n_lag,
        F_stat = g_test$Granger$statistic,
        p_value = g_test$Granger$p.value,
        n_obs = nrow(granger_data) - n_lag,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
}

# Also test aggregate Commodities
if ("Commodities" %in% names(analysis_data)) {
  granger_data_agg <- analysis_data %>%
    dplyr::select(FCI_ENDO_AVG, Commodities) %>%
    na.omit()

  for (n_lag in c(3, 6, 12)) {
    tryCatch({
      var_g <- VAR(granger_data_agg, p = n_lag, type = "const")
      g_test <- causality(var_g, cause = "Commodities")
      granger_results <- rbind(granger_results, data.frame(
        cause = "Commodities_AGG",
        effect = "FCI_ENDO_AVG",
        lag = n_lag,
        F_stat = g_test$Granger$statistic,
        p_value = g_test$Granger$p.value,
        n_obs = nrow(granger_data_agg) - n_lag,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
}

cat("Granger Causality: Commodity → FCI_ENDO\n")
cat(sprintf("%-20s %5s %10s %10s\n", "Cause", "Lag", "F-stat", "p-value"))
cat(strrep("-", 50), "\n")
for (i in seq_len(nrow(granger_results))) {
  stars <- ifelse(granger_results$p_value[i] < 0.01, "***",
                  ifelse(granger_results$p_value[i] < 0.05, "**",
                         ifelse(granger_results$p_value[i] < 0.10, "*", "")))
  cat(sprintf("%-20s %5d %10.2f %10.3f %s\n",
              granger_results$cause[i], granger_results$lag[i],
              granger_results$F_stat[i], granger_results$p_value[i], stars))
}
cat("\n")


################################################################################
# SAVE OUTPUTS
################################################################################

cat("================================================================================\n")
cat("SAVING OUTPUTS\n")
cat("================================================================================\n\n")

write.csv(indiv_results, file.path(COMM_CONFIG$output_dir, "Commodity_Disaggregated_Regressions.csv"),
          row.names = FALSE)
cat("Saved: Commodity_Disaggregated_Regressions.csv\n")

write.csv(model_comparison, file.path(COMM_CONFIG$output_dir, "Commodity_ToT_Model_Comparison.csv"),
          row.names = FALSE)
cat("Saved: Commodity_ToT_Model_Comparison.csv\n")

write.csv(interaction_lp, file.path(COMM_CONFIG$output_dir, "Commodity_Credit_Interaction_LP.csv"),
          row.names = FALSE)
cat("Saved: Commodity_Credit_Interaction_LP.csv\n")

write.csv(granger_results, file.path(COMM_CONFIG$output_dir, "Commodity_Granger_Results.csv"),
          row.names = FALSE)
cat("Saved: Commodity_Granger_Results.csv\n")

# Summary
puzzle_summary <- data.frame(
  Finding = c(
    "Aggregate commodity → FCI",
    "Paraguay ToT → FCI",
    "Disaggregated joint F-test",
    "Export index sign",
    "Oil_Brent sign",
    "Interaction FCI×ToT at h=12"
  ),
  stringsAsFactors = FALSE
)
puzzle_summary$Result <- NA
puzzle_summary$p_value <- NA

# Fill in from results
if (nrow(model_comparison) > 0) {
  m1_row <- model_comparison %>% filter(grepl("Aggregate", Model))
  if (nrow(m1_row) > 0) {
    puzzle_summary$Result[1] <- sprintf("coef=%+.4f", m1_row$Coef[1])
    puzzle_summary$p_value[1] <- m1_row$p_value[1]
  }
  m2_row <- model_comparison %>% filter(grepl("ToT", Model))
  if (nrow(m2_row) > 0) {
    puzzle_summary$Result[2] <- sprintf("coef=%+.4f", m2_row$Coef[1])
    puzzle_summary$p_value[2] <- m2_row$p_value[1]
  }
  m3_row <- model_comparison %>% filter(grepl("Disagg", Model))
  if (nrow(m3_row) > 0) {
    puzzle_summary$Result[3] <- sprintf("R²=%.4f", m3_row$R2[1])
    puzzle_summary$p_value[3] <- m3_row$p_value[1]
  }
}

if (exists("ct_offset")) {
  puzzle_summary$Result[4] <- sprintf("coef=%+.4f", ct_offset["d_Export_index", 1])
  puzzle_summary$p_value[4] <- ct_offset["d_Export_index", 4]
  puzzle_summary$Result[5] <- sprintf("coef=%+.4f", ct_offset["d_Oil_Brent", 1])
  puzzle_summary$p_value[5] <- ct_offset["d_Oil_Brent", 4]
}

int_h12 <- interaction_lp %>% filter(horizon == 12 & variable == "FCI_x_ToT")
if (nrow(int_h12) > 0) {
  puzzle_summary$Result[6] <- sprintf("coef=%+.4f", int_h12$coef[1])
  puzzle_summary$p_value[6] <- int_h12$p_value[1]
}

write.csv(puzzle_summary, file.path(COMM_CONFIG$output_dir, "Commodity_Puzzle_Summary.csv"),
          row.names = FALSE)
cat("Saved: Commodity_Puzzle_Summary.csv\n")


################################################################################
# VISUALIZATIONS
################################################################################

cat("\n--- Generating Visualizations ---\n\n")

# ---- 226: Disaggregated Commodity Coefficients ----
tryCatch({
  if (nrow(indiv_results) > 0) {
    p <- ggplot(indiv_results, aes(x = reorder(commodity, coef), y = coef,
                                    fill = p_value < 0.10)) +
      geom_col(show.legend = FALSE) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3) +
      coord_flip() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Disaggregated Commodity → FCI Coefficients",
           subtitle = "Individual regressions with macro controls | 90% CI",
           x = NULL, y = "Coefficient") +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
      theme_minimal(base_size = 11)

    ggsave(file.path(COMM_CONFIG$output_dir, "226_Disaggregated_Commodity_Coefficients.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 226_Disaggregated_Commodity_Coefficients.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 226 failed:", e$message, "\n"))

# ---- 227: ToT vs Aggregate Comparison ----
tryCatch({
  if (nrow(model_comparison) > 0) {
    p <- ggplot(model_comparison %>% filter(!is.na(Coef)),
                aes(x = Model, y = Coef, fill = p_value < 0.10)) +
      geom_col(show.legend = FALSE) +
      geom_errorbar(aes(ymin = Coef - z_crit * SE, ymax = Coef + z_crit * SE), width = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Commodity Specifications: Coefficient Comparison",
           subtitle = "Aggregate vs ToT vs Disaggregated",
           x = NULL, y = "Coefficient (on FCI)") +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 15, hjust = 1))

    ggsave(file.path(COMM_CONFIG$output_dir, "227_ToT_vs_Aggregate_Comparison.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 227_ToT_vs_Aggregate_Comparison.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 227 failed:", e$message, "\n"))

# ---- 228: Offsetting Effects ----
tryCatch({
  if (exists("ct_offset")) {
    offset_df <- data.frame(
      Variable = c("Export Index\n(Agricultural)", "Oil (Brent)", "Net Effect\n(Export - Oil)"),
      Coef = c(ct_offset["d_Export_index", 1], ct_offset["d_Oil_Brent", 1],
               ct_net["d_Net_commodity", 1]),
      SE = c(ct_offset["d_Export_index", 2], ct_offset["d_Oil_Brent", 2],
             ct_net["d_Net_commodity", 2]),
      p_value = c(ct_offset["d_Export_index", 4], ct_offset["d_Oil_Brent", 4],
                  ct_net["d_Net_commodity", 4])
    )

    p <- ggplot(offset_df, aes(x = Variable, y = Coef, fill = p_value < 0.10)) +
      geom_col(show.legend = FALSE) +
      geom_errorbar(aes(ymin = Coef - z_crit * SE, ymax = Coef + z_crit * SE), width = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Offsetting Commodity Effects on FCI",
           subtitle = "Agricultural exports vs oil imports",
           x = NULL, y = "Coefficient") +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
      theme_minimal(base_size = 11)

    ggsave(file.path(COMM_CONFIG$output_dir, "228_Offsetting_Effects.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 228_Offsetting_Effects.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 228 failed:", e$message, "\n"))

# ---- 229: Lagged Commodity Effects ----
tryCatch({
  if (nrow(lag_results) > 0) {
    p <- ggplot(lag_results, aes(x = max_lag, y = cumulative_effect)) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_point(aes(color = F_p_value < 0.10), size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Cumulative ToT → FCI Effect at Different Lag Structures",
           subtitle = "Blue = jointly significant at 10%",
           x = "Maximum Lag (months)", y = "Cumulative Coefficient") +
      scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60"),
                         name = "Significant") +
      theme_minimal(base_size = 11)

    ggsave(file.path(COMM_CONFIG$output_dir, "229_Lagged_Commodity_Effects.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 229_Lagged_Commodity_Effects.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 229 failed:", e$message, "\n"))

# ---- 230: Nonlinear Effects ----
tryCatch({
  if (exists("ct_thresh") && exists("ct_quad")) {
    nonlin_df <- data.frame(
      Test = c("Large moves (>1SD)", "Small moves (<=1SD)", "Linear (quadratic)", "Squared term"),
      Coef = c(ct_thresh["d_ToT_large", 1], ct_thresh["d_ToT_small", 1],
               ct_quad["d_ToT", 1], ct_quad["d_ToT_sq", 1]),
      p_value = c(ct_thresh["d_ToT_large", 4], ct_thresh["d_ToT_small", 4],
                  ct_quad["d_ToT", 4], ct_quad["d_ToT_sq", 4])
    )

    p <- ggplot(nonlin_df, aes(x = Test, y = Coef, fill = p_value < 0.10)) +
      geom_col(show.legend = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Nonlinear Commodity Effects on FCI",
           subtitle = "Threshold and quadratic specifications",
           x = NULL, y = "Coefficient") +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 15, hjust = 1))

    ggsave(file.path(COMM_CONFIG$output_dir, "230_Nonlinear_Commodity_Effects.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 230_Nonlinear_Commodity_Effects.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 230 failed:", e$message, "\n"))

# ---- 231: Commodity-Credit Interaction LP ----
tryCatch({
  if (nrow(interaction_lp) > 0) {
    int_plot <- interaction_lp %>%
      mutate(variable_label = factor(variable,
                                      levels = c("FCI_exCredit_AVG", "d_ToT", "FCI_x_ToT"),
                                      labels = c("FCI (direct)", "ToT (direct)", "FCI × ToT (interaction)")))

    p <- ggplot(int_plot, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~variable_label, scales = "free_y") +
      labs(title = "Commodity-Credit Interaction LP",
           subtitle = "Credit ~ FCI + ToT + FCI×ToT + controls | 90% CI",
           x = "Horizon (months)", y = "Coefficient") +
      theme_minimal(base_size = 10)

    ggsave(file.path(COMM_CONFIG$output_dir, "231_Commodity_Credit_Interaction_LP.png"),
           p, width = 14, height = 5, dpi = 150)
    cat("Saved: 231_Commodity_Credit_Interaction_LP.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 231 failed:", e$message, "\n"))

# ---- 232: Granger Summary Heatmap ----
tryCatch({
  if (nrow(granger_results) > 0) {
    granger_results$sig_label <- sprintf("%.2f", granger_results$p_value)

    p <- ggplot(granger_results, aes(x = factor(lag), y = cause, fill = -log10(p_value + 0.001))) +
      geom_tile(color = "white") +
      geom_text(aes(label = sig_label), size = 2.5) +
      scale_fill_gradient(low = "white", high = "steelblue", name = "-log10(p)") +
      labs(title = "Granger Causality: Commodity → FCI_ENDO",
           subtitle = "Cell values = p-values",
           x = "Lag Order", y = NULL) +
      theme_minimal(base_size = 10)

    ggsave(file.path(COMM_CONFIG$output_dir, "232_Commodity_Granger_Summary.png"),
           p, width = 10, height = 8, dpi = 150)
    cat("Saved: 232_Commodity_Granger_Summary.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 232 failed:", e$message, "\n"))

# ---- 233: Dashboard ----
tryCatch({
  plots <- list()

  # Individual coefficients
  if (nrow(indiv_results) > 0) {
    plots[[1]] <- ggplot(indiv_results, aes(x = reorder(commodity, coef), y = coef,
                                             fill = p_value < 0.10)) +
      geom_col(show.legend = FALSE) + coord_flip() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Individual β", x = NULL, y = NULL) +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
      theme_minimal(base_size = 8)
  }

  # Interaction FCI×ToT
  fci_x_tot <- interaction_lp %>% filter(variable == "FCI_x_ToT")
  if (nrow(fci_x_tot) > 0) {
    plots[[2]] <- ggplot(fci_x_tot, aes(x = horizon, y = coef)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "FCI×ToT Interaction", x = "Horizon", y = "Coef") +
      theme_minimal(base_size = 8)
  }

  # Lagged effects
  if (nrow(lag_results) > 0) {
    plots[[3]] <- ggplot(lag_results, aes(x = max_lag, y = cumulative_effect)) +
      geom_line(color = "darkgreen") + geom_point(color = "darkgreen") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Cumulative ToT Effect", x = "Max Lag", y = "Cum. β") +
      theme_minimal(base_size = 8)
  }

  # Granger p-values for ToT
  granger_tot <- granger_results %>% filter(cause == "d_ToT")
  if (nrow(granger_tot) > 0) {
    plots[[4]] <- ggplot(granger_tot, aes(x = factor(lag), y = p_value)) +
      geom_col(fill = "steelblue") +
      geom_hline(yintercept = 0.10, linetype = "dashed", color = "red") +
      labs(title = "ToT→FCI Granger", x = "Lag", y = "p-value") +
      theme_minimal(base_size = 8)
  }

  if (length(plots) > 0) {
    combined <- do.call(grid.arrange, c(plots, ncol = 2))
    ggsave(file.path(COMM_CONFIG$output_dir, "233_Commodity_Puzzle_Dashboard.png"),
           combined, width = 12, height = 8, dpi = 150)
    cat("Saved: 233_Commodity_Puzzle_Dashboard.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 233 failed:", e$message, "\n"))


cat("\n################################################################################\n")
cat("SCRIPT 21 COMPLETE\n")
cat("################################################################################\n\n")
cat("Outputs: 9 PNGs (226-234) + 6 CSVs\n\n")
