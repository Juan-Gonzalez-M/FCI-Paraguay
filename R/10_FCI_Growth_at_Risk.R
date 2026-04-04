################################################################################
# FCI GROWTH-AT-RISK ANALYSIS (Phase 3.1)
################################################################################
#
# Project:      Financial Conditions Index - Growth-at-Risk
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Implements Growth-at-Risk (GaR) methodology following
#               Adrian, Boyarchenko, and Giannone (2019).
#               Tests whether FCI affects the left tail of growth distribution
#               more than the median (non-linear/asymmetric effects).
#
# Key Hypothesis:
#   Financial stress disproportionately affects downside risks to growth.
#   |beta_0.05| > |beta_0.50| would indicate FCI matters more for tail risk.
#
# METHODOLOGY:
#   1. Quantile regressions at tau = {0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95}
#   2. Horizons: 1, 3, 6, 12, 24 months ahead
#   3. Test significance of coefficient differences across quantiles
#   4. Construct predictive density "fan charts"
#
# Dependencies: Requires 01_FCI_Complete.R to be run first
#
# Output:
#   - Growth_at_Risk_Results.csv
#   - 100_Growth_at_Risk_Fan.png
#   - 101_GaR_Coefficient_Comparison.png
#
# References:
#   - Adrian, Boyarchenko & Giannone (2019) - Vulnerable Growth
#   - Adrian et al. (2022) - IMF Financial Stability Report methodology
#
# Author:       Departamento de Analisis Macroeconomico
# Created:      2025
#
################################################################################


################################################################################
# 1. SETUP AND CONFIGURATION
################################################################################

# Install quantreg if not available
if (!requireNamespace("quantreg", quietly = TRUE)) {
  cat("Installing quantreg package...\n")
  install.packages("quantreg", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(quantreg)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
})

GAR_CONFIG <- list(
  data_file = "../data/FCI_data_1.xlsx",
  macro_sheet = "Datos_macro",
  output_dir = "../output",

  # Quantiles to estimate
  quantiles = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95),

  # Forecast horizons (months)
  horizons = c(1, 3, 6, 12, 24),

  # Control variables lags
  n_lags = 2,

  # Minimum observations
  min_obs = 50,

  # Confidence level for coefficient comparison
  confidence_level = 0.90,

  verbose = TRUE
)

set.seed(20250126)

cat("\n################################################################################\n")
cat("FCI GROWTH-AT-RISK ANALYSIS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# 2. DATA PREPARATION
################################################################################

cat("================================================================================\n")
cat("DATA PREPARATION\n")
cat("================================================================================\n\n")

# Load FCI results
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  source("01_FCI_Complete.R")
}

# Load macro data
macro_data <- read_excel(GAR_CONFIG$data_file, sheet = GAR_CONFIG$macro_sheet)
fecha_col <- names(macro_data)[grepl("fecha|date", names(macro_data), ignore.case = TRUE)][1]

macro_data <- macro_data %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Calculate growth rates
calc_yoy <- function(x) (x / dplyr::lag(x, 12) - 1) * 100

for (var in c("IMAEP", "IPC")) {
  if (var %in% names(macro_data)) {
    macro_data <- macro_data %>%
      mutate(!!paste0(var, "_yoy") := calc_yoy(.data[[var]]))
  }
}

# Merge with FCI
gar_data <- resultado_fci$all_indices %>%
  dplyr::select(fecha, FCI = FCI_COMP_AVG, FCI_ENDO = FCI_ENDO_AVG) %>%
  inner_join(macro_data %>% dplyr::select(fecha, IMAEP_yoy, IPC_yoy), by = "fecha") %>%
  arrange(fecha) %>%
  na.omit()

cat("GaR data prepared:", nrow(gar_data), "observations\n")
cat("Period:", format(min(gar_data$fecha)), "to", format(max(gar_data$fecha)), "\n\n")


################################################################################
# 3. QUANTILE REGRESSION FUNCTIONS
################################################################################

#' Run quantile regression for a single horizon and quantile
#' @param data Data frame with FCI and target variable
#' @param target_var Name of target variable (e.g., "IMAEP_yoy")
#' @param fci_var Name of FCI variable
#' @param horizon Forecast horizon in months
#' @param tau Quantile (0 to 1)
#' @param n_lags Number of lags to include
#' @return List with coefficient, se, t-stat
run_quantile_regression <- function(data, target_var, fci_var, horizon, tau, n_lags = 2) {

  # Create forward-looking target
  data_h <- data %>%
    mutate(y_fwd = lead(!!sym(target_var), horizon))

  # Add lags of target
  for (i in 1:n_lags) {
    data_h <- data_h %>%
      mutate(!!paste0("y_lag", i) := lag(!!sym(target_var), i))
  }

  # Add lag of FCI
  data_h <- data_h %>%
    mutate(fci_lag = lag(!!sym(fci_var), 1))

  # Construct formula
  lag_vars <- paste0("y_lag", 1:n_lags)
  formula_str <- paste("y_fwd ~", fci_var, "+", paste(lag_vars, collapse = " + "), "+ fci_lag")

  # Subset to complete cases
  reg_data <- data_h %>%
    dplyr::select(y_fwd, !!sym(fci_var), all_of(lag_vars), fci_lag) %>%
    na.omit()

  if (nrow(reg_data) < GAR_CONFIG$min_obs) {
    return(list(coef = NA, se = NA, t_stat = NA, n_obs = nrow(reg_data)))
  }

  # Run quantile regression
  tryCatch({
    qr_model <- rq(as.formula(formula_str), tau = tau, data = reg_data)
    qr_summary <- summary(qr_model, se = "boot", R = 500)

    # Extract FCI coefficient
    coef_idx <- which(rownames(qr_summary$coefficients) == fci_var)

    if (length(coef_idx) > 0) {
      return(list(
        coef = qr_summary$coefficients[coef_idx, "Value"],
        se = qr_summary$coefficients[coef_idx, "Std. Error"],
        t_stat = qr_summary$coefficients[coef_idx, "t value"],
        n_obs = nrow(reg_data)
      ))
    } else {
      return(list(coef = NA, se = NA, t_stat = NA, n_obs = nrow(reg_data)))
    }
  }, error = function(e) {
    return(list(coef = NA, se = NA, t_stat = NA, n_obs = 0))
  })
}


################################################################################
# 4. RUN GROWTH-AT-RISK ANALYSIS
################################################################################

cat("================================================================================\n")
cat("RUNNING GROWTH-AT-RISK QUANTILE REGRESSIONS\n")
cat("================================================================================\n\n")

# Main target variable
target_var <- "IMAEP_yoy"
fci_var <- "FCI"

gar_results <- data.frame()

for (h in GAR_CONFIG$horizons) {
  cat(sprintf("Processing horizon h=%d months...\n", h))

  for (tau in GAR_CONFIG$quantiles) {
    result <- run_quantile_regression(gar_data, target_var, fci_var, h, tau, GAR_CONFIG$n_lags)

    gar_results <- rbind(gar_results, data.frame(
      horizon = h,
      quantile = tau,
      coef = result$coef,
      se = result$se,
      t_stat = result$t_stat,
      n_obs = result$n_obs,
      target = target_var,
      fci = fci_var
    ))
  }
}

# Add significance indicators
gar_results <- gar_results %>%
  mutate(
    p_value = 2 * (1 - pt(abs(t_stat), df = n_obs - 5)),  # Approximate
    significant_10 = abs(t_stat) > 1.65,
    significant_05 = abs(t_stat) > 1.96,
    significant_01 = abs(t_stat) > 2.58
  )

cat("\nQuantile regression results summary:\n")
cat(sprintf("%-10s %-10s %12s %12s %12s\n", "Horizon", "Quantile", "Coefficient", "Std.Error", "t-stat"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (h in GAR_CONFIG$horizons) {
  for (tau in c(0.05, 0.50, 0.95)) {
    r <- gar_results %>% filter(horizon == h, quantile == tau)
    if (nrow(r) > 0 && !is.na(r$coef)) {
      sig_stars <- ifelse(r$significant_01, "***",
                          ifelse(r$significant_05, "**",
                                 ifelse(r$significant_10, "*", "")))
      cat(sprintf("%-10d %-10.2f %+12.3f %12.3f %12.2f%s\n",
                  h, tau, r$coef, r$se, r$t_stat, sig_stars))
    }
  }
  cat("\n")
}


################################################################################
# 5. TEST FOR ASYMMETRIC EFFECTS (KEY HYPOTHESIS)
################################################################################

cat("================================================================================\n")
cat("TESTING FOR ASYMMETRIC EFFECTS (KEY HYPOTHESIS)\n")
cat("================================================================================\n\n")

# For each horizon, test if |beta_0.05| > |beta_0.50|
asymmetry_tests <- data.frame()

for (h in GAR_CONFIG$horizons) {
  coef_05 <- gar_results %>% filter(horizon == h, quantile == 0.05) %>% pull(coef)
  se_05 <- gar_results %>% filter(horizon == h, quantile == 0.05) %>% pull(se)
  coef_50 <- gar_results %>% filter(horizon == h, quantile == 0.50) %>% pull(coef)
  se_50 <- gar_results %>% filter(horizon == h, quantile == 0.50) %>% pull(se)
  coef_95 <- gar_results %>% filter(horizon == h, quantile == 0.95) %>% pull(coef)
  se_95 <- gar_results %>% filter(horizon == h, quantile == 0.95) %>% pull(se)

  if (!any(is.na(c(coef_05, coef_50, coef_95)))) {
    # Difference in absolute effects
    diff_lower <- abs(coef_05) - abs(coef_50)
    diff_upper <- abs(coef_95) - abs(coef_50)

    # Approximate test (assuming independence, which is conservative)
    se_diff_lower <- sqrt(se_05^2 + se_50^2)
    se_diff_upper <- sqrt(se_95^2 + se_50^2)

    z_lower <- diff_lower / se_diff_lower
    z_upper <- diff_upper / se_diff_upper

    asymmetry_tests <- rbind(asymmetry_tests, data.frame(
      horizon = h,
      coef_05 = coef_05,
      coef_50 = coef_50,
      coef_95 = coef_95,
      abs_diff_lower = diff_lower,
      abs_diff_upper = diff_upper,
      z_lower = z_lower,
      z_upper = z_upper,
      p_lower = 2 * (1 - pnorm(abs(z_lower))),
      p_upper = 2 * (1 - pnorm(abs(z_upper))),
      left_tail_stronger = abs(coef_05) > abs(coef_50),
      right_tail_stronger = abs(coef_95) > abs(coef_50)
    ))
  }
}

cat("Asymmetry Test Results:\n")
cat("H0: |beta_tail| = |beta_median| (symmetric effects)\n\n")

cat(sprintf("%-10s %12s %12s %12s %15s %15s\n",
            "Horizon", "Coef(0.05)", "Coef(0.50)", "Coef(0.95)", "Left > Median?", "Right > Median?"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in 1:nrow(asymmetry_tests)) {
  r <- asymmetry_tests[i, ]
  left_sig <- ifelse(r$p_lower < 0.10, paste0("YES (p=", sprintf("%.3f", r$p_lower), ")"), "No")
  right_sig <- ifelse(r$p_upper < 0.10, paste0("YES (p=", sprintf("%.3f", r$p_upper), ")"), "No")
  cat(sprintf("%-10d %+12.3f %+12.3f %+12.3f %15s %15s\n",
              r$horizon, r$coef_05, r$coef_50, r$coef_95, left_sig, right_sig))
}

# Key finding
any_asymmetry <- any(asymmetry_tests$left_tail_stronger & asymmetry_tests$p_lower < 0.10)

cat("\n================================================================================\n")
cat("KEY FINDING:\n")
cat("================================================================================\n\n")

if (any_asymmetry) {
  cat("FCI has ASYMMETRIC effects on growth distribution.\n")
  cat("The left tail (downside risk) is more sensitive to financial conditions\n")
  cat("than the median, supporting the Growth-at-Risk hypothesis.\n")
} else {
  cat("No significant evidence of asymmetric effects at conventional levels.\n")
  cat("FCI affects the growth distribution approximately symmetrically.\n")
  cat("This may reflect:\n")
  cat("  1. Paraguay's economic structure (commodity exporter)\n")
  cat("  2. Sample size limitations\n")
  cat("  3. Monetary policy success in smoothing financial shocks\n")
}


################################################################################
# 6. VISUALIZATIONS
################################################################################

cat("\n================================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Plot 1: Growth-at-Risk Fan Chart (Coefficient across quantiles by horizon)
gar_plot_data <- gar_results %>%
  filter(!is.na(coef)) %>%
  mutate(horizon_label = paste0("h=", horizon, "m"))

p1 <- ggplot(gar_plot_data, aes(x = quantile, y = coef, color = factor(horizon))) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = coef - 1.65 * se, ymax = coef + 1.65 * se, fill = factor(horizon)),
              alpha = 0.1, color = NA) +
  scale_color_brewer(palette = "Set1", name = "Horizon") +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  scale_x_continuous(breaks = GAR_CONFIG$quantiles,
                     labels = scales::percent(GAR_CONFIG$quantiles)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Growth-at-Risk: FCI Effect Across Growth Distribution",
    subtitle = "Quantile regression coefficients | Shaded = 90% CI",
    x = "Quantile of Growth Distribution",
    y = "FCI Coefficient on Future Growth"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

ggsave(file.path(GAR_CONFIG$output_dir, "100_Growth_at_Risk_Fan.png"), p1,
       width = 12, height = 8, dpi = 300)
cat("Saved: 100_Growth_at_Risk_Fan.png\n")

# Plot 2: Coefficient Comparison (5th vs 50th vs 95th)
comparison_data <- gar_results %>%
  filter(quantile %in% c(0.05, 0.50, 0.95), !is.na(coef)) %>%
  mutate(
    quantile_label = case_when(
      quantile == 0.05 ~ "5th (Left Tail)",
      quantile == 0.50 ~ "50th (Median)",
      quantile == 0.95 ~ "95th (Right Tail)"
    )
  )

p2 <- ggplot(comparison_data, aes(x = factor(horizon), y = coef, fill = quantile_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = coef - 1.65 * se, ymax = coef + 1.65 * se),
                position = position_dodge(width = 0.8), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("5th (Left Tail)" = "#E74C3C",
                               "50th (Median)" = "#3498DB",
                               "95th (Right Tail)" = "#27AE60"),
                    name = "Quantile") +
  theme_minimal(base_size = 12) +
  labs(
    title = "FCI Effects: Left Tail vs Median vs Right Tail",
    subtitle = "Comparing 5th, 50th, and 95th percentile effects | Error bars = 90% CI",
    x = "Forecast Horizon (months)",
    y = "FCI Coefficient"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

ggsave(file.path(GAR_CONFIG$output_dir, "101_GaR_Coefficient_Comparison.png"), p2,
       width = 12, height = 7, dpi = 300)
cat("Saved: 101_GaR_Coefficient_Comparison.png\n")

# Plot 3: Predictive Density at Different FCI Levels (Fan Chart)
# Simulate density at low (-1), neutral (0), and high (+1) FCI

if (nrow(gar_results %>% filter(horizon == 12, !is.na(coef))) >= 5) {
  # Get coefficients at h=12
  coefs_12 <- gar_results %>%
    filter(horizon == 12, !is.na(coef)) %>%
    arrange(quantile)

  # Simulate for FCI = -1 (loose), 0 (neutral), +1 (tight)
  fci_scenarios <- c(-1, 0, 1)
  scenario_labels <- c("FCI = -1 (Loose)", "FCI = 0 (Neutral)", "FCI = +1 (Tight)")

  fan_data <- data.frame()
  for (i in seq_along(fci_scenarios)) {
    fci_val <- fci_scenarios[i]
    # Predicted values at each quantile
    for (j in 1:nrow(coefs_12)) {
      # Simplified: just show FCI contribution
      predicted <- coefs_12$coef[j] * fci_val
      fan_data <- rbind(fan_data, data.frame(
        scenario = scenario_labels[i],
        fci_value = fci_val,
        quantile = coefs_12$quantile[j],
        predicted = predicted
      ))
    }
  }

  p3 <- ggplot(fan_data, aes(x = quantile, y = predicted, color = scenario, fill = scenario)) +
    geom_ribbon(aes(ymin = 0, ymax = predicted), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("FCI = -1 (Loose)" = "#27AE60",
                                  "FCI = 0 (Neutral)" = "#3498DB",
                                  "FCI = +1 (Tight)" = "#E74C3C")) +
    scale_fill_manual(values = c("FCI = -1 (Loose)" = "#27AE60",
                                 "FCI = 0 (Neutral)" = "#3498DB",
                                 "FCI = +1 (Tight)" = "#E74C3C")) +
    scale_x_continuous(breaks = GAR_CONFIG$quantiles,
                       labels = scales::percent(GAR_CONFIG$quantiles)) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Predictive Density Shift by FCI Scenario (h=12)",
      subtitle = "How financial conditions shift the growth distribution",
      x = "Quantile",
      y = "FCI Contribution to Growth (pp)",
      color = "Scenario", fill = "Scenario"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  ggsave(file.path(GAR_CONFIG$output_dir, "102_GaR_Predictive_Density.png"), p3,
         width = 10, height = 7, dpi = 300)
  cat("Saved: 102_GaR_Predictive_Density.png\n")
}


################################################################################
# 7. EXPORT RESULTS
################################################################################

cat("\n================================================================================\n")
cat("EXPORTING RESULTS\n")
cat("================================================================================\n\n")

# Export main results
write.csv(gar_results,
          file.path(GAR_CONFIG$output_dir, "Growth_at_Risk_Results.csv"),
          row.names = FALSE)
cat("Saved: Growth_at_Risk_Results.csv\n")

# Export asymmetry tests
write.csv(asymmetry_tests,
          file.path(GAR_CONFIG$output_dir, "GaR_Asymmetry_Tests.csv"),
          row.names = FALSE)
cat("Saved: GaR_Asymmetry_Tests.csv\n")


################################################################################
# 8. SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("GROWTH-AT-RISK SUMMARY\n")
cat("================================================================================\n\n")

cat("1. METHODOLOGY:\n")
cat("   - Quantile regressions at tau = {0.05, 0.10, ..., 0.95}\n")
cat("   - Horizons: 1, 3, 6, 12, 24 months ahead\n")
cat("   - Target: IMAEP year-over-year growth\n\n")

cat("2. KEY RESULTS:\n")
# Best horizon
best_h <- asymmetry_tests %>%
  filter(p_lower == min(p_lower)) %>%
  slice(1)

cat(sprintf("   - Strongest asymmetry at h=%d months\n", best_h$horizon))
cat(sprintf("   - 5th percentile coefficient: %.3f\n", best_h$coef_05))
cat(sprintf("   - 50th percentile coefficient: %.3f\n", best_h$coef_50))
cat(sprintf("   - p-value for difference: %.3f\n\n", best_h$p_lower))

cat("3. POLICY IMPLICATIONS:\n")
if (any_asymmetry) {
  cat("   - FCI tightening increases downside risk disproportionately\n")
  cat("   - Financial stability monitoring should focus on tail risks\n")
  cat("   - Early warning indicators should incorporate FCI levels\n")
} else {
  cat("   - FCI affects growth roughly symmetrically\n")
  cat("   - Standard conditional mean forecasts are appropriate\n")
  cat("   - No special focus on tail risks needed based on FCI\n")
}

cat("\n================================================================================\n")
cat("GROWTH-AT-RISK ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
cat("Plots: 100-102\n")
cat("Output:", GAR_CONFIG$output_dir, "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
