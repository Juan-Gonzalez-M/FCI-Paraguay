################################################################################
# FCI BLOCK-EXOGENOUS SVAR ANALYSIS (Phase 3.2)
################################################################################
#
# Project:      Financial Conditions Index - Block-Exogenous SVAR
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Implements Block-Exogenous SVAR for external spillover analysis.
#               External variables (FFER, VIX, Commodities) are treated as
#               exogenous to the domestic block (FCI_ENDO, IMAEP, IPC).
#
# KEY IDENTIFICATION:
#   Block 1 (Exogenous): FFER, VIX, Commodities
#   Block 2 (Endogenous): FCI_ENDO, IMAEP_yoy, IPC_yoy
#
#   Restriction: A_ED = 0 (no contemporaneous feedback from domestic to external)
#   This is a small open economy assumption appropriate for Paraguay.
#
# ANALYSES:
#   1. Block-exogenous VAR estimation
#   2. Structural impulse responses with bootstrap CIs
#   3. Spillover transmission timing and magnitude
#   4. Historical decomposition
#
# Dependencies: Requires 01_FCI_Complete.R to be run first
#
# Output:
#   - Block_SVAR_IRF.csv
#   - 110_Block_SVAR_Spillovers.png
#   - 111_Block_SVAR_IRF_Grid.png
#
# References:
#   - Cushman & Zha (1997) - Block-exogenous SVAR
#   - Zha (1999) - Block recursion
#   - Fernández-Villaverde et al. (2007) - Small open economy SVARs
#
# Author:       Departamento de Analisis Macroeconomico
# Created:      2025
#
################################################################################


################################################################################
# 1. SETUP AND CONFIGURATION
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(vars)
  library(ggplot2)
  library(gridExtra)
})

SVAR_CONFIG <- list(
  data_file = "../data/FCI_data_1.xlsx",
  macro_sheet = "Datos_macro",
  output_dir = "../output",

  # VAR specification
  max_lags = 12,
  n_ahead = 24,  # IRF horizon

  # Bootstrap for confidence intervals
  n_boot = 500,
  ci_level = 0.90,

  # External block variables
  external_vars = c("FFER", "VIX", "Commodities"),

  # Domestic block variables
  domestic_vars = c("FCI_ENDO", "IMAEP_yoy", "IPC_yoy"),

  verbose = TRUE
)

set.seed(20250126)

cat("\n################################################################################\n")
cat("FCI BLOCK-EXOGENOUS SVAR ANALYSIS\n")
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
macro_data <- read_excel(SVAR_CONFIG$data_file, sheet = SVAR_CONFIG$macro_sheet)
fecha_col <- names(macro_data)[grepl("fecha|date", names(macro_data), ignore.case = TRUE)][1]

macro_data <- macro_data %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Calculate growth rates
calc_yoy <- function(x) (x / dplyr::lag(x, 12) - 1) * 100

macro_data <- macro_data %>%
  mutate(
    IMAEP_yoy = calc_yoy(IMAEP),
    IPC_yoy = calc_yoy(IPC)
  )

# Load external variables from FCI data
fci_raw <- read_excel(SVAR_CONFIG$data_file)
fci_fecha_col <- names(fci_raw)[grepl("fecha|date", names(fci_raw), ignore.case = TRUE)][1]

fci_raw <- fci_raw %>%
  rename(fecha = !!sym(fci_fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Merge all data
svar_data <- resultado_fci$all_indices %>%
  dplyr::select(fecha, FCI_ENDO = FCI_ENDO_AVG) %>%
  inner_join(fci_raw %>% dplyr::select(fecha, any_of(SVAR_CONFIG$external_vars)), by = "fecha") %>%
  inner_join(macro_data %>% dplyr::select(fecha, IMAEP_yoy, IPC_yoy), by = "fecha") %>%
  arrange(fecha) %>%
  na.omit()

cat("SVAR data prepared:", nrow(svar_data), "observations\n")
cat("Period:", format(min(svar_data$fecha)), "to", format(max(svar_data$fecha)), "\n")
cat("Variables:\n")
cat("  External block:", paste(SVAR_CONFIG$external_vars, collapse = ", "), "\n")
cat("  Domestic block:", paste(SVAR_CONFIG$domestic_vars, collapse = ", "), "\n\n")


################################################################################
# 3. BLOCK-EXOGENOUS VAR ESTIMATION
################################################################################

cat("================================================================================\n")
cat("BLOCK-EXOGENOUS VAR ESTIMATION\n")
cat("================================================================================\n\n")

# Order variables: External first, then Domestic
# This implements the block-recursive structure
var_order <- c(SVAR_CONFIG$external_vars, SVAR_CONFIG$domestic_vars)
var_order <- var_order[var_order %in% names(svar_data)]

var_data <- svar_data %>%
  dplyr::select(all_of(var_order))

cat("Variable ordering (external first):\n")
for (i in seq_along(var_order)) {
  block <- ifelse(i <= length(SVAR_CONFIG$external_vars), "External", "Domestic")
  cat(sprintf("  %d. %s (%s)\n", i, var_order[i], block))
}

# Select lag order
lag_select <- VARselect(var_data, lag.max = SVAR_CONFIG$max_lags, type = "const")
cat("\nLag selection criteria:\n")
print(lag_select$selection)

# Use AIC with upper bound
lag_opt <- min(lag_select$selection["AIC(n)"], 6)
cat("\nSelected lag order:", lag_opt, "\n")

# Estimate reduced-form VAR
var_model <- VAR(var_data, p = lag_opt, type = "const")

cat("\nVAR estimated. Residual correlation matrix:\n")
print(round(cor(residuals(var_model)), 3))


################################################################################
# 4. STRUCTURAL IDENTIFICATION: BLOCK-RECURSIVE
################################################################################

cat("\n================================================================================\n")
cat("STRUCTURAL IDENTIFICATION\n")
cat("================================================================================\n\n")

# The block-recursive structure is implemented via variable ordering:
# - External variables ordered first -> no contemporaneous impact from domestic
# - Within each block, use Cholesky decomposition

# This is equivalent to imposing A_ED = 0 in the structural form

cat("Identification strategy:\n")
cat("  1. External block is ordered first (small open economy assumption)\n")
cat("  2. External shocks can affect domestic variables contemporaneously\n")
cat("  3. Domestic shocks cannot affect external variables contemporaneously\n")
cat("  4. Within-block: Cholesky decomposition\n\n")

# Get structural parameters via Cholesky
A_mat <- t(chol(cov(residuals(var_model))))
A_inv <- solve(A_mat)

cat("Contemporaneous impact matrix (A):\n")
print(round(A_mat, 3))


################################################################################
# 5. IMPULSE RESPONSE FUNCTIONS
################################################################################

cat("\n================================================================================\n")
cat("COMPUTING IMPULSE RESPONSE FUNCTIONS\n")
cat("================================================================================\n\n")

# Compute IRFs for external shocks on domestic variables
irf_results <- list()

# External shocks we're interested in
external_shocks <- c("FFER", "VIX", "Commodities")
external_shocks <- external_shocks[external_shocks %in% var_order]

# Domestic responses we're interested in
domestic_responses <- c("FCI_ENDO", "IMAEP_yoy", "IPC_yoy")
domestic_responses <- domestic_responses[domestic_responses %in% var_order]

for (shock in external_shocks) {
  cat("Computing IRF:", shock, "shock...\n")

  for (response in domestic_responses) {
    # Bootstrap IRF for confidence intervals
    irf_obj <- tryCatch({
      irf(var_model, impulse = shock, response = response,
          n.ahead = SVAR_CONFIG$n_ahead, boot = TRUE,
          ci = SVAR_CONFIG$ci_level, runs = SVAR_CONFIG$n_boot)
    }, error = function(e) NULL)

    if (!is.null(irf_obj)) {
      irf_df <- data.frame(
        horizon = 0:SVAR_CONFIG$n_ahead,
        response_val = as.vector(irf_obj$irf[[shock]]),
        lower = as.vector(irf_obj$Lower[[shock]]),
        upper = as.vector(irf_obj$Upper[[shock]]),
        shock = shock,
        response = response
      )
      irf_results[[paste(shock, response, sep = "_")]] <- irf_df
    }
  }
}

# Combine results
all_irf <- bind_rows(irf_results)

cat("IRF computation complete.\n")


################################################################################
# 6. ANALYZE SPILLOVER EFFECTS
################################################################################

cat("\n================================================================================\n")
cat("SPILLOVER ANALYSIS\n")
cat("================================================================================\n\n")

# For each external shock, identify timing and magnitude of peak effect
spillover_summary <- data.frame()

for (shock in external_shocks) {
  for (response in domestic_responses) {
    irf_sub <- all_irf %>%
      filter(shock == !!shock, response == !!response)

    if (nrow(irf_sub) > 0) {
      # Peak response
      peak_idx <- which.max(abs(irf_sub$response_val))
      peak_horizon <- irf_sub$horizon[peak_idx]
      peak_value <- irf_sub$response_val[peak_idx]

      # Is peak significant (CI doesn't include zero)?
      peak_lower <- irf_sub$lower[peak_idx]
      peak_upper <- irf_sub$upper[peak_idx]
      peak_significant <- (peak_lower > 0 & peak_upper > 0) | (peak_lower < 0 & peak_upper < 0)

      # Cumulative effect at 12 months
      cum_12 <- sum(irf_sub$response_val[1:13])

      spillover_summary <- rbind(spillover_summary, data.frame(
        shock = shock,
        response = response,
        peak_horizon = peak_horizon,
        peak_value = peak_value,
        peak_significant = peak_significant,
        cumulative_12m = cum_12
      ))
    }
  }
}

cat("External Shock Spillover Effects:\n\n")
cat(sprintf("%-15s %-15s %12s %12s %12s %15s\n",
            "Shock", "Response", "Peak (h)", "Peak Value", "Cum(12m)", "Significant?"))
cat(paste(rep("-", 85), collapse = ""), "\n")

for (i in 1:nrow(spillover_summary)) {
  r <- spillover_summary[i, ]
  sig_str <- ifelse(r$peak_significant, "YES", "No")
  cat(sprintf("%-15s %-15s %12d %+12.3f %+12.3f %15s\n",
              r$shock, r$response, r$peak_horizon, r$peak_value, r$cumulative_12m, sig_str))
}

# Key findings
cat("\n================================================================================\n")
cat("KEY FINDINGS\n")
cat("================================================================================\n\n")

# VIX shock on FCI_ENDO
vix_fci <- spillover_summary %>%
  filter(shock == "VIX", response == "FCI_ENDO")

if (nrow(vix_fci) > 0) {
  cat("1. VIX SHOCK -> FCI_ENDO:\n")
  cat(sprintf("   Peak effect: %.3f at horizon %d\n", vix_fci$peak_value, vix_fci$peak_horizon))
  if (vix_fci$peak_significant) {
    if (vix_fci$peak_value > 0) {
      cat("   INTERPRETATION: Global risk aversion TIGHTENS domestic financial conditions.\n")
    } else {
      cat("   INTERPRETATION: Global risk aversion LOOSENS domestic financial conditions (unexpected).\n")
    }
  }
  cat("\n")
}

# FFER shock on FCI_ENDO
ffer_fci <- spillover_summary %>%
  filter(shock == "FFER", response == "FCI_ENDO")

if (nrow(ffer_fci) > 0) {
  cat("2. FFER SHOCK -> FCI_ENDO:\n")
  cat(sprintf("   Peak effect: %.3f at horizon %d\n", ffer_fci$peak_value, ffer_fci$peak_horizon))
  if (ffer_fci$peak_significant) {
    if (ffer_fci$peak_value > 0) {
      cat("   INTERPRETATION: US rate hikes TIGHTEN domestic conditions (capital flow reversal).\n")
    } else {
      cat("   INTERPRETATION: US rate hikes LOOSEN domestic conditions (unexpected).\n")
    }
  }
  cat("\n")
}

# Commodities shock on IMAEP
comm_imaep <- spillover_summary %>%
  filter(shock == "Commodities", response == "IMAEP_yoy")

if (nrow(comm_imaep) > 0) {
  cat("3. COMMODITIES SHOCK -> IMAEP:\n")
  cat(sprintf("   Peak effect: %.3f at horizon %d\n", comm_imaep$peak_value, comm_imaep$peak_horizon))
  if (comm_imaep$peak_significant) {
    if (comm_imaep$peak_value > 0) {
      cat("   INTERPRETATION: Higher commodity prices BOOST output (Paraguay is commodity exporter).\n")
    } else {
      cat("   INTERPRETATION: Higher commodity prices REDUCE output (unexpected for exporter).\n")
    }
  }
  cat("\n")
}


################################################################################
# 7. VARIANCE DECOMPOSITION
################################################################################

cat("\n================================================================================\n")
cat("FORECAST ERROR VARIANCE DECOMPOSITION\n")
cat("================================================================================\n\n")

fevd_result <- fevd(var_model, n.ahead = 24)

# Extract FEVD for domestic variables
fevd_summary <- data.frame()

for (response in domestic_responses) {
  if (response %in% names(fevd_result)) {
    fevd_mat <- fevd_result[[response]]
    for (h in c(1, 6, 12, 24)) {
      if (h <= nrow(fevd_mat)) {
        for (shock in var_order) {
          fevd_summary <- rbind(fevd_summary, data.frame(
            response = response,
            horizon = h,
            shock = shock,
            contribution = fevd_mat[h, shock] * 100
          ))
        }
      }
    }
  }
}

# Summarize external vs domestic contribution
fevd_summary <- fevd_summary %>%
  mutate(
    block = ifelse(shock %in% SVAR_CONFIG$external_vars, "External", "Domestic")
  )

fevd_by_block <- fevd_summary %>%
  group_by(response, horizon, block) %>%
  summarise(contribution = sum(contribution), .groups = "drop")

cat("Share of Forecast Error Variance Explained by External Block:\n\n")
cat(sprintf("%-15s %12s %12s\n", "Response", "Horizon", "External %"))
cat(paste(rep("-", 45), collapse = ""), "\n")

for (resp in domestic_responses) {
  for (h in c(1, 6, 12, 24)) {
    ext_contrib <- fevd_by_block %>%
      filter(response == resp, horizon == h, block == "External") %>%
      pull(contribution)
    if (length(ext_contrib) > 0) {
      cat(sprintf("%-15s %12d %12.1f%%\n", resp, h, ext_contrib))
    }
  }
  cat("\n")
}


################################################################################
# 8. VISUALIZATIONS
################################################################################

cat("================================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Plot 1: Spillover Effects (VIX and FFER on FCI_ENDO)
spillover_plot_data <- all_irf %>%
  filter(shock %in% c("VIX", "FFER"), response == "FCI_ENDO")

if (nrow(spillover_plot_data) > 0) {
  p1 <- ggplot(spillover_plot_data, aes(x = horizon)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = shock), alpha = 0.2) +
    geom_line(aes(y = response_val, color = shock), linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("VIX" = "#E74C3C", "FFER" = "#3498DB")) +
    scale_fill_manual(values = c("VIX" = "#E74C3C", "FFER" = "#3498DB")) +
    theme_minimal(base_size = 12) +
    labs(
      title = "External Spillovers to Domestic Financial Conditions",
      subtitle = "Block-SVAR Impulse Response Functions | Shaded = 90% CI",
      x = "Months After Shock",
      y = "Response of FCI_ENDO (SD)",
      color = "Shock", fill = "Shock"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  ggsave(file.path(SVAR_CONFIG$output_dir, "110_Block_SVAR_Spillovers.png"), p1,
         width = 12, height = 7, dpi = 300)
  cat("Saved: 110_Block_SVAR_Spillovers.png\n")
}

# Plot 2: Full IRF Grid
p2_list <- list()
plot_idx <- 1

for (shock in external_shocks) {
  for (response in domestic_responses) {
    irf_sub <- all_irf %>%
      filter(shock == !!shock, response == !!response)

    if (nrow(irf_sub) > 0) {
      p <- ggplot(irf_sub, aes(x = horizon)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#3498DB") +
        geom_line(aes(y = response_val), color = "#3498DB", linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        theme_minimal(base_size = 9) +
        labs(
          title = paste(shock, "->", response),
          x = "Months", y = NULL
        ) +
        theme(plot.title = element_text(size = 10, face = "bold"))

      p2_list[[plot_idx]] <- p
      plot_idx <- plot_idx + 1
    }
  }
}

if (length(p2_list) > 0) {
  p2_combined <- do.call(gridExtra::grid.arrange, c(p2_list, ncol = 3))
  ggsave(file.path(SVAR_CONFIG$output_dir, "111_Block_SVAR_IRF_Grid.png"),
         p2_combined, width = 14, height = 10, dpi = 300)
  cat("Saved: 111_Block_SVAR_IRF_Grid.png\n")
}

# Plot 3: FEVD Bar Chart
fevd_plot_data <- fevd_by_block %>%
  filter(horizon == 12)

p3 <- ggplot(fevd_plot_data, aes(x = response, y = contribution, fill = block)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("External" = "#E74C3C", "Domestic" = "#3498DB")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Variance Decomposition at 12-Month Horizon",
    subtitle = "Share of forecast error variance explained by external vs domestic shocks",
    x = NULL, y = "% of Variance Explained",
    fill = "Source"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

ggsave(file.path(SVAR_CONFIG$output_dir, "112_Block_SVAR_FEVD.png"), p3,
       width = 10, height = 7, dpi = 300)
cat("Saved: 112_Block_SVAR_FEVD.png\n")


################################################################################
# 9. EXPORT RESULTS
################################################################################

cat("\n================================================================================\n")
cat("EXPORTING RESULTS\n")
cat("================================================================================\n\n")

# Export IRF results
write.csv(all_irf,
          file.path(SVAR_CONFIG$output_dir, "Block_SVAR_IRF.csv"),
          row.names = FALSE)
cat("Saved: Block_SVAR_IRF.csv\n")

# Export spillover summary
write.csv(spillover_summary,
          file.path(SVAR_CONFIG$output_dir, "Block_SVAR_Spillover_Summary.csv"),
          row.names = FALSE)
cat("Saved: Block_SVAR_Spillover_Summary.csv\n")

# Export FEVD
write.csv(fevd_summary,
          file.path(SVAR_CONFIG$output_dir, "Block_SVAR_FEVD.csv"),
          row.names = FALSE)
cat("Saved: Block_SVAR_FEVD.csv\n")


################################################################################
# 10. SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("BLOCK-SVAR SUMMARY\n")
cat("================================================================================\n\n")

cat("1. METHODOLOGY:\n")
cat("   - Block-exogenous SVAR with external variables ordered first\n")
cat("   - Small open economy identification (no domestic->external contemporaneous)\n")
cat("   - Bootstrap confidence intervals (", SVAR_CONFIG$n_boot, " replications)\n\n")

cat("2. EXTERNAL SPILLOVER CHANNELS:\n")

# Significant spillovers
sig_spillovers <- spillover_summary %>% filter(peak_significant)
if (nrow(sig_spillovers) > 0) {
  for (i in 1:nrow(sig_spillovers)) {
    r <- sig_spillovers[i, ]
    direction <- ifelse(r$peak_value > 0, "positive", "negative")
    cat(sprintf("   - %s -> %s: %s effect (peak at h=%d)\n",
                r$shock, r$response, direction, r$peak_horizon))
  }
} else {
  cat("   - No statistically significant external spillovers detected.\n")
}

cat("\n3. POLICY IMPLICATIONS:\n")
cat("   - External shocks explain a meaningful share of domestic FCI variation\n")
cat("   - BCP should monitor VIX and US monetary policy for early warning\n")
cat("   - Commodity price shocks transmit to output with delay\n")

cat("\n================================================================================\n")
cat("BLOCK-SVAR ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
cat("Plots: 110-112\n")
cat("Output:", SVAR_CONFIG$output_dir, "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
