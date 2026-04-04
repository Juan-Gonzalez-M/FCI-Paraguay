################################################################################
# FCI ROLLING/EXPANDING WINDOW STABILITY ANALYSIS
################################################################################
#
# Project:      Financial Conditions Index - Parameter Stability
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Addresses referee concerns about parameter stability by
#               analyzing time-varying PCA loadings and weight convergence.
#
# ANALYSES INCLUDED:
#   1. Rolling PCA Loadings (windows: 36, 48, 60, 72 months)
#   2. Expanding Window Weight Convergence
#   3. Stability Metrics:
#      - Time-varying loadings vs full-sample
#      - Rolling correlation with benchmark
#      - Sign consistency across periods
#
# Dependencies: Requires 01_FCI_Complete.R to be run first
#
# Output:
#   - 60_Rolling_PCA_Loadings.png
#   - 61_Expanding_Weight_Convergence.png
#   - 62_Weight_Stability_Correlation.png
#   - 63_Weight_Stability_Heatmap.png
#   - FCI_Weight_Stability_Summary.csv
#   - FCI_Rolling_Loadings.csv
#
# Author:       Departamento de Analisis Macroeconomico
# Created:      2025-01-17
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
  library(ggplot2)
  library(gridExtra)
})

STABILITY_CONFIG <- list(
  # Data
  data_file = "../data/FCI_data_1.xlsx",

  # Rolling window sizes (months)
  window_sizes = c(36, 48, 60, 72),

  # Minimum observations for analysis
  min_obs = 36,

  # Output
  output_dir = "../output"
)

set.seed(20250117)

cat("\n################################################################################\n")
cat("FCI ROLLING/EXPANDING WINDOW STABILITY ANALYSIS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# 2. DATA PREPARATION
################################################################################

cat("================================================================================\n")
cat("DATA PREPARATION\n")
cat("================================================================================\n\n")

# Load FCI results (run core script if needed)
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  source("01_FCI_Complete.R")
}

# Get variable definitions from resultado_fci
VARIABLES <- resultado_fci$variables$level3
all_vars <- c(VARIABLES$rates$vars, VARIABLES$banking$vars, VARIABLES$external$vars)
all_signs <- c(VARIABLES$rates$signs, VARIABLES$banking$signs, VARIABLES$external$signs)

# Load and prepare raw data
datos_raw <- read_excel(STABILITY_CONFIG$data_file)
fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

datos <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Calculate credit growth if needed
if ("Creditos_Sector_privado_totales" %in% names(datos)) {
  datos <- datos %>%
    mutate(Crecimiento_creditos = (Creditos_Sector_privado_totales /
                                    lag(Creditos_Sector_privado_totales, 12) - 1) * 100)
}

# Apply signs
apply_signs <- function(data, vars, signs) {
  data_adj <- data
  for (i in seq_along(vars)) {
    if (vars[i] %in% names(data_adj)) {
      data_adj[[vars[i]]] <- data_adj[[vars[i]]] * signs[i]
    }
  }
  return(data_adj)
}

datos_signed <- apply_signs(datos, all_vars, all_signs)

cat("Data prepared:\n")
cat("  Total observations:", nrow(datos), "\n")
cat("  Period:", format(min(datos$fecha)), "to", format(max(datos$fecha)), "\n")
cat("  Variables:", length(all_vars), "\n\n")


################################################################################
# 3. ROLLING PCA LOADINGS ANALYSIS
################################################################################

cat("================================================================================\n")
cat("ROLLING PCA LOADINGS ANALYSIS\n")
cat("================================================================================\n\n")

# Function to calculate PCA loadings for a given window
calc_pca_loadings <- function(data, vars) {
  vars_available <- intersect(vars, names(data))
  if (length(vars_available) < 2) return(NULL)

  datos_subset <- data[, vars_available, drop = FALSE]

  # Standardize
  datos_std <- scale(datos_subset)

  # Remove rows with NAs
  complete_rows <- complete.cases(datos_std)
  if (sum(complete_rows) < 10) return(NULL)

  datos_clean <- datos_std[complete_rows, , drop = FALSE]

  # Impute remaining NAs with column medians
  for (col in 1:ncol(datos_clean)) {
    if (any(is.na(datos_clean[, col]))) {
      datos_clean[is.na(datos_clean[, col]), col] <- median(datos_clean[, col], na.rm = TRUE)
    }
  }

  tryCatch({
    pca_model <- prcomp(datos_clean, center = FALSE, scale. = FALSE)
    loadings <- pca_model$rotation[, 1]  # First PC loadings

    # Ensure consistent sign (positive correlation with average)
    avg_data <- rowMeans(datos_clean)
    scores <- datos_clean %*% loadings
    if (cor(scores, avg_data) < 0) {
      loadings <- -loadings
    }

    return(loadings)
  }, error = function(e) {
    return(NULL)
  })
}

# Calculate rolling loadings
rolling_loadings <- list()

for (window_size in STABILITY_CONFIG$window_sizes) {
  cat(sprintf("Processing window size: %d months...\n", window_size))

  n_obs <- nrow(datos_signed)
  n_windows <- n_obs - window_size + 1

  if (n_windows < 10) {
    cat(sprintf("  Insufficient data for %d-month window\n", window_size))
    next
  }

  for (i in 1:n_windows) {
    start_idx <- i
    end_idx <- i + window_size - 1

    window_data <- datos_signed[start_idx:end_idx, ]
    window_end_date <- datos_signed$fecha[end_idx]

    loadings <- calc_pca_loadings(window_data, all_vars)

    if (!is.null(loadings)) {
      for (var_name in names(loadings)) {
        rolling_loadings[[length(rolling_loadings) + 1]] <- data.frame(
          window_size = window_size,
          fecha = window_end_date,
          variable = var_name,
          loading = loadings[var_name]
        )
      }
    }
  }
}

rolling_loadings_df <- bind_rows(rolling_loadings)
cat("\nRolling loadings calculated:", nrow(rolling_loadings_df), "observations\n\n")


################################################################################
# 4. FULL-SAMPLE BENCHMARK LOADINGS
################################################################################

cat("================================================================================\n")
cat("FULL-SAMPLE BENCHMARK LOADINGS\n")
cat("================================================================================\n\n")

benchmark_loadings <- calc_pca_loadings(datos_signed, all_vars)

if (!is.null(benchmark_loadings)) {
  cat("Full-sample PCA loadings (PC1):\n")
  for (var_name in names(sort(abs(benchmark_loadings), decreasing = TRUE))) {
    cat(sprintf("  %-30s %+.3f\n", var_name, benchmark_loadings[var_name]))
  }

  # Create benchmark data frame
  benchmark_df <- data.frame(
    variable = names(benchmark_loadings),
    benchmark_loading = as.numeric(benchmark_loadings)
  )
}


################################################################################
# 5. EXPANDING WINDOW CONVERGENCE
################################################################################

cat("\n================================================================================\n")
cat("EXPANDING WINDOW CONVERGENCE ANALYSIS\n")
cat("================================================================================\n\n")

expanding_loadings <- list()
min_window <- 36  # Start with minimum 3 years

n_obs <- nrow(datos_signed)

for (end_idx in min_window:n_obs) {
  window_data <- datos_signed[1:end_idx, ]
  window_end_date <- datos_signed$fecha[end_idx]

  loadings <- calc_pca_loadings(window_data, all_vars)

  if (!is.null(loadings)) {
    for (var_name in names(loadings)) {
      expanding_loadings[[length(expanding_loadings) + 1]] <- data.frame(
        fecha = window_end_date,
        n_obs = end_idx,
        variable = var_name,
        loading = loadings[var_name]
      )
    }
  }
}

expanding_loadings_df <- bind_rows(expanding_loadings)
cat("Expanding window loadings calculated:", nrow(expanding_loadings_df), "observations\n\n")


################################################################################
# 6. STABILITY METRICS
################################################################################

cat("================================================================================\n")
cat("STABILITY METRICS\n")
cat("================================================================================\n\n")

# 6.1 Time-varying loadings vs full-sample
if (!is.null(benchmark_loadings)) {
  stability_metrics <- rolling_loadings_df %>%
    filter(window_size == 60) %>%  # Use 60-month window as reference
    left_join(benchmark_df, by = "variable") %>%
    group_by(variable) %>%
    summarise(
      mean_loading = mean(loading, na.rm = TRUE),
      sd_loading = sd(loading, na.rm = TRUE),
      min_loading = min(loading, na.rm = TRUE),
      max_loading = max(loading, na.rm = TRUE),
      range_loading = max_loading - min_loading,
      benchmark = first(benchmark_loading),
      # Stability = correlation with benchmark loading
      deviation_from_benchmark = mean(abs(loading - benchmark_loading), na.rm = TRUE),
      # Sign consistency
      sign_changes = sum(diff(sign(loading)) != 0, na.rm = TRUE),
      pct_positive = mean(loading > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    mutate(
      # Is sign consistent with benchmark?
      sign_consistent = sign(mean_loading) == sign(benchmark),
      # Coefficient of variation
      cv = abs(sd_loading / mean_loading) * 100
    ) %>%
    arrange(desc(abs(benchmark)))

  cat("Loading Stability (60-month rolling window):\n")
  cat(sprintf("%-25s %8s %8s %8s %8s %8s\n",
              "Variable", "Bench", "Mean", "SD", "CV%", "Sign"))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (i in 1:nrow(stability_metrics)) {
    r <- stability_metrics[i, ]
    cat(sprintf("%-25s %+8.3f %+8.3f %8.3f %8.1f %8s\n",
                r$variable, r$benchmark, r$mean_loading, r$sd_loading,
                r$cv, ifelse(r$sign_consistent, "OK", "CHANGE")))
  }
}

# 6.2 Rolling correlation with benchmark FCI
cat("\n\nRolling Correlation with Full-Sample FCI:\n")

# Calculate FCI for each rolling window and correlate with benchmark
fci_benchmark <- resultado_fci$comprehensive$FCI_COMP_AVG

rolling_fci_corr <- list()

for (window_size in STABILITY_CONFIG$window_sizes) {
  n_obs <- nrow(datos_signed)
  n_windows <- n_obs - window_size + 1

  if (n_windows < 10) next

  for (i in 1:n_windows) {
    start_idx <- i
    end_idx <- i + window_size - 1

    window_data <- datos_signed[start_idx:end_idx, ]

    # Calculate rolling FCI
    loadings <- calc_pca_loadings(window_data, all_vars)

    if (!is.null(loadings)) {
      # Calculate FCI scores for this window
      vars_available <- names(loadings)
      window_std <- scale(window_data[, vars_available, drop = FALSE])
      window_std[is.na(window_std)] <- 0
      fci_rolling <- as.vector(window_std %*% loadings)

      # Get corresponding benchmark FCI values
      benchmark_subset <- fci_benchmark[start_idx:end_idx]

      # Correlation
      if (sum(!is.na(fci_rolling) & !is.na(benchmark_subset)) >= 20) {
        corr <- cor(fci_rolling, benchmark_subset, use = "complete.obs")

        rolling_fci_corr[[length(rolling_fci_corr) + 1]] <- data.frame(
          window_size = window_size,
          fecha = datos_signed$fecha[end_idx],
          correlation = corr
        )
      }
    }
  }
}

rolling_fci_corr_df <- bind_rows(rolling_fci_corr)

# Summary by window size
if (nrow(rolling_fci_corr_df) > 0) {
  corr_summary <- rolling_fci_corr_df %>%
    group_by(window_size) %>%
    summarise(
      mean_corr = mean(correlation, na.rm = TRUE),
      sd_corr = sd(correlation, na.rm = TRUE),
      min_corr = min(correlation, na.rm = TRUE),
      n_windows = n(),
      .groups = "drop"
    )

  cat("\nFCI Correlation with Benchmark by Window Size:\n")
  print(corr_summary)
}


################################################################################
# 7. VISUALIZATIONS
################################################################################

cat("\n================================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

if (!dir.exists(STABILITY_CONFIG$output_dir)) {
  dir.create(STABILITY_CONFIG$output_dir, recursive = TRUE)
}

# Define colors for variables by channel
var_colors <- c(
  # Rates (reds)
  "TPM" = "#E41A1C", "Spread_activas_pasivas" = "#FC8D62", "Spread_mercado_TPM" = "#FDAE61",
  # Banking (blues)
  "Crecimiento_creditos" = "#377EB8", "Ratio_Cred_Depo" = "#4292C6",
  "Morosidad" = "#6BAED6", "Rentabilidad" = "#9ECAE1", "Liquidez" = "#C6DBEF",
  # External (greens)
  "TCN" = "#4DAF4A", "Commodities" = "#ADDD8E",
  "FFER" = "#D9F0A3", "VIX" = "#F7FCB4"
)

# Plot 1: Rolling PCA Loadings by Variable (faceted)
if (nrow(rolling_loadings_df) > 0) {
  # Add benchmark line data
  rolling_plot_data <- rolling_loadings_df %>%
    filter(window_size == 60) %>%
    left_join(benchmark_df, by = "variable")

  p1 <- ggplot(rolling_plot_data, aes(x = fecha)) +
    geom_line(aes(y = loading, color = variable), linewidth = 0.8) +
    geom_hline(aes(yintercept = benchmark_loading), linetype = "dashed",
               color = "red", linewidth = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
    facet_wrap(~variable, scales = "free_y", ncol = 4) +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    scale_color_manual(values = var_colors, guide = "none") +
    theme_minimal(base_size = 10) +
    labs(title = "Rolling PCA Loadings by Variable",
         subtitle = "60-month rolling window | Red dashed = full-sample benchmark",
         x = NULL, y = "PC1 Loading") +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(STABILITY_CONFIG$output_dir, "60_Rolling_PCA_Loadings.png"), p1,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 60_Rolling_PCA_Loadings.png\n")
}

# Plot 2: Expanding Window Convergence
if (nrow(expanding_loadings_df) > 0) {
  # Select key variables
  key_vars <- c("TPM", "Crecimiento_creditos", "TCN", "VIX")
  key_vars <- key_vars[key_vars %in% unique(expanding_loadings_df$variable)]

  expanding_plot_data <- expanding_loadings_df %>%
    filter(variable %in% key_vars) %>%
    left_join(benchmark_df, by = "variable")

  p2 <- ggplot(expanding_plot_data, aes(x = n_obs)) +
    geom_line(aes(y = loading, color = variable), linewidth = 1) +
    geom_hline(aes(yintercept = benchmark_loading, color = variable),
               linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0, color = "gray50") +
    scale_color_manual(values = var_colors) +
    theme_minimal(base_size = 12) +
    labs(title = "Expanding Window: Loading Convergence",
         subtitle = "Solid lines = loadings with increasing sample | Dashed = final value",
         x = "Number of Observations", y = "PC1 Loading",
         color = "Variable") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")

  ggsave(file.path(STABILITY_CONFIG$output_dir, "61_Expanding_Weight_Convergence.png"), p2,
         width = 12, height = 7, dpi = 300)
  cat("Saved: 61_Expanding_Weight_Convergence.png\n")
}

# Plot 3: Rolling Correlation with Benchmark
if (nrow(rolling_fci_corr_df) > 0) {
  p3 <- ggplot(rolling_fci_corr_df, aes(x = fecha, y = correlation,
                                         color = factor(window_size))) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray50") +
    scale_color_brewer(palette = "Set1", name = "Window (months)") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 12) +
    labs(title = "FCI Stability: Correlation with Full-Sample Benchmark",
         subtitle = "Rolling window FCI correlated with full-sample FCI | Dashed line = 0.9 threshold",
         x = NULL, y = "Correlation") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(STABILITY_CONFIG$output_dir, "62_Weight_Stability_Correlation.png"), p3,
         width = 12, height = 6, dpi = 300)
  cat("Saved: 62_Weight_Stability_Correlation.png\n")
}

# Plot 4: Stability Heatmap
if (exists("stability_metrics") && nrow(stability_metrics) > 0) {
  # Create stability score (lower is more stable)
  stability_heatmap <- stability_metrics %>%
    mutate(
      stability_score = cv,  # Use CV as stability measure
      variable = factor(variable, levels = variable[order(cv)])
    )

  p4 <- ggplot(stability_heatmap, aes(x = 1, y = variable, fill = stability_score)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.1f%%", cv)), size = 3.5) +
    scale_fill_gradient2(low = "#27AE60", mid = "#F1C40F", high = "#E74C3C",
                         midpoint = median(stability_heatmap$cv),
                         name = "CV (%)") +
    theme_minimal(base_size = 12) +
    labs(title = "Loading Stability by Variable",
         subtitle = "Coefficient of Variation (lower = more stable)",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  ggsave(file.path(STABILITY_CONFIG$output_dir, "63_Weight_Stability_Heatmap.png"), p4,
         width = 8, height = 8, dpi = 300)
  cat("Saved: 63_Weight_Stability_Heatmap.png\n")
}


################################################################################
# 8. EXPORT RESULTS
################################################################################

cat("\n================================================================================\n")
cat("EXPORTING RESULTS\n")
cat("================================================================================\n\n")

# Export stability metrics
if (exists("stability_metrics") && nrow(stability_metrics) > 0) {
  write.csv(stability_metrics,
            file.path(STABILITY_CONFIG$output_dir, "FCI_Weight_Stability_Summary.csv"),
            row.names = FALSE)
  cat("Saved: FCI_Weight_Stability_Summary.csv\n")
}

# Export rolling loadings
if (nrow(rolling_loadings_df) > 0) {
  write.csv(rolling_loadings_df,
            file.path(STABILITY_CONFIG$output_dir, "FCI_Rolling_Loadings.csv"),
            row.names = FALSE)
  cat("Saved: FCI_Rolling_Loadings.csv\n")
}

# Export expanding loadings
if (nrow(expanding_loadings_df) > 0) {
  write.csv(expanding_loadings_df,
            file.path(STABILITY_CONFIG$output_dir, "FCI_Expanding_Loadings.csv"),
            row.names = FALSE)
  cat("Saved: FCI_Expanding_Loadings.csv\n")
}

# Export rolling correlations
if (nrow(rolling_fci_corr_df) > 0) {
  write.csv(rolling_fci_corr_df,
            file.path(STABILITY_CONFIG$output_dir, "FCI_Rolling_Correlations.csv"),
            row.names = FALSE)
  cat("Saved: FCI_Rolling_Correlations.csv\n")
}

# =============================================================================
# WINDOW SENSITIVITY SUMMARY (Phase 2.3)
# Compare stability metrics across different window sizes
# =============================================================================

cat("\n================================================================================\n")
cat("WINDOW SENSITIVITY ANALYSIS\n")
cat("================================================================================\n\n")

window_sensitivity_results <- data.frame()

for (window_size in STABILITY_CONFIG$window_sizes) {
  # Calculate stability metrics for each window size
  window_loadings <- rolling_loadings_df %>%
    filter(window_size == !!window_size)

  if (nrow(window_loadings) > 0) {
    window_metrics <- window_loadings %>%
      left_join(benchmark_df, by = "variable") %>%
      group_by(variable) %>%
      summarise(
        mean_loading = mean(loading, na.rm = TRUE),
        sd_loading = sd(loading, na.rm = TRUE),
        cv = abs(sd_loading / mean_loading) * 100,
        sign_consistent = sign(mean_loading) == sign(first(benchmark_loading)),
        .groups = "drop"
      )

    # Summary for this window
    window_sensitivity_results <- rbind(window_sensitivity_results, data.frame(
      window_size = window_size,
      mean_cv = mean(window_metrics$cv, na.rm = TRUE),
      median_cv = median(window_metrics$cv, na.rm = TRUE),
      max_cv = max(window_metrics$cv, na.rm = TRUE),
      pct_sign_consistent = mean(window_metrics$sign_consistent, na.rm = TRUE) * 100,
      n_windows = length(unique(window_loadings$fecha))
    ))
  }
}

if (nrow(rolling_fci_corr_df) > 0) {
  # Add mean correlation with benchmark
  corr_by_window <- rolling_fci_corr_df %>%
    group_by(window_size) %>%
    summarise(mean_corr_benchmark = mean(correlation, na.rm = TRUE), .groups = "drop")

  window_sensitivity_results <- window_sensitivity_results %>%
    left_join(corr_by_window, by = "window_size")
}

if (nrow(window_sensitivity_results) > 0) {
  cat("Window Sensitivity Summary:\n")
  cat(sprintf("%-12s %10s %10s %10s %12s %12s\n",
              "Window", "Mean CV", "Median CV", "Max CV", "Sign Cons%", "Corr w/Bench"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in 1:nrow(window_sensitivity_results)) {
    r <- window_sensitivity_results[i, ]
    corr_val <- ifelse("mean_corr_benchmark" %in% names(r) && !is.na(r$mean_corr_benchmark),
                       sprintf("%.3f", r$mean_corr_benchmark), "N/A")
    cat(sprintf("%-12d %10.1f%% %10.1f%% %10.1f%% %12.1f%% %12s\n",
                r$window_size, r$mean_cv, r$median_cv, r$max_cv,
                r$pct_sign_consistent, corr_val))
  }

  cat("\nInterpretation:\n")
  cat("  - Shorter windows (36mo) capture time-varying relationships but higher CV\n")
  cat("  - Longer windows (72mo) more stable but may miss structural changes\n")
  cat("  - 60-month window provides balance between stability and responsiveness\n")

  # Export window sensitivity
  write.csv(window_sensitivity_results,
            file.path(STABILITY_CONFIG$output_dir, "FCI_Window_Sensitivity.csv"),
            row.names = FALSE)
  cat("\nSaved: FCI_Window_Sensitivity.csv\n")

  # Create comparison plot
  if (nrow(rolling_loadings_df) > 0) {
    # Select key variables for comparison
    key_vars <- c("TPM", "Crecimiento_creditos", "TCN", "VIX")
    key_vars <- key_vars[key_vars %in% unique(rolling_loadings_df$variable)]

    if (length(key_vars) > 0) {
      window_plot_data <- rolling_loadings_df %>%
        filter(variable %in% key_vars) %>%
        mutate(window_label = paste0(window_size, "m"))

      p_window_compare <- ggplot(window_plot_data,
                                  aes(x = fecha, y = loading, color = factor(window_size))) +
        geom_line(alpha = 0.7, linewidth = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        facet_wrap(~variable, scales = "free_y", ncol = 2) +
        scale_color_brewer(palette = "Set1", name = "Window (months)") +
        scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
        theme_minimal(base_size = 11) +
        labs(
          title = "Rolling Window Sensitivity: PCA Loadings by Window Size",
          subtitle = "Comparing 36, 48, 60, and 72-month rolling windows",
          x = NULL, y = "PC1 Loading"
        ) +
        theme(
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1)
        )

      ggsave(file.path(STABILITY_CONFIG$output_dir, "64_Rolling_Window_Comparison.png"),
             p_window_compare, width = 12, height = 8, dpi = 300)
      cat("Saved: 64_Rolling_Window_Comparison.png\n")
    }
  }
}


################################################################################
# 9. KEY FINDINGS SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("KEY FINDINGS SUMMARY\n")
cat("================================================================================\n\n")

if (exists("stability_metrics") && nrow(stability_metrics) > 0) {
  # Overall stability assessment
  mean_cv <- mean(stability_metrics$cv, na.rm = TRUE)
  max_cv <- max(stability_metrics$cv, na.rm = TRUE)
  n_sign_consistent <- sum(stability_metrics$sign_consistent)
  n_vars <- nrow(stability_metrics)

  cat("1. LOADING STABILITY:\n")
  cat(sprintf("   Mean CV: %.1f%% (lower is more stable)\n", mean_cv))
  cat(sprintf("   Max CV:  %.1f%% (variable: %s)\n", max_cv,
              stability_metrics$variable[which.max(stability_metrics$cv)]))
  cat(sprintf("   Sign consistency: %d/%d variables (%.0f%%)\n\n",
              n_sign_consistent, n_vars, n_sign_consistent/n_vars*100))

  # Most/least stable variables
  cat("2. MOST STABLE VARIABLES (lowest CV):\n")
  most_stable <- stability_metrics %>% arrange(cv) %>% head(3)
  for (i in 1:nrow(most_stable)) {
    cat(sprintf("   - %s: CV=%.1f%%\n", most_stable$variable[i], most_stable$cv[i]))
  }

  cat("\n3. LEAST STABLE VARIABLES (highest CV):\n")
  least_stable <- stability_metrics %>% arrange(desc(cv)) %>% head(3)
  for (i in 1:nrow(least_stable)) {
    cat(sprintf("   - %s: CV=%.1f%%\n", least_stable$variable[i], least_stable$cv[i]))
  }
}

if (nrow(rolling_fci_corr_df) > 0) {
  min_corr <- min(rolling_fci_corr_df$correlation, na.rm = TRUE)
  mean_corr <- mean(rolling_fci_corr_df$correlation, na.rm = TRUE)

  cat("\n4. FCI CORRELATION WITH BENCHMARK:\n")
  cat(sprintf("   Mean correlation: %.3f\n", mean_corr))
  cat(sprintf("   Minimum correlation: %.3f\n", min_corr))

  if (min_corr > 0.8) {
    cat("   >>> HIGH STABILITY: Rolling FCI closely tracks benchmark\n")
  } else if (min_corr > 0.6) {
    cat("   >>> MODERATE STABILITY: Some variation in rolling FCI\n")
  } else {
    cat("   >>> LOW STABILITY: Significant variation in rolling FCI\n")
  }
}

cat("\n5. INTERPRETATION FOR REFEREES:\n")
if (exists("stability_metrics") && mean_cv < 30 && n_sign_consistent/n_vars > 0.8) {
  cat("   The FCI weights are REASONABLY STABLE over time.\n")
  cat("   Historical FCI values are comparable and interpretation is consistent.\n")
} else if (exists("stability_metrics")) {
  cat("   Some weight instability detected.\n")
  cat("   Consider reporting sub-period analyses or time-varying FCI.\n")
}

cat("\n================================================================================\n")
cat("ROLLING/EXPANDING STABILITY ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
cat("Plots: 60-63\n")
cat("Output:", STABILITY_CONFIG$output_dir, "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
