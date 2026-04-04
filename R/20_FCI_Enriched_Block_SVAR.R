################################################################################
# FCI ENRICHED BLOCK-EXOGENOUS SVAR
################################################################################
#
# Project:      Financial Conditions Index - Enriched External Block Analysis
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Extends Block-SVAR from 3-var to 7-var external block with
#               additional variables (DXY, US10Y, SP500, Selic, ToT). Adds
#               historical decomposition to identify shock drivers of FCI
#               during key episodes (Asian crisis, GFC, COVID).
#
#   SECTION A: Enriched Block-SVAR Estimation
#   SECTION B: Expanded IRFs
#   SECTION C: Enriched Variance Decomposition
#   SECTION D: Historical Decomposition
#   SECTION E: Comparison with Original 3-var Block-SVAR
#
# References:
#   - Cushman & Zha (1997) - Identifying Monetary Policy in a Small Open Economy
#   - Lutkepohl (2005) - New Introduction to Multiple Time Series Analysis
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

BSVAR_CONFIG <- list(
  max_var_lags    = 6,
  max_horizon     = 24,
  n_bootstrap     = 500,
  confidence_level = 0.90,
  output_dir      = "../output",
  verbose         = TRUE
)

set.seed(20260211)

cat("\n################################################################################\n")
cat("FCI ENRICHED BLOCK-EXOGENOUS SVAR\n")
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

# ---- Load new external variables ----
ext_file <- file.path(BSVAR_CONFIG$output_dir, "New_External_Variables.csv")
if (!file.exists(ext_file)) {
  cat("New_External_Variables.csv not found. Running script 17...\n")
  source("17_FCI_New_External_Data.R")
}
ext_data <- read.csv(ext_file)
ext_data$fecha <- as.Date(ext_data$fecha)

# ---- Load macro/credit data ----
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
  left_join(macro_data %>% dplyr::select(fecha, IMAEP_yoy, IPC_yoy, Cred_Real_Total),
            by = "fecha", suffix = c("", ".y")) %>%
  arrange(fecha)

# Clean duplicate columns
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


################################################################################
# SECTION A: ENRICHED BLOCK-SVAR ESTIMATION
################################################################################

cat("================================================================================\n")
cat("SECTION A: ENRICHED BLOCK-SVAR ESTIMATION\n")
cat("================================================================================\n\n")

# External block (7 variables)
external_vars <- c("VIX", "SP500", "DXY", "US_10Y", "FFER", "Selic_rate", "d_ToT")
# Domestic block (4 variables)
domestic_vars <- c("FCI_ENDO_AVG", "Cred_Real_Total", "IMAEP_yoy", "IPC_yoy")

# Check availability
available_ext <- intersect(external_vars, names(analysis_data))
available_dom <- intersect(domestic_vars, names(analysis_data))

cat("External block:", length(available_ext), "vars:", paste(available_ext, collapse = ", "), "\n")
cat("Domestic block:", length(available_dom), "vars:", paste(available_dom, collapse = ", "), "\n")

# Order: external first (block-exogenous), then domestic
var_order <- c(available_ext, available_dom)

var_data <- analysis_data %>%
  dplyr::select(fecha, all_of(var_order)) %>%
  na.omit()

cat("Sample:", nrow(var_data), "obs (",
    format(min(var_data$fecha)), "to", format(max(var_data$fecha)), ")\n")

# Estimate VAR
var_ts <- var_data[, var_order]
var_select <- VARselect(var_ts, lag.max = BSVAR_CONFIG$max_var_lags, type = "const")
lag_opt <- max(1, min(var_select$selection["AIC(n)"], BSVAR_CONFIG$max_var_lags))
cat("Optimal lag (AIC):", lag_opt, "\n\n")

var_model <- VAR(var_ts, p = lag_opt, type = "const")

# Structural identification via Cholesky (external first = block-exogenous)
K <- length(var_order)
n_ext <- length(available_ext)
n_dom <- length(available_dom)


################################################################################
# SECTION B: EXPANDED IRFS
################################################################################

cat("================================================================================\n")
cat("SECTION B: EXPANDED IRFs\n")
cat("================================================================================\n\n")

all_irf <- data.frame()

for (shock in available_ext) {
  for (resp in c(available_dom)) {
    tryCatch({
      irf_obj <- irf(var_model, impulse = shock, response = resp,
                      n.ahead = BSVAR_CONFIG$max_horizon, boot = TRUE,
                      ci = BSVAR_CONFIG$confidence_level, runs = BSVAR_CONFIG$n_bootstrap)

      irf_df <- data.frame(
        shock = shock,
        response = resp,
        horizon = 0:BSVAR_CONFIG$max_horizon,
        irf_val = as.vector(irf_obj$irf[[shock]]),
        lower = as.vector(irf_obj$Lower[[shock]]),
        upper = as.vector(irf_obj$Upper[[shock]]),
        stringsAsFactors = FALSE
      )
      all_irf <- rbind(all_irf, irf_df)
    }, error = function(e) {
      if (BSVAR_CONFIG$verbose) cat(sprintf("  IRF %s → %s failed: %s\n", shock, resp, e$message))
    })
  }
}

cat("IRFs computed:", nrow(all_irf), "rows\n")

# Peak effects summary
cat("\nPeak Effects on FCI_ENDO:\n")
cat(sprintf("%-15s %10s %10s %10s\n", "Shock", "Peak", "Horizon", "Significant"))
cat(strrep("-", 50), "\n")

fci_irfs <- all_irf %>% filter(response == "FCI_ENDO_AVG")
for (shock in available_ext) {
  sub <- fci_irfs %>% filter(shock == !!shock)
  if (nrow(sub) > 0) {
    peak_idx <- which.max(abs(sub$irf_val))
    peak_val <- sub$irf_val[peak_idx]
    peak_h <- sub$horizon[peak_idx]
    sig <- (sub$lower[peak_idx] > 0 & sub$upper[peak_idx] > 0) |
           (sub$lower[peak_idx] < 0 & sub$upper[peak_idx] < 0)
    cat(sprintf("%-15s %+10.4f %10d %10s\n", shock, peak_val, peak_h, ifelse(sig, "YES", "No")))
  }
}
cat("\n")


################################################################################
# SECTION C: ENRICHED VARIANCE DECOMPOSITION
################################################################################

cat("================================================================================\n")
cat("SECTION C: ENRICHED VARIANCE DECOMPOSITION\n")
cat("================================================================================\n\n")

fevd_result <- fevd(var_model, n.ahead = BSVAR_CONFIG$max_horizon)

fevd_all <- data.frame()
for (resp in names(fevd_result)) {
  fevd_mat <- as.data.frame(fevd_result[[resp]])
  fevd_mat$horizon <- 1:nrow(fevd_mat)
  fevd_mat$response <- resp
  fevd_long <- fevd_mat %>%
    pivot_longer(-c(horizon, response), names_to = "shock", values_to = "contribution")
  fevd_all <- rbind(fevd_all, fevd_long)
}

# Block decomposition for FCI
fevd_fci <- fevd_all %>%
  filter(response == "FCI_ENDO_AVG") %>%
  mutate(block = ifelse(shock %in% available_ext, "External", "Domestic"))

fevd_by_block <- fevd_fci %>%
  group_by(horizon, block) %>%
  summarise(contribution = sum(contribution), .groups = "drop")

cat("FCI Variance Decomposition by Block:\n")
for (h in c(1, 6, 12, 24)) {
  row <- fevd_by_block %>% filter(horizon == h)
  ext_share <- row %>% filter(block == "External") %>% pull(contribution) * 100
  dom_share <- row %>% filter(block == "Domestic") %>% pull(contribution) * 100
  if (length(ext_share) > 0 && length(dom_share) > 0) {
    cat(sprintf("  h=%2d: External %.1f%%, Domestic %.1f%%\n", h, ext_share, dom_share))
  }
}

# Individual shock contributions at h=12
cat("\nFCI Variance by Individual Shock (h=12):\n")
fevd_fci_12 <- fevd_fci %>% filter(horizon == 12) %>% arrange(desc(contribution))
for (i in seq_len(nrow(fevd_fci_12))) {
  cat(sprintf("  %-15s: %.1f%%\n", fevd_fci_12$shock[i], fevd_fci_12$contribution[i] * 100))
}
cat("\n")


################################################################################
# SECTION D: HISTORICAL DECOMPOSITION
################################################################################

cat("================================================================================\n")
cat("SECTION D: HISTORICAL DECOMPOSITION\n")
cat("================================================================================\n\n")

# Structural MA representation: Y_t = sum_{s=0}^{t} Theta_s * epsilon_{t-s}
# Where Theta_s = Phi_s * A0^{-1} (with Cholesky)

# Get structural residuals
u_hat <- residuals(var_model)
Sigma <- cov(u_hat)
A0_inv <- t(chol(Sigma))  # Lower triangular Cholesky factor
A0 <- solve(A0_inv)
eps_hat <- t(A0 %*% t(u_hat))  # Structural shocks
T_eff <- nrow(u_hat)

# MA coefficients
coef_list <- Acoef(var_model)
p <- var_model$p

# Compute MA coefficient matrices (Phi_s)
compute_Phi <- function(coef_list, K, p, max_s) {
  # Companion form
  A_comp <- matrix(0, K * p, K * p)
  for (i in 1:p) {
    A_comp[1:K, ((i-1)*K + 1):(i*K)] <- coef_list[[i]]
  }
  if (p > 1) {
    A_comp[(K+1):(K*p), 1:(K*(p-1))] <- diag(K*(p-1))
  }

  J <- cbind(diag(K), matrix(0, K, K*(p-1)))
  Phi_list <- list()
  Phi_list[[1]] <- diag(K)
  A_power <- diag(K * p)
  for (s in 1:max_s) {
    A_power <- A_power %*% A_comp
    Phi_list[[s + 1]] <- J %*% A_power %*% t(J)
  }
  return(Phi_list)
}

Phi_list <- compute_Phi(coef_list, K, p, T_eff)

# Historical decomposition: contribution of shock j to variable i at time t
# HD_{i,j,t} = sum_{s=0}^{t-1} (Phi_s * A0_inv)[i,j] * eps_hat[t-s, j]

cat("Computing historical decomposition...\n")
hd_fci <- matrix(0, T_eff, K)  # Contribution of each shock to FCI_ENDO_AVG
fci_idx <- which(var_order == "FCI_ENDO_AVG")

for (t in 1:T_eff) {
  for (s in 0:min(t - 1, T_eff - 1)) {
    Theta_s <- Phi_list[[s + 1]] %*% A0_inv  # Structural MA coefficient
    for (j in 1:K) {
      hd_fci[t, j] <- hd_fci[t, j] + Theta_s[fci_idx, j] * eps_hat[t - s, j]
    }
  }
}

colnames(hd_fci) <- var_order
hd_dates <- var_data$fecha[(p + 1):nrow(var_data)]

hd_df <- as.data.frame(hd_fci)
hd_df$fecha <- hd_dates

# Actual FCI (demeaned, relative to start)
hd_df$FCI_actual <- var_data$FCI_ENDO_AVG[(p + 1):nrow(var_data)]

# Group by block
hd_long <- hd_df %>%
  pivot_longer(-c(fecha, FCI_actual), names_to = "shock", values_to = "contribution")

hd_by_block <- hd_long %>%
  mutate(block = ifelse(shock %in% available_ext, "External", "Domestic")) %>%
  group_by(fecha, FCI_actual, block) %>%
  summarise(contribution = sum(contribution), .groups = "drop")

cat("Historical decomposition computed for", T_eff, "periods\n\n")


################################################################################
# SECTION E: COMPARISON WITH ORIGINAL 3-VAR BLOCK-SVAR
################################################################################

cat("================================================================================\n")
cat("SECTION E: COMPARISON WITH ORIGINAL 3-VAR BLOCK-SVAR\n")
cat("================================================================================\n\n")

# Load original FEVD if available
orig_fevd_file <- file.path(BSVAR_CONFIG$output_dir, "Block_SVAR_FEVD.csv")
comparison_df <- data.frame()

if (file.exists(orig_fevd_file)) {
  orig_fevd <- read.csv(orig_fevd_file)
  cat("Original Block-SVAR FEVD loaded.\n")

  # Extract external share at key horizons
  # Original FEVD stores contributions as percentages (0-100 scale)
  orig_ext <- orig_fevd %>%
    filter(grepl("FCI", response, ignore.case = TRUE) & block == "External") %>%
    group_by(horizon) %>%
    summarise(external_share_pct = sum(contribution, na.rm = TRUE), .groups = "drop")

  # New (enriched) FEVD: contributions are fractions (0-1 scale) — convert to %
  new_ext <- fevd_by_block %>%
    filter(block == "External") %>%
    mutate(external_share_pct = contribution * 100) %>%
    dplyr::select(horizon, external_share_pct)

  for (h in c(1, 6, 12, 24)) {
    orig_val <- orig_ext %>% filter(horizon == h) %>% pull(external_share_pct)
    new_val <- new_ext %>% filter(horizon == h) %>% pull(external_share_pct)
    orig_val <- if (length(orig_val) > 0) orig_val[1] else NA
    new_val <- if (length(new_val) > 0) new_val[1] else NA

    # Sanity check: FEVD percentages must be in [0, 100]
    if (!is.na(orig_val) && (orig_val < 0 | orig_val > 100)) {
      warning(sprintf("FEVD sanity check failed at h=%d: original=%.1f%% (expected 0-100)", h, orig_val))
    }

    comparison_df <- rbind(comparison_df, data.frame(
      horizon = h,
      original_3var_ext_pct = orig_val,
      enriched_7var_ext_pct = new_val,
      stringsAsFactors = FALSE
    ))
  }

  cat("\nFEVD Comparison (External Share for FCI):\n")
  cat(sprintf("%-8s %15s %15s\n", "Horizon", "Original (3-var)", "Enriched (7-var)"))
  cat(strrep("-", 42), "\n")
  for (i in seq_len(nrow(comparison_df))) {
    cat(sprintf("h=%2d    %14.1f%% %14.1f%%\n",
                comparison_df$horizon[i],
                comparison_df$original_3var_ext_pct[i],
                comparison_df$enriched_7var_ext_pct[i]))
  }
} else {
  cat("Original Block-SVAR FEVD not found. Skipping comparison.\n")
}
cat("\n")


################################################################################
# SAVE OUTPUTS
################################################################################

cat("================================================================================\n")
cat("SAVING OUTPUTS\n")
cat("================================================================================\n\n")

write.csv(all_irf, file.path(BSVAR_CONFIG$output_dir, "Enriched_BSVAR_IRF.csv"),
          row.names = FALSE)
cat("Saved: Enriched_BSVAR_IRF.csv\n")

write.csv(fevd_all, file.path(BSVAR_CONFIG$output_dir, "Enriched_BSVAR_FEVD.csv"),
          row.names = FALSE)
cat("Saved: Enriched_BSVAR_FEVD.csv\n")

# Historical decomposition
hd_export <- hd_df %>%
  dplyr::select(fecha, everything())
write.csv(hd_export, file.path(BSVAR_CONFIG$output_dir, "Enriched_BSVAR_Historical_Decomposition.csv"),
          row.names = FALSE)
cat("Saved: Enriched_BSVAR_Historical_Decomposition.csv\n")

if (nrow(comparison_df) > 0) {
  write.csv(comparison_df, file.path(BSVAR_CONFIG$output_dir, "Enriched_BSVAR_Comparison.csv"),
            row.names = FALSE)
  cat("Saved: Enriched_BSVAR_Comparison.csv\n")
} else {
  # Save enriched FEVD summary as comparison
  write.csv(fevd_by_block, file.path(BSVAR_CONFIG$output_dir, "Enriched_BSVAR_Comparison.csv"),
            row.names = FALSE)
  cat("Saved: Enriched_BSVAR_Comparison.csv (block shares only)\n")
}


################################################################################
# VISUALIZATIONS
################################################################################

cat("\n--- Generating Visualizations ---\n\n")

ext_colors <- c("VIX" = "#E41A1C", "SP500" = "#377EB8", "DXY" = "#4DAF4A",
                "US_10Y" = "#984EA3", "FFER" = "#FF7F00", "Selic_rate" = "#A65628",
                "d_ToT" = "#F781BF")

# ---- 216: IRFs from external shocks → FCI ----
tryCatch({
  fci_irfs <- all_irf %>% filter(response == "FCI_ENDO_AVG")

  if (nrow(fci_irfs) > 0) {
    p <- ggplot(fci_irfs, aes(x = horizon, y = irf_val)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "steelblue", linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~shock, scales = "free_y", ncol = 4) +
      labs(title = "Enriched Block-SVAR: External Shocks → FCI_ENDO",
           subtitle = "7-variable external block | 90% bootstrap CI",
           x = "Horizon (months)", y = "FCI response") +
      theme_minimal(base_size = 10)

    ggsave(file.path(BSVAR_CONFIG$output_dir, "216_Enriched_BSVAR_IRF_to_FCI.png"),
           p, width = 16, height = 8, dpi = 150)
    cat("Saved: 216_Enriched_BSVAR_IRF_to_FCI.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 216 failed:", e$message, "\n"))

# ---- 217: IRFs from external shocks → Credit ----
tryCatch({
  credit_irfs <- all_irf %>% filter(response == "Cred_Real_Total")

  if (nrow(credit_irfs) > 0) {
    p <- ggplot(credit_irfs, aes(x = horizon, y = irf_val)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "darkgreen") +
      geom_line(color = "darkgreen", linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~shock, scales = "free_y", ncol = 4) +
      labs(title = "Enriched Block-SVAR: External Shocks → Real Credit",
           x = "Horizon (months)", y = "Credit response (pp)") +
      theme_minimal(base_size = 10)

    ggsave(file.path(BSVAR_CONFIG$output_dir, "217_Enriched_BSVAR_IRF_to_Credit.png"),
           p, width = 16, height = 8, dpi = 150)
    cat("Saved: 217_Enriched_BSVAR_IRF_to_Credit.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 217 failed:", e$message, "\n"))

# ---- 218: FEVD stacked bars ----
tryCatch({
  fevd_fci_key <- fevd_fci %>%
    filter(horizon %in% c(1, 6, 12, 24))

  if (nrow(fevd_fci_key) > 0) {
    p <- ggplot(fevd_fci_key, aes(x = factor(horizon), y = contribution * 100, fill = shock)) +
      geom_col(position = "stack") +
      labs(title = "FCI Variance Decomposition by External Shock",
           subtitle = "Enriched 7-variable external block",
           x = "Horizon (months)", y = "Variance share (%)", fill = "Shock") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "right")

    ggsave(file.path(BSVAR_CONFIG$output_dir, "218_Enriched_BSVAR_FEVD.png"),
           p, width = 12, height = 6, dpi = 150)
    cat("Saved: 218_Enriched_BSVAR_FEVD.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 218 failed:", e$message, "\n"))

# ---- 219: FEVD Comparison (3-var vs 7-var) ----
tryCatch({
  if (nrow(comparison_df) > 0) {
    comp_long <- comparison_df %>%
      pivot_longer(-horizon, names_to = "specification", values_to = "external_pct")

    p <- ggplot(comp_long, aes(x = factor(horizon), y = external_pct, fill = specification)) +
      geom_col(position = "dodge") +
      labs(title = "FEVD Comparison: External Share for FCI",
           subtitle = "Original 3-var vs Enriched 7-var external block",
           x = "Horizon (months)", y = "External share (%)") +
      scale_fill_manual(values = c("original_3var_ext_pct" = "gray60",
                                    "enriched_7var_ext_pct" = "steelblue"),
                        labels = c("Original (3-var)", "Enriched (7-var)")) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")

    ggsave(file.path(BSVAR_CONFIG$output_dir, "219_Enriched_BSVAR_FEVD_Comparison.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 219_Enriched_BSVAR_FEVD_Comparison.png\n")
  } else {
    # Plot only enriched FEVD
    p <- ggplot(fevd_by_block %>% filter(horizon %in% c(1, 6, 12, 24)),
                aes(x = factor(horizon), y = contribution * 100, fill = block)) +
      geom_col(position = "stack") +
      labs(title = "FCI Variance Decomposition by Block",
           x = "Horizon", y = "Share (%)") +
      theme_minimal(base_size = 11)
    ggsave(file.path(BSVAR_CONFIG$output_dir, "219_Enriched_BSVAR_FEVD_Comparison.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 219_Enriched_BSVAR_FEVD_Comparison.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 219 failed:", e$message, "\n"))

# ---- 220: Historical Decomposition (full sample) ----
tryCatch({
  p <- ggplot(hd_by_block, aes(x = fecha, y = contribution, fill = block)) +
    geom_area(alpha = 0.6) +
    geom_line(aes(y = FCI_actual, fill = NULL), color = "black", linewidth = 0.6) +
    labs(title = "Historical Decomposition of FCI_ENDO",
         subtitle = "Contribution of external vs domestic shocks | Black line = actual FCI",
         x = NULL, y = "FCI contribution") +
    scale_fill_manual(values = c("External" = "steelblue", "Domestic" = "firebrick")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(BSVAR_CONFIG$output_dir, "220_Historical_Decomposition_FCI.png"),
         p, width = 14, height = 6, dpi = 150)
  cat("Saved: 220_Historical_Decomposition_FCI.png\n")
}, error = function(e) cat("  WARNING: Plot 220 failed:", e$message, "\n"))

# ---- 221: Historical Decomposition - Key Episodes ----
tryCatch({
  episodes <- list(
    "Asian/Russian Crisis\n(1997-1998)" = c(as.Date("1997-01-01"), as.Date("1999-06-01")),
    "Regional Contagion\n(2002-2003)" = c(as.Date("2001-06-01"), as.Date("2004-01-01")),
    "Global Financial Crisis\n(2008-2009)" = c(as.Date("2007-06-01"), as.Date("2010-01-01")),
    "COVID-19\n(2020-2021)" = c(as.Date("2019-06-01"), as.Date("2022-01-01"))
  )

  plots <- list()
  for (ep_name in names(episodes)) {
    ep_range <- episodes[[ep_name]]
    ep_data <- hd_by_block %>%
      filter(fecha >= ep_range[1] & fecha <= ep_range[2])

    if (nrow(ep_data) > 0) {
      plots[[ep_name]] <- ggplot(ep_data, aes(x = fecha, y = contribution, fill = block)) +
        geom_area(alpha = 0.6) +
        geom_line(aes(y = FCI_actual, fill = NULL), color = "black", linewidth = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        labs(title = ep_name, x = NULL, y = "FCI") +
        scale_fill_manual(values = c("External" = "steelblue", "Domestic" = "firebrick")) +
        theme_minimal(base_size = 9) +
        theme(legend.position = "none",
              plot.title = element_text(size = 9, face = "bold"))
    }
  }

  if (length(plots) > 0) {
    combined <- do.call(grid.arrange, c(plots, ncol = 2))
    ggsave(file.path(BSVAR_CONFIG$output_dir, "221_Historical_Decomposition_Episodes.png"),
           combined, width = 14, height = 8, dpi = 150)
    cat("Saved: 221_Historical_Decomposition_Episodes.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 221 failed:", e$message, "\n"))

# ---- 222: Dashboard ----
tryCatch({
  plots <- list()

  # Peak effects bar
  peak_data <- fci_irfs %>%
    group_by(shock) %>%
    summarise(peak_effect = irf_val[which.max(abs(irf_val))],
              peak_h = horizon[which.max(abs(irf_val))],
              .groups = "drop")
  plots[[1]] <- ggplot(peak_data, aes(x = reorder(shock, abs(peak_effect)),
                                       y = peak_effect, fill = peak_effect > 0)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    labs(title = "Peak FCI Effect by Shock", x = NULL, y = "Peak IRF") +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue")) +
    theme_minimal(base_size = 9)

  # FEVD at h=12
  fevd12 <- fevd_fci %>% filter(horizon == 12) %>% arrange(desc(contribution))
  plots[[2]] <- ggplot(fevd12, aes(x = reorder(shock, contribution), y = contribution * 100)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "FCI FEVD at h=12", x = NULL, y = "Share (%)") +
    theme_minimal(base_size = 9)

  # Block shares over time
  plots[[3]] <- ggplot(fevd_by_block, aes(x = horizon, y = contribution * 100, color = block)) +
    geom_line(linewidth = 0.8) +
    labs(title = "Block Shares Over Horizon", x = "Horizon", y = "Share (%)") +
    scale_color_manual(values = c("External" = "steelblue", "Domestic" = "firebrick")) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom")

  # HD summary for key dates
  plots[[4]] <- ggplot(hd_by_block, aes(x = fecha, y = contribution, fill = block)) +
    geom_area(alpha = 0.5) +
    labs(title = "Historical Decomposition", x = NULL, y = "FCI") +
    scale_fill_manual(values = c("External" = "steelblue", "Domestic" = "firebrick")) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom")

  combined <- do.call(grid.arrange, c(plots, ncol = 2))
  ggsave(file.path(BSVAR_CONFIG$output_dir, "222_Enriched_BSVAR_Dashboard.png"),
         combined, width = 14, height = 10, dpi = 150)
  cat("Saved: 222_Enriched_BSVAR_Dashboard.png\n")
}, error = function(e) cat("  WARNING: Plot 222 failed:", e$message, "\n"))


cat("\n################################################################################\n")
cat("SCRIPT 20 COMPLETE\n")
cat("################################################################################\n\n")
cat("Outputs: 7 PNGs (216-222) + 4 CSVs\n\n")
