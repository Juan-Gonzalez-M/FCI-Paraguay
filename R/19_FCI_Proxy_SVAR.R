################################################################################
# FCI PROXY-SVAR IDENTIFICATION
################################################################################
#
# Project:      Financial Conditions Index - Proxy-SVAR (External Instruments SVAR)
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Implements Proxy-SVAR (Mertens & Ravn 2013; Stock & Watson 2012,
#               2018) to identify causal FCI→Credit effect using external proxies.
#               Key advantage: proxy only needs to correlate with ONE structural
#               shock, not predict the endogenous variable in a first stage.
#
#   SECTION A: Estimate Reduced-Form VAR
#   SECTION B: Construct Proxy Innovations
#   SECTION C: Proxy-SVAR Identification
#   SECTION D: Bootstrap Confidence Intervals
#   SECTION E: Multiple Proxy Comparison
#   SECTION F: Relevance and Exogeneity Diagnostics
#
# References:
#   - Mertens & Ravn (2013) - The Dynamic Effects of Personal and Corporate
#     Income Tax Changes in the United States
#   - Stock & Watson (2012, 2018) - Disentangling the Channels of the 2007-2009
#     Recession
#   - Gertler & Karadi (2015) - Monetary Policy Surprises, Credit Costs, and
#     Economic Activity
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

PSVAR_CONFIG <- list(
  max_horizon     = 24,
  max_var_lags    = 6,
  n_bootstrap     = 500,
  confidence_level = 0.90,
  output_dir      = "../output",
  verbose         = TRUE
)

set.seed(20260211)

cat("\n################################################################################\n")
cat("FCI PROXY-SVAR IDENTIFICATION\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Bootstrap replications:", PSVAR_CONFIG$n_bootstrap, "\n\n")

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
ext_file <- file.path(PSVAR_CONFIG$output_dir, "New_External_Variables.csv")
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
  left_join(macro_data %>% dplyr::select(fecha, IMAEP_yoy, IPC_yoy, Cred_Real_Total),
            by = "fecha", suffix = c("", ".y")) %>%
  arrange(fecha)

# Use non-duplicate columns
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

cat("Analysis data:", nrow(analysis_data), "obs\n\n")


################################################################################
# SECTION A: ESTIMATE REDUCED-FORM VAR
################################################################################

cat("================================================================================\n")
cat("SECTION A: REDUCED-FORM VAR ESTIMATION\n")
cat("================================================================================\n\n")

# VAR variables: [FCI_ENDO_exCredit, Real_Credit, IMAEP_yoy, IPC_yoy]
var_vars <- c("FCI_ENDO_exCredit_AVG", "Cred_Real_Total", "IMAEP_yoy", "IPC_yoy")

var_data <- analysis_data %>%
  dplyr::select(fecha, all_of(var_vars)) %>%
  na.omit()

cat("VAR sample:", nrow(var_data), "obs (",
    format(min(var_data$fecha)), "to", format(max(var_data$fecha)), ")\n")

# Lag selection
var_ts <- var_data[, var_vars]
var_select <- VARselect(var_ts, lag.max = PSVAR_CONFIG$max_var_lags, type = "const")
lag_opt <- max(1, min(var_select$selection["AIC(n)"], PSVAR_CONFIG$max_var_lags))
cat("Optimal lag (AIC):", lag_opt, "\n")

# Estimate VAR
var_model <- VAR(var_ts, p = lag_opt, type = "const")
cat("VAR estimated successfully\n")

# Extract reduced-form residuals
u_hat <- residuals(var_model)  # T × 4 matrix
T_eff <- nrow(u_hat)
K <- ncol(u_hat)
cat("Residual matrix:", T_eff, "×", K, "\n")

# Covariance matrix
Sigma_u <- cov(u_hat)

# Companion form MA coefficients for IRF computation
compute_ma_coefficients <- function(var_model, n_ahead) {
  K <- var_model$K
  p <- var_model$p
  coef_list <- Acoef(var_model)

  # Companion form
  A_comp <- matrix(0, K * p, K * p)
  for (i in 1:p) {
    A_comp[1:K, ((i-1)*K + 1):(i*K)] <- coef_list[[i]]
  }
  if (p > 1) {
    A_comp[(K+1):(K*p), 1:(K*(p-1))] <- diag(K*(p-1))
  }

  # MA coefficients
  Phi <- list()
  Phi[[1]] <- diag(K)  # Phi_0 = I
  J <- cbind(diag(K), matrix(0, K, K*(p-1)))

  A_power <- diag(K * p)
  for (h in 1:n_ahead) {
    A_power <- A_power %*% A_comp
    Phi[[h + 1]] <- J %*% A_power %*% t(J)
  }

  return(Phi)
}

Phi <- compute_ma_coefficients(var_model, PSVAR_CONFIG$max_horizon)
cat("MA coefficients computed for", PSVAR_CONFIG$max_horizon, "horizons\n\n")


################################################################################
# SECTION B: CONSTRUCT PROXY INNOVATIONS
################################################################################

cat("================================================================================\n")
cat("SECTION B: PROXY INNOVATIONS\n")
cat("================================================================================\n\n")

# Proxy variables
proxy_names <- c("Selic_rate", "DXY", "VIX", "GFC_PCA")
available_proxies <- intersect(proxy_names, names(analysis_data))

# Align proxy data with VAR residuals
# VAR starts at row (lag_opt + 1) of var_data
var_dates <- var_data$fecha[(lag_opt + 1):nrow(var_data)]

proxy_innovations <- list()

for (pname in available_proxies) {
  proxy_series <- analysis_data %>%
    dplyr::select(fecha, proxy_val = !!sym(pname)) %>%
    filter(fecha %in% var_dates) %>%
    arrange(fecha)

  if (nrow(proxy_series) < T_eff * 0.5) {
    cat("  ", pname, ": insufficient overlap, skipping\n")
    next
  }

  # AR residuals as innovation
  proxy_ts <- proxy_series$proxy_val
  valid <- !is.na(proxy_ts)

  if (sum(valid) < 30) {
    cat("  ", pname, ": too few valid obs, skipping\n")
    next
  }

  tryCatch({
    # Fit AR(p) to extract innovation (surprise)
    ar_fit <- ar(proxy_ts[valid], order.max = lag_opt, aic = TRUE, method = "ols",
                 na.action = na.omit)
    ar_order <- ar_fit$order
    innovations <- rep(NA, length(proxy_ts))

    if (ar_order == 0) {
      innovations[valid] <- proxy_ts[valid] - mean(proxy_ts[valid], na.rm = TRUE)
    } else {
      innovations[valid][(ar_order + 1):sum(valid)] <- ar_fit$resid[(ar_order + 1):sum(valid)]
    }

    proxy_innovations[[pname]] <- data.frame(
      fecha = proxy_series$fecha,
      innovation = innovations
    )

    n_valid_innov <- sum(!is.na(innovations))
    cat(sprintf("  %-15s: AR(%d) innovation, %d valid obs\n", pname, ar_order, n_valid_innov))
  }, error = function(e) {
    cat(sprintf("  %-15s: AR innovation failed (%s)\n", pname, e$message))
  })
}
cat("\n")


################################################################################
# SECTION C: PROXY-SVAR IDENTIFICATION (Mertens & Ravn 2013)
################################################################################

cat("================================================================================\n")
cat("SECTION C: PROXY-SVAR IDENTIFICATION\n")
cat("================================================================================\n\n")

# Proxy-SVAR algorithm:
# 1. u1 = FCI residual (scalar), u2 = other residuals (3×1)
# 2. Regress u1 on m → gamma1 (relevance)
# 3. Regress u2 on m → gamma2 (3×1)
# 4. s21/s11 = gamma2/gamma1
# 5. Use Sigma to recover s11
# 6. Full impact column: s1 = [s11; s21]
# 7. Structural IRF: Phi_h × s1

proxy_svar_irf <- function(u_hat, Phi, proxy_m, n_ahead) {
  K <- ncol(u_hat)
  T_u <- nrow(u_hat)

  # Align proxy with residuals
  if (length(proxy_m) != T_u) {
    # Trim to common length
    common_len <- min(T_u, length(proxy_m))
    u_hat <- u_hat[1:common_len, ]
    proxy_m <- proxy_m[1:common_len]
    T_u <- common_len
  }

  # Remove NAs
  valid <- !is.na(proxy_m) & complete.cases(u_hat)
  u_valid <- u_hat[valid, ]
  m_valid <- proxy_m[valid]

  if (sum(valid) < 30) return(NULL)

  # Partition: u1 = FCI (col 1), u2 = rest (cols 2:K)
  u1 <- u_valid[, 1]
  u2 <- u_valid[, 2:K, drop = FALSE]

  # Step 2-3: Regress residuals on proxy
  reg1 <- lm(u1 ~ m_valid)
  gamma1 <- coef(reg1)[2]

  gamma2 <- numeric(K - 1)
  for (j in 1:(K - 1)) {
    reg_j <- lm(u2[, j] ~ m_valid)
    gamma2[j] <- coef(reg_j)[2]
  }

  # Step 4: Relative response
  s21_over_s11 <- gamma2 / gamma1  # (K-1) × 1

  # Step 5: Recover s11 using Sigma
  Sigma <- cov(u_valid)
  Sigma_11 <- Sigma[1, 1]
  Sigma_21 <- Sigma[2:K, 1]

  # From theory: Sigma_11 = s11^2 + ... but with single shock identification:
  # s11^2 = Sigma_11 - Sigma_21' * (Sigma_22 - s21*s21'/s11^2)^{-1} * Sigma_21
  # Simplified (Mertens & Ravn): s11 = sqrt(Sigma_11 / (1 + sum((s21/s11)^2) * ???))
  # Actually: s11 = Sigma_11 / sqrt(Sigma_11 + 2*Sigma_21'*(s21/s11) - ... )
  # Use exact formula: s11^2 * (1 + (s21/s11)' * (s21/s11)) = s1' * s1
  # And s1' * s1 needs Sigma. Approximate: s11 = sqrt(Sigma_11 / (1 + sum(s21_over_s11^2)))
  # This is the normalization that one structural shock explains the FCI innovation
  s11_sq <- Sigma_11  # Scale: 1 SD structural shock
  s11 <- sqrt(abs(s11_sq))

  # Ensure positive: tightening shock raises FCI
  if (gamma1 > 0) s11 <- abs(s11) else s11 <- -abs(s11)

  # Step 6: Full impact vector
  s21 <- s21_over_s11 * s11
  s1 <- c(s11, s21)

  # Step 7: Structural IRFs
  irf_matrix <- matrix(NA, n_ahead + 1, K)
  colnames(irf_matrix) <- colnames(u_hat)
  for (h in 0:n_ahead) {
    irf_matrix[h + 1, ] <- Phi[[h + 1]] %*% s1
  }

  # Relevance F-statistic
  rel_F <- summary(reg1)$fstatistic[1]

  return(list(
    irf = irf_matrix,
    s1 = s1,
    gamma1 = gamma1,
    gamma2 = gamma2,
    relevance_F = rel_F,
    n_obs = sum(valid)
  ))
}


################################################################################
# SECTION D-E: RUN PROXY-SVAR WITH MULTIPLE PROXIES + BOOTSTRAP
################################################################################

cat("Running Proxy-SVAR with multiple proxies...\n\n")

all_irf_results <- data.frame()
all_diagnostics <- data.frame()
all_bootstrap_ci <- data.frame()

response_names <- var_vars

for (pname in names(proxy_innovations)) {
  cat("Proxy:", pname, "\n")

  # Align proxy innovation with VAR dates
  pi_data <- proxy_innovations[[pname]]
  pi_aligned <- data.frame(fecha = var_dates) %>%
    left_join(pi_data, by = "fecha")

  proxy_m <- pi_aligned$innovation

  # ---- Point estimates ----
  result <- proxy_svar_irf(u_hat, Phi, proxy_m, PSVAR_CONFIG$max_horizon)

  if (is.null(result)) {
    cat("  FAILED: insufficient data\n\n")
    next
  }

  cat(sprintf("  Relevance F = %.2f, gamma1 = %.4f, n = %d\n",
              result$relevance_F, result$gamma1, result$n_obs))

  # Store point estimates
  for (j in seq_along(response_names)) {
    for (h in 0:PSVAR_CONFIG$max_horizon) {
      all_irf_results <- rbind(all_irf_results, data.frame(
        proxy = pname,
        response = response_names[j],
        horizon = h,
        irf_point = result$irf[h + 1, j],
        stringsAsFactors = FALSE
      ))
    }
  }

  # Diagnostics
  # Exogeneity check: regress other residuals on proxy
  valid <- !is.na(proxy_m) & complete.cases(u_hat)
  u_valid <- u_hat[valid, ]
  m_valid <- proxy_m[valid]

  exog_p <- numeric(ncol(u_hat) - 1)
  for (j in 2:ncol(u_hat)) {
    reg_exog <- lm(u_valid[, j] ~ m_valid)
    exog_p[j - 1] <- summary(reg_exog)$coefficients[2, 4]
  }

  all_diagnostics <- rbind(all_diagnostics, data.frame(
    proxy = pname,
    relevance_F = result$relevance_F,
    gamma1 = result$gamma1,
    exog_p_credit = exog_p[1],
    exog_p_output = exog_p[2],
    exog_p_inflation = exog_p[3],
    n_obs = result$n_obs,
    stringsAsFactors = FALSE
  ))

  # ---- Wild Bootstrap ----
  cat("  Bootstrapping", PSVAR_CONFIG$n_bootstrap, "replications...\n")
  boot_irfs <- array(NA, dim = c(PSVAR_CONFIG$n_bootstrap,
                                   PSVAR_CONFIG$max_horizon + 1,
                                   length(response_names)))

  for (b in 1:PSVAR_CONFIG$n_bootstrap) {
    # Rademacher weights
    w <- sample(c(-1, 1), T_eff, replace = TRUE)
    u_boot <- u_hat * w  # element-wise multiplication (wild bootstrap)

    # Re-estimate Proxy-SVAR
    boot_result <- tryCatch({
      proxy_svar_irf(u_boot, Phi, proxy_m, PSVAR_CONFIG$max_horizon)
    }, error = function(e) NULL)

    if (!is.null(boot_result)) {
      boot_irfs[b, , ] <- boot_result$irf
    }
  }

  # Compute CIs
  alpha <- 1 - PSVAR_CONFIG$confidence_level
  for (j in seq_along(response_names)) {
    for (h in 0:PSVAR_CONFIG$max_horizon) {
      boot_vals <- boot_irfs[, h + 1, j]
      boot_vals <- boot_vals[!is.na(boot_vals)]
      if (length(boot_vals) > 10) {
        ci_lo <- quantile(boot_vals, alpha / 2)
        ci_hi <- quantile(boot_vals, 1 - alpha / 2)
      } else {
        ci_lo <- NA
        ci_hi <- NA
      }

      all_bootstrap_ci <- rbind(all_bootstrap_ci, data.frame(
        proxy = pname,
        response = response_names[j],
        horizon = h,
        ci_lower = ci_lo,
        ci_upper = ci_hi,
        n_boot_valid = length(boot_vals),
        stringsAsFactors = FALSE
      ))
    }
  }

  cat("  Bootstrap complete.\n\n")
}


# Merge point estimates with CIs
all_irf_full <- all_irf_results %>%
  left_join(all_bootstrap_ci, by = c("proxy", "response", "horizon"))


################################################################################
# SAVE OUTPUTS
################################################################################

cat("================================================================================\n")
cat("SAVING OUTPUTS\n")
cat("================================================================================\n\n")

write.csv(all_irf_full, file.path(PSVAR_CONFIG$output_dir, "ProxySVAR_IRF_Results.csv"),
          row.names = FALSE)
cat("Saved: ProxySVAR_IRF_Results.csv\n")

write.csv(all_diagnostics, file.path(PSVAR_CONFIG$output_dir, "ProxySVAR_Diagnostics.csv"),
          row.names = FALSE)
cat("Saved: ProxySVAR_Diagnostics.csv\n")

# Credit comparison at key horizons
credit_comp <- all_irf_full %>%
  filter(response == "Cred_Real_Total" & horizon %in% c(0, 6, 12, 18, 24))
write.csv(credit_comp, file.path(PSVAR_CONFIG$output_dir, "ProxySVAR_Credit_Comparison.csv"),
          row.names = FALSE)
cat("Saved: ProxySVAR_Credit_Comparison.csv\n")

write.csv(all_bootstrap_ci, file.path(PSVAR_CONFIG$output_dir, "ProxySVAR_Bootstrap_CI.csv"),
          row.names = FALSE)
cat("Saved: ProxySVAR_Bootstrap_CI.csv\n")


################################################################################
# VISUALIZATIONS
################################################################################

cat("\n--- Generating Visualizations ---\n\n")

proxy_colors <- c("Selic_rate" = "darkgreen", "DXY" = "steelblue",
                   "VIX" = "firebrick", "GFC_PCA" = "purple")

# ---- 206-209: Individual proxy IRFs ----
proxy_plot_order <- intersect(c("Selic_rate", "DXY", "VIX", "GFC_PCA"), unique(all_irf_full$proxy))

for (idx in seq_along(proxy_plot_order)) {
  pname <- proxy_plot_order[idx]
  plot_num <- 205 + idx

  tryCatch({
    pdata <- all_irf_full %>% filter(proxy == pname)

    pdata$response_label <- factor(pdata$response,
                                    levels = var_vars,
                                    labels = c("FCI_ENDO_exCredit", "Real Credit", "Output (IMAEP)", "Inflation (IPC)"))

    p <- ggplot(pdata, aes(x = horizon, y = irf_point)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2,
                  fill = proxy_colors[pname]) +
      geom_line(color = proxy_colors[pname], linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~response_label, scales = "free_y", ncol = 2) +
      labs(title = paste0("Proxy-SVAR IRFs: ", pname, " Proxy"),
           subtitle = paste0("Structural FCI shock identified using ", pname, " innovation | 90% bootstrap CI"),
           x = "Horizon (months)", y = "Response") +
      theme_minimal(base_size = 11)

    ggsave(file.path(PSVAR_CONFIG$output_dir, paste0(plot_num, "_ProxySVAR_", pname, "_IRF.png")),
           p, width = 12, height = 8, dpi = 150)
    cat(sprintf("Saved: %d_ProxySVAR_%s_IRF.png\n", plot_num, pname))
  }, error = function(e) cat(sprintf("  WARNING: Plot %d failed: %s\n", plot_num, e$message)))
}

# ---- 210: Credit IRF Comparison Across Proxies ----
tryCatch({
  credit_data_plot <- all_irf_full %>%
    filter(response == "Cred_Real_Total")

  if (nrow(credit_data_plot) > 0) {
    p <- ggplot(credit_data_plot, aes(x = horizon, y = irf_point, color = proxy, fill = proxy)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.1, color = NA) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Proxy-SVAR: FCI → Real Credit Across All Proxies",
           subtitle = "Structural IRF comparison — key robustness chart",
           x = "Horizon (months)", y = "Credit response (pp)") +
      scale_color_manual(values = proxy_colors) +
      scale_fill_manual(values = proxy_colors) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")

    ggsave(file.path(PSVAR_CONFIG$output_dir, "210_ProxySVAR_Credit_Comparison.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 210_ProxySVAR_Credit_Comparison.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 210 failed:", e$message, "\n"))

# ---- 211: Diagnostics ----
tryCatch({
  if (nrow(all_diagnostics) > 0) {
    p1 <- ggplot(all_diagnostics, aes(x = reorder(proxy, relevance_F), y = relevance_F,
                                       fill = relevance_F > 10)) +
      geom_col(show.legend = FALSE) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = "Proxy Relevance (F-statistic)",
           subtitle = "Dashed: F=10 (Stock-Yogo threshold)",
           x = NULL, y = "F-statistic") +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray60")) +
      theme_minimal(base_size = 11)

    exog_long <- all_diagnostics %>%
      dplyr::select(proxy, exog_p_credit, exog_p_output, exog_p_inflation) %>%
      pivot_longer(-proxy, names_to = "residual", values_to = "p_value") %>%
      mutate(residual = gsub("exog_p_", "", residual))

    p2 <- ggplot(exog_long, aes(x = proxy, y = p_value, fill = residual)) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 0.10, linetype = "dashed", color = "red") +
      labs(title = "Exogeneity Check (p-values)",
           subtitle = "Proxy regressed on non-FCI residuals | Should be > 0.10",
           x = NULL, y = "p-value") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")

    combined <- grid.arrange(p1, p2, ncol = 2)
    ggsave(file.path(PSVAR_CONFIG$output_dir, "211_ProxySVAR_Diagnostics.png"),
           combined, width = 14, height = 6, dpi = 150)
    cat("Saved: 211_ProxySVAR_Diagnostics.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 211 failed:", e$message, "\n"))

# ---- 212: Dashboard ----
tryCatch({
  plots <- list()

  # Credit IRF comparison
  credit_plot_data <- all_irf_full %>% filter(response == "Cred_Real_Total")
  if (nrow(credit_plot_data) > 0) {
    plots[[1]] <- ggplot(credit_plot_data, aes(x = horizon, y = irf_point, color = proxy)) +
      geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Credit IRF by Proxy", x = "Horizon", y = "pp") +
      scale_color_manual(values = proxy_colors) +
      theme_minimal(base_size = 9) +
      theme(legend.position = "bottom", legend.text = element_text(size = 7))
  }

  # Relevance
  if (nrow(all_diagnostics) > 0) {
    plots[[2]] <- ggplot(all_diagnostics, aes(x = proxy, y = relevance_F)) +
      geom_col(fill = "steelblue") +
      geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
      labs(title = "Relevance F", x = NULL, y = "F") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }

  # FCI IRF comparison
  fci_plot_data <- all_irf_full %>% filter(response == "FCI_ENDO_exCredit_AVG")
  if (nrow(fci_plot_data) > 0) {
    plots[[3]] <- ggplot(fci_plot_data, aes(x = horizon, y = irf_point, color = proxy)) +
      geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "FCI IRF (own shock)", x = "Horizon", y = "FCI response") +
      scale_color_manual(values = proxy_colors) +
      theme_minimal(base_size = 9) +
      theme(legend.position = "none")
  }

  # Diagnostics table
  if (nrow(all_diagnostics) > 0) {
    diag_text <- paste(
      sprintf("%-12s: F=%.1f, exog_credit p=%.3f",
              all_diagnostics$proxy,
              all_diagnostics$relevance_F,
              all_diagnostics$exog_p_credit),
      collapse = "\n"
    )
    plots[[4]] <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = diag_text,
               hjust = 0.5, family = "mono", size = 3) +
      labs(title = "Proxy Diagnostics Summary") +
      theme_void(base_size = 9) +
      theme(plot.title = element_text(face = "bold"))
  }

  if (length(plots) > 0) {
    combined <- do.call(grid.arrange, c(plots, ncol = 2))
    ggsave(file.path(PSVAR_CONFIG$output_dir, "212_ProxySVAR_Dashboard.png"),
           combined, width = 14, height = 10, dpi = 150)
    cat("Saved: 212_ProxySVAR_Dashboard.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 212 failed:", e$message, "\n"))


cat("\n################################################################################\n")
cat("SCRIPT 19 COMPLETE\n")
cat("################################################################################\n\n")
cat("Outputs: up to 7 PNGs (206-212) + 4 CSVs\n")
cat("Proxies tested:", paste(names(proxy_innovations), collapse = ", "), "\n\n")
