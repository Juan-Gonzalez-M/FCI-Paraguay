################################################################################
# FCI QUARTERLY LP WITH PUBLISHED BCP DATA
################################################################################
#
# Project:      Financial Conditions Index - Published Data Validation
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Replaces Section 5.4 and Appendix F with publicly available
#               quarterly national accounts data from FCI_data_1.xlsx.
#
#   PART A: Data Preparation (quarterly aggregation)
#   PART B: Standard Quarterly LP -- Supply Side
#   PART C: Standard Quarterly LP -- Demand Side
#   PART D: Endo vs Exo Decomposition
#   PART E: Asymmetric LP (Tightening vs Easing)
#   PART F: Post-IT Subsample
#   PART G: TOST Equivalence Test for Agriculture
#   PART H: VECM with Quarterly Data
#   PART I: Comparison Dashboard
#
# Data:
#   - FCI_data_1.xlsx: Quarterly_SA (quarterly national accounts)
#   - FCI_data_1.xlsx: FCI variables + macro controls
#
# References:
#   - Jordà (2005) AER - Local Projections
#   - Johansen (1988, 1991) - Cointegration Tests
#   - Hatzius et al. (2010) - FCI motivation
#
################################################################################


################################################################################
# 1. SETUP AND CONFIGURATION
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(sandwich)
  library(lmtest)
  library(urca)
  library(vars)
  library(zoo)
  library(lubridate)
})

CONFIG <- list(
  # Data files
  published_file = "../data/FCI_data_1.xlsx",
  fci_file       = "../data/FCI_data_1.xlsx",
  macro_sheet    = "Datos_macro",

  # LP parameters
  max_horizon_q = 8,        # 8 quarters = 24 months
  horizons_report = c(2, 4, 6, 8),  # key reporting horizons (quarters)
  n_lags = 2,
  confidence_level = 0.90,
  min_obs = 25,             # minimum observations for quarterly LP

  # VECM parameters
  vecm_lag_max = 6,

  # Output
  output_dir = "../output"
)

# Supply-side sectors (from Quarterly_SA)
SUPPLY_SECTORS <- c(
  "Agricultura",
  "Ganaderia_fp",
  "Manufactura",
  "Electricidad_y_agua",
  "Construccion",
  "Servicios",
  "PIB",
  "PIB_exAgri"
)

SUPPLY_LABELS <- c(
  "Agricultura"         = "Agriculture",
  "Ganaderia_fp"        = "Livestock+",
  "Manufactura"         = "Manufacturing",
  "Electricidad_y_agua" = "Electricity & Water",
  "Construccion"        = "Construction",
  "Servicios"           = "Services",
  "PIB"                 = "GDP",
  "PIB_exAgri"          = "GDP ex-Agriculture"
)

# Demand-side components (from Quarterly_SA)
DEMAND_SECTORS <- c(
  "Consumo_Privado",
  "Consumo_Publico",
  "FBCF",
  "Exportaciones",
  "Importaciones"
)

DEMAND_LABELS <- c(
  "Consumo_Privado" = "Private Consumption",
  "Consumo_Publico" = "Public Consumption",
  "FBCF"            = "Investment (GFCF)",
  "Exportaciones"   = "Exports",
  "Importaciones"   = "Imports"
)

ALL_LABELS <- c(SUPPLY_LABELS, DEMAND_LABELS)

set.seed(20260320)

cat("\n################################################################################\n")
cat("FCI QUARTERLY LP WITH PUBLISHED BCP DATA\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# 2. LP FUNCTIONS
################################################################################

#' Standard Local Projection (quarterly)
run_lp_quarterly <- function(data, y_var, fci_var, max_h, n_lags = 2,
                             control_vars = NULL, conf_level = 0.90,
                             min_obs = 25) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - conf_level) / 2)

  for (h in 1:max_h) {
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        fci_lag1 = lag(!!sym(fci_var), 1)
      )

    for (i in 2:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", i) := lag(!!sym(y_var), i))
    }

    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")

    control_cols <- c()
    if (!is.null(control_vars)) {
      for (cv in control_vars) {
        if (cv %in% names(data)) {
          data_h <- data_h %>%
            mutate(!!cv := !!sym(cv))
          control_cols <- c(control_cols, cv)
          for (j in 1:n_lags) {
            cv_lag_name <- paste0(cv, "_lag", j)
            data_h <- data_h %>%
              mutate(!!cv_lag_name := lag(!!sym(cv), j))
            control_cols <- c(control_cols, cv_lag_name)
          }
        }
      }
    }

    formula_str <- paste("y_fwd ~", fci_var, "+",
                         paste(c(lag_vars, control_cols), collapse = " + "))

    reg_data <- data_h %>%
      dplyr::select(y_fwd, !!sym(fci_var), all_of(lag_vars),
                    all_of(control_cols)) %>%
      na.omit()

    if (nrow(reg_data) < min_obs) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    idx <- which(rownames(coef_test) == fci_var)

    results <- rbind(results, data.frame(
      horizon = h,
      coef = coef_test[idx, 1],
      se = coef_test[idx, 2],
      t_stat = coef_test[idx, 3],
      p_value = coef_test[idx, 4],
      ci_lower = coef_test[idx, 1] - z_crit * coef_test[idx, 2],
      ci_upper = coef_test[idx, 1] + z_crit * coef_test[idx, 2],
      n_obs = nrow(reg_data),
      r_squared = summary(model)$r.squared
    ))
  }

  return(results)
}

#' Asymmetric (Non-linear) Local Projection (quarterly)
run_lp_asymmetric_q <- function(data, y_var, fci_var, max_h, n_lags = 2,
                                control_vars = NULL, conf_level = 0.90,
                                min_obs = 25) {

  results <- data.frame()
  z_crit <- qnorm(1 - (1 - conf_level) / 2)

  fci_pos <- paste0(fci_var, "_pos")
  fci_neg <- paste0(fci_var, "_neg")

  for (h in 1:max_h) {
    data_h <- data %>%
      mutate(
        y_fwd = lead(!!sym(y_var), h),
        y_lag1 = lag(!!sym(y_var), 1),
        fci_lag1 = lag(!!sym(fci_var), 1)
      )

    for (i in 2:n_lags) {
      data_h <- data_h %>%
        mutate(!!paste0("y_lag", i) := lag(!!sym(y_var), i))
    }

    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")

    control_cols <- c()
    if (!is.null(control_vars)) {
      for (cv in control_vars) {
        if (cv %in% names(data)) {
          data_h <- data_h %>%
            mutate(!!cv := !!sym(cv))
          control_cols <- c(control_cols, cv)
          for (j in 1:n_lags) {
            cv_lag_name <- paste0(cv, "_lag", j)
            data_h <- data_h %>%
              mutate(!!cv_lag_name := lag(!!sym(cv), j))
            control_cols <- c(control_cols, cv_lag_name)
          }
        }
      }
    }

    formula_str <- paste("y_fwd ~", fci_pos, "+", fci_neg, "+",
                         paste(c(lag_vars, control_cols), collapse = " + "))

    reg_data <- data_h %>%
      dplyr::select(y_fwd, !!sym(fci_pos), !!sym(fci_neg),
                    all_of(lag_vars), all_of(control_cols)) %>%
      na.omit()

    if (nrow(reg_data) < min_obs) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    idx_pos <- which(rownames(coef_test) == fci_pos)
    idx_neg <- which(rownames(coef_test) == fci_neg)

    coef_diff <- coef_test[idx_pos, 1] - coef_test[idx_neg, 1]
    var_diff <- vcov_hac[idx_pos, idx_pos] + vcov_hac[idx_neg, idx_neg] -
                2 * vcov_hac[idx_pos, idx_neg]
    wald_p <- 1 - pchisq(coef_diff^2 / var_diff, df = 1)

    results <- rbind(results, data.frame(
      horizon = h,
      coef_tight = coef_test[idx_pos, 1],
      se_tight = coef_test[idx_pos, 2],
      p_tight = coef_test[idx_pos, 4],
      ci_tight_lo = coef_test[idx_pos, 1] - z_crit * coef_test[idx_pos, 2],
      ci_tight_hi = coef_test[idx_pos, 1] + z_crit * coef_test[idx_pos, 2],
      coef_ease = coef_test[idx_neg, 1],
      se_ease = coef_test[idx_neg, 2],
      p_ease = coef_test[idx_neg, 4],
      ci_ease_lo = coef_test[idx_neg, 1] - z_crit * coef_test[idx_neg, 2],
      ci_ease_hi = coef_test[idx_neg, 1] + z_crit * coef_test[idx_neg, 2],
      asym_p = wald_p,
      n_obs = nrow(reg_data)
    ))
  }

  return(results)
}


################################################################################
# 3. PART A: DATA PREPARATION
################################################################################

cat("================================================================================\n")
cat("PART A: DATA PREPARATION\n")
cat("================================================================================\n\n")

if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)

# --- Load FCI ---
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  CONFIG_27 <- CONFIG  # Save our CONFIG before source overwrites it
  source("01_FCI_Complete.R")
  CONFIG <- CONFIG_27  # Restore our CONFIG
}

fci_monthly <- resultado_fci$all_indices %>%
  dplyr::select(fecha,
                FCI_COMP = FCI_COMP_AVG,
                FCI_ENDO = FCI_ENDO_AVG,
                FCI_EXO = FCI_EXO_AVG,
                FCI_exCredit = FCI_exCredit_AVG,
                FCI_ENDO_exCredit = FCI_ENDO_exCredit_AVG)

# --- Load macro data (monthly controls: IPC, credit) ---
cat("Loading macro data for controls...\n")
macro_raw <- read_excel(CONFIG$fci_file, sheet = CONFIG$macro_sheet)
fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

macro_monthly <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  mutate(
    IPC_yoy   = (IPC / lag(IPC, 12) - 1) * 100,
    Cred_Real_yoy = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100
  ) %>%
  dplyr::select(fecha, IPC_yoy, Cred_Real_yoy)

# --- Aggregate monthly FCI to quarterly (end-of-quarter: last month) ---
cat("Aggregating monthly FCI to quarterly (end-of-quarter)...\n")

fci_quarterly <- fci_monthly %>%
  mutate(
    year = year(fecha),
    quarter = quarter(fecha),
    month_in_q = month(fecha) %% 3
  ) %>%
  # End-of-quarter = months 3, 6, 9, 12 (i.e., month_in_q == 0)
  filter(month_in_q == 0) %>%
  mutate(
    # Create quarterly date as first day of the quarter-ending month
    q_date = fecha
  ) %>%
  dplyr::select(q_date, FCI_COMP, FCI_ENDO, FCI_EXO, FCI_exCredit, FCI_ENDO_exCredit)

cat(sprintf("  FCI quarterly: %d observations, %s to %s\n",
            nrow(fci_quarterly), format(min(fci_quarterly$q_date)),
            format(max(fci_quarterly$q_date))))

# --- Aggregate monthly controls to quarterly (end-of-quarter) ---
controls_quarterly <- macro_monthly %>%
  mutate(month_in_q = month(fecha) %% 3) %>%
  filter(month_in_q == 0) %>%
  mutate(q_date = fecha) %>%
  dplyr::select(q_date, IPC_yoy, Cred_Real_yoy)

# --- Load quarterly national accounts ---
cat("Loading FCI_data_1.xlsx (Quarterly_SA)...\n")

quarterly_raw <- read_excel(CONFIG$published_file, sheet = "Quarterly_SA")
names(quarterly_raw)[1] <- "fecha"

quarterly_data <- quarterly_raw %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Rename columns to safe names
col_mapping <- c(
  "Agricultura"         = "Agricultura",
  "Ganadería forestal,  pesca y minería" = "Ganaderia_fp",
  "Manufactura"         = "Manufactura",
  "Electricidad y agua" = "Electricidad_y_agua",
  "Construcción"        = "Construccion",
  "Servicios"           = "Servicios",
  "PIB"                 = "PIB",
  "Consumo Privado"     = "Consumo_Privado",
  "Consumo Público"     = "Consumo_Publico",
  "Formación bruta de capital fijo" = "FBCF",
  "Exportaciones"       = "Exportaciones",
  "Importaciones"       = "Importaciones"
)

for (old_name in names(col_mapping)) {
  if (old_name %in% names(quarterly_data)) {
    names(quarterly_data)[names(quarterly_data) == old_name] <- col_mapping[old_name]
  }
}

# Compute PIB_exAgri = PIB - Agricultura
quarterly_data <- quarterly_data %>%
  mutate(PIB_exAgri = PIB - Agricultura)

# Compute YoY growth rates (Y_t / Y_{t-4} - 1) * 100
all_sectors <- c(SUPPLY_SECTORS, DEMAND_SECTORS)

for (s in all_sectors) {
  if (s %in% names(quarterly_data)) {
    yoy_name <- paste0(s, "_yoy")
    quarterly_data <- quarterly_data %>%
      mutate(!!yoy_name := (!!sym(s) / lag(!!sym(s), 4) - 1) * 100)
  }
}

# Use end-of-quarter dates for merging
quarterly_data <- quarterly_data %>%
  mutate(q_date = fecha)

cat(sprintf("  Quarterly national accounts: %d observations, %s to %s\n",
            nrow(quarterly_data), format(min(quarterly_data$fecha)),
            format(max(quarterly_data$fecha))))
cat(sprintf("  Available supply sectors: %s\n",
            paste(intersect(SUPPLY_SECTORS, names(quarterly_data)), collapse = ", ")))
cat(sprintf("  Available demand components: %s\n",
            paste(intersect(DEMAND_SECTORS, names(quarterly_data)), collapse = ", ")))

# --- Merge all quarterly data ---
analysis_q <- fci_quarterly %>%
  inner_join(quarterly_data, by = "q_date") %>%
  left_join(controls_quarterly, by = "q_date") %>%
  mutate(
    FCI_COMP_pos = pmax(FCI_COMP, 0),
    FCI_COMP_neg = abs(pmin(FCI_COMP, 0)),
    FCI_exCredit_pos = pmax(FCI_exCredit, 0),
    FCI_exCredit_neg = abs(pmin(FCI_exCredit, 0)),
    FCI_ENDO_pos = pmax(FCI_ENDO, 0),
    FCI_ENDO_neg = abs(pmin(FCI_ENDO, 0)),
    FCI_EXO_pos = pmax(FCI_EXO, 0),
    FCI_EXO_neg = abs(pmin(FCI_EXO, 0)),
    FCI_ENDO_exCredit_pos = pmax(FCI_ENDO_exCredit, 0),
    FCI_ENDO_exCredit_neg = abs(pmin(FCI_ENDO_exCredit, 0))
  ) %>%
  arrange(q_date)

cat(sprintf("\nMerged quarterly dataset: %d observations, %s to %s\n",
            nrow(analysis_q), format(min(analysis_q$q_date)),
            format(max(analysis_q$q_date))))

# Report available YoY growth columns
avail_yoy <- names(analysis_q)[grepl("_yoy$", names(analysis_q))]
cat(sprintf("  YoY growth variables: %d\n", length(avail_yoy)))
for (v in avail_yoy) {
  n_ok <- sum(!is.na(analysis_q[[v]]))
  cat(sprintf("    %s: %d non-NA\n", v, n_ok))
}


################################################################################
# 4. PART B: STANDARD QUARTERLY LP -- SUPPLY SIDE
################################################################################

cat("\n================================================================================\n")
cat("PART B: SUPPLY-SIDE QUARTERLY LOCAL PROJECTIONS\n")
cat("================================================================================\n\n")

supply_results <- list()

for (sector in SUPPLY_SECTORS) {
  y_var <- paste0(sector, "_yoy")
  if (!y_var %in% names(analysis_q)) {
    cat(sprintf("  SKIP: %s not found\n", sector))
    next
  }

  data_sub <- analysis_q %>% filter(!is.na(!!sym(y_var)))
  cat(sprintf("  %s (n=%d) ... ", SUPPLY_LABELS[sector], nrow(data_sub)))

  res <- run_lp_quarterly(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon_q, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy"),
    min_obs = CONFIG$min_obs
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- SUPPLY_LABELS[sector]
    res$side <- "Supply"
    supply_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_supply <- bind_rows(supply_results)
cat(sprintf("\n  Supply results: %d sector-horizon pairs\n\n", nrow(all_supply)))

# Summary table
cat("SUPPLY-SIDE SUMMARY (quarterly):\n")
cat(sprintf("%-25s %10s %8s  %10s %8s\n", "Sector", "h=2Q", "p", "h=4Q", "p"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (sector in SUPPLY_SECTORS) {
  v2 <- all_supply %>% filter(sector == !!sector, horizon == 2)
  v4 <- all_supply %>% filter(sector == !!sector, horizon == 4)

  fmt <- function(v) {
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      sprintf("%+.2f%s", v$coef, stars)
    } else "NA"
  }
  p_fmt <- function(v) if (nrow(v) > 0) sprintf("%.3f", v$p_value) else "NA"

  cat(sprintf("%-25s %10s %8s  %10s %8s\n",
              SUPPLY_LABELS[sector], fmt(v2), p_fmt(v2), fmt(v4), p_fmt(v4)))
}

# Plot 260: Supply-side IRFs
if (nrow(all_supply) > 0) {
  p260 <- ggplot(all_supply, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#2C3E50", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#2C3E50", linewidth = 1) +
    geom_point(aes(y = coef), color = "#2C3E50", size = 2) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    theme_minimal(base_size = 11) +
    labs(title = "Sectoral Output Response to FCI Tightening (Quarterly, Full Sample)",
         subtitle = "FCI_COMP \u2192 Output growth (YoY) | Published BCP data | 90% CI",
         x = "Quarters", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(CONFIG$output_dir, "260_PubQ_LP_Supply_Full.png"), p260,
         width = 14, height = 8, dpi = 300)
  cat("Saved: 260_PubQ_LP_Supply_Full.png\n")
}

write.csv(all_supply, file.path(CONFIG$output_dir, "csv", "PubQ_LP_Supply_Full.csv"),
          row.names = FALSE)
cat("Saved: PubQ_LP_Supply_Full.csv\n")


################################################################################
# 5. PART C: STANDARD QUARTERLY LP -- DEMAND SIDE
################################################################################

cat("\n================================================================================\n")
cat("PART C: DEMAND-SIDE QUARTERLY LOCAL PROJECTIONS\n")
cat("================================================================================\n\n")

demand_results <- list()

for (sector in DEMAND_SECTORS) {
  y_var <- paste0(sector, "_yoy")
  if (!y_var %in% names(analysis_q)) {
    cat(sprintf("  SKIP: %s not found\n", sector))
    next
  }

  data_sub <- analysis_q %>% filter(!is.na(!!sym(y_var)))
  cat(sprintf("  %s (n=%d) ... ", DEMAND_LABELS[sector], nrow(data_sub)))

  res <- run_lp_quarterly(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon_q, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy"),
    min_obs = CONFIG$min_obs
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- DEMAND_LABELS[sector]
    res$side <- "Demand"
    demand_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_demand <- bind_rows(demand_results)
cat(sprintf("\n  Demand results: %d sector-horizon pairs\n\n", nrow(all_demand)))

# Summary table
cat("DEMAND-SIDE SUMMARY (quarterly):\n")
cat(sprintf("%-25s %10s %8s  %10s %8s\n", "Sector", "h=2Q", "p", "h=4Q", "p"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (sector in DEMAND_SECTORS) {
  v2 <- all_demand %>% filter(sector == !!sector, horizon == 2)
  v4 <- all_demand %>% filter(sector == !!sector, horizon == 4)

  fmt <- function(v) {
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      sprintf("%+.2f%s", v$coef, stars)
    } else "NA"
  }
  p_fmt <- function(v) if (nrow(v) > 0) sprintf("%.3f", v$p_value) else "NA"

  cat(sprintf("%-25s %10s %8s  %10s %8s\n",
              DEMAND_LABELS[sector], fmt(v2), p_fmt(v2), fmt(v4), p_fmt(v4)))
}

# Plot 262: Demand-side IRFs
if (nrow(all_demand) > 0) {
  p262 <- ggplot(all_demand, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#8E44AD", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#8E44AD", linewidth = 1) +
    geom_point(aes(y = coef), color = "#8E44AD", size = 2) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    theme_minimal(base_size = 11) +
    labs(title = "Demand-Side Response to FCI Tightening (Quarterly, Full Sample)",
         subtitle = "FCI_COMP \u2192 Component growth (YoY) | Published BCP data | 90% CI",
         x = "Quarters", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(CONFIG$output_dir, "262_PubQ_LP_Demand_Full.png"), p262,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 262_PubQ_LP_Demand_Full.png\n")
}

write.csv(all_demand, file.path(CONFIG$output_dir, "csv", "PubQ_LP_Demand_Full.csv"),
          row.names = FALSE)
cat("Saved: PubQ_LP_Demand_Full.csv\n")

# Combined heatmap at h=2Q and h=4Q
all_both <- bind_rows(all_supply, all_demand)

heatmap_q <- all_both %>%
  filter(horizon %in% c(2, 4)) %>%
  mutate(
    horizon_label = paste0("h=", horizon, "Q"),
    label = sprintf("%.1f%s", coef,
                    ifelse(p_value < 0.01, "***",
                           ifelse(p_value < 0.05, "**",
                                  ifelse(p_value < 0.10, "*", ""))))
  )

if (nrow(heatmap_q) > 0) {
  p263 <- ggplot(heatmap_q, aes(x = horizon_label, y = sector_label, fill = coef)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "Effect (pp)") +
    facet_wrap(~side, scales = "free_y", ncol = 2) +
    theme_minimal(base_size = 12) +
    labs(title = "Quarterly FCI Effect on Output Components",
         subtitle = "Effect in pp | *** p<0.01, ** p<0.05, * p<0.10",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size = 10))

  ggsave(file.path(CONFIG$output_dir, "263_PubQ_LP_Heatmap.png"), p263,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 263_PubQ_LP_Heatmap.png\n")
}


################################################################################
# 6. PART D: ENDO vs EXO DECOMPOSITION
################################################################################

cat("\n================================================================================\n")
cat("PART D: ENDO vs EXO DECOMPOSITION\n")
cat("================================================================================\n\n")

fci_output_types <- c("FCI_ENDO" = "Endogenous (Domestic)",
                       "FCI_EXO" = "Exogenous (External)")
colors_endo_exo <- c("Endogenous (Domestic)" = "#E74C3C",
                      "Exogenous (External)" = "#3498DB")

endoexo_results <- list()

for (sector in c(SUPPLY_SECTORS, DEMAND_SECTORS)) {
  y_var <- paste0(sector, "_yoy")
  if (!y_var %in% names(analysis_q)) next

  data_sub <- analysis_q %>% filter(!is.na(!!sym(y_var)))

  for (fci_var in names(fci_output_types)) {
    cat(sprintf("  %s x %s (n=%d) ... ",
                ALL_LABELS[sector], fci_output_types[fci_var], nrow(data_sub)))

    res <- run_lp_quarterly(
      data_sub, y_var, fci_var, CONFIG$max_horizon_q, CONFIG$n_lags,
      control_vars = c("Cred_Real_yoy", "IPC_yoy"),
      min_obs = CONFIG$min_obs
    )

    if (nrow(res) > 0) {
      res$sector <- sector
      res$sector_label <- ALL_LABELS[sector]
      res$side <- ifelse(sector %in% SUPPLY_SECTORS, "Supply", "Demand")
      res$fci_var <- fci_var
      res$fci_label <- fci_output_types[fci_var]
      endoexo_results[[paste(sector, fci_var)]] <- res
      cat("done\n")
    } else {
      cat("insufficient data\n")
    }
  }
}

all_endoexo <- bind_rows(endoexo_results)

# Plot 264: Endo vs Exo -- Supply
endoexo_supply <- all_endoexo %>% filter(side == "Supply")
if (nrow(endoexo_supply) > 0) {
  p264 <- ggplot(endoexo_supply,
                 aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    scale_color_manual(values = colors_endo_exo) +
    scale_fill_manual(values = colors_endo_exo) +
    theme_minimal(base_size = 10) +
    labs(title = "Supply Side: Endogenous vs Exogenous FCI (Quarterly)",
         subtitle = "FCI_ENDO vs FCI_EXO \u2192 Output growth (YoY) | 90% CI",
         x = "Quarters", y = "Effect (pp)",
         color = "FCI Component", fill = "FCI Component") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "264_PubQ_LP_EndoExo_Supply.png"), p264,
         width = 14, height = 8, dpi = 300)
  cat("Saved: 264_PubQ_LP_EndoExo_Supply.png\n")
}

# Plot 265: Endo vs Exo -- Demand
endoexo_demand <- all_endoexo %>% filter(side == "Demand")
if (nrow(endoexo_demand) > 0) {
  p265 <- ggplot(endoexo_demand,
                 aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    scale_color_manual(values = colors_endo_exo) +
    scale_fill_manual(values = colors_endo_exo) +
    theme_minimal(base_size = 10) +
    labs(title = "Demand Side: Endogenous vs Exogenous FCI (Quarterly)",
         subtitle = "FCI_ENDO vs FCI_EXO \u2192 Component growth (YoY) | 90% CI",
         x = "Quarters", y = "Effect (pp)",
         color = "FCI Component", fill = "FCI Component") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "265_PubQ_LP_EndoExo_Demand.png"), p265,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 265_PubQ_LP_EndoExo_Demand.png\n")
}

write.csv(all_endoexo %>% filter(side == "Supply"),
          file.path(CONFIG$output_dir, "csv", "PubQ_LP_EndoExo_Supply.csv"), row.names = FALSE)
write.csv(all_endoexo %>% filter(side == "Demand"),
          file.path(CONFIG$output_dir, "csv", "PubQ_LP_EndoExo_Demand.csv"), row.names = FALSE)
cat("Saved: PubQ_LP_EndoExo_Supply.csv, PubQ_LP_EndoExo_Demand.csv\n")


################################################################################
# 7. PART E: ASYMMETRIC LP (TIGHTENING vs EASING)
################################################################################

cat("\n================================================================================\n")
cat("PART E: ASYMMETRIC LP (Tightening vs Easing)\n")
cat("================================================================================\n\n")

colors_asym <- c("Tightening" = "#E74C3C", "Easing" = "#3498DB")

# --- Supply-side asymmetric ---
cat("--- Supply-Side Asymmetric ---\n\n")
asym_supply_results <- list()

for (sector in SUPPLY_SECTORS) {
  y_var <- paste0(sector, "_yoy")
  if (!y_var %in% names(analysis_q)) next

  data_sub <- analysis_q %>% filter(!is.na(!!sym(y_var)))
  cat(sprintf("  %s (n=%d) ... ", SUPPLY_LABELS[sector], nrow(data_sub)))

  res <- run_lp_asymmetric_q(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon_q, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy"),
    min_obs = CONFIG$min_obs
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- SUPPLY_LABELS[sector]
    res$side <- "Supply"
    asym_supply_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_asym_supply <- bind_rows(asym_supply_results)

# --- Demand-side asymmetric ---
cat("\n--- Demand-Side Asymmetric ---\n\n")
asym_demand_results <- list()

for (sector in DEMAND_SECTORS) {
  y_var <- paste0(sector, "_yoy")
  if (!y_var %in% names(analysis_q)) next

  data_sub <- analysis_q %>% filter(!is.na(!!sym(y_var)))
  cat(sprintf("  %s (n=%d) ... ", DEMAND_LABELS[sector], nrow(data_sub)))

  res <- run_lp_asymmetric_q(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon_q, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy"),
    min_obs = CONFIG$min_obs
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- DEMAND_LABELS[sector]
    res$side <- "Demand"
    asym_demand_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_asym_demand <- bind_rows(asym_demand_results)
all_asym <- bind_rows(all_asym_supply, all_asym_demand)

# Plot 266: Asymmetric Supply
if (nrow(all_asym_supply) > 0) {
  asym_s_plot <- all_asym_supply %>%
    pivot_longer(cols = c(coef_tight, coef_ease),
                 names_to = "effect_type", values_to = "coef") %>%
    mutate(
      effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
      ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
      ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
    )

  p266 <- ggplot(asym_s_plot, aes(x = horizon, y = coef,
                                    color = effect_label, fill = effect_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    scale_color_manual(values = colors_asym) +
    scale_fill_manual(values = colors_asym) +
    theme_minimal(base_size = 10) +
    labs(title = "Asymmetric Supply-Side Response (Quarterly)",
         subtitle = "Tightening (FCI>0) vs Easing (|FCI<0|) | FCI_COMP | 90% CI",
         x = "Quarters", y = "Effect (pp)",
         color = "FCI State", fill = "FCI State") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "266_PubQ_LP_Asym_Supply.png"), p266,
         width = 14, height = 8, dpi = 300)
  cat("Saved: 266_PubQ_LP_Asym_Supply.png\n")
}

# Plot 267: Asymmetric Demand
if (nrow(all_asym_demand) > 0) {
  asym_d_plot <- all_asym_demand %>%
    pivot_longer(cols = c(coef_tight, coef_ease),
                 names_to = "effect_type", values_to = "coef") %>%
    mutate(
      effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
      ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
      ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
    )

  p267 <- ggplot(asym_d_plot, aes(x = horizon, y = coef,
                                    color = effect_label, fill = effect_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    scale_color_manual(values = colors_asym) +
    scale_fill_manual(values = colors_asym) +
    theme_minimal(base_size = 10) +
    labs(title = "Asymmetric Demand-Side Response (Quarterly)",
         subtitle = "Tightening (FCI>0) vs Easing (|FCI<0|) | FCI_COMP | 90% CI",
         x = "Quarters", y = "Effect (pp)",
         color = "FCI State", fill = "FCI State") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "267_PubQ_LP_Asym_Demand.png"), p267,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 267_PubQ_LP_Asym_Demand.png\n")
}

write.csv(all_asym,
          file.path(CONFIG$output_dir, "csv", "PubQ_LP_Asymmetric.csv"), row.names = FALSE)
cat("Saved: PubQ_LP_Asymmetric.csv\n")


################################################################################
# 8. PART F: POST-IT SUBSAMPLE
################################################################################

cat("\n================================================================================\n")
cat("PART F: POST-IT SUBSAMPLE (Q2-2011+)\n")
cat("================================================================================\n\n")

# IT adoption: May 2011 -> Q2-2011 starts June 2011
analysis_postIT <- analysis_q %>%
  filter(q_date >= as.Date("2011-06-01"))

cat(sprintf("Post-IT sample: %d observations, %s to %s\n",
            nrow(analysis_postIT), format(min(analysis_postIT$q_date)),
            format(max(analysis_postIT$q_date))))

# Reduce max horizon for smaller sample
max_h_postIT <- min(CONFIG$max_horizon_q, 6)  # 6Q = 18M max

# --- Post-IT Supply ---
cat("\n--- Post-IT Supply-Side ---\n\n")
postIT_supply_results <- list()

for (sector in SUPPLY_SECTORS) {
  y_var <- paste0(sector, "_yoy")
  if (!y_var %in% names(analysis_postIT)) next

  data_sub <- analysis_postIT %>% filter(!is.na(!!sym(y_var)))
  cat(sprintf("  %s (n=%d) ... ", SUPPLY_LABELS[sector], nrow(data_sub)))

  res <- run_lp_quarterly(
    data_sub, y_var, "FCI_COMP", max_h_postIT, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy"),
    min_obs = 20  # relaxed for smaller sample
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- SUPPLY_LABELS[sector]
    res$side <- "Supply"
    postIT_supply_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_postIT_supply <- bind_rows(postIT_supply_results)

# --- Post-IT Demand ---
cat("\n--- Post-IT Demand-Side ---\n\n")
postIT_demand_results <- list()

for (sector in DEMAND_SECTORS) {
  y_var <- paste0(sector, "_yoy")
  if (!y_var %in% names(analysis_postIT)) next

  data_sub <- analysis_postIT %>% filter(!is.na(!!sym(y_var)))
  cat(sprintf("  %s (n=%d) ... ", DEMAND_LABELS[sector], nrow(data_sub)))

  res <- run_lp_quarterly(
    data_sub, y_var, "FCI_COMP", max_h_postIT, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy"),
    min_obs = 20
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- DEMAND_LABELS[sector]
    res$side <- "Demand"
    postIT_demand_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_postIT_demand <- bind_rows(postIT_demand_results)
all_postIT <- bind_rows(all_postIT_supply, all_postIT_demand)

# Plot 268: Post-IT Supply
if (nrow(all_postIT_supply) > 0) {
  p268 <- ggplot(all_postIT_supply, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#27AE60", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#27AE60", linewidth = 1) +
    geom_point(aes(y = coef), color = "#27AE60", size = 2) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_x_continuous(breaks = 1:max_h_postIT) +
    theme_minimal(base_size = 10) +
    labs(title = "Supply-Side Response (Post-IT, Quarterly)",
         subtitle = paste0("FCI_COMP \u2192 Output growth (YoY) | Q2-2011+ (n=",
                           nrow(analysis_postIT), ") | 90% CI"),
         x = "Quarters", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(CONFIG$output_dir, "268_PubQ_LP_PostIT_Supply.png"), p268,
         width = 14, height = 8, dpi = 300)
  cat("Saved: 268_PubQ_LP_PostIT_Supply.png\n")
}

# Plot 269: Post-IT Demand
if (nrow(all_postIT_demand) > 0) {
  p269 <- ggplot(all_postIT_demand, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#F39C12", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#F39C12", linewidth = 1) +
    geom_point(aes(y = coef), color = "#F39C12", size = 2) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = 1:max_h_postIT) +
    theme_minimal(base_size = 10) +
    labs(title = "Demand-Side Response (Post-IT, Quarterly)",
         subtitle = paste0("FCI_COMP \u2192 Component growth (YoY) | Q2-2011+ (n=",
                           nrow(analysis_postIT), ") | 90% CI"),
         x = "Quarters", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(CONFIG$output_dir, "269_PubQ_LP_PostIT_Demand.png"), p269,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 269_PubQ_LP_PostIT_Demand.png\n")
}

write.csv(all_postIT,
          file.path(CONFIG$output_dir, "csv", "PubQ_LP_PostIT.csv"), row.names = FALSE)
cat("Saved: PubQ_LP_PostIT.csv\n")


################################################################################
# 9. PART G: TOST EQUIVALENCE TEST FOR AGRICULTURE
################################################################################

cat("\n================================================================================\n")
cat("PART G: TOST EQUIVALENCE TEST FOR AGRICULTURE\n")
cat("================================================================================\n\n")

tost_results <- data.frame()
tost_bound <- 3  # +/- 3 pp equivalence bound (quarterly)

# Test at h=4Q (approx 12M)
agri_4q <- all_supply %>% filter(sector == "Agricultura", horizon == 4)

if (nrow(agri_4q) > 0) {
  beta_hat <- agri_4q$coef
  se_hat <- agri_4q$se
  n_obs <- agri_4q$n_obs

  # TOST: Two one-sided tests
  # H0: |beta| >= bound  vs  H1: |beta| < bound
  t_upper <- (beta_hat - tost_bound) / se_hat
  t_lower <- (beta_hat + tost_bound) / se_hat

  # One-sided p-values
  p_upper <- pt(t_upper, df = n_obs - 1)             # P(beta < +bound)
  p_lower <- 1 - pt(t_lower, df = n_obs - 1)         # P(beta > -bound)

  tost_p <- max(p_upper, p_lower)  # TOST p-value

  # Power analysis
  ncp <- tost_bound / se_hat  # non-centrality parameter under exact boundary
  power <- pt(qt(0.05, df = n_obs - 1), df = n_obs - 1, ncp = ncp)

  tost_results <- data.frame(
    sector = "Agricultura",
    horizon = 4,
    beta = beta_hat,
    se = se_hat,
    n_obs = n_obs,
    tost_bound = tost_bound,
    t_upper = t_upper,
    t_lower = t_lower,
    p_upper = p_upper,
    p_lower = p_lower,
    tost_p = tost_p,
    equivalence = ifelse(tost_p < 0.05, "EQUIVALENT", "INCONCLUSIVE"),
    power = power
  )

  cat("TOST Results for Agriculture at h=4Q:\n")
  cat(sprintf("  Beta = %+.3f (SE = %.3f)\n", beta_hat, se_hat))
  cat(sprintf("  Equivalence bound: +/- %.1f pp\n", tost_bound))
  cat(sprintf("  TOST p-value: %.4f\n", tost_p))
  cat(sprintf("  Conclusion: %s\n", tost_results$equivalence))
  cat(sprintf("  Power (at boundary): %.3f\n\n", power))

  cat("NOTE: Quarterly frequency reduces sample size relative to monthly analysis.\n")
  cat("  Power may be limited with ~", n_obs, " usable observations.\n\n")
} else {
  cat("  Agriculture results not available at h=4Q.\n\n")
}

write.csv(tost_results,
          file.path(CONFIG$output_dir, "csv", "PubQ_TOST_Agriculture.csv"), row.names = FALSE)
cat("Saved: PubQ_TOST_Agriculture.csv\n")


################################################################################
# 10. PART H: VECM WITH QUARTERLY DATA
################################################################################

cat("\n================================================================================\n")
cat("PART H: VECM WITH QUARTERLY DATA\n")
cat("================================================================================\n\n")

# Prepare quarterly credit growth for VECM
# Aggregate monthly credit to quarterly (end-of-quarter)
datos_raw <- read_excel(CONFIG$fci_file)
fecha_col_raw <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]

credit_monthly <- datos_raw %>%
  rename(fecha = !!sym(fecha_col_raw)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  mutate(
    Credit_yoy = (Creditos_Sector_privado_totales /
                   lag(Creditos_Sector_privado_totales, 12) - 1) * 100
  ) %>%
  dplyr::select(fecha, Credit_yoy)

credit_quarterly <- credit_monthly %>%
  mutate(month_in_q = month(fecha) %% 3) %>%
  filter(month_in_q == 0) %>%
  mutate(q_date = fecha) %>%
  dplyr::select(q_date, Credit_yoy)

# Merge with quarterly data
vecm_data <- analysis_q %>%
  left_join(credit_quarterly, by = "q_date") %>%
  dplyr::select(q_date, FCI_COMP, Credit_yoy, PIB_yoy) %>%
  na.omit()

cat(sprintf("VECM data: %d observations, %s to %s\n\n",
            nrow(vecm_data), format(min(vecm_data$q_date)),
            format(max(vecm_data$q_date))))

vecm_ts <- as.matrix(vecm_data[, c("FCI_COMP", "Credit_yoy", "PIB_yoy")])

vecm_summary <- data.frame()

tryCatch({
  # Lag selection
  lag_select <- vars::VARselect(vecm_ts, lag.max = CONFIG$vecm_lag_max, type = "const")
  cat("Information Criteria:\n")
  print(lag_select$selection)

  bic_lag <- lag_select$selection["SC(n)"]
  optimal_lag <- max(bic_lag, 2)

  # Johansen cointegration test
  K_johansen <- max(min(optimal_lag, 4), 2)
  cat(sprintf("\nUsing K = %d lags for Johansen test\n\n", K_johansen))

  johansen_trace <- urca::ca.jo(vecm_ts, type = "trace", ecdet = "const", K = K_johansen)
  cat("Johansen Trace Test:\n")
  print(summary(johansen_trace))

  # Determine rank
  trace_stats <- johansen_trace@teststat
  trace_crit <- johansen_trace@cval[, "5pct"]
  coint_rank <- sum(trace_stats > trace_crit)

  cat(sprintf("\nCointegration rank (5%%): r = %d\n", coint_rank))

  if (coint_rank > 0 && coint_rank < ncol(vecm_ts)) {
    cat(">>> COINTEGRATION DETECTED\n\n")

    vecm_model <- urca::cajorls(johansen_trace, r = coint_rank)
    cat("Cointegrating Vector:\n")
    print(vecm_model$beta)

    cat("\nError Correction Coefficients:\n")
    print(vecm_model$rlm$coefficients[1:coint_rank, ])

    vecm_summary <- data.frame(
      test = "Johansen Trace",
      K = K_johansen,
      rank = coint_rank,
      frequency = "Quarterly",
      n_obs = nrow(vecm_data)
    )
  } else if (coint_rank >= ncol(vecm_ts)) {
    cat(">>> FULL RANK: Variables appear stationary. No VECM needed.\n\n")
    vecm_summary <- data.frame(
      test = "Johansen Trace", K = K_johansen,
      rank = coint_rank, frequency = "Quarterly", n_obs = nrow(vecm_data)
    )
  } else {
    cat(">>> NO COINTEGRATION DETECTED\n\n")
    vecm_summary <- data.frame(
      test = "Johansen Trace", K = K_johansen,
      rank = 0, frequency = "Quarterly", n_obs = nrow(vecm_data)
    )
  }
}, error = function(e) {
  cat("Error in VECM analysis:", e$message, "\n")
  vecm_summary <<- data.frame(
    test = "Johansen Trace", K = NA, rank = NA,
    frequency = "Quarterly", n_obs = nrow(vecm_data)
  )
})

# --- LP: FCI -> FBCF and FCI -> Consumo Privado (replacing monthly FBKf/Consumo) ---
cat("\n--- Quarterly LP: FCI -> FBCF and Consumo Privado ---\n\n")

transmission_results <- list()

for (sector in c("FBCF", "Consumo_Privado")) {
  y_var <- paste0(sector, "_yoy")
  if (!y_var %in% names(analysis_q)) next

  data_sub <- analysis_q %>% filter(!is.na(!!sym(y_var)))
  cat(sprintf("  %s (n=%d) ... ", ALL_LABELS[sector], nrow(data_sub)))

  res <- run_lp_quarterly(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon_q, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy"),
    min_obs = CONFIG$min_obs
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- ALL_LABELS[sector]
    transmission_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

# --- Credit mediation analysis ---
cat("\n--- Credit Mediation Analysis (quarterly) ---\n\n")

mediation_results <- data.frame()

# Step 1: FCI -> Credit
cat("  FCI -> Credit_yoy ... ")
mediation_fci_credit <- NULL
if ("Credit_yoy" %in% names(analysis_q %>% left_join(credit_quarterly, by = "q_date"))) {
  med_data <- analysis_q %>%
    left_join(credit_quarterly, by = "q_date") %>%
    filter(!is.na(Credit_yoy), !is.na(PIB_yoy))

  mediation_fci_credit <- run_lp_quarterly(
    med_data, "Credit_yoy", "FCI_COMP", CONFIG$max_horizon_q, CONFIG$n_lags,
    control_vars = c("IPC_yoy"), min_obs = CONFIG$min_obs
  )
  if (!is.null(mediation_fci_credit) && nrow(mediation_fci_credit) > 0) {
    mediation_fci_credit$step <- "FCI->Credit"
    cat("done\n")
  }
}

# Step 2: FCI -> PIB without credit control
cat("  FCI -> PIB (without credit) ... ")
mediation_fci_pib_nocred <- run_lp_quarterly(
  analysis_q %>% filter(!is.na(PIB_yoy)),
  "PIB_yoy", "FCI_COMP", CONFIG$max_horizon_q, CONFIG$n_lags,
  control_vars = c("IPC_yoy"), min_obs = CONFIG$min_obs
)
if (nrow(mediation_fci_pib_nocred) > 0) {
  mediation_fci_pib_nocred$step <- "FCI->PIB (no credit)"
  cat("done\n")
}

# Step 3: FCI -> PIB with credit control (mediation test)
cat("  FCI -> PIB (with credit control) ... ")
med_data2 <- analysis_q %>%
  left_join(credit_quarterly, by = "q_date") %>%
  filter(!is.na(PIB_yoy), !is.na(Credit_yoy))

mediation_fci_pib_cred <- run_lp_quarterly(
  med_data2,
  "PIB_yoy", "FCI_COMP", CONFIG$max_horizon_q, CONFIG$n_lags,
  control_vars = c("IPC_yoy", "Credit_yoy"), min_obs = CONFIG$min_obs
)
if (nrow(mediation_fci_pib_cred) > 0) {
  mediation_fci_pib_cred$step <- "FCI->PIB (credit controlled)"
  cat("done\n")
}

mediation_all <- bind_rows(
  mediation_fci_credit,
  mediation_fci_pib_nocred,
  mediation_fci_pib_cred
)

if (nrow(mediation_all) > 0) {
  cat("\nMEDIATION RESULTS at h=4Q:\n")
  for (s in unique(mediation_all$step)) {
    v <- mediation_all %>% filter(step == s, horizon == 4)
    if (nrow(v) > 0) {
      cat(sprintf("  %s: coef = %+.3f, p = %.3f\n", s, v$coef, v$p_value))
    }
  }
}

# Plot 271: VECM and transmission
if (length(transmission_results) > 0) {
  trans_df <- bind_rows(transmission_results)

  p271 <- ggplot(trans_df, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#D35400", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#D35400", linewidth = 1) +
    geom_point(aes(y = coef), color = "#D35400", size = 2) +
    facet_wrap(~sector_label, scales = "free_y") +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    theme_minimal(base_size = 11) +
    labs(title = "Transmission Channels: FCI -> Investment & Consumption (Quarterly)",
         subtitle = "FCI_COMP \u2192 YoY growth | Published BCP data | 90% CI",
         x = "Quarters", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(CONFIG$output_dir, "271_PubQ_Transmission_FBCF_Consumo.png"), p271,
         width = 10, height = 5, dpi = 300)
  cat("Saved: 271_PubQ_Transmission_FBCF_Consumo.png\n")
}

write.csv(vecm_summary,
          file.path(CONFIG$output_dir, "csv", "PubQ_VECM_Summary.csv"), row.names = FALSE)
write.csv(mediation_all,
          file.path(CONFIG$output_dir, "csv", "PubQ_Mediation.csv"), row.names = FALSE)
cat("Saved: PubQ_VECM_Summary.csv, PubQ_Mediation.csv\n")


################################################################################
# 11. PART I: COMPARISON DASHBOARD
################################################################################

cat("\n================================================================================\n")
cat("PART I: COMPARISON DASHBOARD\n")
cat("================================================================================\n\n")

# Build comparison: quarterly h=2Q vs monthly h=6M, h=4Q vs h=12M, h=6Q vs h=18M
# Load original monthly results if available
orig_monthly_file <- file.path(CONFIG$output_dir, "csv", "LP_Sectoral_Output_Full.csv")
if (!file.exists(orig_monthly_file)) {
  orig_monthly_file <- file.path(CONFIG$output_dir, "LP_Sectoral_Output_Full.csv")
}

comparison_data <- data.frame()

if (file.exists(orig_monthly_file)) {
  cat("Loading original monthly results for comparison...\n")
  orig_monthly <- read.csv(orig_monthly_file, stringsAsFactors = FALSE)

  # Map quarterly sectors to monthly sectors
  sector_map <- data.frame(
    q_sector = c("Agricultura", "Manufactura", "Construccion", "Servicios", "PIB", "PIB_exAgri"),
    m_sector = c("Agricultura", "Manufactura", "Construcción", "Otros servicios", "PIB", "PIB SANB"),
    stringsAsFactors = FALSE
  )

  # Horizon mapping
  horizon_map <- data.frame(
    q_h = c(2, 4, 6),
    m_h = c(6, 12, 18),
    label = c("2Q vs 6M", "4Q vs 12M", "6Q vs 18M"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(sector_map)) {
    for (j in 1:nrow(horizon_map)) {
      q_res <- all_supply %>%
        filter(sector == sector_map$q_sector[i], horizon == horizon_map$q_h[j])
      m_res <- orig_monthly %>%
        filter(sector == sector_map$m_sector[i], horizon == horizon_map$m_h[j])

      if (nrow(q_res) > 0 && nrow(m_res) > 0) {
        # z-test for equality of coefficients
        z_stat <- (q_res$coef - m_res$coef) / sqrt(q_res$se^2 + m_res$se^2)
        z_p <- 2 * (1 - pnorm(abs(z_stat)))

        comparison_data <- rbind(comparison_data, data.frame(
          sector = sector_map$q_sector[i],
          sector_label = SUPPLY_LABELS[sector_map$q_sector[i]],
          horizon_match = horizon_map$label[j],
          q_coef = q_res$coef,
          q_se = q_res$se,
          q_p = q_res$p_value,
          q_n = q_res$n_obs,
          m_coef = m_res$coef,
          m_se = m_res$se,
          m_p = m_res$p_value,
          m_n = m_res$n_obs,
          z_stat = z_stat,
          z_p = z_p,
          equal = ifelse(z_p > 0.10, "Yes", "No")
        ))
      }
    }
  }

  if (nrow(comparison_data) > 0) {
    cat("\nCOMPARISON: Quarterly Published vs Original Monthly\n")
    cat(sprintf("%-20s %-12s %10s %10s %8s %8s\n",
                "Sector", "Horizons", "Q.Coef", "M.Coef", "z-stat", "Equal?"))
    cat(paste(rep("-", 75), collapse = ""), "\n")

    for (i in 1:nrow(comparison_data)) {
      r <- comparison_data[i, ]
      cat(sprintf("%-20s %-12s %+10.2f %+10.2f %8.2f %8s\n",
                  r$sector_label, r$horizon_match,
                  r$q_coef, r$m_coef, r$z_stat, r$equal))
    }
    cat("\nEqual? = z-test p > 0.10 (coefficients not statistically different)\n")
  }
} else {
  cat("Original monthly results not found. Skipping comparison.\n")
  cat("  Expected file: ", orig_monthly_file, "\n")
}

# --- Dual-engine test summary ---
cat("\n\nDUAL-ENGINE ECONOMY TEST RESULTS (Quarterly):\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

agri_4q <- all_supply %>% filter(sector == "Agricultura", horizon == 4)
pib_exagri_4q <- all_supply %>% filter(sector == "PIB_exAgri", horizon == 4)
fbcf_4q <- all_demand %>% filter(sector == "FBCF", horizon == 4)
constr_4q <- all_supply %>% filter(sector == "Construccion", horizon == 4)
exports_4q <- all_demand %>% filter(sector == "Exportaciones", horizon == 4)

cat("Key predictions and results at h=4Q (approx 12M):\n\n")

print_test <- function(name, result, expected_sig) {
  if (nrow(result) > 0) {
    sig <- result$p_value < 0.10
    match <- sig == expected_sig
    cat(sprintf("  %-25s coef=%+7.2f  p=%.3f  %s  %s\n",
                name, result$coef, result$p_value,
                ifelse(sig, "SIG", "n.s."),
                ifelse(match, "[CONFIRMED]", "[UNEXPECTED]")))
  } else {
    cat(sprintf("  %-25s  NOT AVAILABLE\n", name))
  }
}

print_test("Agriculture (expect n.s.)", agri_4q, FALSE)
print_test("GDP ex-Agri (expect sig)", pib_exagri_4q, TRUE)
print_test("Investment/FBCF (expect sig)", fbcf_4q, TRUE)
print_test("Construction (expect sig)", constr_4q, TRUE)
print_test("Exports (expect n.s.)", exports_4q, FALSE)

# --- Dashboard plot 272 ---
if (nrow(all_supply) > 0 || nrow(all_demand) > 0) {

  # Panel 1: Key supply sectors overlay
  key_supply <- all_supply %>%
    filter(sector %in% c("Agricultura", "PIB_exAgri", "Construccion", "Manufactura"))

  p1 <- ggplot(key_supply, aes(x = horizon, y = coef, color = sector_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    theme_minimal(base_size = 10) +
    labs(title = "A. Key Supply Sectors", x = "Quarters", y = "Effect (pp)",
         color = NULL) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  # Panel 2: Key demand components overlay
  key_demand <- all_demand %>%
    filter(sector %in% c("FBCF", "Consumo_Privado", "Exportaciones"))

  p2 <- ggplot(key_demand, aes(x = horizon, y = coef, color = sector_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    theme_minimal(base_size = 10) +
    labs(title = "B. Key Demand Components", x = "Quarters", y = "Effect (pp)",
         color = NULL) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  # Panel 3: Agriculture vs GDP ex-Agri with CIs
  agri_vs_gdp <- all_supply %>%
    filter(sector %in% c("Agricultura", "PIB_exAgri"))

  p3 <- ggplot(agri_vs_gdp, aes(x = horizon, color = sector_label, fill = sector_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(aes(y = coef), linewidth = 1) +
    geom_point(aes(y = coef), size = 2) +
    scale_x_continuous(breaks = 1:CONFIG$max_horizon_q) +
    scale_color_manual(values = c("Agriculture" = "#27AE60", "GDP ex-Agriculture" = "#2C3E50")) +
    scale_fill_manual(values = c("Agriculture" = "#27AE60", "GDP ex-Agriculture" = "#2C3E50")) +
    theme_minimal(base_size = 10) +
    labs(title = "C. Dual-Engine: Agriculture vs GDP ex-Agri",
         x = "Quarters", y = "Effect (pp)", color = NULL, fill = NULL) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  # Panel 4: Comparison heatmap
  hm_data <- all_both %>%
    filter(horizon == 4) %>%
    mutate(label = sprintf("%.1f%s", coef,
                           ifelse(p_value < 0.01, "***",
                                  ifelse(p_value < 0.05, "**",
                                         ifelse(p_value < 0.10, "*", "")))))

  p4 <- ggplot(hm_data, aes(x = side, y = sector_label, fill = coef)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = label), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "pp") +
    theme_minimal(base_size = 10) +
    labs(title = "D. FCI Effect at h=4Q", x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"))

  p272 <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2,
    top = grid::textGrob("Published Quarterly Data: Dual-Engine Economy Test",
                          gp = grid::gpar(fontsize = 14, fontface = "bold")))

  ggsave(file.path(CONFIG$output_dir, "272_PubQ_Comparison_Dashboard.png"), p272,
         width = 16, height = 12, dpi = 300)
  cat("Saved: 272_PubQ_Comparison_Dashboard.png\n")
}

write.csv(comparison_data,
          file.path(CONFIG$output_dir, "csv", "PubQ_Comparison_Summary.csv"), row.names = FALSE)
cat("Saved: PubQ_Comparison_Summary.csv\n")


################################################################################
# 12. FINAL SUMMARY
################################################################################

cat("\n################################################################################\n")
cat("QUARTERLY LP WITH PUBLISHED DATA - COMPLETE\n")
cat("################################################################################\n\n")

cat("Outputs generated:\n")
cat("  PNGs:  260-272 series\n")
cat("  CSVs:  PubQ_LP_Supply_Full, PubQ_LP_Demand_Full, PubQ_LP_EndoExo_Supply/Demand,\n")
cat("         PubQ_LP_Asymmetric, PubQ_LP_PostIT, PubQ_TOST_Agriculture,\n")
cat("         PubQ_VECM_Summary, PubQ_Mediation, PubQ_Comparison_Summary\n\n")

cat("Key verification checks:\n")
cat("  1. Agriculture insensitive at h=2Q, h=4Q? -> Check 260 and summary above\n")
cat("  2. FBCF and Construction significant? -> Check 262 and summary above\n")
cat("  3. GDP ex-Agri more significant than GDP? -> Check supply summary\n")
cat("  4. Exports insensitive (falsification)? -> Check demand summary\n")
cat("  5. Quarterly h=4Q ~ monthly h=12? -> Check comparison dashboard\n\n")

cat("Output directory:", CONFIG$output_dir, "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
