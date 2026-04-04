################################################################################
# FCI SECTORAL OUTPUT PUZZLE INVESTIGATION & EXTENDED LP
################################################################################
#
# Project:      Financial Conditions Index - Sectoral Analysis
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Disaggregated sectoral analysis to resolve the "output puzzle"
#               (FCI predicts credit but not aggregate output).
#
#   PART A: SECTORAL CREDIT LP (FCI → Sectoral Credit Growth)
#     A1 - Long sample (Dec 2001+): 5 sectors
#     A2 - Post-2016 sample: 12 sectors
#     NOTE: Skipped if sectoral credit data not available in consolidated file
#
#   PART B: SECTORAL OUTPUT LP (FCI → Sectoral Output Growth)
#     B1 - Full sample (1994Q1+): all sectors (quarterly data from Quarterly_SA)
#     B2 - Post-IT sample (2011+): all sectors
#
#   PART C: EXTENDED HORIZON REPLICATION (h = 1..36)
#     Replicates Chapter 5 LP plots with extended horizon
#
#   PART D: SUMMARY REPORT
#     Generates markdown report with key findings
#
# Data:
#   - FCI_data_1.xlsx: Quarterly_SA (sectoral GDP), Datos_macro (macro + IPC)
#
# References:
#   - Jordà (2005) AER - Local Projections
#   - Hatzius et al. (2010) - FCI motivation
#   - Adrian, Boyarchenko & Giannone (2019) - Vulnerable Growth
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
  library(lubridate)
})

CONFIG <- list(
  # Data files
  puzzle_file = "../data/FCI_data_1.xlsx",
  fci_file    = "../data/FCI_data_1.xlsx",
  macro_sheet = "Datos_macro",

  # LP parameters (quarterly: 8 quarters = 24 months)
  max_horizon = 8,
  extended_horizon = 36,
  n_lags = 2,
  confidence_level = 0.90,
  min_obs = 25,

  # Sectoral credit definitions (only used if Credits sheet available)
  credit_long = c("AGRICULTURA", "COMERCIO AL POR MAYOR", "CONSUMO",
                   "GANADERIA", "INDUSTRIA"),
  credit_short = c("AGRICULTURA", "COMERCIO AL POR MAYOR", "CONSUMO",
                    "GANADERIA", "INDUSTRIA", "CONSTRUCCION",
                    "COMERCIO AL POR MENOR", "ACTIVIDADES INMOBILIARIAS",
                    "SERVICIOS", "SECTOR FINANCIERO", "OTROS", "VIVIENDA"),

  # Activity sectors (from Quarterly_SA sheet)
  activity_sectors = c("Agricultura", "Ganaderia_fp", "Manufactura",
                        "Electricidad_y_agua", "Construccion", "Servicios",
                        "PIB", "PIB_exAgri", "Consumo_Privado", "FBCF"),

  # Output
  output_dir = "../output"
)

# Short labels for plotting
CREDIT_SHORT_LABELS <- c(
  "AGRICULTURA" = "Agriculture",
  "COMERCIO AL POR MAYOR" = "Wholesale",
  "CONSUMO" = "Consumer",
  "GANADERIA" = "Livestock",
  "INDUSTRIA" = "Industry",
  "CONSTRUCCION" = "Construction",
  "COMERCIO AL POR MENOR" = "Retail",
  "ACTIVIDADES INMOBILIARIAS" = "Real Estate",
  "SERVICIOS" = "Services",

  "SECTOR FINANCIERO" = "Financial",
  "OTROS" = "Other",
  "VIVIENDA" = "Housing"
)

ACTIVITY_LABELS <- c(
  "Agricultura"         = "Agriculture",
  "Ganaderia_fp"        = "Livestock+",
  "Manufactura"         = "Manufacturing",
  "Electricidad_y_agua" = "Electricity & Water",
  "Construccion"        = "Construction",
  "Servicios"           = "Services",
  "PIB"                 = "GDP",
  "PIB_exAgri"          = "GDP ex-Agriculture",
  "Consumo_Privado"     = "Private Consumption",
  "FBCF"                = "Investment (GFCF)"
)

set.seed(20250209)

cat("\n################################################################################\n")
cat("FCI SECTORAL OUTPUT PUZZLE INVESTIGATION & EXTENDED LP\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# 2. LP FUNCTIONS (self-contained, copied from script 05)
################################################################################

#' Standard Local Projection
run_lp_standard <- function(data, y_var, fci_var, max_h, n_lags = 2,
                            control_vars = NULL, conf_level = 0.90) {

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

    if (nrow(reg_data) < 30) next

    model <- lm(as.formula(formula_str), data = reg_data)
    vcov_hac <- sandwich::NeweyWest(model, lag = h + 1, prewhite = FALSE)
    coef_test <- lmtest::coeftest(model, vcov = vcov_hac)

    idx <- which(rownames(coef_test) == fci_var)

    results <- rbind(results, data.frame(
      horizon = h,
      coef = coef_test[idx, 1],
      se = coef_test[idx, 2],
      p_value = coef_test[idx, 4],
      ci_lower = coef_test[idx, 1] - z_crit * coef_test[idx, 2],
      ci_upper = coef_test[idx, 1] + z_crit * coef_test[idx, 2],
      n_obs = nrow(reg_data)
    ))
  }

  return(results)
}

#' Asymmetric (Non-linear) Local Projection
run_lp_asymmetric <- function(data, y_var, fci_var, max_h, n_lags = 2,
                              control_vars = NULL, conf_level = 0.90) {

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

    if (nrow(reg_data) < 30) next

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
# 3. DATA LOADING AND PREPARATION
################################################################################

cat("================================================================================\n")
cat("DATA LOADING AND PREPARATION\n")
cat("================================================================================\n\n")

if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)

# --- Load FCI ---
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  CONFIG_14 <- CONFIG  # Save our CONFIG before source overwrites it
  source("01_FCI_Complete.R")
  CONFIG <- CONFIG_14  # Restore our CONFIG
}

fci_monthly <- resultado_fci$all_indices %>%
  dplyr::select(fecha,
                FCI_COMP = FCI_COMP_AVG,
                FCI_ENDO = FCI_ENDO_AVG,
                FCI_EXO = FCI_EXO_AVG,
                FCI_exCredit = FCI_exCredit_AVG,
                FCI_ENDO_exCredit = FCI_ENDO_exCredit_AVG)

# --- Aggregate monthly FCI to quarterly (end-of-quarter: months 3, 6, 9, 12) ---
cat("Aggregating monthly FCI to quarterly (end-of-quarter)...\n")

fci_quarterly <- fci_monthly %>%
  mutate(
    year = year(fecha),
    quarter = quarter(fecha),
    month_in_q = month(fecha) %% 3
  ) %>%
  filter(month_in_q == 0) %>%
  mutate(q_date = fecha) %>%
  dplyr::select(q_date, FCI_COMP, FCI_ENDO, FCI_EXO, FCI_exCredit, FCI_ENDO_exCredit)

cat(sprintf("  FCI quarterly: %d observations, %s to %s\n",
            nrow(fci_quarterly), format(min(fci_quarterly$q_date)),
            format(max(fci_quarterly$q_date))))

# --- Load macro data (monthly, then aggregate to quarterly) ---
macro_raw <- read_excel(CONFIG$fci_file, sheet = CONFIG$macro_sheet)
fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

macro_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  mutate(
    IMAEP_yoy = (IMAEP / lag(IMAEP, 12) - 1) * 100,
    IPC_yoy   = (IPC / lag(IPC, 12) - 1) * 100,
    Cred_Real_yoy = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100
  )

# Aggregate monthly controls to quarterly (end-of-quarter)
controls_quarterly <- macro_data %>%
  dplyr::select(fecha, IPC_yoy, Cred_Real_yoy) %>%
  mutate(month_in_q = month(fecha) %% 3) %>%
  filter(month_in_q == 0) %>%
  mutate(q_date = fecha) %>%
  dplyr::select(q_date, IPC_yoy, Cred_Real_yoy)

# Get IPC level for credit deflation (monthly, for Part A if available)
ipc_series <- macro_data %>% dplyr::select(fecha, IPC)

# --- Check for sectoral credit data ---
credit_data_available <- FALSE
credit_growth <- NULL

tryCatch({
  cat("Checking for sectoral credit data (Credits sheet)...\n")
  credit_raw <- read_excel(CONFIG$puzzle_file, sheet = "Credits")
  names(credit_raw)[1] <- "fecha"
  credit_raw <- credit_raw %>%
    mutate(fecha = as.Date(fecha)) %>%
    arrange(fecha)

  # Deflate credit stocks by IPC and compute YoY real growth
  credit_sectors <- credit_raw %>%
    left_join(ipc_series, by = "fecha")

  credit_cols <- setdiff(names(credit_raw), c("fecha", "Total", "Total general",
                                               "ADMINISTRACION PUBLICA"))

  credit_growth <- credit_sectors %>% dplyr::select(fecha)

  for (col in credit_cols) {
    if (col %in% names(credit_sectors) && col != "IPC") {
      real_stock <- credit_sectors[[col]] / credit_sectors$IPC
      yoy_growth <- (real_stock / lag(real_stock, 12) - 1) * 100
      safe_name <- paste0("CG_", gsub(" ", "_", col))
      credit_growth[[safe_name]] <- yoy_growth
    }
  }

  cat(sprintf("  Credit sectors: %d columns, %d observations\n",
              length(credit_cols), nrow(credit_growth)))
  credit_data_available <- TRUE
}, error = function(e) {
  cat("NOTE: Sectoral credit data not available in consolidated file. Skipping Part A.\n\n")
})

# --- Load sectoral activity data from Quarterly_SA ---
cat("Loading sectoral activity data from Quarterly_SA...\n")
quarterly_raw <- read_excel(CONFIG$puzzle_file, sheet = "Quarterly_SA")
names(quarterly_raw)[1] <- "fecha"

quarterly_data <- quarterly_raw %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Rename columns to safe names (same mapping as script 27)
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

# Compute YoY growth rates (Y_t / Y_{t-4} - 1) * 100 for quarterly data
for (s in CONFIG$activity_sectors) {
  if (s %in% names(quarterly_data)) {
    yoy_name <- paste0("AG_", s)
    quarterly_data <- quarterly_data %>%
      mutate(!!yoy_name := (!!sym(s) / lag(!!sym(s), 4) - 1) * 100)
  }
}

quarterly_data <- quarterly_data %>%
  mutate(q_date = fecha)

cat(sprintf("  Quarterly activity: %d observations, %s to %s\n",
            nrow(quarterly_data), format(min(quarterly_data$fecha)),
            format(max(quarterly_data$fecha))))

# --- Merge quarterly datasets ---
analysis_full <- fci_quarterly %>%
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
  rename(fecha = q_date) %>%
  arrange(fecha)

cat(sprintf("\nMerged quarterly dataset: %d observations, %s to %s\n\n",
            nrow(analysis_full),
            format(min(analysis_full$fecha)),
            format(max(analysis_full$fecha))))

# Color palette for many sectors
sector_colors <- c(
  "#2C3E50", "#E74C3C", "#3498DB", "#27AE60", "#F39C12",
  "#8E44AD", "#1ABC9C", "#D35400", "#34495E", "#E91E63",
  "#009688", "#795548"
)


################################################################################
# 4. PART A: SECTORAL CREDIT LP (FCI → Sectoral Credit Growth)
#    NOTE: Only runs if sectoral credit data is available (Credits sheet)
################################################################################

if (credit_data_available) {

# --- Need to merge credit data with monthly FCI for Part A ---
# Part A uses monthly data, so we need a monthly analysis dataset
fci_data_monthly <- fci_monthly
macro_controls_monthly <- macro_data %>%
  dplyr::select(fecha, IMAEP_yoy, IPC_yoy, Cred_Real_yoy)

analysis_full_monthly <- fci_data_monthly %>%
  inner_join(macro_controls_monthly, by = "fecha") %>%
  inner_join(credit_growth, by = "fecha") %>%
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
  arrange(fecha)

# Part A uses monthly horizons (up to 24 months)
max_horizon_monthly <- 24

cat("================================================================================\n")
cat("PART A: SECTORAL CREDIT LOCAL PROJECTIONS\n")
cat("================================================================================\n\n")

# --- A1: Long sample (Dec 2001+, 5 sectors) ---
cat("--- A1: Long Sample (5 sectors, Dec 2001+) ---\n\n")

credit_long_results <- list()

for (sector in CONFIG$credit_long) {
  y_var <- paste0("CG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full_monthly)) {
    cat(sprintf("  SKIP: %s not found\n", sector))
    next
  }

  # Filter to long-sample available data
  data_sub <- analysis_full_monthly %>% filter(!is.na(!!sym(y_var)))

  cat(sprintf("  %s (n=%d) ... ", CREDIT_SHORT_LABELS[sector], nrow(data_sub)))

  res <- run_lp_standard(
    data_sub, y_var, "FCI_exCredit", max_horizon_monthly, CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- CREDIT_SHORT_LABELS[sector]
    credit_long_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_credit_long <- bind_rows(credit_long_results)

cat(sprintf("\n  Long sample results: %d sector-horizon pairs\n\n", nrow(all_credit_long)))

# --- A2: Short sample (Jan 2016+, 12 sectors) ---
cat("--- A2: Short Sample (12 sectors, Jan 2016+) ---\n\n")

credit_short_results <- list()

for (sector in CONFIG$credit_short) {
  y_var <- paste0("CG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full_monthly)) {
    cat(sprintf("  SKIP: %s not found\n", sector))
    next
  }

  # Filter to post-2016 data
  data_sub <- analysis_full_monthly %>%
    filter(fecha >= as.Date("2016-01-01"), !is.na(!!sym(y_var)))

  cat(sprintf("  %s (n=%d) ... ", CREDIT_SHORT_LABELS[sector], nrow(data_sub)))

  res <- run_lp_standard(
    data_sub, y_var, "FCI_exCredit", max_horizon_monthly, CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- CREDIT_SHORT_LABELS[sector]
    credit_short_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_credit_short <- bind_rows(credit_short_results)

cat(sprintf("\n  Short sample results: %d sector-horizon pairs\n\n", nrow(all_credit_short)))

# --- A Summary ---
cat("PART A SUMMARY at h=12:\n")
cat(sprintf("%-25s %10s %10s %8s  %s\n", "Sector", "Coef", "SE", "p-value", "Sample"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (sector in CONFIG$credit_long) {
  v <- all_credit_long %>% filter(sector == !!sector, horizon == 12)
  if (nrow(v) > 0) {
    stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                    ifelse(v$p_value < 0.10, "*", "")))
    cat(sprintf("%-25s %+10.2f %10.2f %8.3f%s  Long\n",
                CREDIT_SHORT_LABELS[sector], v$coef, v$se, v$p_value, stars))
  }
}
cat("\n")
for (sector in CONFIG$credit_short) {
  v <- all_credit_short %>% filter(sector == !!sector, horizon == 12)
  if (nrow(v) > 0) {
    stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                    ifelse(v$p_value < 0.10, "*", "")))
    cat(sprintf("%-25s %+10.2f %10.2f %8.3f%s  Short\n",
                CREDIT_SHORT_LABELS[sector], v$coef, v$se, v$p_value, stars))
  }
}


################################################################################
# 5. PART A VISUALIZATIONS
################################################################################

cat("\n\n================================================================================\n")
cat("PART A VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Plot 140: Long sample sectoral credit IRFs
if (nrow(all_credit_long) > 0) {
  p140 <- ggplot(all_credit_long, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#377EB8", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#377EB8", linewidth = 1) +
    geom_point(aes(y = coef), color = "#377EB8", size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 3) +
    theme_minimal(base_size = 11) +
    labs(title = "Sectoral Credit Response to FCI Tightening (Long Sample)",
         subtitle = "FCI_exCredit → Real credit growth (YoY) | Dec 2001+ | 90% CI",
         x = "Months", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(CONFIG$output_dir, "140_LP_Sectoral_Credit_Long.png"), p140,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 140_LP_Sectoral_Credit_Long.png\n")
}

# Plot 141: Short sample sectoral credit IRFs
if (nrow(all_credit_short) > 0) {
  p141 <- ggplot(all_credit_short, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#E74C3C", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#E74C3C", linewidth = 1) +
    geom_point(aes(y = coef), color = "#E74C3C", size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    theme_minimal(base_size = 10) +
    labs(title = "Sectoral Credit Response to FCI Tightening (Post-2016)",
         subtitle = "FCI_exCredit → Real credit growth (YoY) | Jan 2016+ | 90% CI",
         x = "Months", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9))

  ggsave(file.path(CONFIG$output_dir, "141_LP_Sectoral_Credit_Short.png"), p141,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 141_LP_Sectoral_Credit_Short.png\n")
}

# Plot 142: Heatmap at h=12 (both samples)
heatmap_long <- all_credit_long %>%
  filter(horizon == 12) %>%
  mutate(sample = "Long (2001+)")

heatmap_short <- all_credit_short %>%
  filter(horizon == 12) %>%
  mutate(sample = "Short (2016+)")

heatmap_credit <- bind_rows(heatmap_long, heatmap_short) %>%
  mutate(label = sprintf("%.1f%s", coef,
                         ifelse(p_value < 0.01, "***",
                                ifelse(p_value < 0.05, "**",
                                       ifelse(p_value < 0.10, "*", "")))))

if (nrow(heatmap_credit) > 0) {
  p142 <- ggplot(heatmap_credit, aes(x = sample, y = sector_label, fill = coef)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "Effect (pp)") +
    theme_minimal(base_size = 12) +
    labs(title = "Sectoral Credit: FCI Effect at 12-Month Horizon",
         subtitle = "Effect in pp | *** p<0.01, ** p<0.05, * p<0.10",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size = 10))

  ggsave(file.path(CONFIG$output_dir, "142_LP_Sectoral_Credit_Heatmap.png"), p142,
         width = 8, height = 8, dpi = 300)
  cat("Saved: 142_LP_Sectoral_Credit_Heatmap.png\n")
}

# Export Part A CSVs
write.csv(all_credit_long, file.path(CONFIG$output_dir, "LP_Sectoral_Credit_Long.csv"),
          row.names = FALSE)
write.csv(all_credit_short, file.path(CONFIG$output_dir, "LP_Sectoral_Credit_Short.csv"),
          row.names = FALSE)
cat("Saved: LP_Sectoral_Credit_Long.csv, LP_Sectoral_Credit_Short.csv\n")


################################################################################
# 5b. PART A ASYMMETRIC: SECTORAL CREDIT (Tightening vs Easing)
################################################################################

cat("\n================================================================================\n")
cat("PART A ASYMMETRIC: SECTORAL CREDIT (Tightening vs Easing)\n")
cat("================================================================================\n\n")

# --- A1 Asymmetric: Long sample ---
cat("--- A1 Asymmetric: Long Sample (5 sectors) ---\n\n")

credit_long_asym <- list()

for (sector in CONFIG$credit_long) {
  y_var <- paste0("CG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full_monthly)) next

  data_sub <- analysis_full_monthly %>% filter(!is.na(!!sym(y_var)))

  cat(sprintf("  %s (n=%d) ... ", CREDIT_SHORT_LABELS[sector], nrow(data_sub)))

  res <- run_lp_asymmetric(
    data_sub, y_var, "FCI_exCredit", max_horizon_monthly, CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- CREDIT_SHORT_LABELS[sector]
    credit_long_asym[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_credit_long_asym <- bind_rows(credit_long_asym)

# --- A2 Asymmetric: Short sample ---
cat("\n--- A2 Asymmetric: Short Sample (12 sectors) ---\n\n")

credit_short_asym <- list()

for (sector in CONFIG$credit_short) {
  y_var <- paste0("CG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full_monthly)) next

  data_sub <- analysis_full_monthly %>%
    filter(fecha >= as.Date("2016-01-01"), !is.na(!!sym(y_var)))

  cat(sprintf("  %s (n=%d) ... ", CREDIT_SHORT_LABELS[sector], nrow(data_sub)))

  res <- run_lp_asymmetric(
    data_sub, y_var, "FCI_exCredit", max_horizon_monthly, CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- CREDIT_SHORT_LABELS[sector]
    credit_short_asym[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_credit_short_asym <- bind_rows(credit_short_asym)

# --- Asymmetric credit summary ---
cat("\nASYMMETRIC CREDIT SUMMARY at h=12:\n")
cat(sprintf("%-20s %10s %10s %10s %8s  %s\n",
            "Sector", "Tight", "Ease", "Diff", "Asym.p", "Sample"))
cat(paste(rep("-", 75), collapse = ""), "\n")

for (sector in CONFIG$credit_long) {
  v <- all_credit_long_asym %>% filter(sector == !!sector, horizon == 12)
  if (nrow(v) > 0) {
    asym_star <- ifelse(v$asym_p < 0.05, "**", ifelse(v$asym_p < 0.10, "*", ""))
    cat(sprintf("%-20s %+10.2f %+10.2f %+10.2f %8.3f%s  Long\n",
                CREDIT_SHORT_LABELS[sector], v$coef_tight, v$coef_ease,
                v$coef_tight - v$coef_ease, v$asym_p, asym_star))
  }
}
cat("\n")
for (sector in CONFIG$credit_short) {
  v <- all_credit_short_asym %>% filter(sector == !!sector, horizon == 12)
  if (nrow(v) > 0) {
    asym_star <- ifelse(v$asym_p < 0.05, "**", ifelse(v$asym_p < 0.10, "*", ""))
    cat(sprintf("%-20s %+10.2f %+10.2f %+10.2f %8.3f%s  Short\n",
                CREDIT_SHORT_LABELS[sector], v$coef_tight, v$coef_ease,
                v$coef_tight - v$coef_ease, v$asym_p, asym_star))
  }
}

# --- Asymmetric credit visualizations ---
cat("\n\nGenerating asymmetric credit visualizations...\n")

# Plot 154: Asymmetric credit -- long sample
if (nrow(all_credit_long_asym) > 0) {
  asym_cr_long_plot <- all_credit_long_asym %>%
    pivot_longer(cols = c(coef_tight, coef_ease),
                 names_to = "effect_type", values_to = "coef") %>%
    mutate(
      effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
      ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
      ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
    )

  p154 <- ggplot(asym_cr_long_plot, aes(x = horizon, y = coef,
                                         color = effect_label, fill = effect_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("Tightening" = "#E74C3C", "Easing" = "#3498DB")) +
    scale_fill_manual(values = c("Tightening" = "#E74C3C", "Easing" = "#3498DB")) +
    theme_minimal(base_size = 11) +
    labs(title = "Asymmetric Sectoral Credit Response (Long Sample)",
         subtitle = "Tightening (FCI>0) vs Easing (|FCI<0|) | FCI_exCredit | Dec 2001+ | 90% CI",
         x = "Months", y = "Effect on credit growth (pp)",
         color = "FCI State", fill = "FCI State") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "154_LP_Sectoral_Credit_Asym_Long.png"), p154,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 154_LP_Sectoral_Credit_Asym_Long.png\n")
}

# Plot 155: Asymmetric credit -- short sample
if (nrow(all_credit_short_asym) > 0) {
  asym_cr_short_plot <- all_credit_short_asym %>%
    pivot_longer(cols = c(coef_tight, coef_ease),
                 names_to = "effect_type", values_to = "coef") %>%
    mutate(
      effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
      ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
      ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
    )

  p155 <- ggplot(asym_cr_short_plot, aes(x = horizon, y = coef,
                                          color = effect_label, fill = effect_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c("Tightening" = "#E74C3C", "Easing" = "#3498DB")) +
    scale_fill_manual(values = c("Tightening" = "#E74C3C", "Easing" = "#3498DB")) +
    theme_minimal(base_size = 10) +
    labs(title = "Asymmetric Sectoral Credit Response (Post-2016)",
         subtitle = "Tightening (FCI>0) vs Easing (|FCI<0|) | FCI_exCredit | Jan 2016+ | 90% CI",
         x = "Months", y = "Effect on credit growth (pp)",
         color = "FCI State", fill = "FCI State") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "155_LP_Sectoral_Credit_Asym_Short.png"), p155,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 155_LP_Sectoral_Credit_Asym_Short.png\n")
}

# Plot 156: Asymmetric credit heatmap at h=12
hm_asym_cr_long <- all_credit_long_asym %>%
  filter(horizon == 12) %>% mutate(sample = "Long (2001+)")
hm_asym_cr_short <- all_credit_short_asym %>%
  filter(horizon == 12) %>% mutate(sample = "Short (2016+)")

hm_asym_credit <- bind_rows(hm_asym_cr_long, hm_asym_cr_short) %>%
  pivot_longer(cols = c(coef_tight, coef_ease),
               names_to = "effect_type", values_to = "coef_val") %>%
  mutate(
    p_val = ifelse(effect_type == "coef_tight", p_tight, p_ease),
    effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
    label = sprintf("%.1f%s", coef_val,
                    ifelse(p_val < 0.01, "***",
                           ifelse(p_val < 0.05, "**",
                                  ifelse(p_val < 0.10, "*", ""))))
  )

if (nrow(hm_asym_credit) > 0) {
  p156 <- ggplot(hm_asym_credit,
                 aes(x = interaction(effect_label, sample, sep = "\n"),
                     y = sector_label, fill = coef_val)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "Effect (pp)") +
    theme_minimal(base_size = 11) +
    labs(title = "Asymmetric Sectoral Credit Effects at h=12",
         subtitle = "Tightening vs Easing | *** p<0.01, ** p<0.05, * p<0.10",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(size = 9))

  ggsave(file.path(CONFIG$output_dir, "156_LP_Sectoral_Credit_Asym_Heatmap.png"), p156,
         width = 10, height = 8, dpi = 300)
  cat("Saved: 156_LP_Sectoral_Credit_Asym_Heatmap.png\n")
}

# Export asymmetric credit CSVs
write.csv(all_credit_long_asym,
          file.path(CONFIG$output_dir, "LP_Sectoral_Credit_Asym_Long.csv"), row.names = FALSE)
write.csv(all_credit_short_asym,
          file.path(CONFIG$output_dir, "LP_Sectoral_Credit_Asym_Short.csv"), row.names = FALSE)
cat("Saved: LP_Sectoral_Credit_Asym_Long.csv, LP_Sectoral_Credit_Asym_Short.csv\n")


################################################################################
# 5c. PART A ENDO/EXO: SECTORAL CREDIT (Domestic vs External FCI)
################################################################################

cat("\n================================================================================\n")
cat("PART A ENDO/EXO: SECTORAL CREDIT (Domestic vs External FCI)\n")
cat("================================================================================\n\n")

# For credit-LHS: use FCI_ENDO_exCredit (endogeneity-corrected) and FCI_EXO (no correction needed)
fci_credit_types <- c("FCI_ENDO_exCredit" = "Endogenous (Domestic)",
                       "FCI_EXO" = "Exogenous (External)")
colors_endo_exo <- c("Endogenous (Domestic)" = "#E74C3C",
                      "Exogenous (External)" = "#3498DB")

# --- A1 Endo/Exo: Long sample ---
cat("--- A1 Endo/Exo: Long Sample (5 sectors) ---\n\n")

credit_long_endoexo <- list()

for (sector in CONFIG$credit_long) {
  y_var <- paste0("CG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full_monthly)) next

  data_sub <- analysis_full_monthly %>% filter(!is.na(!!sym(y_var)))

  for (fci_var in names(fci_credit_types)) {
    cat(sprintf("  %s x %s (n=%d) ... ",
                CREDIT_SHORT_LABELS[sector], fci_credit_types[fci_var], nrow(data_sub)))

    res <- run_lp_standard(
      data_sub, y_var, fci_var, max_horizon_monthly, CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    )

    if (nrow(res) > 0) {
      res$sector <- sector
      res$sector_label <- CREDIT_SHORT_LABELS[sector]
      res$fci_var <- fci_var
      res$fci_label <- fci_credit_types[fci_var]
      credit_long_endoexo[[paste(sector, fci_var)]] <- res
      cat("done\n")
    } else {
      cat("insufficient data\n")
    }
  }
}

all_credit_long_endoexo <- bind_rows(credit_long_endoexo)

# --- A2 Endo/Exo: Short sample ---
cat("\n--- A2 Endo/Exo: Short Sample (12 sectors) ---\n\n")

credit_short_endoexo <- list()

for (sector in CONFIG$credit_short) {
  y_var <- paste0("CG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full_monthly)) next

  data_sub <- analysis_full_monthly %>%
    filter(fecha >= as.Date("2016-01-01"), !is.na(!!sym(y_var)))

  for (fci_var in names(fci_credit_types)) {
    cat(sprintf("  %s x %s (n=%d) ... ",
                CREDIT_SHORT_LABELS[sector], fci_credit_types[fci_var], nrow(data_sub)))

    res <- run_lp_standard(
      data_sub, y_var, fci_var, max_horizon_monthly, CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    )

    if (nrow(res) > 0) {
      res$sector <- sector
      res$sector_label <- CREDIT_SHORT_LABELS[sector]
      res$fci_var <- fci_var
      res$fci_label <- fci_credit_types[fci_var]
      credit_short_endoexo[[paste(sector, fci_var)]] <- res
      cat("done\n")
    } else {
      cat("insufficient data\n")
    }
  }
}

all_credit_short_endoexo <- bind_rows(credit_short_endoexo)

# --- Endo/Exo credit summary ---
cat("\nENDO/EXO CREDIT SUMMARY at h=6, h=12:\n")
cat(sprintf("%-20s %10s %10s %10s %10s  %s\n",
            "Sector", "Endo h=6", "Exo h=6", "Endo h=12", "Exo h=12", "Sample"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (sector in CONFIG$credit_long) {
  ve6 <- all_credit_long_endoexo %>% filter(sector == !!sector, horizon == 6, fci_var == "FCI_ENDO_exCredit")
  vx6 <- all_credit_long_endoexo %>% filter(sector == !!sector, horizon == 6, fci_var == "FCI_EXO")
  ve12 <- all_credit_long_endoexo %>% filter(sector == !!sector, horizon == 12, fci_var == "FCI_ENDO_exCredit")
  vx12 <- all_credit_long_endoexo %>% filter(sector == !!sector, horizon == 12, fci_var == "FCI_EXO")

  fmt <- function(v) {
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      sprintf("%+.2f%s", v$coef, stars)
    } else "NA"
  }

  cat(sprintf("%-20s %10s %10s %10s %10s  Long\n",
              CREDIT_SHORT_LABELS[sector], fmt(ve6), fmt(vx6), fmt(ve12), fmt(vx12)))
}
cat("\n")

# --- Endo/Exo credit visualizations ---
cat("\nGenerating endo/exo credit visualizations...\n")

# Plot 160: Endo vs Exo credit -- long sample
if (nrow(all_credit_long_endoexo) > 0) {
  p160 <- ggplot(all_credit_long_endoexo,
                 aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 3) +
    scale_color_manual(values = colors_endo_exo) +
    scale_fill_manual(values = colors_endo_exo) +
    theme_minimal(base_size = 11) +
    labs(title = "Sectoral Credit: Endogenous vs Exogenous FCI (Long Sample)",
         subtitle = "FCI_ENDO_exCredit vs FCI_EXO → Real credit growth (YoY) | Dec 2001+ | 90% CI",
         x = "Months", y = "Effect on credit growth (pp)",
         color = "FCI Component", fill = "FCI Component") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "160_LP_Sectoral_Credit_EndoExo_Long.png"), p160,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 160_LP_Sectoral_Credit_EndoExo_Long.png\n")
}

# Plot 161: Endo vs Exo credit -- short sample
if (nrow(all_credit_short_endoexo) > 0) {
  p161 <- ggplot(all_credit_short_endoexo,
                 aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = colors_endo_exo) +
    scale_fill_manual(values = colors_endo_exo) +
    theme_minimal(base_size = 10) +
    labs(title = "Sectoral Credit: Endogenous vs Exogenous FCI (Post-2016)",
         subtitle = "FCI_ENDO_exCredit vs FCI_EXO → Real credit growth (YoY) | Jan 2016+ | 90% CI",
         x = "Months", y = "Effect on credit growth (pp)",
         color = "FCI Component", fill = "FCI Component") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "161_LP_Sectoral_Credit_EndoExo_Short.png"), p161,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 161_LP_Sectoral_Credit_EndoExo_Short.png\n")
}

# Plot 162: Endo/Exo credit heatmap at h=6 and h=12
hm_endoexo_cr_long <- all_credit_long_endoexo %>%
  filter(horizon %in% c(6, 12)) %>% mutate(sample = "Long (2001+)")
hm_endoexo_cr_short <- all_credit_short_endoexo %>%
  filter(horizon %in% c(6, 12)) %>% mutate(sample = "Short (2016+)")

hm_endoexo_credit <- bind_rows(hm_endoexo_cr_long, hm_endoexo_cr_short) %>%
  mutate(
    horizon_label = paste0("h=", horizon),
    label = sprintf("%.1f%s", coef,
                    ifelse(p_value < 0.01, "***",
                           ifelse(p_value < 0.05, "**",
                                  ifelse(p_value < 0.10, "*", ""))))
  )

if (nrow(hm_endoexo_credit) > 0) {
  p162 <- ggplot(hm_endoexo_credit,
                 aes(x = interaction(fci_label, horizon_label, sample, sep = "\n"),
                     y = sector_label, fill = coef)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 2.8, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "Effect (pp)") +
    theme_minimal(base_size = 11) +
    labs(title = "Sectoral Credit: Endo vs Exo FCI Effects at h=6 and h=12",
         subtitle = "Effect in pp | *** p<0.01, ** p<0.05, * p<0.10",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(size = 7, angle = 0))

  ggsave(file.path(CONFIG$output_dir, "162_LP_Sectoral_Credit_EndoExo_Heatmap.png"), p162,
         width = 16, height = 8, dpi = 300)
  cat("Saved: 162_LP_Sectoral_Credit_EndoExo_Heatmap.png\n")
}

# Export Endo/Exo credit CSVs
write.csv(all_credit_long_endoexo,
          file.path(CONFIG$output_dir, "LP_Sectoral_Credit_EndoExo_Long.csv"), row.names = FALSE)
write.csv(all_credit_short_endoexo,
          file.path(CONFIG$output_dir, "LP_Sectoral_Credit_EndoExo_Short.csv"), row.names = FALSE)
cat("Saved: LP_Sectoral_Credit_EndoExo_Long.csv, LP_Sectoral_Credit_EndoExo_Short.csv\n")

} else {
  cat("\n================================================================================\n")
  cat("PART A: SKIPPED (sectoral credit data not available)\n")
  cat("================================================================================\n\n")

  # Initialize empty data frames so Part D summary does not error
  all_credit_long <- data.frame()
  all_credit_short <- data.frame()
  all_credit_long_asym <- data.frame()
  all_credit_short_asym <- data.frame()
  all_credit_long_endoexo <- data.frame()
  all_credit_short_endoexo <- data.frame()
}
# END of credit_data_available conditional


################################################################################
# 6. PART B: SECTORAL OUTPUT LP (FCI → Sectoral Output Growth)
################################################################################

cat("\n================================================================================\n")
cat("PART B: SECTORAL OUTPUT LOCAL PROJECTIONS (Quarterly)\n")
cat("================================================================================\n\n")

# --- B1: Full sample ---
cat("--- B1: Full Sample (1994Q1+, quarterly) ---\n\n")

output_full_results <- list()

for (sector in CONFIG$activity_sectors) {
  y_var <- paste0("AG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full)) {
    cat(sprintf("  SKIP: %s not found\n", sector))
    next
  }

  data_sub <- analysis_full %>% filter(!is.na(!!sym(y_var)))

  cat(sprintf("  %s (n=%d) ... ", ACTIVITY_LABELS[sector], nrow(data_sub)))

  res <- run_lp_standard(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy")
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- ACTIVITY_LABELS[sector]
    output_full_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_output_full <- bind_rows(output_full_results)

cat(sprintf("\n  Full sample results: %d sector-horizon pairs\n\n", nrow(all_output_full)))

# --- B2: Post-IT sample (2011Q2+) ---
cat("--- B2: Post-IT Sample (2011Q2+, quarterly) ---\n\n")

output_postIT_results <- list()

for (sector in CONFIG$activity_sectors) {
  y_var <- paste0("AG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full)) {
    cat(sprintf("  SKIP: %s not found\n", sector))
    next
  }

  data_sub <- analysis_full %>%
    filter(fecha >= as.Date("2011-06-01"), !is.na(!!sym(y_var)))

  cat(sprintf("  %s (n=%d) ... ", ACTIVITY_LABELS[sector], nrow(data_sub)))

  res <- run_lp_standard(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy")
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- ACTIVITY_LABELS[sector]
    output_postIT_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_output_postIT <- bind_rows(output_postIT_results)

cat(sprintf("\n  Post-IT results: %d sector-horizon pairs\n\n", nrow(all_output_postIT)))

# --- B Summary ---
cat("PART B SUMMARY at h=4 (4 quarters = 12 months):\n")
cat(sprintf("%-25s %10s %8s  %10s %8s\n", "Sector", "Full Coef", "p", "Post-IT Coef", "p"))
cat(paste(rep("-", 70), collapse = ""), "\n")

for (sector in CONFIG$activity_sectors) {
  vf <- all_output_full %>% filter(sector == !!sector, horizon == 4)
  vp <- all_output_postIT %>% filter(sector == !!sector, horizon == 4)

  coef_full <- if(nrow(vf) > 0) sprintf("%+.2f%s", vf$coef,
                   ifelse(vf$p_value < 0.01, "***", ifelse(vf$p_value < 0.05, "**",
                   ifelse(vf$p_value < 0.10, "*", "")))) else "NA"
  p_full <- if(nrow(vf) > 0) sprintf("%.3f", vf$p_value) else "NA"
  coef_pit <- if(nrow(vp) > 0) sprintf("%+.2f%s", vp$coef,
                   ifelse(vp$p_value < 0.01, "***", ifelse(vp$p_value < 0.05, "**",
                   ifelse(vp$p_value < 0.10, "*", "")))) else "NA"
  p_pit <- if(nrow(vp) > 0) sprintf("%.3f", vp$p_value) else "NA"

  cat(sprintf("%-25s %10s %8s  %12s %8s\n",
              ACTIVITY_LABELS[sector], coef_full, p_full, coef_pit, p_pit))
}


################################################################################
# 7. PART B VISUALIZATIONS
################################################################################

cat("\n\n================================================================================\n")
cat("PART B VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Plot 143: Full sample sectoral output IRFs
if (nrow(all_output_full) > 0) {
  p143 <- ggplot(all_output_full, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#2C3E50", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#2C3E50", linewidth = 1) +
    geom_point(aes(y = coef), color = "#2C3E50", size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    theme_minimal(base_size = 10) +
    labs(title = "Sectoral Output Response to FCI Tightening (Full Sample, Quarterly)",
         subtitle = "FCI_COMP → Output growth (YoY) | 1994Q1+ | 90% CI",
         x = "Quarters", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9))

  ggsave(file.path(CONFIG$output_dir, "143_LP_Sectoral_Output_Full.png"), p143,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 143_LP_Sectoral_Output_Full.png\n")
}

# Plot 144: Post-IT sectoral output IRFs
if (nrow(all_output_postIT) > 0) {
  p144 <- ggplot(all_output_postIT, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#8E44AD", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#8E44AD", linewidth = 1) +
    geom_point(aes(y = coef), color = "#8E44AD", size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    theme_minimal(base_size = 10) +
    labs(title = "Sectoral Output Response to FCI Tightening (Post-IT, Quarterly)",
         subtitle = "FCI_COMP → Output growth (YoY) | 2011Q2+ | 90% CI",
         x = "Quarters", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9))

  ggsave(file.path(CONFIG$output_dir, "144_LP_Sectoral_Output_PostIT.png"), p144,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 144_LP_Sectoral_Output_PostIT.png\n")
}

# Plot 145: Output heatmap at h=4 (4 quarters = 12 months)
hm_full <- all_output_full %>%
  filter(horizon == 4) %>%
  mutate(sample = "Full (1994Q1+)")

hm_pit <- all_output_postIT %>%
  filter(horizon == 4) %>%
  mutate(sample = "Post-IT (2011Q2+)")

heatmap_output <- bind_rows(hm_full, hm_pit) %>%
  mutate(label = sprintf("%.1f%s", coef,
                         ifelse(p_value < 0.01, "***",
                                ifelse(p_value < 0.05, "**",
                                       ifelse(p_value < 0.10, "*", "")))))

if (nrow(heatmap_output) > 0) {
  p145 <- ggplot(heatmap_output, aes(x = sample, y = sector_label, fill = coef)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "Effect (pp)") +
    theme_minimal(base_size = 12) +
    labs(title = "Sectoral Output: FCI Effect at h=4 (4 Quarters)",
         subtitle = "Effect in pp | *** p<0.01, ** p<0.05, * p<0.10",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size = 10))

  ggsave(file.path(CONFIG$output_dir, "145_LP_Sectoral_Output_Heatmap.png"), p145,
         width = 8, height = 8, dpi = 300)
  cat("Saved: 145_LP_Sectoral_Output_Heatmap.png\n")
}

# Export Part B CSVs
write.csv(all_output_full, file.path(CONFIG$output_dir, "LP_Sectoral_Output_Full.csv"),
          row.names = FALSE)
write.csv(all_output_postIT, file.path(CONFIG$output_dir, "LP_Sectoral_Output_PostIT.csv"),
          row.names = FALSE)
cat("Saved: LP_Sectoral_Output_Full.csv, LP_Sectoral_Output_PostIT.csv\n")


################################################################################
# 7b. PART B ASYMMETRIC: SECTORAL OUTPUT (Tightening vs Easing)
################################################################################

cat("\n================================================================================\n")
cat("PART B ASYMMETRIC: SECTORAL OUTPUT (Tightening vs Easing)\n")
cat("================================================================================\n\n")

# --- B1 Asymmetric: Full sample ---
cat("--- B1 Asymmetric: Full Sample (quarterly) ---\n\n")

output_full_asym <- list()

for (sector in CONFIG$activity_sectors) {
  y_var <- paste0("AG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full)) next

  data_sub <- analysis_full %>% filter(!is.na(!!sym(y_var)))

  cat(sprintf("  %s (n=%d) ... ", ACTIVITY_LABELS[sector], nrow(data_sub)))

  res <- run_lp_asymmetric(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy")
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- ACTIVITY_LABELS[sector]
    output_full_asym[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_output_full_asym <- bind_rows(output_full_asym)

# --- B2 Asymmetric: Post-IT sample ---
cat("\n--- B2 Asymmetric: Post-IT Sample (2011Q2+, quarterly) ---\n\n")

output_postIT_asym <- list()

for (sector in CONFIG$activity_sectors) {
  y_var <- paste0("AG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full)) next

  data_sub <- analysis_full %>%
    filter(fecha >= as.Date("2011-06-01"), !is.na(!!sym(y_var)))

  cat(sprintf("  %s (n=%d) ... ", ACTIVITY_LABELS[sector], nrow(data_sub)))

  res <- run_lp_asymmetric(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy")
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- ACTIVITY_LABELS[sector]
    output_postIT_asym[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_output_postIT_asym <- bind_rows(output_postIT_asym)

# --- Asymmetric output summary ---
cat("\nASYMMETRIC OUTPUT SUMMARY at h=4 (4 quarters):\n")
cat(sprintf("%-25s %10s %10s %8s  %10s %10s %8s\n",
            "Sector", "Full:Tight", "Full:Ease", "Asym.p",
            "PIT:Tight", "PIT:Ease", "Asym.p"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (sector in CONFIG$activity_sectors) {
  vf <- all_output_full_asym %>% filter(sector == !!sector, horizon == 4)
  vp <- all_output_postIT_asym %>% filter(sector == !!sector, horizon == 4)

  ft <- if(nrow(vf)>0) sprintf("%+.2f", vf$coef_tight) else "NA"
  fe <- if(nrow(vf)>0) sprintf("%+.2f", vf$coef_ease) else "NA"
  fp <- if(nrow(vf)>0) {
    asym_star <- ifelse(vf$asym_p < 0.05, "**", ifelse(vf$asym_p < 0.10, "*", ""))
    sprintf("%.3f%s", vf$asym_p, asym_star)
  } else "NA"
  pt <- if(nrow(vp)>0) sprintf("%+.2f", vp$coef_tight) else "NA"
  pe <- if(nrow(vp)>0) sprintf("%+.2f", vp$coef_ease) else "NA"
  pp <- if(nrow(vp)>0) {
    asym_star <- ifelse(vp$asym_p < 0.05, "**", ifelse(vp$asym_p < 0.10, "*", ""))
    sprintf("%.3f%s", vp$asym_p, asym_star)
  } else "NA"

  cat(sprintf("%-25s %10s %10s %8s  %10s %10s %8s\n",
              ACTIVITY_LABELS[sector], ft, fe, fp, pt, pe, pp))
}

# --- Asymmetric output visualizations ---
cat("\n\nGenerating asymmetric output visualizations...\n")

colors_asym_out <- c("Tightening" = "#E74C3C", "Easing" = "#3498DB")

# Plot 157: Asymmetric output -- full sample
if (nrow(all_output_full_asym) > 0) {
  asym_out_full_plot <- all_output_full_asym %>%
    pivot_longer(cols = c(coef_tight, coef_ease),
                 names_to = "effect_type", values_to = "coef") %>%
    mutate(
      effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
      ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
      ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
    )

  p157 <- ggplot(asym_out_full_plot, aes(x = horizon, y = coef,
                                          color = effect_label, fill = effect_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = colors_asym_out) +
    scale_fill_manual(values = colors_asym_out) +
    theme_minimal(base_size = 10) +
    labs(title = "Asymmetric Sectoral Output Response (Full Sample, Quarterly)",
         subtitle = "Tightening (FCI>0) vs Easing (|FCI<0|) | FCI_COMP | 1994Q1+ | 90% CI",
         x = "Quarters", y = "Effect on output growth (pp)",
         color = "FCI State", fill = "FCI State") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "157_LP_Sectoral_Output_Asym_Full.png"), p157,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 157_LP_Sectoral_Output_Asym_Full.png\n")
}

# Plot 158: Asymmetric output -- post-IT
if (nrow(all_output_postIT_asym) > 0) {
  asym_out_pit_plot <- all_output_postIT_asym %>%
    pivot_longer(cols = c(coef_tight, coef_ease),
                 names_to = "effect_type", values_to = "coef") %>%
    mutate(
      effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
      ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
      ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
    )

  p158 <- ggplot(asym_out_pit_plot, aes(x = horizon, y = coef,
                                          color = effect_label, fill = effect_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = colors_asym_out) +
    scale_fill_manual(values = colors_asym_out) +
    theme_minimal(base_size = 10) +
    labs(title = "Asymmetric Sectoral Output Response (Post-IT, Quarterly)",
         subtitle = "Tightening (FCI>0) vs Easing (|FCI<0|) | FCI_COMP | 2011Q2+ | 90% CI",
         x = "Quarters", y = "Effect on output growth (pp)",
         color = "FCI State", fill = "FCI State") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "158_LP_Sectoral_Output_Asym_PostIT.png"), p158,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 158_LP_Sectoral_Output_Asym_PostIT.png\n")
}

# Plot 159: Asymmetric output heatmap at h=4 (4 quarters)
hm_asym_out_full <- all_output_full_asym %>%
  filter(horizon == 4) %>% mutate(sample = "Full (1994Q1+)")
hm_asym_out_pit <- all_output_postIT_asym %>%
  filter(horizon == 4) %>% mutate(sample = "Post-IT (2011Q2+)")

hm_asym_output <- bind_rows(hm_asym_out_full, hm_asym_out_pit) %>%
  pivot_longer(cols = c(coef_tight, coef_ease),
               names_to = "effect_type", values_to = "coef_val") %>%
  mutate(
    p_val = ifelse(effect_type == "coef_tight", p_tight, p_ease),
    effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
    label = sprintf("%.1f%s", coef_val,
                    ifelse(p_val < 0.01, "***",
                           ifelse(p_val < 0.05, "**",
                                  ifelse(p_val < 0.10, "*", ""))))
  )

if (nrow(hm_asym_output) > 0) {
  p159 <- ggplot(hm_asym_output,
                 aes(x = interaction(effect_label, sample, sep = "\n"),
                     y = sector_label, fill = coef_val)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "Effect (pp)") +
    theme_minimal(base_size = 11) +
    labs(title = "Asymmetric Sectoral Output Effects at h=4 (4 Quarters)",
         subtitle = "Tightening vs Easing | *** p<0.01, ** p<0.05, * p<0.10",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(size = 8))

  ggsave(file.path(CONFIG$output_dir, "159_LP_Sectoral_Output_Asym_Heatmap.png"), p159,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 159_LP_Sectoral_Output_Asym_Heatmap.png\n")
}

# Export asymmetric output CSVs
write.csv(all_output_full_asym,
          file.path(CONFIG$output_dir, "LP_Sectoral_Output_Asym_Full.csv"), row.names = FALSE)
write.csv(all_output_postIT_asym,
          file.path(CONFIG$output_dir, "LP_Sectoral_Output_Asym_PostIT.csv"), row.names = FALSE)
cat("Saved: LP_Sectoral_Output_Asym_Full.csv, LP_Sectoral_Output_Asym_PostIT.csv\n")


################################################################################
# 7c. PART B ENDO/EXO: SECTORAL OUTPUT (Domestic vs External FCI)
################################################################################

cat("\n================================================================================\n")
cat("PART B ENDO/EXO: SECTORAL OUTPUT (Domestic vs External FCI)\n")
cat("================================================================================\n\n")

# For output-LHS: use FCI_ENDO and FCI_EXO directly (no endogeneity issue)
fci_output_types <- c("FCI_ENDO" = "Endogenous (Domestic)",
                       "FCI_EXO" = "Exogenous (External)")

# --- B1 Endo/Exo: Full sample ---
cat("--- B1 Endo/Exo: Full Sample (quarterly) ---\n\n")

output_full_endoexo <- list()

for (sector in CONFIG$activity_sectors) {
  y_var <- paste0("AG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full)) next

  data_sub <- analysis_full %>% filter(!is.na(!!sym(y_var)))

  for (fci_var in names(fci_output_types)) {
    cat(sprintf("  %s x %s (n=%d) ... ",
                ACTIVITY_LABELS[sector], fci_output_types[fci_var], nrow(data_sub)))

    res <- run_lp_standard(
      data_sub, y_var, fci_var, CONFIG$max_horizon, CONFIG$n_lags,
      control_vars = c("Cred_Real_yoy", "IPC_yoy")
    )

    if (nrow(res) > 0) {
      res$sector <- sector
      res$sector_label <- ACTIVITY_LABELS[sector]
      res$fci_var <- fci_var
      res$fci_label <- fci_output_types[fci_var]
      output_full_endoexo[[paste(sector, fci_var)]] <- res
      cat("done\n")
    } else {
      cat("insufficient data\n")
    }
  }
}

all_output_full_endoexo <- bind_rows(output_full_endoexo)

# --- B2 Endo/Exo: Post-IT sample ---
cat("\n--- B2 Endo/Exo: Post-IT Sample (2011Q2+, quarterly) ---\n\n")

output_postIT_endoexo <- list()

for (sector in CONFIG$activity_sectors) {
  y_var <- paste0("AG_", gsub(" ", "_", sector))
  if (!y_var %in% names(analysis_full)) next

  data_sub <- analysis_full %>%
    filter(fecha >= as.Date("2011-06-01"), !is.na(!!sym(y_var)))

  for (fci_var in names(fci_output_types)) {
    cat(sprintf("  %s x %s (n=%d) ... ",
                ACTIVITY_LABELS[sector], fci_output_types[fci_var], nrow(data_sub)))

    res <- run_lp_standard(
      data_sub, y_var, fci_var, CONFIG$max_horizon, CONFIG$n_lags,
      control_vars = c("Cred_Real_yoy", "IPC_yoy")
    )

    if (nrow(res) > 0) {
      res$sector <- sector
      res$sector_label <- ACTIVITY_LABELS[sector]
      res$fci_var <- fci_var
      res$fci_label <- fci_output_types[fci_var]
      output_postIT_endoexo[[paste(sector, fci_var)]] <- res
      cat("done\n")
    } else {
      cat("insufficient data\n")
    }
  }
}

all_output_postIT_endoexo <- bind_rows(output_postIT_endoexo)

# --- Endo/Exo output summary ---
cat("\nENDO/EXO OUTPUT SUMMARY at h=2, h=4 (quarters):\n")
cat(sprintf("%-25s %10s %10s %10s %10s  %s\n",
            "Sector", "Endo h=2", "Exo h=2", "Endo h=4", "Exo h=4", "Sample"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (sector in CONFIG$activity_sectors) {
  ve2 <- all_output_full_endoexo %>% filter(sector == !!sector, horizon == 2, fci_var == "FCI_ENDO")
  vx2 <- all_output_full_endoexo %>% filter(sector == !!sector, horizon == 2, fci_var == "FCI_EXO")
  ve4 <- all_output_full_endoexo %>% filter(sector == !!sector, horizon == 4, fci_var == "FCI_ENDO")
  vx4 <- all_output_full_endoexo %>% filter(sector == !!sector, horizon == 4, fci_var == "FCI_EXO")

  fmt <- function(v) {
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      sprintf("%+.2f%s", v$coef, stars)
    } else "NA"
  }

  cat(sprintf("%-25s %10s %10s %10s %10s  Full\n",
              ACTIVITY_LABELS[sector], fmt(ve2), fmt(vx2), fmt(ve4), fmt(vx4)))
}
cat("\n")
for (sector in CONFIG$activity_sectors) {
  ve2 <- all_output_postIT_endoexo %>% filter(sector == !!sector, horizon == 2, fci_var == "FCI_ENDO")
  vx2 <- all_output_postIT_endoexo %>% filter(sector == !!sector, horizon == 2, fci_var == "FCI_EXO")
  ve4 <- all_output_postIT_endoexo %>% filter(sector == !!sector, horizon == 4, fci_var == "FCI_ENDO")
  vx4 <- all_output_postIT_endoexo %>% filter(sector == !!sector, horizon == 4, fci_var == "FCI_EXO")

  fmt <- function(v) {
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      sprintf("%+.2f%s", v$coef, stars)
    } else "NA"
  }

  cat(sprintf("%-25s %10s %10s %10s %10s  Post-IT\n",
              ACTIVITY_LABELS[sector], fmt(ve2), fmt(vx2), fmt(ve4), fmt(vx4)))
}

# --- Endo/Exo output visualizations ---
cat("\n\nGenerating endo/exo output visualizations...\n")

# Plot 163: Endo vs Exo output -- full sample
if (nrow(all_output_full_endoexo) > 0) {
  p163 <- ggplot(all_output_full_endoexo,
                 aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = colors_endo_exo) +
    scale_fill_manual(values = colors_endo_exo) +
    theme_minimal(base_size = 10) +
    labs(title = "Sectoral Output: Endogenous vs Exogenous FCI (Full Sample, Quarterly)",
         subtitle = "FCI_ENDO vs FCI_EXO → Output growth (YoY) | 1994Q1+ | 90% CI",
         x = "Quarters", y = "Effect on output growth (pp)",
         color = "FCI Component", fill = "FCI Component") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "163_LP_Sectoral_Output_EndoExo_Full.png"), p163,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 163_LP_Sectoral_Output_EndoExo_Full.png\n")
}

# Plot 164: Endo vs Exo output -- post-IT
if (nrow(all_output_postIT_endoexo) > 0) {
  p164 <- ggplot(all_output_postIT_endoexo,
                 aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = colors_endo_exo) +
    scale_fill_manual(values = colors_endo_exo) +
    theme_minimal(base_size = 10) +
    labs(title = "Sectoral Output: Endogenous vs Exogenous FCI (Post-IT, Quarterly)",
         subtitle = "FCI_ENDO vs FCI_EXO → Output growth (YoY) | 2011Q2+ | 90% CI",
         x = "Quarters", y = "Effect on output growth (pp)",
         color = "FCI Component", fill = "FCI Component") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "164_LP_Sectoral_Output_EndoExo_PostIT.png"), p164,
         width = 14, height = 10, dpi = 300)
  cat("Saved: 164_LP_Sectoral_Output_EndoExo_PostIT.png\n")
}

# Plot 165: Endo/Exo output heatmap at h=2 and h=4 (quarters)
hm_endoexo_out_full <- all_output_full_endoexo %>%
  filter(horizon %in% c(2, 4)) %>% mutate(sample = "Full (1994Q1+)")
hm_endoexo_out_pit <- all_output_postIT_endoexo %>%
  filter(horizon %in% c(2, 4)) %>% mutate(sample = "Post-IT (2011Q2+)")

hm_endoexo_output <- bind_rows(hm_endoexo_out_full, hm_endoexo_out_pit) %>%
  mutate(
    horizon_label = paste0("h=", horizon),
    label = sprintf("%.1f%s", coef,
                    ifelse(p_value < 0.01, "***",
                           ifelse(p_value < 0.05, "**",
                                  ifelse(p_value < 0.10, "*", ""))))
  )

if (nrow(hm_endoexo_output) > 0) {
  p165 <- ggplot(hm_endoexo_output,
                 aes(x = interaction(fci_label, horizon_label, sample, sep = "\n"),
                     y = sector_label, fill = coef)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 2.8, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "Effect (pp)") +
    theme_minimal(base_size = 11) +
    labs(title = "Sectoral Output: Endo vs Exo FCI Effects at h=2 and h=4 (Quarters)",
         subtitle = "Effect in pp | *** p<0.01, ** p<0.05, * p<0.10",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(size = 7, angle = 0))

  ggsave(file.path(CONFIG$output_dir, "165_LP_Sectoral_Output_EndoExo_Heatmap.png"), p165,
         width = 16, height = 8, dpi = 300)
  cat("Saved: 165_LP_Sectoral_Output_EndoExo_Heatmap.png\n")
}

# Export Endo/Exo output CSVs
write.csv(all_output_full_endoexo,
          file.path(CONFIG$output_dir, "LP_Sectoral_Output_EndoExo_Full.csv"), row.names = FALSE)
write.csv(all_output_postIT_endoexo,
          file.path(CONFIG$output_dir, "LP_Sectoral_Output_EndoExo_PostIT.csv"), row.names = FALSE)
cat("Saved: LP_Sectoral_Output_EndoExo_Full.csv, LP_Sectoral_Output_EndoExo_PostIT.csv\n")


################################################################################
# 8. PART C: EXTENDED HORIZON REPLICATION (h = 1..36)
################################################################################

cat("\n================================================================================\n")
cat("PART C: EXTENDED HORIZON LP (h = 1..36)\n")
cat("================================================================================\n\n")

# Prepare the same analysis_data as script 05 for consistency
# Load credit data for extended LP
datos_raw_05 <- read_excel(CONFIG$fci_file)
fecha_col_05 <- names(datos_raw_05)[grepl("fecha|date", names(datos_raw_05), ignore.case = TRUE)][1]

ipc_data_05 <- macro_data %>%
  dplyr::select(fecha, IPC, Creditos_deflactados)

credit_data_05 <- datos_raw_05 %>%
  rename(fecha = !!sym(fecha_col_05)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  left_join(ipc_data_05, by = "fecha") %>%
  mutate(
    Cred_Total = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100,
    Cred_MN = (Creditos_Sector_privado_MN / lag(Creditos_Sector_privado_MN, 12) - 1) * 100,
    Cred_USD = (Creditos_Sector_privado_USD_equivalente /
                lag(Creditos_Sector_privado_USD_equivalente, 12) - 1) * 100,
    Cred_Real_MN = ((Creditos_Sector_privado_MN / IPC) /
                     lag(Creditos_Sector_privado_MN / IPC, 12) - 1) * 100
  ) %>%
  dplyr::select(fecha, Cred_Total, Cred_MN, Cred_USD, Cred_Real_MN)

macro_data_05 <- macro_data %>%
  left_join(credit_data_05 %>% dplyr::select(fecha, Cred_Real_yoy_05 = Cred_Total), by = "fecha") %>%
  dplyr::select(fecha, IMAEP_yoy, IPC_yoy, Cred_Real_yoy = Cred_Real_yoy_05)

ext_data <- fci_monthly %>%
  inner_join(credit_data_05, by = "fecha") %>%
  inner_join(macro_data_05, by = "fecha") %>%
  arrange(fecha) %>%
  mutate(
    FCI_COMP_pos = pmax(FCI_COMP, 0),
    FCI_COMP_neg = abs(pmin(FCI_COMP, 0)),
    FCI_ENDO_pos = pmax(FCI_ENDO, 0),
    FCI_ENDO_neg = abs(pmin(FCI_ENDO, 0)),
    FCI_EXO_pos = pmax(FCI_EXO, 0),
    FCI_EXO_neg = abs(pmin(FCI_EXO, 0)),
    FCI_exCredit_pos = pmax(FCI_exCredit, 0),
    FCI_exCredit_neg = abs(pmin(FCI_exCredit, 0)),
    FCI_ENDO_exCredit_pos = pmax(FCI_ENDO_exCredit, 0),
    FCI_ENDO_exCredit_neg = abs(pmin(FCI_ENDO_exCredit, 0))
  ) %>%
  na.omit()

cat(sprintf("Extended LP data: %d observations, %s to %s\n\n",
            nrow(ext_data), format(min(ext_data$fecha)), format(max(ext_data$fecha))))

H <- CONFIG$extended_horizon

# --- C1: Macro LP (IMAEP, IPC, Credit) ---
cat("C1: Macro LP (h=1..36) ...\n")

ext_macro_std <- list()

# IMAEP
ext_macro_std[["IMAEP"]] <- run_lp_standard(
  ext_data, "IMAEP_yoy", "FCI_COMP", H, CONFIG$n_lags,
  control_vars = c("IPC_yoy", "Cred_Real_yoy")
) %>% mutate(variable = "IMAEP")

# IPC
ext_macro_std[["IPC"]] <- run_lp_standard(
  ext_data, "IPC_yoy", "FCI_COMP", H, CONFIG$n_lags,
  control_vars = c("IMAEP_yoy", "Cred_Real_yoy")
) %>% mutate(variable = "IPC")

# Credit (with FCI_exCredit)
ext_macro_std[["Creditos"]] <- run_lp_standard(
  ext_data, "Cred_Total", "FCI_exCredit", H, CONFIG$n_lags,
  control_vars = c("IMAEP_yoy", "IPC_yoy")
) %>% mutate(variable = "Creditos")

all_ext_macro <- bind_rows(ext_macro_std)
cat("  done\n")

# Plot 146
p146 <- all_ext_macro %>%
  ggplot(aes(x = horizon)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#377EB8", alpha = 0.2) +
  geom_line(aes(y = coef), color = "#377EB8", linewidth = 1) +
  geom_point(aes(y = coef), color = "#377EB8", size = 1.5) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  theme_minimal(base_size = 11) +
  labs(title = "Extended LP: FCI Effect on Macro Variables (h=1..36)",
       subtitle = "Response to 1 SD FCI shock | 90% CI",
       x = "Months", y = "Effect (pp)") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(CONFIG$output_dir, "146_LP_Macro_Extended36.png"), p146,
       width = 14, height = 5, dpi = 300)
cat("Saved: 146_LP_Macro_Extended36.png\n")

# --- C2: Credit channel decomposition ---
cat("C2: Credit channel decomposition (h=1..36) ...\n")

FCI_LABELS <- c("FCI_COMP" = "Comprehensive",
                "FCI_ENDO" = "Endogenous (Domestic)",
                "FCI_EXO" = "Exogenous (External)")

CREDIT_LABELS_05 <- c("Total" = "Real Total Credit",
                       "USD" = "USD (Dollarized)",
                       "Real_MN" = "Real MN Credit")

fci_types <- c("FCI_COMP", "FCI_ENDO", "FCI_EXO")
credit_types <- c("Total", "USD", "Real_MN")

fci_credit_swap <- c("FCI_COMP" = "FCI_exCredit",
                     "FCI_ENDO" = "FCI_ENDO_exCredit",
                     "FCI_EXO" = "FCI_EXO")

ext_credit_std <- list()

for (fci in fci_types) {
  for (cred in credit_types) {
    y_var <- paste0("Cred_", cred)
    key <- paste(fci, cred, sep = "_")
    fci_actual <- fci_credit_swap[fci]

    ext_credit_std[[key]] <- run_lp_standard(
      ext_data, y_var, fci_actual, H, CONFIG$n_lags,
      control_vars = c("IMAEP_yoy", "IPC_yoy")
    ) %>% mutate(fci_type = fci, credit_type = cred,
                 fci_label = FCI_LABELS[fci], credit_label = CREDIT_LABELS_05[cred])
  }
}

all_ext_credit <- bind_rows(ext_credit_std)
cat("  done\n")

colors_fci <- c("Comprehensive" = "#2C3E50",
                "Endogenous (Domestic)" = "#E74C3C",
                "Exogenous (External)" = "#3498DB")

colors_credit <- c("Real Total Credit" = "#2C3E50",
                   "USD (Dollarized)" = "#8E44AD",
                   "Real MN Credit" = "#27AE60")

colors_asym <- c("Tightening" = "#E74C3C", "Easing" = "#3498DB")

# Plot 147: Credit by type (facet by credit, color by FCI)
p147 <- all_ext_credit %>%
  ggplot(aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~credit_label, ncol = 3) +
  scale_color_manual(values = colors_fci) +
  scale_fill_manual(values = colors_fci) +
  theme_minimal(base_size = 11) +
  labs(title = "Extended Credit Channel: FCI Types (h=1..36)",
       subtitle = "Response of credit growth (YoY) to 1 SD FCI shock | 90% CI",
       x = "Months", y = "Effect on credit growth (pp)",
       color = "FCI Type", fill = "FCI Type") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(CONFIG$output_dir, "147_LP_Credit_byType_Extended36.png"), p147,
       width = 14, height = 5, dpi = 300)
cat("Saved: 147_LP_Credit_byType_Extended36.png\n")

# Plot 148: Credit by FCI (facet by FCI, color by credit)
p148 <- all_ext_credit %>%
  ggplot(aes(x = horizon, y = coef, color = credit_label, fill = credit_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~fci_label, ncol = 3) +
  scale_color_manual(values = colors_credit) +
  scale_fill_manual(values = colors_credit) +
  theme_minimal(base_size = 11) +
  labs(title = "Extended Credit Channel: Credit Types (h=1..36)",
       subtitle = "Comparing Real Total, USD, and Real MN credit | 90% CI",
       x = "Months", y = "Effect on credit growth (pp)",
       color = "Credit Type", fill = "Credit Type") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(CONFIG$output_dir, "148_LP_Credit_byFCI_Extended36.png"), p148,
       width = 14, height = 5, dpi = 300)
cat("Saved: 148_LP_Credit_byFCI_Extended36.png\n")

# Plot 149: Heatmap at h=12, 24, 36
heatmap_ext <- all_ext_credit %>%
  filter(horizon %in% c(12, 24, 36)) %>%
  mutate(
    horizon_label = paste0("h=", horizon),
    label = sprintf("%.1f%s", coef,
                    ifelse(p_value < 0.01, "***",
                           ifelse(p_value < 0.05, "**",
                                  ifelse(p_value < 0.10, "*", ""))))
  )

if (nrow(heatmap_ext) > 0) {
  p149 <- ggplot(heatmap_ext,
                 aes(x = interaction(credit_label, horizon_label, sep = "\n"),
                     y = fci_label, fill = coef)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "Effect (pp)") +
    theme_minimal(base_size = 11) +
    labs(title = "Extended Credit Heatmap (h=12, 24, 36)",
         subtitle = "Effect in pp | *** p<0.01, ** p<0.05, * p<0.10",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(size = 8))

  ggsave(file.path(CONFIG$output_dir, "149_LP_Credit_Heatmap_Extended36.png"), p149,
         width = 14, height = 6, dpi = 300)
  cat("Saved: 149_LP_Credit_Heatmap_Extended36.png\n")
}

# --- C3: Asymmetric LP ---
cat("C3: Asymmetric LP (h=1..36) ...\n")

# Endogenous asymmetric
ext_asym_endo <- list()
for (cred in credit_types) {
  y_var <- paste0("Cred_", cred)
  ext_asym_endo[[cred]] <- run_lp_asymmetric(
    ext_data, y_var, "FCI_ENDO_exCredit", H, CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  ) %>% mutate(credit_type = cred, credit_label = CREDIT_LABELS_05[cred])
}

all_ext_asym_endo <- bind_rows(ext_asym_endo)

# Exogenous asymmetric
ext_asym_exo <- list()
for (cred in credit_types) {
  y_var <- paste0("Cred_", cred)
  ext_asym_exo[[cred]] <- run_lp_asymmetric(
    ext_data, y_var, "FCI_EXO", H, CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  ) %>% mutate(credit_type = cred, credit_label = CREDIT_LABELS_05[cred])
}

all_ext_asym_exo <- bind_rows(ext_asym_exo)
cat("  done\n")

# Plot 150: Asymmetric Endo
asym_endo_plot <- all_ext_asym_endo %>%
  pivot_longer(cols = c(coef_tight, coef_ease),
               names_to = "effect_type", values_to = "coef") %>%
  mutate(
    effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
    ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
    ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
  )

p150 <- ggplot(asym_endo_plot, aes(x = horizon, y = coef,
                                    color = effect_label, fill = effect_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~credit_label, ncol = 3) +
  scale_color_manual(values = colors_asym) +
  scale_fill_manual(values = colors_asym) +
  theme_minimal(base_size = 11) +
  labs(title = "Extended Asymmetric: Endogenous FCI on Credit (h=1..36)",
       subtitle = "Tightening vs easing | 90% CI",
       x = "Months", y = "Effect on credit growth (pp)",
       color = "FCI State", fill = "FCI State") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(CONFIG$output_dir, "150_LP_Asymmetric_Endo_Extended36.png"), p150,
       width = 14, height = 5, dpi = 300)
cat("Saved: 150_LP_Asymmetric_Endo_Extended36.png\n")

# Plot 151: Asymmetric Exo
asym_exo_plot <- all_ext_asym_exo %>%
  pivot_longer(cols = c(coef_tight, coef_ease),
               names_to = "effect_type", values_to = "coef") %>%
  mutate(
    effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
    ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
    ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
  )

p151 <- ggplot(asym_exo_plot, aes(x = horizon, y = coef,
                                    color = effect_label, fill = effect_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~credit_label, ncol = 3) +
  scale_color_manual(values = colors_asym) +
  scale_fill_manual(values = colors_asym) +
  theme_minimal(base_size = 11) +
  labs(title = "Extended Asymmetric: Exogenous FCI on Credit (h=1..36)",
       subtitle = "Tightening vs easing of external conditions | 90% CI",
       x = "Months", y = "Effect on credit growth (pp)",
       color = "FCI State", fill = "FCI State") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(CONFIG$output_dir, "151_LP_Asymmetric_Exo_Extended36.png"), p151,
       width = 14, height = 5, dpi = 300)
cat("Saved: 151_LP_Asymmetric_Exo_Extended36.png\n")

# --- C4: NPL extended ---
cat("C4: NPL extended LP (h=1..36) ...\n")

# Load NPL data
datos_npl <- read_excel(CONFIG$fci_file)
fecha_col_npl <- names(datos_npl)[grepl("fecha|date", names(datos_npl), ignore.case = TRUE)][1]

npl_data <- datos_npl %>%
  rename(fecha = !!sym(fecha_col_npl)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  dplyr::select(fecha, Morosidad) %>%
  mutate(
    NPL_level = Morosidad,
    NPL_change = Morosidad - lag(Morosidad, 12)
  ) %>%
  na.omit()

npl_ext <- fci_monthly %>%
  inner_join(npl_data, by = "fecha") %>%
  inner_join(macro_data_05, by = "fecha") %>%
  mutate(
    FCI_COMP_pos = pmax(FCI_COMP, 0),
    FCI_COMP_neg = abs(pmin(FCI_COMP, 0)),
    FCI_ENDO_pos = pmax(FCI_ENDO, 0),
    FCI_ENDO_neg = abs(pmin(FCI_ENDO, 0)),
    FCI_EXO_pos = pmax(FCI_EXO, 0),
    FCI_EXO_neg = abs(pmin(FCI_EXO, 0))
  ) %>%
  na.omit()

npl_ext_results <- list()
npl_ext_asym <- list()

for (fci in fci_types) {
  npl_ext_results[[fci]] <- run_lp_standard(
    npl_ext, "NPL_level", fci, H, CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  ) %>% mutate(fci_type = fci, fci_label = FCI_LABELS[fci])

  npl_ext_asym[[fci]] <- run_lp_asymmetric(
    npl_ext, "NPL_level", fci, H, CONFIG$n_lags,
    control_vars = c("IMAEP_yoy", "IPC_yoy")
  ) %>% mutate(fci_type = fci, fci_label = FCI_LABELS[fci])
}

all_npl_ext <- bind_rows(npl_ext_results)
all_npl_ext_asym <- bind_rows(npl_ext_asym)
cat("  done\n")

# Plot 152: NPL by FCI type
p152 <- all_npl_ext %>%
  ggplot(aes(x = horizon, y = coef, color = fci_label, fill = fci_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = colors_fci) +
  scale_fill_manual(values = colors_fci) +
  theme_minimal(base_size = 11) +
  labs(title = "Extended NPL Response to FCI (h=1..36)",
       subtitle = "Response of NPL ratio to 1 SD FCI shock | 90% CI",
       x = "Months", y = "Effect on NPL (pp)",
       color = "FCI Type", fill = "FCI Type") +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(CONFIG$output_dir, "152_LP_NPL_Extended36.png"), p152,
       width = 12, height = 6, dpi = 300)
cat("Saved: 152_LP_NPL_Extended36.png\n")

# Plot 153: NPL asymmetric
npl_asym_plot <- all_npl_ext_asym %>%
  pivot_longer(cols = c(coef_tight, coef_ease),
               names_to = "effect_type", values_to = "coef") %>%
  mutate(
    effect_label = ifelse(effect_type == "coef_tight", "Tightening", "Easing"),
    ci_lo = ifelse(effect_type == "coef_tight", ci_tight_lo, ci_ease_lo),
    ci_hi = ifelse(effect_type == "coef_tight", ci_tight_hi, ci_ease_hi)
  )

p153 <- ggplot(npl_asym_plot, aes(x = horizon, y = coef,
                                    color = effect_label, fill = effect_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~fci_label, ncol = 3) +
  scale_color_manual(values = colors_asym) +
  scale_fill_manual(values = colors_asym) +
  theme_minimal(base_size = 11) +
  labs(title = "Extended Asymmetric FCI Effects on NPL (h=1..36)",
       subtitle = "Tightening vs easing | 90% CI",
       x = "Months", y = "Effect on NPL (pp)",
       color = "FCI State", fill = "FCI State") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(CONFIG$output_dir, "153_LP_NPL_Asymmetric_Extended36.png"), p153,
       width = 14, height = 5, dpi = 300)
cat("Saved: 153_LP_NPL_Asymmetric_Extended36.png\n")

# Export Part C CSVs
write.csv(all_ext_macro, file.path(CONFIG$output_dir, "LP_Macro_Extended36.csv"),
          row.names = FALSE)
write.csv(all_ext_credit, file.path(CONFIG$output_dir, "LP_Credit_Extended36.csv"),
          row.names = FALSE)
write.csv(all_ext_asym_endo, file.path(CONFIG$output_dir, "LP_Asym_Endo_Extended36.csv"),
          row.names = FALSE)
write.csv(all_ext_asym_exo, file.path(CONFIG$output_dir, "LP_Asym_Exo_Extended36.csv"),
          row.names = FALSE)
write.csv(all_npl_ext, file.path(CONFIG$output_dir, "LP_NPL_Extended36.csv"),
          row.names = FALSE)
write.csv(all_npl_ext_asym, file.path(CONFIG$output_dir, "LP_NPL_Asym_Extended36.csv"),
          row.names = FALSE)
cat("Saved: Extended LP CSVs\n")


################################################################################
# 9. PART D: SUMMARY REPORT
################################################################################

cat("\n================================================================================\n")
cat("PART D: GENERATING SUMMARY REPORT\n")
cat("================================================================================\n\n")

# Collect key results for report
report_lines <- c(
  "# Sectoral Output Puzzle Investigation",
  "",
  paste0("**Generated**: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. Motivation",
  "",
  "The FCI strongly predicts credit growth but not aggregate output (IMAEP) --- the",
  "\"output puzzle.\" This analysis uses disaggregated sectoral data to test whether",
  "specific sectors respond to financial conditions, potentially resolving the puzzle.",
  "",
  "## 2. Part A: Sectoral Credit Response",
  "",
  "### A1: Long Sample (Dec 2001+, 5 sectors)",
  "",
  "| Sector | h=6 | h=12 | h=18 | h=24 |",
  "|--------|-----|------|------|------|"
)

for (sector in CONFIG$credit_long) {
  vals <- c()
  for (h in c(6, 12, 18, 24)) {
    v <- all_credit_long %>% filter(sector == !!sector, horizon == h)
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      vals <- c(vals, sprintf("%+.1f%s", v$coef, stars))
    } else {
      vals <- c(vals, "---")
    }
  }
  report_lines <- c(report_lines,
    sprintf("| %s | %s | %s | %s | %s |",
            CREDIT_SHORT_LABELS[sector], vals[1], vals[2], vals[3], vals[4]))
}

report_lines <- c(report_lines, "",
  "*FCI_exCredit used to address endogeneity. Controls: IMAEP_yoy, IPC_yoy. Newey-West HAC SE.*",
  "",
  "### A2: Short Sample (Jan 2016+, 12 sectors)",
  "",
  "| Sector | h=6 | h=12 | h=18 | h=24 |",
  "|--------|-----|------|------|------|"
)

for (sector in CONFIG$credit_short) {
  vals <- c()
  for (h in c(6, 12, 18, 24)) {
    v <- all_credit_short %>% filter(sector == !!sector, horizon == h)
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      vals <- c(vals, sprintf("%+.1f%s", v$coef, stars))
    } else {
      vals <- c(vals, "---")
    }
  }
  report_lines <- c(report_lines,
    sprintf("| %s | %s | %s | %s | %s |",
            CREDIT_SHORT_LABELS[sector], vals[1], vals[2], vals[3], vals[4]))
}

report_lines <- c(report_lines, "",
  "**Figures**: 140 (long sample IRFs), 141 (short sample IRFs), 142 (heatmap)",
  "",
  "## 3. Part B: Sectoral Output Response (Quarterly LP)",
  "",
  "### B1: Full Sample (1994Q1+)",
  "",
  "| Sector | h=2q | h=4q | h=6q | h=8q |",
  "|--------|------|------|------|------|"
)

for (sector in CONFIG$activity_sectors) {
  vals <- c()
  for (h in c(2, 4, 6, 8)) {
    v <- all_output_full %>% filter(sector == !!sector, horizon == h)
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      vals <- c(vals, sprintf("%+.1f%s", v$coef, stars))
    } else {
      vals <- c(vals, "---")
    }
  }
  report_lines <- c(report_lines,
    sprintf("| %s | %s | %s | %s | %s |",
            ACTIVITY_LABELS[sector], vals[1], vals[2], vals[3], vals[4]))
}

report_lines <- c(report_lines, "",
  "### B2: Post-IT Sample (2011Q2+)",
  "",
  "| Sector | h=2q | h=4q | h=6q | h=8q |",
  "|--------|------|------|------|------|"
)

for (sector in CONFIG$activity_sectors) {
  vals <- c()
  for (h in c(2, 4, 6, 8)) {
    v <- all_output_postIT %>% filter(sector == !!sector, horizon == h)
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      vals <- c(vals, sprintf("%+.1f%s", v$coef, stars))
    } else {
      vals <- c(vals, "---")
    }
  }
  report_lines <- c(report_lines,
    sprintf("| %s | %s | %s | %s | %s |",
            ACTIVITY_LABELS[sector], vals[1], vals[2], vals[3], vals[4]))
}

report_lines <- c(report_lines, "",
  "*FCI_COMP used. Controls: Cred_Real_yoy, IPC_yoy. Newey-West HAC SE.*",
  "",
  "**Figures**: 143 (full sample IRFs), 144 (post-IT IRFs), 145 (heatmap)",
  "",
  "## 4. Part C: Extended Horizon Results (h=1..36)",
  "",
  "### Macro variables at h=12, 24, 36",
  "",
  "| Variable | h=12 | h=24 | h=36 |",
  "|----------|------|------|------|"
)

for (var in c("IMAEP", "IPC", "Creditos")) {
  vals <- c()
  for (h in c(12, 24, 36)) {
    v <- all_ext_macro %>% filter(variable == var, horizon == h)
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      vals <- c(vals, sprintf("%+.2f%s", v$coef, stars))
    } else {
      vals <- c(vals, "---")
    }
  }
  report_lines <- c(report_lines,
    sprintf("| %s | %s | %s | %s |", var, vals[1], vals[2], vals[3]))
}

report_lines <- c(report_lines, "",
  "### Credit channel decomposition at h=12, 24, 36",
  "",
  "| FCI Type | Credit | h=12 | h=24 | h=36 |",
  "|----------|--------|------|------|------|"
)

for (fci in fci_types) {
  for (cred in credit_types) {
    vals <- c()
    for (h in c(12, 24, 36)) {
      v <- all_ext_credit %>% filter(fci_type == fci, credit_type == cred, horizon == h)
      if (nrow(v) > 0) {
        stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                        ifelse(v$p_value < 0.10, "*", "")))
        vals <- c(vals, sprintf("%+.2f%s", v$coef, stars))
      } else {
        vals <- c(vals, "---")
      }
    }
    report_lines <- c(report_lines,
      sprintf("| %s | %s | %s | %s | %s |",
              FCI_LABELS[fci], CREDIT_LABELS_05[cred], vals[1], vals[2], vals[3]))
  }
}

report_lines <- c(report_lines, "",
  "**Figures**: 146-153 (extended LP plots and CSVs)",
  "",
  "## 5. Conclusions",
  ""
)

# Identify which credit sectors are significantly affected
sig_credit <- all_credit_long %>%
  filter(horizon == 12, p_value < 0.10)

sig_output_full <- all_output_full %>%
  filter(horizon == 4, p_value < 0.10)

sig_output_pit <- all_output_postIT %>%
  filter(horizon == 4, p_value < 0.10)

report_lines <- c(report_lines,
  sprintf("**Credit sectors significantly affected at h=12 months (long sample)**: %s",
    if(nrow(sig_credit) > 0)
      paste(sig_credit$sector_label, collapse = ", ")
    else "None (or Part A skipped)"),
  "",
  sprintf("**Output sectors significantly affected at h=4 quarters (full sample)**: %s",
    if(nrow(sig_output_full) > 0)
      paste(sig_output_full$sector_label, collapse = ", ")
    else "None"),
  "",
  sprintf("**Output sectors significantly affected at h=4 quarters (post-IT)**: %s",
    if(nrow(sig_output_pit) > 0)
      paste(sig_output_pit$sector_label, collapse = ", ")
    else "None"),
  "",
  "### Key Findings",
  "",
  "1. **Credit channel**: Sectoral disaggregation reveals which credit segments",
  "   are most responsive to financial conditions tightening.",
  "",
  "2. **Output puzzle**: The sectoral analysis identifies whether specific sectors",
  "   (e.g., construction, manufacturing) respond to FCI even when aggregate",
  "   output does not.",
  "",
  "3. **Post-IT regime**: The post-IT sample tests whether the output puzzle",
  "   resolves once the monetary regime stabilizes.",
  "",
  "4. **Extended horizons**: The h=36 analysis confirms whether credit effects",
  "   persist or revert at longer horizons.",
  "",
  "### Policy Implications",
  "",
  "- Financial conditions monitoring should focus on sectors with significant FCI sensitivity",
  "- The output puzzle may reflect agricultural dominance masking the financial channel",
  "- Sectoral heterogeneity suggests differentiated macroprudential policy responses",
  "",
  "---",
  "",
  "*** p<0.01, ** p<0.05, * p<0.10 | All LP use Newey-West HAC(h+1) standard errors | 90% CI",
  ""
)

writeLines(report_lines, file.path(CONFIG$output_dir, "Output_Puzzle_Sectoral_Report.md"))
cat("Saved: Output_Puzzle_Sectoral_Report.md\n")


################################################################################
# 10. FINAL SUMMARY
################################################################################

cat("\n################################################################################\n")
cat("SECTORAL OUTPUT PUZZLE ANALYSIS COMPLETE\n")
cat("################################################################################\n\n")

cat("Outputs generated:\n")
cat("  PNGs:  140-165 series (Part A credit PNGs only if credit data available)\n")
cat("  CSVs:  LP_Sectoral_Output_Full/PostIT (quarterly),\n")
cat("         LP_Sectoral_Output_Asym_Full/PostIT, LP_Sectoral_Output_EndoExo_Full/PostIT,\n")
cat("         LP_Macro_Extended36, LP_Credit_Extended36, LP_Asym_Endo/Exo_Extended36,\n")
cat("         LP_NPL_Extended36, LP_NPL_Asym_Extended36\n")
if (credit_data_available) {
  cat("         LP_Sectoral_Credit_Long/Short, LP_Sectoral_Credit_Asym_Long/Short,\n")
  cat("         LP_Sectoral_Credit_EndoExo_Long/Short\n")
}
cat("  Report: Output_Puzzle_Sectoral_Report.md\n\n")
cat("Output directory:", CONFIG$output_dir, "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
