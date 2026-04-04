################################################################################
# FCI MONTHLY LP WITH PUBLISHED BCP DATA (SUPPLEMENTARY)
################################################################################
#
# Project:      Financial Conditions Index - Published Data Validation
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Monthly-frequency validation using broad sectors from
#               FCI_data_1.xlsx (Monthly_SA). Serves as robustness for
#               quarterly results and provides a purely post-IT check.
#
#   PART A: Data Preparation (merge monthly published with FCI)
#   PART B: Standard Monthly LP (h = 1...24)
#   PART C: Core Dual-Engine Test
#   PART D: Comparison with Original
#
# Data:
#   - FCI_data_1.xlsx: Monthly_SA (broad sectors, Jan 2014+)
#   - FCI_data_1.xlsx: FCI variables + macro controls
#
# References:
#   - Jordà (2005) AER - Local Projections
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
  library(zoo)
  library(lubridate)
})

CONFIG <- list(
  # Data files
  published_file = "../data/FCI_data_1.xlsx",
  fci_file       = "../data/FCI_data_1.xlsx",
  macro_sheet    = "Datos_macro",

  # LP parameters
  max_horizon = 24,
  horizons_report = c(3, 6, 12, 18, 24),
  n_lags = 2,
  confidence_level = 0.90,

  # Output
  output_dir = "../output"
)

# Monthly published sectors (from Monthly_SA sheet)
MONTHLY_SECTORS <- c(
  "Primario",
  "Secundario",
  "Manufacturero",
  "Servicios",
  "IMAEP",
  "IMAEP_SANB"
)

MONTHLY_LABELS <- c(
  "Primario"      = "Primary Sector",
  "Secundario"    = "Secondary Sector",
  "Manufacturero" = "Manufacturing",
  "Servicios"     = "Services",
  "IMAEP"         = "Aggregate GDP (IMAEP)",
  "IMAEP_SANB"    = "Non-Agro GDP (IMAEP ex-Agri)"
)

set.seed(20260320)

cat("\n################################################################################\n")
cat("FCI MONTHLY LP WITH PUBLISHED BCP DATA (SUPPLEMENTARY)\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


################################################################################
# 2. LP FUNCTIONS (standard and asymmetric, same as codebase)
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
  CONFIG_28 <- CONFIG
  source("01_FCI_Complete.R")
  CONFIG <- CONFIG_28
}

fci_monthly <- resultado_fci$all_indices %>%
  dplyr::select(fecha,
                FCI_COMP = FCI_COMP_AVG,
                FCI_ENDO = FCI_ENDO_AVG,
                FCI_EXO = FCI_EXO_AVG,
                FCI_exCredit = FCI_exCredit_AVG,
                FCI_ENDO_exCredit = FCI_ENDO_exCredit_AVG)

# --- Load macro controls ---
cat("Loading macro controls...\n")
macro_raw <- read_excel(CONFIG$fci_file, sheet = CONFIG$macro_sheet)
fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]

macro_controls <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  mutate(
    IPC_yoy = (IPC / lag(IPC, 12) - 1) * 100,
    Cred_Real_yoy = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100
  ) %>%
  dplyr::select(fecha, IPC_yoy, Cred_Real_yoy)

# --- Load monthly published data ---
cat("Loading FCI_data_1.xlsx (Monthly_SA)...\n")

monthly_raw <- read_excel(CONFIG$published_file, sheet = "Monthly_SA")
names(monthly_raw)[1] <- "fecha"

monthly_data <- monthly_raw %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Rename columns to safe names
monthly_col_map <- c(
  "Sector Primario"                      = "Primario",
  "Sector Secundario"                    = "Secundario",
  "Sector Manufacturero"                 = "Manufacturero",
  "Servicios"                            = "Servicios",
  "IMAEP"                               = "IMAEP",
  "IMAEP sin agricultura ni binacionales" = "IMAEP_SANB"
)

for (old_name in names(monthly_col_map)) {
  if (old_name %in% names(monthly_data)) {
    names(monthly_data)[names(monthly_data) == old_name] <- monthly_col_map[old_name]
  }
}

# Compute YoY growth: (Y_t / Y_{t-12} - 1) * 100
for (s in MONTHLY_SECTORS) {
  if (s %in% names(monthly_data)) {
    yoy_name <- paste0(s, "_yoy")
    monthly_data <- monthly_data %>%
      mutate(!!yoy_name := (!!sym(s) / lag(!!sym(s), 12) - 1) * 100)
  }
}

cat(sprintf("  Monthly published: %d observations, %s to %s\n",
            nrow(monthly_data), format(min(monthly_data$fecha)),
            format(max(monthly_data$fecha))))

# --- Merge all monthly data ---
analysis_m <- fci_monthly %>%
  inner_join(monthly_data, by = "fecha") %>%
  left_join(macro_controls, by = "fecha") %>%
  mutate(
    FCI_COMP_pos = pmax(FCI_COMP, 0),
    FCI_COMP_neg = abs(pmin(FCI_COMP, 0))
  ) %>%
  arrange(fecha)

cat(sprintf("\nMerged monthly dataset: %d observations, %s to %s\n",
            nrow(analysis_m), format(min(analysis_m$fecha)),
            format(max(analysis_m$fecha))))

# Report available YoY variables
avail_yoy <- names(analysis_m)[grepl("_yoy$", names(analysis_m))]
cat(sprintf("  YoY growth variables: %d\n", length(avail_yoy)))
for (v in avail_yoy) {
  n_ok <- sum(!is.na(analysis_m[[v]]))
  cat(sprintf("    %s: %d non-NA\n", v, n_ok))
}


################################################################################
# 4. PART B: STANDARD MONTHLY LP (h = 1...24)
################################################################################

cat("\n================================================================================\n")
cat("PART B: STANDARD MONTHLY LP (h = 1...24)\n")
cat("================================================================================\n\n")

monthly_results <- list()

for (sector in MONTHLY_SECTORS) {
  y_var <- paste0(sector, "_yoy")
  if (!y_var %in% names(analysis_m)) {
    cat(sprintf("  SKIP: %s not found\n", sector))
    next
  }

  data_sub <- analysis_m %>% filter(!is.na(!!sym(y_var)))
  cat(sprintf("  %s (n=%d) ... ", MONTHLY_LABELS[sector], nrow(data_sub)))

  res <- run_lp_standard(
    data_sub, y_var, "FCI_COMP", CONFIG$max_horizon, CONFIG$n_lags,
    control_vars = c("Cred_Real_yoy", "IPC_yoy")
  )

  if (nrow(res) > 0) {
    res$sector <- sector
    res$sector_label <- MONTHLY_LABELS[sector]
    monthly_results[[sector]] <- res
    cat("done\n")
  } else {
    cat("insufficient data\n")
  }
}

all_monthly <- bind_rows(monthly_results)
cat(sprintf("\n  Monthly results: %d sector-horizon pairs\n\n", nrow(all_monthly)))

# Summary table
cat("MONTHLY LP SUMMARY (Published data, Jan 2014+):\n")
cat(sprintf("%-30s %10s %8s  %10s %8s  %10s %8s\n",
            "Sector", "h=6", "p", "h=12", "p", "h=18", "p"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (sector in MONTHLY_SECTORS) {
  v6 <- all_monthly %>% filter(sector == !!sector, horizon == 6)
  v12 <- all_monthly %>% filter(sector == !!sector, horizon == 12)
  v18 <- all_monthly %>% filter(sector == !!sector, horizon == 18)

  fmt <- function(v) {
    if (nrow(v) > 0) {
      stars <- ifelse(v$p_value < 0.01, "***", ifelse(v$p_value < 0.05, "**",
                      ifelse(v$p_value < 0.10, "*", "")))
      sprintf("%+.2f%s", v$coef, stars)
    } else "NA"
  }
  p_fmt <- function(v) if (nrow(v) > 0) sprintf("%.3f", v$p_value) else "NA"

  cat(sprintf("%-30s %10s %8s  %10s %8s  %10s %8s\n",
              MONTHLY_LABELS[sector],
              fmt(v6), p_fmt(v6), fmt(v12), p_fmt(v12), fmt(v18), p_fmt(v18)))
}

# Plot 280: Monthly published -- all sectors
if (nrow(all_monthly) > 0) {
  p280 <- ggplot(all_monthly, aes(x = horizon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#3498DB", alpha = 0.2) +
    geom_line(aes(y = coef), color = "#3498DB", linewidth = 1) +
    geom_point(aes(y = coef), color = "#3498DB", size = 1.5) +
    facet_wrap(~sector_label, scales = "free_y", ncol = 3) +
    theme_minimal(base_size = 11) +
    labs(title = "Sectoral Output Response to FCI Tightening (Monthly, Published Data)",
         subtitle = paste0("FCI_COMP \u2192 Output growth (YoY) | Jan 2014+ (n=",
                           nrow(analysis_m), ") | 90% CI"),
         x = "Months", y = "Effect (pp)") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(CONFIG$output_dir, "280_PubM_LP_Monthly_Sectors.png"), p280,
         width = 12, height = 8, dpi = 300)
  cat("Saved: 280_PubM_LP_Monthly_Sectors.png\n")
}

# Plot 281: Heatmap at h=6, h=12
heatmap_m <- all_monthly %>%
  filter(horizon %in% c(6, 12)) %>%
  mutate(
    horizon_label = paste0("h=", horizon),
    label = sprintf("%.1f%s", coef,
                    ifelse(p_value < 0.01, "***",
                           ifelse(p_value < 0.05, "**",
                                  ifelse(p_value < 0.10, "*", ""))))
  )

if (nrow(heatmap_m) > 0) {
  p281 <- ggplot(heatmap_m, aes(x = horizon_label, y = sector_label, fill = coef)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#27AE60",
                         midpoint = 0, name = "Effect (pp)") +
    theme_minimal(base_size = 12) +
    labs(title = "Monthly FCI Effect on Published Sectors",
         subtitle = "Effect in pp | Published data Jan 2014+ | *** p<0.01, ** p<0.05, * p<0.10",
         x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(CONFIG$output_dir, "281_PubM_LP_Heatmap.png"), p281,
         width = 8, height = 6, dpi = 300)
  cat("Saved: 281_PubM_LP_Heatmap.png\n")
}

write.csv(all_monthly,
          file.path(CONFIG$output_dir, "csv", "PubM_LP_Monthly_Sectors.csv"), row.names = FALSE)
cat("Saved: PubM_LP_Monthly_Sectors.csv\n")


################################################################################
# 5. PART C: CORE DUAL-ENGINE TEST
################################################################################

cat("\n================================================================================\n")
cat("PART C: CORE DUAL-ENGINE TEST (Monthly Published)\n")
cat("================================================================================\n\n")

cat("Key test: Sector Primario non-responsive vs IMAEP_SANB responsive?\n")
cat("Caveat: Primario includes Ganaderia (partially finance-sensitive)\n\n")

primario_12 <- all_monthly %>% filter(sector == "Primario", horizon == 12)
sanb_12 <- all_monthly %>% filter(sector == "IMAEP_SANB", horizon == 12)
imaep_12 <- all_monthly %>% filter(sector == "IMAEP", horizon == 12)

if (nrow(primario_12) > 0) {
  cat(sprintf("  Sector Primario (h=12): coef = %+.3f, p = %.3f\n",
              primario_12$coef, primario_12$p_value))
}
if (nrow(sanb_12) > 0) {
  cat(sprintf("  IMAEP_SANB (h=12):      coef = %+.3f, p = %.3f\n",
              sanb_12$coef, sanb_12$p_value))
}
if (nrow(imaep_12) > 0) {
  cat(sprintf("  IMAEP (h=12):           coef = %+.3f, p = %.3f\n",
              imaep_12$coef, imaep_12$p_value))
}

cat("\nINTERPRETATION:\n")
if (nrow(primario_12) > 0 && nrow(sanb_12) > 0) {
  if (primario_12$p_value >= 0.10 && sanb_12$p_value < 0.10) {
    cat("  >>> DUAL-ENGINE CONFIRMED: Primary sector unresponsive, non-agro GDP responsive\n")
  } else if (primario_12$p_value >= 0.10 && sanb_12$p_value >= 0.10) {
    cat("  >>> Primary sector unresponsive (as expected), but IMAEP_SANB also non-significant\n")
    cat("      May reflect limited sample size (post-2014 only)\n")
  } else {
    cat("  >>> Mixed result. See coefficients above for interpretation.\n")
  }
}

# Plot 283: Dual-engine test
dual_sectors <- all_monthly %>%
  filter(sector %in% c("Primario", "IMAEP_SANB"))

if (nrow(dual_sectors) > 0) {
  p283 <- ggplot(dual_sectors, aes(x = horizon, color = sector_label, fill = sector_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(aes(y = coef), linewidth = 1.2) +
    geom_point(aes(y = coef), size = 2) +
    scale_color_manual(values = c("Primary Sector" = "#27AE60",
                                   "Non-Agro GDP (IMAEP ex-Agri)" = "#2C3E50")) +
    scale_fill_manual(values = c("Primary Sector" = "#27AE60",
                                  "Non-Agro GDP (IMAEP ex-Agri)" = "#2C3E50")) +
    theme_minimal(base_size = 12) +
    labs(title = "Dual-Engine Test: Primary Sector vs Non-Agro GDP (Monthly Published)",
         subtitle = paste0("FCI_COMP \u2192 Output growth (YoY) | Jan 2014+ (n=",
                           nrow(analysis_m), ") | 90% CI"),
         x = "Months", y = "Effect (pp)",
         color = NULL, fill = NULL) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(CONFIG$output_dir, "283_PubM_DualEngine_Test.png"), p283,
         width = 10, height = 6, dpi = 300)
  cat("Saved: 283_PubM_DualEngine_Test.png\n")
}


################################################################################
# 6. PART D: COMPARISON WITH ORIGINAL
################################################################################

cat("\n================================================================================\n")
cat("PART D: COMPARISON WITH ORIGINAL MONTHLY RESULTS\n")
cat("================================================================================\n\n")

# Load original monthly results (from script 14 Part B or script 12)
orig_file <- file.path(CONFIG$output_dir, "csv", "LP_Sectoral_Output_Full.csv")
if (!file.exists(orig_file)) {
  orig_file <- file.path(CONFIG$output_dir, "LP_Sectoral_Output_Full.csv")
}

comparison_m <- data.frame()

if (file.exists(orig_file)) {
  cat("Loading original monthly results...\n")
  orig_data <- read.csv(orig_file, stringsAsFactors = FALSE)

  # Match: Published IMAEP vs Original PIB; Published IMAEP_SANB vs Original PIB SANB
  matches <- data.frame(
    pub_sector = c("IMAEP", "IMAEP_SANB"),
    orig_sector = c("PIB", "PIB SANB"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(matches)) {
    for (h in CONFIG$horizons_report) {
      pub_res <- all_monthly %>%
        filter(sector == matches$pub_sector[i], horizon == h)
      orig_res <- orig_data %>%
        filter(sector == matches$orig_sector[i], horizon == h)

      if (nrow(pub_res) > 0 && nrow(orig_res) > 0) {
        z_stat <- (pub_res$coef - orig_res$coef) / sqrt(pub_res$se^2 + orig_res$se^2)
        z_p <- 2 * (1 - pnorm(abs(z_stat)))

        comparison_m <- rbind(comparison_m, data.frame(
          variable = matches$pub_sector[i],
          horizon = h,
          pub_coef = pub_res$coef,
          pub_se = pub_res$se,
          pub_p = pub_res$p_value,
          pub_n = pub_res$n_obs,
          orig_coef = orig_res$coef,
          orig_se = orig_res$se,
          orig_p = orig_res$p_value,
          orig_n = orig_res$n_obs,
          z_stat = z_stat,
          z_p = z_p,
          equal = ifelse(z_p > 0.10, "Yes", "No")
        ))
      }
    }
  }

  if (nrow(comparison_m) > 0) {
    cat("\nCOMPARISON: Published Monthly vs Internal Monthly\n")
    cat(sprintf("%-15s %4s %10s %10s %8s %8s\n",
                "Variable", "h", "Pub.Coef", "Orig.Coef", "z-stat", "Equal?"))
    cat(paste(rep("-", 65), collapse = ""), "\n")

    for (i in 1:nrow(comparison_m)) {
      r <- comparison_m[i, ]
      cat(sprintf("%-15s %4d %+10.2f %+10.2f %8.2f %8s\n",
                  r$variable, r$horizon, r$pub_coef, r$orig_coef,
                  r$z_stat, r$equal))
    }
    cat("\nEqual? = z-test p > 0.10\n")
    cat("Note: Published data is post-2014 only; original uses full 1996+ sample\n")
  }

  # Plot 284: Overlay comparison
  overlay_data <- data.frame()

  for (i in 1:nrow(matches)) {
    pub_res <- all_monthly %>%
      filter(sector == matches$pub_sector[i]) %>%
      mutate(source = "Published (2014+)",
             variable = matches$pub_sector[i])

    orig_res <- orig_data %>%
      filter(sector == matches$orig_sector[i]) %>%
      mutate(source = "Internal (1996+)",
             variable = matches$pub_sector[i])

    overlay_data <- bind_rows(overlay_data, pub_res, orig_res)
  }

  if (nrow(overlay_data) > 0) {
    p284 <- ggplot(overlay_data, aes(x = horizon, y = coef,
                                      color = source, fill = source)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.10, color = NA) +
      geom_line(linewidth = 1) +
      geom_point(size = 1.5) +
      facet_wrap(~variable, scales = "free_y") +
      scale_color_manual(values = c("Published (2014+)" = "#3498DB",
                                     "Internal (1996+)" = "#E74C3C")) +
      scale_fill_manual(values = c("Published (2014+)" = "#3498DB",
                                    "Internal (1996+)" = "#E74C3C")) +
      theme_minimal(base_size = 11) +
      labs(title = "Comparison: Published vs Internal Monthly Data",
           subtitle = "FCI_COMP \u2192 Output growth (YoY) | 90% CI",
           x = "Months", y = "Effect (pp)",
           color = "Data Source", fill = "Data Source") +
      theme(plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom")

    ggsave(file.path(CONFIG$output_dir, "284_PubM_Comparison_Overlay.png"), p284,
           width = 12, height = 5, dpi = 300)
    cat("Saved: 284_PubM_Comparison_Overlay.png\n")
  }
} else {
  cat("Original monthly results not found. Skipping comparison.\n")
  cat("  Expected file:", orig_file, "\n")
}

write.csv(comparison_m,
          file.path(CONFIG$output_dir, "csv", "PubM_Comparison.csv"), row.names = FALSE)
cat("Saved: PubM_Comparison.csv\n")


################################################################################
# 7. FINAL SUMMARY
################################################################################

cat("\n################################################################################\n")
cat("MONTHLY LP WITH PUBLISHED DATA - COMPLETE\n")
cat("################################################################################\n\n")

cat("Outputs generated:\n")
cat("  PNGs:  280-284 series\n")
cat("  CSVs:  PubM_LP_Monthly_Sectors, PubM_Comparison\n\n")

cat("Key findings:\n")
cat("  1. Dual-engine test (Primary vs IMAEP_SANB)? -> Check Part C above\n")
cat("  2. Monthly published ~ original internal? -> Check Part D comparison\n")
cat("  3. Sample covers Jan 2014+: purely post-IT check\n\n")

cat("Output directory:", CONFIG$output_dir, "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
