################################################################################
# FCI NEW EXTERNAL DATA CONSTRUCTION & PRELIMINARY ANALYSIS
################################################################################
#
# Project:      Financial Conditions Index - Expanded External Variables
# Institution:  Banco Central del Paraguay
# Department:   Departamento de Analisis Macroeconomico
#
# Description:  Loads new external variables from FCI_data_1.xlsx, constructs
#               a Global Financial Conditions (GFC) index, a regional financial
#               conditions proxy (Brazil Selic), disaggregated commodity YoY
#               changes, and runs preliminary correlation analysis.
#
#   SECTION A: Setup & Data Loading
#   SECTION B: Construct Global Financial Conditions Index (GFC)
#   SECTION C: Construct Regional Financial Conditions Proxy
#   SECTION D: Compute YoY Changes for Commodity Prices
#   SECTION E: Preliminary Correlation Analysis
#   SECTION F: Save Constructed Indices
#
# References:
#   - Brave & Butters (2011) - Chicago Fed NFCI
#   - Hatzius et al. (2010) - Goldman Sachs FCI
#   - Miranda-Agrippino & Rey (2020) - Global Financial Cycle
#
################################################################################


################################################################################
# SECTION A: SETUP & DATA LOADING
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(FactoMineR)
  library(ggplot2)
  library(gridExtra)
  library(sandwich)
  library(lmtest)
})

EXT_CONFIG <- list(
  fci_data_file   = "../data/FCI_data_1.xlsx",
  exo_data_file   = "../data/FCI_data_1.xlsx",
  macro_sheet     = "Datos_macro",
  gfc_sheet       = "Global_Financial_Conditions",
  comm_sheet      = "Main_Commodities_Prices",
  rolling_window  = 60,
  output_dir      = "../output",
  verbose         = TRUE
)

set.seed(20260211)

cat("\n################################################################################\n")
cat("FCI NEW EXTERNAL DATA CONSTRUCTION & PRELIMINARY ANALYSIS\n")
cat("################################################################################\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---- Source script 01 if needed ----
if (!exists("resultado_fci")) {
  cat("Loading FCI from 01_FCI_Complete.R...\n")
  CONFIG_SAVE <- if (exists("CONFIG")) CONFIG else NULL
  source("01_FCI_Complete.R")
  if (!is.null(CONFIG_SAVE)) CONFIG <- CONFIG_SAVE
  rm(CONFIG_SAVE)
  suppressPackageStartupMessages(library(dplyr))
}

# ---- Load Global Financial Conditions sheet ----
cat("Loading Global Financial Conditions data...\n")
gfc_raw <- read_excel(EXT_CONFIG$exo_data_file, sheet = EXT_CONFIG$gfc_sheet)
names(gfc_raw)[1] <- "fecha"
gfc_raw <- gfc_raw %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

# Rename S&P 500 to valid R name
if ("S&P 500" %in% names(gfc_raw)) {
  gfc_raw <- gfc_raw %>% rename(SP500 = `S&P 500`)
}

cat("  GFC sheet:", nrow(gfc_raw), "obs,", format(min(gfc_raw$fecha)), "to", format(max(gfc_raw$fecha)), "\n")
cat("  Variables:", paste(setdiff(names(gfc_raw), "fecha"), collapse = ", "), "\n\n")

# ---- Load Main Commodities Prices sheet ----
cat("Loading Main Commodities Prices data...\n")
comm_raw <- read_excel(EXT_CONFIG$exo_data_file, sheet = EXT_CONFIG$comm_sheet)
names(comm_raw)[1] <- "fecha"
comm_raw <- comm_raw %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha)

cat("  Commodities sheet:", nrow(comm_raw), "obs,", format(min(comm_raw$fecha)), "to", format(max(comm_raw$fecha)), "\n")
cat("  Commodities:", paste(setdiff(names(comm_raw), "fecha"), collapse = ", "), "\n\n")

# ---- Merge both sheets ----
new_ext_data <- gfc_raw %>%
  full_join(comm_raw, by = "fecha") %>%
  arrange(fecha)

cat("Merged external data:", nrow(new_ext_data), "obs,",
    format(min(new_ext_data$fecha)), "to", format(max(new_ext_data$fecha)), "\n")

# ---- Merge with FCI indices ----
fci_indices <- resultado_fci$all_indices %>%
  dplyr::select(fecha,
                FCI_COMP_AVG, FCI_ENDO_AVG, FCI_EXO_AVG,
                FCI_exCredit_AVG, FCI_ENDO_exCredit_AVG)

new_ext_data <- new_ext_data %>%
  left_join(fci_indices, by = "fecha")

# ---- Load macro data ----
macro_raw <- read_excel(EXT_CONFIG$fci_data_file, sheet = EXT_CONFIG$macro_sheet)
fecha_col_macro <- names(macro_raw)[grepl("fecha|date", names(macro_raw), ignore.case = TRUE)][1]
macro_data <- macro_raw %>%
  rename(fecha = !!sym(fecha_col_macro)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  mutate(
    IMAEP_yoy = (IMAEP / lag(IMAEP, 12) - 1) * 100,
    IPC_yoy   = (IPC / lag(IPC, 12) - 1) * 100
  )

# Real credit growth
if ("Creditos_deflactados" %in% names(macro_data)) {
  macro_data <- macro_data %>%
    mutate(Cred_Real_Total = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100)
}

new_ext_data <- new_ext_data %>%
  left_join(macro_data %>% dplyr::select(fecha, IMAEP_yoy, IPC_yoy,
                                          any_of("Cred_Real_Total")),
            by = "fecha")

# ---- Also load original FCI variables (FFER, VIX, Commodities) ----
datos_raw <- read_excel(EXT_CONFIG$fci_data_file)
fecha_col <- names(datos_raw)[grepl("fecha|date", names(datos_raw), ignore.case = TRUE)][1]
fci_vars_orig <- datos_raw %>%
  rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  dplyr::select(fecha, any_of(c("FFER", "VIX", "Commodities")))

new_ext_data <- new_ext_data %>%
  left_join(fci_vars_orig, by = "fecha")

# ---- Missing data summary ----
cat("\n--- Missing Data Summary ---\n")
vars_to_check <- c("US_10Y", "SP500", "DXY", "Selic_rate", "IPE", "Paraguay_ToT",
                    "Soybean", "Oil_Brent", "FFER", "VIX", "Commodities",
                    "FCI_ENDO_AVG", "Cred_Real_Total")
for (v in vars_to_check) {
  if (v %in% names(new_ext_data)) {
    n_na <- sum(is.na(new_ext_data[[v]]))
    first_valid <- min(which(!is.na(new_ext_data[[v]])))
    last_valid  <- max(which(!is.na(new_ext_data[[v]])))
    cat(sprintf("  %-25s: %3d NAs, valid from %s to %s\n",
                v, n_na,
                format(new_ext_data$fecha[first_valid]),
                format(new_ext_data$fecha[last_valid])))
  }
}

# Determine common sample
sample_vars <- c("US_10Y", "SP500", "DXY", "FCI_ENDO_AVG")
available_vars <- intersect(sample_vars, names(new_ext_data))
common_sample <- new_ext_data %>%
  dplyr::select(fecha, all_of(available_vars)) %>%
  na.omit()
cat("\nCommon sample (GFC + FCI):", nrow(common_sample), "obs,",
    format(min(common_sample$fecha)), "to", format(max(common_sample$fecha)), "\n\n")


################################################################################
# SECTION B: CONSTRUCT GLOBAL FINANCIAL CONDITIONS INDEX (GFC)
################################################################################

cat("================================================================================\n")
cat("SECTION B: GLOBAL FINANCIAL CONDITIONS INDEX (GFC)\n")
cat("================================================================================\n\n")

# Rolling standardization helper
rolling_z <- function(x, window = EXT_CONFIG$rolling_window) {
  rm <- zoo::rollapply(x, width = window, FUN = mean, na.rm = TRUE,
                       fill = NA, align = "right", partial = TRUE)
  rsd <- zoo::rollapply(x, width = window, FUN = sd, na.rm = TRUE,
                        fill = NA, align = "right", partial = TRUE)
  rsd <- pmax(rsd, 1e-10)
  (x - rm) / rsd
}

# ---- Option A: Simple Average GFC ----
# Signs: VIX(+), DXY(+), US10Y(+), SP500(-), FFER(+) → higher = tighter
gfc_vars <- c("VIX", "DXY", "US_10Y", "SP500", "FFER")
gfc_signs <- c(+1, +1, +1, -1, +1)
gfc_names <- c("z_VIX", "z_DXY", "z_US10Y", "z_SP500_neg", "z_FFER")

for (i in seq_along(gfc_vars)) {
  v <- gfc_vars[i]
  if (v %in% names(new_ext_data)) {
    new_ext_data[[gfc_names[i]]] <- rolling_z(new_ext_data[[v]] * gfc_signs[i])
  }
}

available_gfc <- gfc_names[gfc_names %in% names(new_ext_data)]
if (length(available_gfc) > 0) {
  new_ext_data$GFC_simple <- rowMeans(new_ext_data[, available_gfc, drop = FALSE], na.rm = TRUE)
  new_ext_data$GFC_simple[rowSums(!is.na(new_ext_data[, available_gfc])) < 3] <- NA
  cat("GFC_simple constructed from", length(available_gfc), "z-scored variables\n")
} else {
  new_ext_data$GFC_simple <- NA
  cat("WARNING: Could not construct GFC_simple\n")
}

# ---- Option B: PCA-Based GFC ----
cat("Computing PCA-based GFC...\n")
gfc_pca_data <- new_ext_data %>%
  dplyr::select(fecha, all_of(available_gfc)) %>%
  na.omit()

gfc_pca_result <- NULL
if (nrow(gfc_pca_data) > 30) {
  tryCatch({
    pca_input <- gfc_pca_data[, available_gfc, drop = FALSE]
    pca_model <- PCA(pca_input, ncp = 1, graph = FALSE)

    gfc_pca_scores <- pca_model$ind$coord[, 1]
    gfc_pca_loadings <- pca_model$var$coord[, 1]
    gfc_pca_var_explained <- pca_model$eig[1, 2]

    # Ensure positive = tighter (positive correlation with z_VIX)
    if (cor(gfc_pca_scores, pca_input$z_VIX, use = "complete.obs") < 0) {
      gfc_pca_scores <- -gfc_pca_scores
      gfc_pca_loadings <- -gfc_pca_loadings
    }

    gfc_pca_data$GFC_PCA <- gfc_pca_scores

    # Merge back
    new_ext_data <- new_ext_data %>%
      left_join(gfc_pca_data %>% dplyr::select(fecha, GFC_PCA), by = "fecha")

    gfc_pca_result <- list(
      loadings = gfc_pca_loadings,
      var_explained = gfc_pca_var_explained,
      n_obs = nrow(gfc_pca_data)
    )

    cat("  PC1 variance explained:", round(gfc_pca_var_explained, 1), "%\n")
    cat("  Loadings:\n")
    for (j in seq_along(available_gfc)) {
      cat(sprintf("    %-20s: %+.3f\n", available_gfc[j], gfc_pca_loadings[j]))
    }

    # Correlation between simple and PCA
    corr_gfc <- cor(new_ext_data$GFC_simple, new_ext_data$GFC_PCA, use = "complete.obs")
    cat("  Correlation GFC_simple vs GFC_PCA:", round(corr_gfc, 3), "\n\n")
  }, error = function(e) {
    cat("  WARNING: PCA failed:", e$message, "\n")
    new_ext_data$GFC_PCA <<- NA
  })
} else {
  new_ext_data$GFC_PCA <- NA
  cat("  WARNING: Insufficient observations for PCA\n")
}


################################################################################
# SECTION C: CONSTRUCT REGIONAL FINANCIAL CONDITIONS PROXY
################################################################################

cat("================================================================================\n")
cat("SECTION C: REGIONAL FINANCIAL CONDITIONS PROXY\n")
cat("================================================================================\n\n")

if ("Selic_rate" %in% names(new_ext_data)) {
  new_ext_data$Regional_FC <- rolling_z(new_ext_data$Selic_rate)
  n_valid_selic <- sum(!is.na(new_ext_data$Regional_FC))
  cat("Regional_FC (z-scored Selic) constructed:", n_valid_selic, "valid obs\n")

  # Correlation with FCI_ENDO
  if ("FCI_ENDO_AVG" %in% names(new_ext_data)) {
    corr_selic_fci <- cor(new_ext_data$Regional_FC, new_ext_data$FCI_ENDO_AVG,
                          use = "complete.obs")
    cat("  Correlation Regional_FC vs FCI_ENDO:", round(corr_selic_fci, 3), "\n\n")
  }
} else {
  new_ext_data$Regional_FC <- NA
  cat("WARNING: Selic_rate not found\n\n")
}


################################################################################
# SECTION D: COMPUTE YoY CHANGES FOR COMMODITY PRICES
################################################################################

cat("================================================================================\n")
cat("SECTION D: COMMODITY YoY CHANGES\n")
cat("================================================================================\n\n")

commodity_vars <- c("Soybean", "Soybean_oil", "Soybean_flour", "Cotton",
                    "Corn", "Wheat", "Meat", "Sugar", "Oil_Brent")

for (cv in commodity_vars) {
  yoy_name <- paste0("d_", cv)
  if (cv %in% names(new_ext_data)) {
    new_ext_data[[yoy_name]] <- (new_ext_data[[cv]] / lag(new_ext_data[[cv]], 12) - 1) * 100
    n_valid <- sum(!is.na(new_ext_data[[yoy_name]]))
    cat(sprintf("  %-20s: %d valid obs\n", yoy_name, n_valid))
  }
}

# ToT YoY change
if ("Paraguay_ToT" %in% names(new_ext_data)) {
  new_ext_data$d_ToT <- (new_ext_data$Paraguay_ToT / lag(new_ext_data$Paraguay_ToT, 12) - 1) * 100
  cat(sprintf("  %-20s: %d valid obs\n", "d_ToT", sum(!is.na(new_ext_data$d_ToT))))
}

cat("\n")


################################################################################
# SECTION E: PRELIMINARY CORRELATION ANALYSIS
################################################################################

cat("================================================================================\n")
cat("SECTION E: PRELIMINARY CORRELATION ANALYSIS\n")
cat("================================================================================\n\n")

# ---- E1: New variables vs FCI variants ----
new_level_vars <- c("US_10Y", "SP500", "DXY", "Selic_rate", "IPE", "Paraguay_ToT",
                    "GFC_simple", "GFC_PCA", "Regional_FC")
new_yoy_vars <- paste0("d_", c(commodity_vars, "ToT"))
fci_vars_corr <- c("FCI_COMP_AVG", "FCI_ENDO_AVG", "FCI_EXO_AVG",
                    "FCI_exCredit_AVG", "FCI_ENDO_exCredit_AVG")

# Build correlation matrix: new vars vs FCI
all_corr_vars <- c(new_level_vars, new_yoy_vars)
available_corr <- intersect(all_corr_vars, names(new_ext_data))
available_fci  <- intersect(fci_vars_corr, names(new_ext_data))

corr_matrix_fci <- matrix(NA, nrow = length(available_corr), ncol = length(available_fci),
                           dimnames = list(available_corr, available_fci))

for (v in available_corr) {
  for (f in available_fci) {
    corr_matrix_fci[v, f] <- cor(new_ext_data[[v]], new_ext_data[[f]], use = "complete.obs")
  }
}

cat("Correlations with FCI_ENDO_AVG:\n")
if ("FCI_ENDO_AVG" %in% colnames(corr_matrix_fci)) {
  sorted_corr <- sort(corr_matrix_fci[, "FCI_ENDO_AVG"], decreasing = TRUE, na.last = TRUE)
  for (nm in names(sorted_corr)) {
    cat(sprintf("  %-25s: %+.3f\n", nm, sorted_corr[nm]))
  }
}

# ---- E2: New variables vs credit growth ----
cat("\nCorrelations with Real Credit Growth:\n")
if ("Cred_Real_Total" %in% names(new_ext_data)) {
  for (v in available_corr) {
    cr <- cor(new_ext_data[[v]], new_ext_data$Cred_Real_Total, use = "complete.obs")
    if (!is.na(cr)) {
      cat(sprintf("  %-25s: %+.3f\n", v, cr))
    }
  }
}

# ---- E3: Cross-correlation among new variables (multicollinearity check) ----
cat("\nMulticollinearity check (cross-correlations among new variables):\n")
multicol_vars <- intersect(c("US_10Y", "SP500", "DXY", "Selic_rate",
                              "VIX", "FFER", "GFC_simple"), names(new_ext_data))
if (length(multicol_vars) > 1) {
  multicol_matrix <- cor(new_ext_data[, multicol_vars], use = "pairwise.complete.obs")
  high_corr <- which(abs(multicol_matrix) > 0.7 & upper.tri(multicol_matrix), arr.ind = TRUE)
  if (nrow(high_corr) > 0) {
    cat("  HIGH correlations (|r| > 0.7):\n")
    for (k in seq_len(nrow(high_corr))) {
      i <- high_corr[k, 1]
      j <- high_corr[k, 2]
      cat(sprintf("    %s × %s = %.3f\n",
                  rownames(multicol_matrix)[i], colnames(multicol_matrix)[j],
                  multicol_matrix[i, j]))
    }
  } else {
    cat("  No pairwise correlations exceed |0.7|\n")
  }
}

# ---- E4: Compare disaggregated commodities vs aggregate ----
cat("\nDisaggregated commodities vs FCI_ENDO_AVG (vs aggregate Commodities):\n")
if ("FCI_ENDO_AVG" %in% names(new_ext_data) && "Commodities" %in% names(new_ext_data)) {
  agg_corr <- cor(new_ext_data$Commodities, new_ext_data$FCI_ENDO_AVG, use = "complete.obs")
  cat(sprintf("  %-25s: %+.3f (AGGREGATE)\n", "Commodities", agg_corr))
}
for (cv in commodity_vars) {
  yoy_name <- paste0("d_", cv)
  if (yoy_name %in% names(new_ext_data) && "FCI_ENDO_AVG" %in% names(new_ext_data)) {
    cr <- cor(new_ext_data[[yoy_name]], new_ext_data$FCI_ENDO_AVG, use = "complete.obs")
    if (!is.na(cr)) cat(sprintf("  %-25s: %+.3f\n", yoy_name, cr))
  }
}
cat("\n")


################################################################################
# SECTION F: SAVE CONSTRUCTED INDICES & VISUALIZATIONS
################################################################################

cat("================================================================================\n")
cat("SECTION F: SAVE OUTPUTS\n")
cat("================================================================================\n\n")

# ---- F1: Save main dataset for downstream scripts ----
export_vars <- c("fecha",
                 # GFC indices
                 "GFC_simple", "GFC_PCA", "Regional_FC",
                 # Level variables
                 "US_10Y", "SP500", "DXY", "Selic_rate", "IPE", "Paraguay_ToT",
                 # Original external vars
                 "FFER", "VIX", "Commodities",
                 # FCI indices
                 "FCI_COMP_AVG", "FCI_ENDO_AVG", "FCI_EXO_AVG",
                 "FCI_exCredit_AVG", "FCI_ENDO_exCredit_AVG",
                 # Macro
                 "IMAEP_yoy", "IPC_yoy", "Cred_Real_Total",
                 # Commodity YoY
                 paste0("d_", commodity_vars), "d_ToT",
                 # Raw commodity levels
                 commodity_vars, "Oil_Brent")

available_export <- intersect(export_vars, names(new_ext_data))
export_df <- new_ext_data[, available_export]

write.csv(export_df, file.path(EXT_CONFIG$output_dir, "New_External_Variables.csv"),
          row.names = FALSE)
cat("Saved: New_External_Variables.csv (", nrow(export_df), " obs, ",
    length(available_export), " vars)\n", sep = "")

# ---- F2: GFC Diagnostics CSV ----
gfc_diag <- data.frame(
  Variable = available_gfc,
  PCA_Loading = if (!is.null(gfc_pca_result)) gfc_pca_result$loadings else rep(NA, length(available_gfc)),
  stringsAsFactors = FALSE
)
gfc_diag$Method <- "GFC_PCA"
gfc_diag$Variance_Explained <- if (!is.null(gfc_pca_result)) gfc_pca_result$var_explained else NA
gfc_diag$GFC_Simple_PCA_Correlation <- if ("GFC_PCA" %in% names(new_ext_data)) {
  cor(new_ext_data$GFC_simple, new_ext_data$GFC_PCA, use = "complete.obs")
} else NA

write.csv(gfc_diag, file.path(EXT_CONFIG$output_dir, "GFC_Index_Diagnostics.csv"),
          row.names = FALSE)
cat("Saved: GFC_Index_Diagnostics.csv\n")

# ---- F3: Full correlation matrix CSV ----
full_corr_df <- as.data.frame(corr_matrix_fci)
full_corr_df$Variable <- rownames(corr_matrix_fci)
full_corr_df <- full_corr_df %>% dplyr::select(Variable, everything())
write.csv(full_corr_df, file.path(EXT_CONFIG$output_dir, "New_Variables_Correlation_Matrix.csv"),
          row.names = FALSE)
cat("Saved: New_Variables_Correlation_Matrix.csv\n")

# ---- F4: Disaggregated commodity correlations CSV ----
comm_corr_list <- list()
for (cv in commodity_vars) {
  yoy_name <- paste0("d_", cv)
  if (yoy_name %in% names(new_ext_data)) {
    fci_cor <- if ("FCI_ENDO_AVG" %in% names(new_ext_data)) {
      cor(new_ext_data[[yoy_name]], new_ext_data$FCI_ENDO_AVG, use = "complete.obs")
    } else NA
    cred_cor <- if ("Cred_Real_Total" %in% names(new_ext_data)) {
      cor(new_ext_data[[yoy_name]], new_ext_data$Cred_Real_Total, use = "complete.obs")
    } else NA
    comm_corr_list[[cv]] <- data.frame(
      Commodity = cv,
      Corr_FCI_ENDO = fci_cor,
      Corr_Credit = cred_cor,
      stringsAsFactors = FALSE
    )
  }
}
if (length(comm_corr_list) > 0) {
  comm_corr_df <- do.call(rbind, comm_corr_list)

  # Add aggregate
  if ("Commodities" %in% names(new_ext_data)) {
    agg_row <- data.frame(
      Commodity = "Commodities_AGGREGATE",
      Corr_FCI_ENDO = cor(new_ext_data$Commodities, new_ext_data$FCI_ENDO_AVG, use = "complete.obs"),
      Corr_Credit = if ("Cred_Real_Total" %in% names(new_ext_data)) {
        cor(new_ext_data$Commodities, new_ext_data$Cred_Real_Total, use = "complete.obs")
      } else NA,
      stringsAsFactors = FALSE
    )
    comm_corr_df <- rbind(comm_corr_df, agg_row)
  }
  write.csv(comm_corr_df, file.path(EXT_CONFIG$output_dir, "Commodity_Disaggregated_Correlations.csv"),
            row.names = FALSE)
  cat("Saved: Commodity_Disaggregated_Correlations.csv\n")
}

# ---- F5: Sample info CSV ----
sample_info <- data.frame(
  Variable = vars_to_check,
  stringsAsFactors = FALSE
)
sample_info$N_Total <- nrow(new_ext_data)
sample_info$N_Missing <- sapply(vars_to_check, function(v) {
  if (v %in% names(new_ext_data)) sum(is.na(new_ext_data[[v]])) else NA
})
sample_info$First_Valid <- sapply(vars_to_check, function(v) {
  if (v %in% names(new_ext_data) && any(!is.na(new_ext_data[[v]]))) {
    format(new_ext_data$fecha[min(which(!is.na(new_ext_data[[v]])))])
  } else NA
})
sample_info$Last_Valid <- sapply(vars_to_check, function(v) {
  if (v %in% names(new_ext_data) && any(!is.na(new_ext_data[[v]]))) {
    format(new_ext_data$fecha[max(which(!is.na(new_ext_data[[v]])))])
  } else NA
})
write.csv(sample_info, file.path(EXT_CONFIG$output_dir, "New_External_Sample_Info.csv"),
          row.names = FALSE)
cat("Saved: New_External_Sample_Info.csv\n")


################################################################################
# VISUALIZATIONS
################################################################################

cat("\n--- Generating Visualizations ---\n\n")

# ---- 186: GFC Index Construction ----
tryCatch({
  plot_data <- new_ext_data %>%
    dplyr::select(fecha, GFC_simple, any_of("GFC_PCA")) %>%
    filter(!is.na(GFC_simple)) %>%
    pivot_longer(-fecha, names_to = "Method", values_to = "Value")

  p1 <- ggplot(plot_data, aes(x = fecha, y = Value, color = Method)) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Global Financial Conditions Index (GFC)",
         subtitle = "Simple average vs PCA-based | Higher = tighter global conditions",
         x = NULL, y = "GFC Index") +
    scale_color_manual(values = c("GFC_simple" = "steelblue", "GFC_PCA" = "firebrick")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  # PCA loadings bar chart
  if (!is.null(gfc_pca_result)) {
    loadings_df <- data.frame(
      Variable = available_gfc,
      Loading = gfc_pca_result$loadings
    )
    p2 <- ggplot(loadings_df, aes(x = reorder(Variable, Loading), y = Loading, fill = Loading > 0)) +
      geom_col(show.legend = FALSE) +
      coord_flip() +
      labs(title = paste0("PCA Loadings (Var. Explained: ", round(gfc_pca_result$var_explained, 1), "%)"),
           x = NULL, y = "PC1 Loading") +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick")) +
      theme_minimal(base_size = 11)

    combined <- grid.arrange(p1, p2, ncol = 2, widths = c(2, 1))
    ggsave(file.path(EXT_CONFIG$output_dir, "186_GFC_Index_Construction.png"),
           combined, width = 14, height = 6, dpi = 150)
  } else {
    ggsave(file.path(EXT_CONFIG$output_dir, "186_GFC_Index_Construction.png"),
           p1, width = 10, height = 6, dpi = 150)
  }
  cat("Saved: 186_GFC_Index_Construction.png\n")
}, error = function(e) cat("  WARNING: Plot 186 failed:", e$message, "\n"))

# ---- 187: Regional FC Proxy ----
tryCatch({
  plot_data <- new_ext_data %>%
    dplyr::select(fecha, Regional_FC, FCI_ENDO_AVG) %>%
    filter(!is.na(Regional_FC) & !is.na(FCI_ENDO_AVG)) %>%
    pivot_longer(-fecha, names_to = "Variable", values_to = "Value")

  p <- ggplot(plot_data, aes(x = fecha, y = Value, color = Variable)) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Regional Financial Conditions: Brazil Selic vs Paraguay FCI_ENDO",
         subtitle = "Higher = tighter conditions",
         x = NULL, y = "Standardized Index") +
    scale_color_manual(values = c("Regional_FC" = "darkgreen", "FCI_ENDO_AVG" = "steelblue"),
                       labels = c("Regional_FC" = "Brazil Selic (z-scored)",
                                  "FCI_ENDO_AVG" = "Paraguay FCI_ENDO")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(EXT_CONFIG$output_dir, "187_Regional_FC_Proxy.png"),
         p, width = 10, height = 6, dpi = 150)
  cat("Saved: 187_Regional_FC_Proxy.png\n")
}, error = function(e) cat("  WARNING: Plot 187 failed:", e$message, "\n"))

# ---- 188: Correlation Heatmap ----
tryCatch({
  # Build full correlation matrix for heatmap
  heatmap_row_vars <- intersect(c("US_10Y", "SP500", "DXY", "Selic_rate", "Paraguay_ToT",
                                   "GFC_simple", "GFC_PCA", "Regional_FC",
                                   "d_Soybean", "d_Oil_Brent", "d_ToT",
                                   "VIX", "FFER", "Commodities"),
                                names(new_ext_data))
  heatmap_col_vars <- intersect(c("FCI_COMP_AVG", "FCI_ENDO_AVG", "FCI_EXO_AVG",
                                   "FCI_exCredit_AVG", "Cred_Real_Total"),
                                names(new_ext_data))

  heatmap_mat <- matrix(NA, length(heatmap_row_vars), length(heatmap_col_vars),
                         dimnames = list(heatmap_row_vars, heatmap_col_vars))
  for (rv in heatmap_row_vars) {
    for (cv in heatmap_col_vars) {
      heatmap_mat[rv, cv] <- cor(new_ext_data[[rv]], new_ext_data[[cv]], use = "complete.obs")
    }
  }

  heatmap_df <- as.data.frame(as.table(heatmap_mat))
  names(heatmap_df) <- c("Row", "Col", "Correlation")

  p <- ggplot(heatmap_df, aes(x = Col, y = Row, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2.5) +
    scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "New External Variables: Correlation with FCI & Credit",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(EXT_CONFIG$output_dir, "188_New_Variables_Correlation_Heatmap.png"),
         p, width = 10, height = 8, dpi = 150)
  cat("Saved: 188_New_Variables_Correlation_Heatmap.png\n")
}, error = function(e) cat("  WARNING: Plot 188 failed:", e$message, "\n"))

# ---- 189: Disaggregated Commodities vs Aggregate ----
tryCatch({
  if (exists("comm_corr_df") && nrow(comm_corr_df) > 0) {
    p <- ggplot(comm_corr_df, aes(x = reorder(Commodity, Corr_FCI_ENDO), y = Corr_FCI_ENDO,
                                   fill = Commodity == "Commodities_AGGREGATE")) +
      geom_col(show.legend = FALSE) +
      coord_flip() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Commodity Correlations with FCI_ENDO",
           subtitle = "Disaggregated (YoY changes) vs Aggregate Index",
           x = NULL, y = "Correlation") +
      scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue")) +
      theme_minimal(base_size = 11)

    ggsave(file.path(EXT_CONFIG$output_dir, "189_Disaggregated_Commodities_vs_Aggregate.png"),
           p, width = 10, height = 6, dpi = 150)
    cat("Saved: 189_Disaggregated_Commodities_vs_Aggregate.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 189 failed:", e$message, "\n"))

# ---- 190: Commodity Prices Time Series ----
tryCatch({
  comm_plot_vars <- paste0("d_", commodity_vars)
  available_comm_plot <- intersect(comm_plot_vars, names(new_ext_data))

  if (length(available_comm_plot) > 0) {
    plot_data <- new_ext_data %>%
      dplyr::select(fecha, all_of(available_comm_plot)) %>%
      filter(rowSums(!is.na(.[, -1])) > 0) %>%
      pivot_longer(-fecha, names_to = "Commodity", values_to = "YoY_Change")

    p <- ggplot(plot_data, aes(x = fecha, y = YoY_Change, color = Commodity)) +
      geom_line(linewidth = 0.4, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(title = "Disaggregated Commodity Prices (YoY % Change)",
           x = NULL, y = "YoY % Change") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 7))

    ggsave(file.path(EXT_CONFIG$output_dir, "190_Commodity_Prices_TimeSeries.png"),
           p, width = 12, height = 6, dpi = 150)
    cat("Saved: 190_Commodity_Prices_TimeSeries.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 190 failed:", e$message, "\n"))

# ---- 191: New External Variables Summary Dashboard ----
tryCatch({
  dash_vars <- list(
    list(name = "US_10Y", label = "US 10Y Treasury"),
    list(name = "SP500", label = "S&P 500"),
    list(name = "DXY", label = "Dollar Index (DXY)"),
    list(name = "Selic_rate", label = "Brazil Selic Rate"),
    list(name = "Paraguay_ToT", label = "Paraguay ToT"),
    list(name = "GFC_simple", label = "GFC Index (Simple)")
  )

  plots <- list()
  for (dv in dash_vars) {
    if (dv$name %in% names(new_ext_data)) {
      pd <- new_ext_data %>%
        dplyr::select(fecha, value = !!sym(dv$name)) %>%
        filter(!is.na(value))

      plots[[dv$name]] <- ggplot(pd, aes(x = fecha, y = value)) +
        geom_line(color = "steelblue", linewidth = 0.5) +
        labs(title = dv$label, x = NULL, y = NULL) +
        theme_minimal(base_size = 8) +
        theme(plot.title = element_text(size = 9, face = "bold"))
    }
  }

  if (length(plots) > 0) {
    combined <- do.call(grid.arrange, c(plots, ncol = 3))
    ggsave(file.path(EXT_CONFIG$output_dir, "191_New_External_Variables_Summary.png"),
           combined, width = 14, height = 8, dpi = 150)
    cat("Saved: 191_New_External_Variables_Summary.png\n")
  }
}, error = function(e) cat("  WARNING: Plot 191 failed:", e$message, "\n"))


cat("\n################################################################################\n")
cat("SCRIPT 17 COMPLETE\n")
cat("################################################################################\n\n")
cat("Outputs: 6 PNGs (186-191) + 5 CSVs\n")
cat("Key dataset: New_External_Variables.csv (loaded by scripts 18-22)\n\n")
