################################################################################
# 72_Regen_FigureA2.R
#
# Regenerates ONE figure: output/png/09_FCI_Channel_Contributions_Integral.png
# (Online Appendix Figure A.2, "Channel Contributions to Comprehensive FCI").
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# Two defects were found in the version of this figure produced by
# 01_FCI_Complete.R (Plot 8, lines ~1304-1339):
#
#   (1) AXIS.  scale_x_date(date_breaks = "2 years") with ggplot2's default 5%
#       continuous expansion pushes the panel about 1.5 years past the data,
#       so a "2027" tick is drawn even though the series ends December 2025.
#
#   (2) DIVISOR.  Plot 7's shared setup (line 1204) computes
#           n_vars <- length(intersect(all_vars, names(var_contributions))) - 1
#       with the comment "# exclude fecha".  But `all_vars` does not contain
#       "fecha", so the intersect is already 12 and the "- 1" makes the divisor
#       11.  Each variable's contribution is therefore z_v/11 rather than
#       z_v/12, so the stacked bars and the black line labelled "Total FCI" are
#       (12/11) = 9.1% too large, and that line is NOT the FCI.  This script
#       uses 12 and GATES on the identity sum_v z_v/12 == FCI_COMP_ZSCORE.
#
# WHAT THIS SCRIPT DOES NOT DO
# ----------------------------
# It performs no estimation and rewrites no CSV.  It evaluates only the
# deterministic PREFIX of 01_FCI_Complete.R (through the construction of
# `datos_std`: data load, sign application, rolling standardization -- no PCA,
# VAR or MARSS, all of which are defined later and never called here), so the
# standardized inputs are byte-identical to the published pipeline's rather
# than reimplemented.  `write.csv` is stubbed out during that evaluation so the
# prefix's one incidental write (FCI_Sample_Verification.csv) does not fire.
#
# Run from R/:  Rscript 72_Regen_FigureA2.R
################################################################################

SRC        <- "01_FCI_Complete.R"
PREFIX_END <- 650L   # last line of the `datos_std <- ...` statement
TARGET_PNG <- "../output/png/09_FCI_Channel_Contributions_Integral.png"
ARCHIVE    <- "../output/csv/FCI_Complete_Results.csv"

cat("=== Regenerating Appendix Figure A.2 ===\n\n")

## ---- 1. Evaluate the deterministic prefix of script 01 ---------------------
stopifnot(file.exists(SRC))
src_lines <- readLines(SRC, warn = FALSE)

anchor <- src_lines[PREFIX_END]
if (!grepl("rolling_standardize\\(CONFIG\\$rolling_window\\)", anchor)) {
  stop("PREFIX_END no longer points at the end of the datos_std statement; ",
       "01_FCI_Complete.R has changed. Line ", PREFIX_END, " is:\n  ", anchor)
}

prefix <- paste(src_lines[seq_len(PREFIX_END)], collapse = "\n")
if (grepl("ggsave", prefix)) stop("prefix unexpectedly contains ggsave")

write.csv <- function(...) invisible(NULL)   # neutralize the one prefix write
cat("Evaluating", SRC, "lines 1-", PREFIX_END, "(no estimation, no writes)...\n")
eval(parse(text = prefix), envir = globalenv())
rm(write.csv)

stopifnot(exists("datos_std"), exists("all_vars"), exists("VARIABLES"),
          exists("COLORS"))
cat("  datos_std:", nrow(datos_std), "rows x", ncol(datos_std) - 1, "variables\n\n")

## ---- 2. Contributions, with the divisor corrected to 12 --------------------
var_contributions <- datos_std |>
  dplyr::select(fecha, dplyr::any_of(all_vars)) |>
  na.omit()

n_vars_correct <- length(intersect(all_vars, names(var_contributions)))
n_vars_asshipped <- n_vars_correct - 1L
cat("Divisor: correct =", n_vars_correct,
    "| as shipped in script 01 =", n_vars_asshipped, "\n")

contrib_long <- var_contributions |>
  tidyr::pivot_longer(-fecha, names_to = "Variable", values_to = "Contribution") |>
  dplyr::mutate(Contribution = Contribution / n_vars_correct)

channel_map <- data.frame(
  Variable = c(VARIABLES$rates$vars, VARIABLES$banking$vars, VARIABLES$external$vars),
  Channel  = c(rep("Rates",    length(VARIABLES$rates$vars)),
               rep("Banking",  length(VARIABLES$banking$vars)),
               rep("External", length(VARIABLES$external$vars)))
)

contrib_long <- contrib_long |>
  dplyr::left_join(channel_map, by = "Variable") |>
  dplyr::filter(!is.na(Channel))

total_fci <- contrib_long |>
  dplyr::group_by(fecha) |>
  dplyr::summarise(FCI_Total = sum(Contribution, na.rm = TRUE), .groups = "drop")

## ---- 3. GATE: the corrected total must equal the archived Z-Score ----------
arch <- read.csv(ARCHIVE, stringsAsFactors = FALSE)
arch$fecha <- as.Date(arch$fecha)
cmp <- merge(total_fci, arch[, c("fecha", "FCI_COMP_ZSCORE")], by = "fecha")
cmp <- cmp[!is.na(cmp$FCI_COMP_ZSCORE), ]

d12 <- max(abs(cmp$FCI_Total - cmp$FCI_COMP_ZSCORE))
d11 <- max(abs(cmp$FCI_Total * n_vars_correct / n_vars_asshipped - cmp$FCI_COMP_ZSCORE))

cat("\nGate (", nrow(cmp), "common months):\n", sep = "")
cat("  max |sum(z)/12  - FCI_COMP_ZSCORE| =", format(d12, scientific = TRUE), "\n")
cat("  max |sum(z)/11  - FCI_COMP_ZSCORE| =", format(d11, scientific = TRUE), "\n")
if (d12 > 1e-8) {
  stop("GATE FAILED: reconstructed contributions do not sum to the archived ",
       "Z-Score index. Refusing to overwrite the figure.")
}
cat("  -> PASS: divisor 12 reproduces the published index exactly.\n")
cat("     Peak of plotted line: corrected", round(max(cmp$FCI_Total), 4),
    "vs as-shipped", round(max(cmp$FCI_Total) * 12 / 11, 4), "\n\n")

## ---- 4. Re-plot: identical aesthetics, axis constrained to the data --------
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

channel_contrib <- contrib_long |>
  dplyr::group_by(fecha, Channel) |>
  dplyr::summarise(Contribution = sum(Contribution, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(Channel = factor(Channel, levels = c("Rates", "Banking", "External")))

# Bars are 28 days wide, so a half-bar of margin is the minimum.  Explicit
# breaks (rather than date_breaks) keep the published 1997/1999/.../2025 labels
# and prevent a break being drawn outside the data range.  The lower limit is
# pulled back just far enough to keep the leading 1997 tick inside the panel
# (the series starts 1997-02, since na.omit drops months where any of the
# twelve rolling z-scores is undefined); the upper limit is a half-bar, which
# is what removes the spurious 2027 tick.
x_lo <- min(channel_contrib$fecha) - 45
x_hi <- max(channel_contrib$fecha) + 16
x_breaks <- seq(as.Date("1997-01-01"), as.Date("2025-01-01"), by = "2 years")

p8 <- ggplot(channel_contrib, aes(x = fecha, y = Contribution, fill = Channel)) +
  geom_col(position = "stack", width = 28) +
  geom_line(data = total_fci, aes(x = fecha, y = FCI_Total, fill = NULL),
            color = "black", linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5) +
  scale_fill_manual(values = COLORS$level3) +
  scale_x_date(breaks = x_breaks, date_labels = "%Y",
               limits = c(x_lo, x_hi), expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Channel Contributions to Comprehensive FCI",
    subtitle = "Stacked bars = Channel contributions (sum of variables) | Black line = Total FCI",
    x = NULL, y = "FCI / Contribution", fill = "Channel"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(TARGET_PNG, p8, width = 12, height = 6, dpi = 300)
cat("Axis limits:", format(x_lo), "to", format(x_hi),
    "| breaks:", format(min(x_breaks), "%Y"), "-", format(max(x_breaks), "%Y"), "\n")
cat("Saved:", TARGET_PNG, "\n")
cat("\nNOTE: 01_FCI_Complete.R still carries both defects. Fixing line 1204's\n")
cat("divisor there would also change Plot 7a/7b (variable-level contributions);\n")
cat("this script deliberately regenerates only Appendix Figure A.2.\n")
