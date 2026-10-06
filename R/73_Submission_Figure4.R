# ============================================================================
# 73_Submission_Figure4.R — Round-23 colleague item.
#
# PROBLEM.  Submission Figure 4 (\label{fig:loo}) is Figure_10.png, a byte-for-
# byte copy of output/revision/png/286_FCI_LOO_Forest.png produced by
# 42_FCI_Component_Exclusion.R L116-127.  That script draws an embedded title
# "Credit effect at h = 12 across FCI component-exclusion variants" and a
# subtitle duplicating the caption.  The manuscript caption was recalibrated to
# "h = 12 Credit Coefficient Across FCI Component-Exclusion Variants", so the
# shipped image contradicts both the caption and the response letter's
# terminology paragraph.
#
# Re-running script 42 is the wrong tool: it re-estimates the whole leave-one-out
# battery from the Excel source.  This script is PLOTTING ONLY — it reads the
# archived battery CSV and redraws the identical forest plot with no embedded
# title or subtitle, following the convention set in 62_Submission_Figure1.R:
# journals typeset the caption, and an internal title duplicates it.
#
# Source: output/revision/csv/Rev_FCI_LOO_Battery.csv (written by script 42).
# Output: ../output/submission/Figure_10.png (manuscript Figure 4)
# ============================================================================

suppressMessages(library(ggplot2))

Z90 <- qnorm(0.95)   # revision_helpers.R L31 — 90% intervals, as the caption says

tab <- read.csv("../output/revision/csv/Rev_FCI_LOO_Battery.csv",
                stringsAsFactors = FALSE)
f12 <- tab[tab$horizon == 12, ]

# Labels and display order verbatim from 42_FCI_Component_Exclusion.R L107-115
labs_map <- c(
  zs_exCredit            = "exCredit (z-score baseline)",
  zs_exCredit_exNPL      = "- NPL",
  zs_exCredit_exROE      = "- ROE",
  zs_exCredit_exLiquidez = "- Liquidity",
  zs_exCredit_exRatioCD  = "- Credit/deposit ratio",
  zs_exCredit_exTCN      = "- Exchange rate",
  zs_exCredit_QtyLag1    = "Banking quantities at t-1",
  FCI_exCredit_AVG       = "exCredit (published 4-method)",
  FCI_RATES_AVG          = "Rates-only FCI (co-baseline)")

stopifnot(setequal(f12$variant, names(labs_map)))
f12 <- f12[match(names(labs_map), f12$variant), ]
f12$vlab <- factor(labs_map[f12$variant], levels = rev(labs_map))

# ---- Gate: the plotted values must equal the archived estimates ------------
# Guards against a silently changed CSV. Values are transcribed, never computed.
expect <- c(zs_exCredit = -10.2202288660727, zs_exCredit_exNPL = -10.7108232577062,
            zs_exCredit_exROE = -11.7317839671018, zs_exCredit_exLiquidez = -10.0460317315402,
            zs_exCredit_exRatioCD = -9.16213774918637, zs_exCredit_exTCN = -7.21733860718106,
            zs_exCredit_QtyLag1 = -12.140791366752, FCI_exCredit_AVG = -7.65628048038673,
            FCI_RATES_AVG = -5.43944476708963)
stopifnot(max(abs(f12$coef - expect[f12$variant])) < 1e-10)
cat(sprintf("Gate passed: 9 h=12 coefficients match the archive (max dev %.2e).\n",
            max(abs(f12$coef - expect[f12$variant]))))

baseline <- f12$coef[f12$variant == "zs_exCredit"]

p <- ggplot(f12, aes(coef, vlab)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = baseline, colour = "steelblue", alpha = 0.5) +
  geom_errorbarh(aes(xmin = coef - Z90 * se, xmax = coef + Z90 * se),
                 height = 0.25, colour = "grey40") +
  geom_point(size = 2.4, colour = "steelblue4") +
  labs(x = "Coefficient (pp per FCI unit)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_blank(), plot.subtitle = element_blank())

ggsave("../output/submission/Figure_10.png", p,
       width = 10, height = 6.5, dpi = 150, bg = "white")
cat("Figure_10.png rewritten: 9 variants, 90% NW intervals, no embedded title.\n")
