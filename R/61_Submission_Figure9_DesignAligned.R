# ============================================================================
# 61_Submission_Figure9_DesignAligned.R — Round-11 S6: rebuilds submission
# Figure_9.png (bank-level currency-composition margin) with DESIGN-ALIGNED
# 90% intervals: Driscoll-Kraay (time-dependence-robust) bands instead of
# the bank-clustered sensitivity convention.  Panel (a): FX-credit-exposure
# gradient (full profile + ex-COVID point at h=12); panel (b): split-sample
# shock x USD coefficients by sector orientation.  Data: regenerated
# symmetric-deflation outputs (scripts 33, 48).
# ============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
Z90 <- qnorm(0.95)

dc <- read.csv("../output/micro/csv/Micro_DesignC_Hedging_Results.csv")

grad <- dc %>%
  filter(variant == "gradient_fxshare", param == "shk_usd_fxsh", h <= 18) %>%
  transmute(horizon = h, coef = b, se = se_dk,
            panel = "(a) FX-credit-exposure gradient (per SD of exposure)")

# Round-29: panel (b) is drawn from the SPLIT-SAMPLE regressions
# (Micro_DesignC_Split_DK.csv), which are the estimates main-text Table 7 and
# Section 5.3 print (+2.06 / +0.43 at h = 12).  It previously used the
# "sharp" variant of Micro_DesignC_Hedging_Results.csv, whose shk_usd
# coefficient is a different object (+2.07 at h = 12) and whose second
# parameter is a differential rather than the locally-oriented level -- so
# the panel neither matched Table 7 nor plotted the comparison arm the
# caption promises.  Both groups are now shown.
# The binary split is an institutional preselection of export-oriented
# sectors (cattle ranks sixth of 11 by standardized 2016 exposure), so
# exposure language stays with the continuous gradient in panel (a).
split <- read.csv("../output/micro/csv/Micro_DesignC_Split_DK.csv") %>%
  filter(h <= 18) %>%
  transmute(horizon = h, coef = b, se = se_dk,
            series = ifelse(group == "hedged",
                            "Export-oriented (agriculture, cattle)",
                            "Remaining estimation sample"),
            panel = "(b) Split by institutional export orientation: shock x USD")

grad <- grad %>% mutate(series = "FX-credit-exposure gradient")

both <- bind_rows(grad, split) %>%
  mutate(lo = coef - Z90 * se, hi = coef + Z90 * se,
         series = factor(series,
                         levels = c("FX-credit-exposure gradient",
                                    "Export-oriented (agriculture, cattle)",
                                    "Remaining estimation sample")))

pal <- c("FX-credit-exposure gradient" = "darkorange4",
         "Export-oriented (agriculture, cattle)" = "darkorange4",
         "Remaining estimation sample" = "steelblue4")

p <- ggplot(both, aes(horizon, coef, colour = series, fill = series)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.4) +
  facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 18, 3)) +
  scale_colour_manual(values = pal, breaks = levels(both$series)[2:3]) +
  scale_fill_manual(values = pal, breaks = levels(both$series)[2:3]) +
  labs(x = "Horizon (months)",
       y = "Within-cell USD-PYG differential (pp per 1 SD shock)",
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave("../output/submission/Figure_9.png", p, width = 11, height = 4.8, dpi = 300)
cat("Figure_9.png rebuilt with Driscoll-Kraay (design-aligned) 90% bands.\n")
print(both %>% filter(horizon %in% c(6, 12, 18)) %>% mutate(across(where(is.numeric), ~round(., 3))))
