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

split <- dc %>%
  filter(variant == "sharp", h <= 18) %>%
  mutate(group = ifelse(param == "shk_usd", "High-FX-credit",
                        "Differential: locally oriented - high-FX-credit")) %>%
  filter(param == "shk_usd") %>%
  # Round-14: sectors are classified by their 2016 FX-credit share, not by
  # observed export income, so the label names the classification variable.
  transmute(horizon = h, coef = b, se = se_dk,
            panel = "(b) High-FX-credit sectors (agriculture, cattle): shock x USD")

both <- bind_rows(grad, split) %>%
  mutate(lo = coef - Z90 * se, hi = coef + Z90 * se)

p <- ggplot(both, aes(horizon, coef)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "darkorange3", alpha = 0.22) +
  geom_line(color = "darkorange4", linewidth = 0.9) +
  geom_point(color = "darkorange4", size = 1.4) +
  facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 18, 3)) +
  labs(x = "Horizon (months)",
       y = "Within-cell USD-PYG differential (pp per 1 SD shock)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

ggsave("../output/submission/Figure_9.png", p, width = 11, height = 4.8, dpi = 300)
cat("Figure_9.png rebuilt with Driscoll-Kraay (design-aligned) 90% bands.\n")
print(both %>% filter(horizon %in% c(6, 12, 18)) %>% mutate(across(where(is.numeric), ~round(., 3))))
