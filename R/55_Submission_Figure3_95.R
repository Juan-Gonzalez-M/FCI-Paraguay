# ============================================================================
# 55_Submission_Figure3_95.R — Regenerates submission Figure_3.png (main-text
# credit-channel IRF, post-IT re-estimated FCI_ENDO_exCredit -> total real
# credit) with 95% Newey-West bands (round-9 convention: principal results at
# 95%; 90% retained for AR sets and exploratory profiles).
# Source data: output/csv/PostIT_LP_Credit.csv (script 15).
# ============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

lp <- read.csv("../output/csv/PostIT_LP_Credit.csv") %>%
  filter(fci_type == "FCI_ENDO_exCredit", credit_type == "Cred_Total",
         horizon <= 18) %>%
  mutate(lo95 = coef - qnorm(0.975) * se, hi95 = coef + qnorm(0.975) * se)

p <- ggplot(lp, aes(horizon, coef)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue4", linewidth = 0.9) +
  geom_point(color = "steelblue4", size = 1.6) +
  scale_x_continuous(breaks = seq(0, 18, 3)) +
  labs(x = "Horizon (months)",
       y = "Real credit growth response (pp per index unit)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("../output/submission/Figure_3.png", p, width = 8, height = 5, dpi = 300)
cat("Figure_3.png regenerated with 95% Newey-West bands.\n")
