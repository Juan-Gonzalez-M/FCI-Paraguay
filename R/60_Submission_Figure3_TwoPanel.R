# ============================================================================
# 60_Submission_Figure3_TwoPanel.R — Round-11 S6: rebuilds submission
# Figure_3.png as a two-panel display: (a) PRIMARY one-sided
# FCI_exCredit z-score credit IRF (post-IT, the declared primary
# specification of Table 4 Panel A; script 52 pipeline conventions:
# NW(h+1), lagged controls); (b) the post-IT regime-specific composite
# FCI_ENDO_exCredit IRF, labeled sensitivity (previous Figure_3).
# 95% Newey-West bands (round-9 convention for principal results).
# ============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(sandwich); library(lmtest) })

d <- read.csv("../output/csv/New_External_Variables.csv")
d$fecha <- as.Date(d$fecha)
rv <- read.csv("../output/csv/FCI_Robustness_Versions.csv")
rv$fecha <- as.Date(rv$fecha)
d <- d %>%
  left_join(rv %>% dplyr::select(fecha, FCI_exCredit_ZSCORE), by = "fecha") %>%
  arrange(fecha) %>%
  mutate(IMAEP_yoy_L1 = lag(IMAEP_yoy, 1), IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
         IPC_yoy_L1   = lag(IPC_yoy, 1),   IPC_yoy_L2   = lag(IPC_yoy, 2)) %>%
  filter(fecha >= as.Date("2011-05-01"))

lp_h <- function(data, xvar, h) {
  dh <- data %>% mutate(y_fwd = lead(Cred_Real_Total, h),
                        y_lag1 = lag(Cred_Real_Total, 1),
                        y_lag2 = lag(Cred_Real_Total, 2),
                        x_lag1 = lag(.data[[xvar]], 1))
  rhs <- c(xvar, "y_lag1", "y_lag2", "x_lag1",
           "IMAEP_yoy","IMAEP_yoy_L1","IMAEP_yoy_L2",
           "IPC_yoy","IPC_yoy_L1","IPC_yoy_L2")
  rd <- dh %>% dplyr::select(y_fwd, all_of(rhs)) %>% na.omit()
  if (nrow(rd) < 50) return(NULL)
  m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
  ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE))
  data.frame(horizon = h, coef = ct[xvar, 1], se = ct[xvar, 2])
}
prim <- bind_rows(lapply(1:18, function(h) lp_h(d, "FCI_exCredit_ZSCORE", h))) %>%
  mutate(panel = "(a) Primary: one-sided FCI (post-IT)")

comp <- read.csv("../output/csv/PostIT_LP_Credit.csv") %>%
  filter(fci_type == "FCI_ENDO_exCredit", credit_type == "Cred_Total",
         horizon <= 18) %>%
  transmute(horizon, coef, se,
            panel = "(b) Sensitivity: regime-specific composite")

both <- bind_rows(prim, comp) %>%
  mutate(lo95 = coef - qnorm(0.975) * se, hi95 = coef + qnorm(0.975) * se)

p <- ggplot(both, aes(horizon, coef)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue4", linewidth = 0.9) +
  geom_point(color = "steelblue4", size = 1.4) +
  facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 18, 3)) +
  labs(x = "Horizon (months)",
       y = "Real credit growth response (pp per index unit)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

ggsave("../output/submission/Figure_3.png", p, width = 11, height = 4.8, dpi = 300)
cat("Figure_3.png rebuilt: primary one-sided IRF + composite sensitivity, 95% NW bands.\n")
print(prim %>% filter(horizon %in% c(6, 12, 18)) %>% mutate(across(where(is.numeric), ~round(., 3))))
