# ============================================================================
# 38_CHR_Plausibly_Exogenous.R — Conley, Hansen & Rossi (2012, REStat)
# "plausibly exogenous" bounds for the DXY-only IV-LP (Coworker 1, WP1-Step 3).
#
# Instead of defending a zero direct effect of DXY on credit (trade/commodity
# channels), allow a direct effect gamma and trace the 2SLS credit-channel
# coefficient as gamma ranges over [gamma_min, 0]. Report the breakdown value
# gamma* (smallest |gamma| at which the 90% CI first includes 0) and compare
# it with a data-calibrated direct effect:
#   gamma_cal = (DXY -> ToT growth) x (ToT growth -> credit at h)   [channel 1]
#   plus a quarterly (DXY -> export volume growth) gauge            [channel 2]
#
# Spec mirrors script 23's DXY-only IV-LP exactly (endo FCI_ENDO_exCredit_AVG,
# instrument DXY level, outcome Cred_Real_Total, macro + ER controls).
# ADDITIVE: writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
cat("=== 38: Conley-Hansen-Rossi plausibly-exogenous bounds ===\n")

d <- load_ext_data()
H_SET <- c(6, 12, 18)
ENDO  <- "FCI_ENDO_exCredit_AVG"
YVAR  <- "Cred_Real_Total"
INSTR <- "DXY"

# ---------------------------------------------------------------------------
# 1. Reduced-form scale: total effect of DXY on y_{t+h} (per DXY point)
#    gamma grid runs from 0 to the full reduced form (attributing everything
#    to the direct channel)
# ---------------------------------------------------------------------------
rf <- list()
for (h in H_SET) {
  dh <- as.data.frame(d)
  dh$y_fwd <- dplyr::lead(dh[[YVAR]], h)
  dh$y_lag1 <- dplyr::lag(dh[[YVAR]], 1); dh$y_lag2 <- dplyr::lag(dh[[YVAR]], 2)
  dh$fci_lag1 <- dplyr::lag(dh[[ENDO]], 1)
  exog <- c("y_lag1", "y_lag2", "fci_lag1", EXOG_MACRO, ER_CTRLS)
  reg <- na.omit(dh[, c("y_fwd", exog, INSTR)])
  m <- lm(as.formula(paste("y_fwd ~", INSTR, "+", paste(exog, collapse = "+"))), reg)
  ct <- coeftest(m, vcov. = NeweyWest(m, lag = h + 1, prewhite = FALSE))
  rf[[length(rf) + 1]] <- data.table(horizon = h, rf_coef = ct[INSTR, 1],
                                     rf_se = ct[INSTR, 2])
}
rf <- rbindlist(rf)
cat("Reduced form (credit on DXY level, per DXY point):\n"); print(rf)

# ---------------------------------------------------------------------------
# 2. CHR bounds: 2SLS on (y - gamma*DXY) over the gamma grid
# ---------------------------------------------------------------------------
res <- list()
for (h in H_SET) {
  g_end <- rf[horizon == h, rf_coef]           # full reduced form
  grid  <- seq(0, g_end, length.out = 41)      # 0 -> total attribution
  for (g in grid) {
    r <- iv_lp_h(d, YVAR, ENDO, INSTR, h, gamma = g,
                 ar_grid = seq(-300, 150, by = 1))
    if (!is.null(r)) res[[length(res) + 1]] <- r
  }
  cat(sprintf("  h=%d done (gamma grid 0 .. %.3f)\n", h, g_end))
}
res <- rbindlist(res)

# Breakdown gamma*: smallest |gamma| where the 90% CI includes 0
bd <- res[, {
  s <- .SD[order(abs(gamma))]
  cross <- s[ci_lower <= 0 & ci_upper >= 0]
  .(beta_gamma0 = s$coef[1], ci_lo_gamma0 = s$ci_lower[1], ci_hi_gamma0 = s$ci_upper[1],
    gamma_breakdown = if (nrow(cross)) cross$gamma[1] else NA_real_,
    gamma_full_rf = max(abs(gamma)) * sign(gamma[which.max(abs(gamma))]),
    first_stage_F = s$first_stage_F[1], eff_F = s$eff_F[1])
}, by = horizon]

# ---------------------------------------------------------------------------
# 3. Calibrate the plausible direct effect from own data
# ---------------------------------------------------------------------------
# Channel 1 (monthly): DXY -> ToT growth, times ToT growth -> credit at h
m1 <- lm(d_ToT ~ DXY + IMAEP_yoy + IPC_yoy, data = d)
b_dxy_tot <- coeftest(m1, vcov. = NeweyWest(m1, lag = 13, prewhite = FALSE))["DXY", 1]
cal <- list()
for (h in H_SET) {
  lp_tot <- run_lp(d, YVAR, "d_ToT", max_h = h, controls = c("IMAEP_yoy", "IPC_yoy"))
  b_tot_cred <- lp_tot[horizon == h, coef]
  g1 <- b_dxy_tot * b_tot_cred
  cal[[length(cal) + 1]] <- data.table(horizon = h, channel = "DXY->ToT->credit",
                                       gamma_calibrated = g1)
}
# Channel 2 (quarterly gauge): DXY -> export volume growth (Quarterly_SA)
qsa <- as.data.table(suppressMessages(read_excel(micro_paths$fci_xlsx,
                                                 sheet = "Quarterly_SA")))
qsa[, fecha := as.Date(Date)]
setorder(qsa, fecha)
qsa[, exp_yoy := 100 * (Exportaciones / shift(Exportaciones, 4) - 1)]
dq <- d[, .(qtr = as.Date(cut(fecha, "quarter")), DXY)][
  , .(DXY = mean(DXY, na.rm = TRUE)), by = qtr]
qm <- merge(qsa[, .(qtr = as.Date(cut(fecha, "quarter")), exp_yoy)], dq, by = "qtr")
m2 <- lm(exp_yoy ~ DXY, qm)
b_dxy_exp <- coeftest(m2, vcov. = NeweyWest(m2, lag = 5, prewhite = FALSE))["DXY", ]
cal <- rbindlist(cal)

cat("\nCHR breakdown values vs calibrated direct effects:\n")
comp <- merge(bd, cal, by = "horizon")
comp[, ratio_cal_to_breakdown := abs(gamma_calibrated) / abs(gamma_breakdown)]
print(comp[, .(horizon, beta_gamma0 = round(beta_gamma0, 2),
               gamma_breakdown = round(gamma_breakdown, 4),
               gamma_calibrated = round(gamma_calibrated, 4),
               ratio = round(ratio_cal_to_breakdown, 3))])
cat(sprintf("\nQuarterly gauge: exports_yoy on DXY: b=%.3f (p=%.3f)\n",
            b_dxy_exp[1], b_dxy_exp[4]))

write_rev_csv(res,  "Rev_CHR_Bounds_Full.csv")
write_rev_csv(comp, "Rev_CHR_Bounds.csv")

# ---------------------------------------------------------------------------
# 4. Figure: beta bounds as a function of assumed direct effect (h = 12)
# ---------------------------------------------------------------------------
p12 <- res[horizon == 12]
gcal12 <- cal[horizon == 12, gamma_calibrated]
pC <- ggplot(p12, aes(gamma, coef)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue4", linewidth = 1) +
  geom_vline(xintercept = gcal12, color = "firebrick", linewidth = 0.8) +
  geom_vline(xintercept = bd[horizon == 12, gamma_breakdown],
             color = "grey30", linetype = 3, linewidth = 0.8) +
  annotate("text", x = gcal12, y = max(p12$ci_upper), hjust = -0.05, vjust = 1,
           label = "calibrated trade channel", color = "firebrick", size = 3.2) +
  annotate("text", x = bd[horizon == 12, gamma_breakdown], y = min(p12$ci_lower),
           hjust = -0.05, vjust = 0, label = "breakdown gamma*",
           color = "grey30", size = 3.2) +
  labs(title = "Plausibly-exogenous bounds (Conley-Hansen-Rossi 2012), h = 12",
       subtitle = "2SLS credit coefficient allowing DXY a direct effect gamma on credit; 90% CI (Newey-West)",
       x = "Assumed direct effect of DXY on credit growth (pp per DXY point)",
       y = "Credit-channel coefficient (pp per FCI unit)") +
  theme_micro()
save_rev_png(pC, "282_CHR_Bounds.png")

cat(sprintf("=== 38 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
