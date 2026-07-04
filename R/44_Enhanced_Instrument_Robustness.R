# ============================================================================
# 44_Enhanced_Instrument_Robustness.R — Instrument robustness battery
# (Coworker 1 WP1-Steps 1-2 + WP6):
#   1. Enhanced purged-DXY: re-purge Delta-ln DXY adding Selic changes and BRL
#      returns (Wu-Xia-substituted FFER when available) to the existing set.
#   2. Broad-dollar robustness: first stage + IV-LP with the Fed broad dollar
#      index instead of DXY (the "DXY is 58% euro" objection).
#   3. EUR placebo instrument: expected NO first stage for Paraguay's FCI.
#   4. Dollar-funding-stress horse race: TED spread vs DXY in the first stage
#      and in the credit LP.
# Consumes Rev_External_Series.csv from script 43.
# ADDITIVE: writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
cat("=== 44: Enhanced instrument robustness ===\n")

d   <- load_ext_data()
ext <- fread(file.path(rev_paths$out_csv, "Rev_External_Series.csv"))
ext[, fecha := as.Date(fecha)]
d <- merge(d, ext, by = "fecha", all.x = TRUE)
setorder(d, fecha)

ENDO <- "FCI_ENDO_exCredit_AVG"
YVAR <- "Cred_Real_Total"

# ---------------------------------------------------------------------------
# 1. Shock construction (monthly innovations) and purge variants
# ---------------------------------------------------------------------------
d[, `:=`(dln_dxy = 100 * c(NA, diff(log(DXY))),
         dln_bd  = 100 * c(NA, diff(log(BroadDollar))),
         dln_brl = 100 * c(NA, diff(log(BRL))),
         dln_eur = 100 * c(NA, diff(log(EUR))),
         dln_eurgbp = 100 * c(NA, diff(log(EUR / GBP))),  # dollar-free cross
         dln_sp  = 100 * c(NA, diff(log(SP500))),
         dln_ipe = 100 * c(NA, diff(log(IPE))),
         d_vix   = c(NA, diff(VIX)),
         d_us10  = c(NA, diff(US_10Y)),
         d_ffer  = c(NA, diff(FFER)),
         d_selic = c(NA, diff(Selic_rate)),
         d_ted   = c(NA, diff(TED)))]
if ("WuXia" %in% names(d)) {
  d[, FFER_wx := fifelse(!is.na(WuXia) & fecha >= as.Date("2009-01-01") &
                           fecha <= as.Date("2015-12-01"), WuXia, FFER)]
  d[, d_ffer_wx := c(NA, diff(FFER_wx))]
} else d[, d_ffer_wx := d_ffer]

ar_res <- function(y, X = NULL, nlag = 3) {
  L <- sapply(1:nlag, function(l) shift(y, l))
  df <- data.frame(y = y, L)
  if (!is.null(X)) df <- cbind(df, X)
  out <- rep(NA_real_, length(y))
  ok <- complete.cases(df)
  if (sum(ok) > 30) out[ok] <- resid(lm(y ~ ., df[ok, ]))
  out
}

purge_v1 <- d[, .(d_vix, dln_ipe, d_ToT, d_ffer, dln_sp, d_us10)]        # script-30 set
purge_v2 <- d[, .(d_vix, dln_ipe, d_ToT, d_ffer_wx, dln_sp, d_us10,
                  d_selic, dln_brl)]                                     # + Brazil, Wu-Xia
d[, `:=`(shk_dxy_raw    = ar_res(dln_dxy),
         shk_dxy_purged = ar_res(dln_dxy, purge_v1),
         shk_dxy_purge2 = ar_res(dln_dxy, purge_v2),
         shk_bd_raw     = ar_res(dln_bd),
         shk_eur        = ar_res(dln_eur),
         shk_eurgbp     = ar_res(dln_eurgbp),
         shk_ted        = ar_res(d_ted))]
cat(sprintf("  purged v1 vs v2 corr: %.3f | DXY vs BroadDollar innovation corr: %.3f\n",
            d[, cor(shk_dxy_purged, shk_dxy_purge2, use = "complete.obs")],
            d[, cor(shk_dxy_raw, shk_bd_raw, use = "complete.obs")]))

# ---------------------------------------------------------------------------
# 2. First-stage battery (levels and innovations), conventional + effective F
# ---------------------------------------------------------------------------
fs_battery <- function(instr, label, controls = TRUE) {
  dd <- as.data.frame(d)
  dd$fci_lag1 <- dplyr::lag(dd[[ENDO]], 1)
  exog <- if (controls) c("fci_lag1", EXOG_MACRO, ER_CTRLS) else NULL
  cols <- c(ENDO, instr, exog)
  reg <- na.omit(dd[, cols])
  if (nrow(reg) < 60) return(NULL)
  f <- as.formula(paste(ENDO, "~", paste(c(instr, exog), collapse = "+")))
  m <- lm(f, reg)
  ct  <- coeftest(m, vcov. = NeweyWest(m, lag = 13, prewhite = FALSE))
  m0  <- if (!is.null(exog)) lm(as.formula(paste(ENDO, "~", paste(exog, collapse = "+"))), reg)
         else lm(as.formula(paste(ENDO, "~ 1")), reg)
  pR2 <- 1 - sum(resid(m)^2) / sum(resid(m0)^2)
  data.table(instrument = label, controls = controls, n = nrow(reg),
             coef = ct[instr, 1], t_nw = ct[instr, 3],
             F_conv = coef(summary(m))[instr, 3]^2, F_eff = ct[instr, 3]^2,
             partial_R2_pct = 100 * pR2)
}

fs <- rbindlist(list(
  fs_battery("DXY",            "DXY level (baseline)"),
  fs_battery("BroadDollar",    "Broad dollar level"),
  fs_battery("shk_dxy_raw",    "DXY innovation (raw)"),
  fs_battery("shk_dxy_purged", "DXY innovation (purged v1)"),
  fs_battery("shk_dxy_purge2", "DXY innovation (purged v2: +Selic,BRL,WuXia)"),
  fs_battery("shk_bd_raw",     "Broad dollar innovation"),
  fs_battery("shk_eur",        "EUR/USD innovation (mirror of the dollar - NOT a placebo)"),
  fs_battery("shk_eurgbp",     "EUR/GBP cross innovation (dollar-free PLACEBO)"),
  fs_battery("shk_ted",        "TED spread innovation (funding stress)")),
  fill = TRUE)
cat("\nFirst-stage battery (endo = FCI_ENDO_exCredit, with controls):\n")
print(fs[, .(instrument, n, coef = round(coef, 4), t_nw = round(t_nw, 2),
             F_eff = round(F_eff, 1), pR2 = round(partial_R2_pct, 1))])
write_rev_csv(fs, "Rev_Instrument_Comparison.csv")

# ---------------------------------------------------------------------------
# 3. IV-LP: broad dollar vs DXY (levels convention of script 23) + purged v2
# ---------------------------------------------------------------------------
iv_specs <- list(
  c(instr = "DXY",            label = "DXY (baseline)"),
  c(instr = "BroadDollar",    label = "Broad dollar"),
  c(instr = "shk_dxy_raw",    label = "DXY innovation (raw)"),
  c(instr = "shk_dxy_purged", label = "Purged DXY v1"),
  c(instr = "shk_dxy_purge2", label = "Purged DXY v2"))
iv_out <- list()
for (s in iv_specs) {
  for (h in c(6, 12, 18)) {
    r <- iv_lp_h(d, YVAR, ENDO, s["instr"], h, ar_grid = seq(-300, 150, by = 1))
    if (!is.null(r)) iv_out[[length(iv_out) + 1]] <- cbind(instrument = s["label"], r)
  }
}
iv_out <- rbindlist(iv_out)
cat("\nIV-LP across instruments:\n")
print(iv_out[, .(instrument, h = horizon, b = round(coef, 2), p = round(p_value, 3),
                 F_eff = round(eff_F, 1), ar_lo = round(ar_lo, 1),
                 ar_hi = round(ar_hi, 1))])
write_rev_csv(iv_out, "Rev_BroadDollar_IV_LP.csv")

# ---------------------------------------------------------------------------
# 4. Funding-stress horse race: DXY vs TED in the credit LP
# ---------------------------------------------------------------------------
hr <- list()
for (h in c(6, 12, 18)) {
  dh <- as.data.frame(d)
  dh$y_fwd <- dplyr::lead(dh[[YVAR]], h)
  dh$y_l1 <- dplyr::lag(dh[[YVAR]], 1); dh$y_l2 <- dplyr::lag(dh[[YVAR]], 2)
  reg <- na.omit(dh[, c("y_fwd", "shk_dxy_purged", "shk_ted", "y_l1", "y_l2",
                        EXOG_MACRO)])
  m <- lm(y_fwd ~ ., reg)
  ct <- coeftest(m, vcov. = NeweyWest(m, lag = h + 1, prewhite = FALSE))
  for (v in c("shk_dxy_purged", "shk_ted"))
    hr[[length(hr) + 1]] <- data.table(horizon = h, shock = v, coef = ct[v, 1],
                                       se = ct[v, 2], p_value = ct[v, 4],
                                       n = nrow(reg))
}
hr <- rbindlist(hr)
cat("\nHorse race (credit LP, joint):\n")
print(dcast(hr, horizon ~ shock, value.var = c("coef", "p_value")))
write_rev_csv(hr, "Rev_FundingStress_HorseRace.csv")

# ---------------------------------------------------------------------------
# 5. Figures
# ---------------------------------------------------------------------------
fsp <- fs[!is.na(F_eff)]
fsp[, ilab := factor(instrument, levels = rev(instrument))]
p1 <- ggplot(fsp, aes(F_eff, ilab, fill = grepl("PLACEBO", instrument))) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c(`FALSE` = "steelblue4", `TRUE` = "grey60")) +
  geom_vline(xintercept = 10, linetype = 2, color = "firebrick") +
  labs(title = "First-stage strength across dollar instruments (effective F)",
       subtitle = "Endogenous: FCI_ENDO_exCredit; dashed line = conventional F = 10 threshold; grey = EUR placebo",
       x = "Effective (robust) first-stage F", y = NULL) +
  theme_micro()
save_rev_png(p1, "287_Instrument_Comparison.png", w = 10, h = 5.5)

p2 <- ggplot(iv_out, aes(factor(horizon), coef, color = instrument)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.5)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "IV-LP credit response: DXY vs broad dollar vs enhanced purged instrument",
       subtitle = "2SLS with script-23 controls; 90% CIs (Newey-West)",
       x = "Horizon (months)", y = "Credit response (pp)", color = NULL) +
  theme_micro()
save_rev_png(p2, "288_BroadDollar_vs_DXY_IVLP.png")

cat(sprintf("=== 44 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
