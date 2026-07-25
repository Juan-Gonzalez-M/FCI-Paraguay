# ============================================================================
# 58_Small_Robustness_Items.R — Round-11 colleague items S2-S5 + E2 timing.
#
#   Part 1 (S2): log-level / cumulative real-credit LP.  The aggregate
#     outcome elsewhere is YoY growth at t+h; for h < 12 the annual window
#     straddles the shock, so the h ~ 12 peak may partly reflect the growth
#     filter.  Here: y_h = 100*(ln L_{t+h} - ln L_t) (cumulative log
#     response), same controls/inference.
#   Part 2 (S3): predetermined terms-of-trade state.  The paper's headline
#     interaction uses contemporaneous d_ToT; here the state is the
#     12-month trailing MA of d_ToT lagged one month (known before t),
#     interacted with the composite ex-credit index and the one-sided
#     z-score index.
#   Part 3 (S4): matched ex-TCN comparison.  The paper's full-sample
#     comparison removed four variables (TCN + all three external inputs).
#     Here: 8-variable domestic ex-credit vs 7-variable domestic
#     ex-credit/ex-TCN, SAME z-score construction (rolling and expanding),
#     same observations: credit LP + conditional DXY first stage (classical
#     conditional F convention, stated).
#   Part 4 (S5): effective contemporaneous weights by method and sample
#     (regression of each method index on the 12 signed rolling-z inputs),
#     plus targeted sign-reversal sensitivity for the two genuinely
#     ambiguous inputs (Liquidez, Ratio_Cred_Depo) in the one-sided index.
#   Part 5 (E2 timing): exposure-timing checks for the bank gradient —
#     (a) leave-2016-out (2016-average exposure, shocks from Jan 2017);
#     (b) Jan-2016-only exposure (first month), shocks from Feb 2016.
#     Requires the symmetric-deflation panel from the script 30-31 rerun.
#
# ADDITIVE: writes only output/revision/csv/Rev_S*.csv.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
suppressMessages(library(zoo))
set.seed(20260716)
cat("=== 58: small robustness items (S2-S5, E2 timing) ===\n")

POSTIT <- as.Date("2011-05-01")
d <- load_ext_data()

macro <- as.data.table(readxl::read_excel(micro_paths$fci_xlsx, sheet = "Datos_macro"))
macro[, fecha := as.Date(Fecha)]
d <- merge(d, macro[, .(fecha, Creditos_deflactados)], by = "fecha", all.x = TRUE)
d[, lcred := log(Creditos_deflactados)]

rv <- fread("../output/csv/FCI_Robustness_Versions.csv")
rv[, fecha := as.Date(fecha)]
d <- merge(d, rv[, .(fecha, FCI_exCredit_ZSCORE)], by = "fecha", all.x = TRUE)
setorder(d, fecha)
samples <- list(full = d, postIT = d[fecha >= POSTIT])

# ---------------------------------------------------------------------------
# Part 1 (S2): cumulative log-level credit LP
# ---------------------------------------------------------------------------
cat("\n[1] Cumulative log-level credit LP...\n")
lp_cum <- function(data, xvar, h, sample) {
  dh <- as.data.frame(data)
  dh$y_fwd  <- 100 * (dplyr::lead(dh$lcred, h) - dh$lcred)
  dh$y_lag1 <- dplyr::lag(dh$Cred_Real_Total, 1)
  dh$y_lag2 <- dplyr::lag(dh$Cred_Real_Total, 2)
  dh$x_lag1 <- dplyr::lag(dh[[xvar]], 1)
  rhs <- c(xvar, "y_lag1", "y_lag2", "x_lag1", EXOG_MACRO)
  rd <- na.omit(dh[, c("y_fwd", rhs)])
  if (nrow(rd) < 50) return(NULL)
  m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
  ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE))
  data.table(sample = sample, index = xvar, horizon = h, n = nrow(rd),
             coef = ct[xvar, 1], se = ct[xvar, 2], p_value = ct[xvar, 4])
}
cum <- list()
for (s in names(samples)) for (v in c("FCI_exCredit_AVG", "FCI_exCredit_ZSCORE"))
  for (h in c(3, 6, 9, 12, 15, 18, 24))
    cum[[paste(s, v, h)]] <- lp_cum(samples[[s]], v, h, s)
cum <- rbindlist(cum)
write_rev_csv(cum, "Rev_S2_CumLogLevel_LP.csv")
print(cum[horizon %in% c(6, 12, 18)][, .(sample, index, horizon, coef = round(coef, 2),
                                          p = round(p_value, 4), n)])

# ---------------------------------------------------------------------------
# Part 2 (S3): predetermined ToT state interaction
# ---------------------------------------------------------------------------
cat("\n[2] Predetermined ToT-state interaction...\n")
d[, ToT_MA12 := shift(frollmean(d_ToT, 12, align = "right"), 1)]
samples <- list(full = d, postIT = d[fecha >= POSTIT])
int_lp <- function(data, xvar, state, h, sample, statelab) {
  dh <- as.data.frame(data)
  dh$y_fwd  <- dplyr::lead(dh$Cred_Real_Total, h)
  dh$y_lag1 <- dplyr::lag(dh$Cred_Real_Total, 1)
  dh$y_lag2 <- dplyr::lag(dh$Cred_Real_Total, 2)
  dh$x_lag1 <- dplyr::lag(dh[[xvar]], 1)
  dh$FCIxS  <- dh[[xvar]] * dh[[state]]
  rhs <- c(xvar, state, "FCIxS", "y_lag1", "y_lag2", "x_lag1", EXOG_MACRO)
  rd <- na.omit(dh[, c("y_fwd", rhs)])
  if (nrow(rd) < 50) return(NULL)
  m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
  ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE))
  data.table(sample = sample, index = xvar, state = statelab, horizon = h, n = nrow(rd),
             b_fci = ct[xvar, 1], p_fci = ct[xvar, 4],
             b_inter = ct["FCIxS", 1], se_inter = ct["FCIxS", 2],
             p_inter = ct["FCIxS", 4])
}
it <- list()
for (s in names(samples)) for (v in c("FCI_exCredit_AVG", "FCI_exCredit_ZSCORE"))
  for (h in c(6, 12, 18)) {
    it[[paste(s, v, h, "MA")]]  <- int_lp(samples[[s]], v, "ToT_MA12", h, s,
                                          "MA12 d_ToT, lagged (predetermined)")
    it[[paste(s, v, h, "ct")]]  <- int_lp(samples[[s]], v, "d_ToT", h, s,
                                          "contemporaneous d_ToT (paper baseline)")
  }
it <- rbindlist(it)
# crude joint check across the three horizons per spec (Bonferroni on min p)
it[, p_joint_bonf3 := pmin(1, 3 * min(p_inter)), by = .(sample, index, state)]
write_rev_csv(it, "Rev_S3_ToT_Predetermined.csv")
print(it[sample == "full" & horizon == 12][, .(index, state, b_inter = round(b_inter, 3),
                                               p_inter = round(p_inter, 4),
                                               p_joint_bonf3 = round(p_joint_bonf3, 4))])

# ---------------------------------------------------------------------------
# Part 3 (S4): matched ex-TCN comparison (8 vs 7 domestic vars, same pipeline)
# ---------------------------------------------------------------------------
cat("\n[3] Matched ex-TCN comparison...\n")
inputs <- load_fci_inputs()
V8 <- c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM",
        "Ratio_Cred_Depo", "Morosidad", "Rentabilidad", "Liquidez", "TCN")
V7 <- setdiff(V8, "TCN")
aux <- data.table(fecha = inputs$fecha,
  DOM8_ROLL = build_zscore_fci(inputs, V8, std_fn = roll_z),
  DOM7_ROLL = build_zscore_fci(inputs, V7, std_fn = roll_z),
  DOM8_EXP  = build_zscore_fci(inputs, V8, std_fn = expanding_z),
  DOM7_EXP  = build_zscore_fci(inputs, V7, std_fn = expanding_z))
d <- merge(d, aux, by = "fecha", all.x = TRUE)
samples <- list(full = d, postIT = d[fecha >= POSTIT])

fs_conditional <- function(data, yvar, h, sample) {
  dh <- as.data.frame(data)
  dh$y_fwd  <- dplyr::lead(dh$Cred_Real_Total, h)
  dh$y_lag1 <- dplyr::lag(dh$Cred_Real_Total, 1)
  dh$y_lag2 <- dplyr::lag(dh$Cred_Real_Total, 2)
  dh$fci_lag1 <- dplyr::lag(dh[[yvar]], 1)
  ex <- c("y_lag1", "y_lag2", "fci_lag1", EXOG_MACRO, intersect(ER_CTRLS, names(dh)))
  rd <- na.omit(dh[, c("y_fwd", yvar, ex, "DXY")])
  if (nrow(rd) < 50) return(NULL)
  xs <- paste(ex, collapse = " + ")
  fs   <- lm(as.formula(paste(yvar, "~ DXY +", xs)), data = rd)
  fs_r <- lm(as.formula(paste(yvar, "~", xs)), data = rd)
  ct_hac <- lmtest::coeftest(fs, vcov = sandwich::NeweyWest(fs, lag = h + 1, prewhite = FALSE))
  data.table(sample = sample, index = yvar, horizon = h, n = nrow(rd),
             F_classical_conditional = round(anova(fs_r, fs)$F[2], 2),
             t_hac_nw = round(ct_hac["DXY", 3], 2))
}
mx_lp <- list(); mx_fs <- list()
for (s in names(samples)) for (v in c("DOM8_ROLL", "DOM7_ROLL", "DOM8_EXP", "DOM7_EXP")) {
  # matched sample: complete cases on BOTH the 8- and 7-var indices
  base <- copy(samples[[s]])[!is.na(DOM8_ROLL) & !is.na(DOM7_ROLL) &
                             !is.na(DOM8_EXP) & !is.na(DOM7_EXP)]
  r <- run_lp(base, "Cred_Real_Total", v, max_h = 18,
              controls = c("IMAEP_yoy", "IPC_yoy"))
  if (nrow(r)) {
    sdx <- sd(na.omit(base[[v]]))
    mx_lp[[paste(s, v)]] <- data.table(sample = s, index = v,
                                       r[horizon %in% c(6, 12, 18)],
                                       coef_per_sd = r[horizon %in% c(6, 12, 18)]$coef * sdx)
  }
  for (h in c(12)) mx_fs[[paste(s, v, h)]] <- fs_conditional(base, v, h, s)
}
mx_lp <- rbindlist(mx_lp); mx_fs <- rbindlist(mx_fs)
write_rev_csv(mx_lp, "Rev_S4_ExTCN_Matched_LP.csv")
write_rev_csv(mx_fs, "Rev_S4_ExTCN_Matched_FirstStage.csv")
print(mx_lp[horizon == 12][, .(sample, index, coef = round(coef, 2),
                               p = round(p_value, 4), coef_per_sd = round(coef_per_sd, 2))])
print(mx_fs)

# ---------------------------------------------------------------------------
# Part 4 (S5): effective weights + targeted sign-reversal sensitivity
# ---------------------------------------------------------------------------
cat("\n[4] Effective contemporaneous weights by method and sample...\n")
comp <- fread("../output/csv/FCI_Complete_Results.csv")
comp[, fecha := as.Date(fecha)]
Zroll <- as.data.table(sapply(names(FCI_SIGNS),
                              function(v) roll_z(inputs[[v]]) * FCI_SIGNS[v]))
Zroll[, fecha := inputs$fecha]
wd <- merge(comp[, .(fecha, FCI_COMP_ZSCORE, FCI_COMP_PCA, FCI_COMP_VAR,
                     FCI_COMP_DFM, FCI_COMP_AVG)], Zroll, by = "fecha")
eff_w <- list()
for (s in c("full", "postIT")) {
  ws <- if (s == "postIT") wd[fecha >= POSTIT] else wd
  for (m in c("FCI_COMP_ZSCORE", "FCI_COMP_PCA", "FCI_COMP_VAR",
              "FCI_COMP_DFM", "FCI_COMP_AVG")) {
    rd <- na.omit(ws[, c(m, names(FCI_SIGNS)), with = FALSE])
    if (nrow(rd) < 60) next
    mm <- lm(as.formula(paste(m, "~", paste(names(FCI_SIGNS), collapse = "+"))),
             data = rd)
    b <- coef(mm)[names(FCI_SIGNS)]
    eff_w[[paste(s, m)]] <- data.table(
      sample = s, method = m, variable = names(FCI_SIGNS),
      effective_weight = round(unname(b), 4),
      weight_share_abs = round(unname(abs(b) / sum(abs(b))), 4),
      R2 = round(summary(mm)$r.squared, 3), n = nrow(rd))
  }
}
eff_w <- rbindlist(eff_w)
write_rev_csv(eff_w, "Rev_S5_Effective_Weights.csv")
cat("  R2 of weight regressions:\n")
print(unique(eff_w[, .(sample, method, R2)]))

cat("\n[4b] Sign-reversal sensitivity (Liquidez, Ratio_Cred_Depo)...\n")
V11 <- setdiff(names(FCI_SIGNS), "Crecimiento_creditos")   # one-sided ex-credit
flip_variant <- function(flips) {
  Z <- sapply(V11, function(v) {
    sgn <- FCI_SIGNS[v] * ifelse(v %in% flips, -1, 1)
    roll_z(inputs[[v]]) * sgn
  })
  rowMeans(Z, na.rm = FALSE)
}
sv <- data.table(fecha = inputs$fecha,
                 OS_base      = flip_variant(character(0)),
                 OS_flipLiq   = flip_variant("Liquidez"),
                 OS_flipLDR   = flip_variant("Ratio_Cred_Depo"),
                 OS_flipBoth  = flip_variant(c("Liquidez", "Ratio_Cred_Depo")))
d <- merge(d, sv, by = "fecha", all.x = TRUE)
samples <- list(full = d, postIT = d[fecha >= POSTIT])
sr <- list()
for (s in names(samples)) for (v in c("OS_base", "OS_flipLiq", "OS_flipLDR", "OS_flipBoth")) {
  r <- run_lp(samples[[s]], "Cred_Real_Total", v, max_h = 18,
              controls = c("IMAEP_yoy", "IPC_yoy"))
  if (nrow(r)) sr[[paste(s, v)]] <- data.table(sample = s, variant = v,
                                               r[horizon %in% c(6, 12, 18)])
}
sr <- rbindlist(sr)
write_rev_csv(sr, "Rev_S5_SignReversal_LP.csv")
print(sr[horizon == 12][, .(sample, variant, coef = round(coef, 2), p = round(p_value, 4))])

# ---------------------------------------------------------------------------
# Part 5 (E2): exposure-timing checks for the bank FX-credit-exposure gradient
# ---------------------------------------------------------------------------
cat("\n[5] Exposure-timing checks (bank gradient)...\n")
res5 <- tryCatch({
  suppressMessages(library(fixest))
  dA <- read_rds_micro("micro_designA_panel.rds")
  p2 <- read_rds_micro("micro_p2_sector.rds")
  CTRLS <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
  FE_A  <- "bst + bsc"
  grad_run <- function(dd, fx_cut_ym, shock_from_ym, label, B = 999) {
    fx16 <- p2[ym <= fx_cut_ym, .(tot = sum(total)), by = .(sector, cur)]
    fx16 <- dcast(fx16, sector ~ cur, value.var = "tot")
    setnames(fx16, c("6200", "6900"), c("fc", "pyg"))
    fx16[, fx_share := fc / (fc + pyg)]
    dd <- merge(dd, fx16[, .(sector, fx_share)], by = "sector", all.x = TRUE)
    dd[, fxsh_z := zstd(fx_share)]
    dd[, shk_usd_fxsh := shk_usd * fxsh_z]
    ds <- dd[ok12 == TRUE & ym >= shock_from_ym]
    rhs <- c("shk_usd", "shk_usd_fxsh", CTRLS)
    m <- feols(as.formula(paste("g12 ~", paste(rhs, collapse = "+"), "|", FE_A)),
               data = ds, cluster = ~bank, notes = FALSE)
    ct <- coeftable(m)
    p_wcb_sec <- wcb_pval(ds, "g12", rhs, FE_A, "shk_usd_fxsh",
                          cluster = "sector", B = B, weights = "webb")
    data.table(variant = label, n = nobs(m),
               b_gradient = ct["shk_usd_fxsh", 1], se_bank = ct["shk_usd_fxsh", 2],
               p_bank = ct["shk_usd_fxsh", 4], p_wcb = p_wcb_sec)
  }
  et <- rbindlist(list(
    grad_run(copy(dA), to_ym("2016-12-01"), to_ym("2016-01-01"),
             "baseline: 2016-avg exposure, shocks from 2016m1"),
    grad_run(copy(dA), to_ym("2016-12-01"), to_ym("2017-01-01"),
             "leave-2016-out: 2016-avg exposure, shocks from 2017m1"),
    grad_run(copy(dA), to_ym("2016-01-01"), to_ym("2016-02-01"),
             "Jan-2016-only exposure, shocks from 2016m2")))
  write_rev_csv(et, "Rev_S6_Exposure_Timing.csv")
  print(et[, .(variant, n, b_gradient = round(b_gradient, 3),
               p_bank = round(p_bank, 4), p_wcb = round(p_wcb, 4))])
  TRUE
}, error = function(e) { cat("  SKIPPED (panel not ready):", conditionMessage(e), "\n"); FALSE })

cat(sprintf("\n=== 58 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
