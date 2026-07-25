# ============================================================================
# 56_Aligned_Identification_Chain.R — Round-11 colleague item E1.
#
# Both round-11 colleague audits identify the same core publication risk:
# the paper's persistence evidence (K.1d), headline OLS projections
# (Table 4 Panel A), and IV/AR/Conley machinery (Tables 5-7) sit on THREE
# different index objects.  Required correction: declare ONE causal index
# ex ante and run the entire identification sequence on that same object.
#
# Declared index (colleague-recommended): FCI_OS — the strictly one-sided,
# credit-excluded, TCN-excluded domestic z-score index.  Normalization:
# expanding window (strictly backward-looking AND persistence-preserving,
# so the levels chain is well-posed; the trailing-60m rolling variant is
# reported as a DIFFERENT ESTIMAND, not an alternative normalization).
# Identical construction to script 50 Block 6 (FCI_EXP_ENDO_exC_exTCN).
#
#   Part 0  construction + correlations with existing objects
#   Part 1  unit-root battery for the index (both samples)
#   Part 2  ARDL/PSS bounds persistence assessment (BIC + whitened refit,
#           exact finite-sample CVs) — same-index levels relation
#   Part 3  conditional DXY first stage (classical + HAC conventions),
#           h = 6/12/18, full and post-IT
#   Part 4  reduced form, 2SLS = RF/FS, AR 90% sets (iv_lp_h), h=6/12/18/24
#   Part 5  Conley/CHR breakdown gamma* on the aligned index; compared with
#           the calibrated trade-channel gamma (index-invariant a*b legs
#           from script 38/49, read from Rev_CHR_Calibration_Interval.csv)
#   Part 6  credit LP h=1..18 on FCI_OS (expanding) and the rolling
#           one-sided variant, full + post-IT, per-unit and per-SD
#   Part 7  S1: null-imposed moving-block bootstrap p-value for the DXY
#           reduced form (persistence-preserving: DXY* resampled as a
#           random walk from block-resampled increments; residuals
#           block-resampled with L=18 to respect the overlap structure)
#   Part 8  verdict table (decision point: causal vs predictive framing)
#
# ADDITIVE: writes only output/revision/csv/Rev_Aligned_*.csv.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
suppressMessages({library(tseries); library(ARDL)})
set.seed(20260716)
cat("=== 56: Aligned identification chain (one-sided ex-credit ex-TCN domestic index) ===\n")

POSTIT <- as.Date("2011-05-01")

# ---------------------------------------------------------------------------
# Part 0 — construction
# ---------------------------------------------------------------------------
cat("\n[0] Constructing FCI_OS (expanding + rolling one-sided, 7 domestic vars)...\n")
d <- load_ext_data()
inputs <- load_fci_inputs()
VARS_OS <- c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM",
             "Ratio_Cred_Depo", "Morosidad", "Rentabilidad", "Liquidez")
aux <- data.table(
  fecha        = inputs$fecha,
  FCI_OS_EXP   = build_zscore_fci(inputs, VARS_OS, std_fn = expanding_z),
  FCI_OS_ROLL  = build_zscore_fci(inputs, VARS_OS, std_fn = roll_z))
d <- merge(d, aux, by = "fecha", all.x = TRUE)

exp_fci <- fread(file.path(rev_paths$out_csv, "Rev_ExpandingFCI_Series.csv"))
exp_fci[, fecha := as.Date(fecha)]
d <- merge(d, exp_fci[, .(fecha, FCI_EXP_ENDO_exCredit)], by = "fecha", all.x = TRUE)
rv <- fread("../output/csv/FCI_Robustness_Versions.csv")
rv[, fecha := as.Date(fecha)]
d <- merge(d, rv[, .(fecha, FCI_exCredit_ZSCORE)], by = "fecha", all.x = TRUE)
setorder(d, fecha)
samples <- list(full = d, postIT = d[fecha >= POSTIT])

cons <- data.table(
  pair = c("FCI_OS_EXP vs FCI_OS_ROLL",
           "FCI_OS_EXP vs composite rolling FCI_ENDO_exCredit_AVG",
           "FCI_OS_EXP vs expanding composite-domestic (script 39)",
           "FCI_OS_EXP vs paper one-sided comprehensive FCI_exCredit_ZSCORE",
           "FCI_OS_ROLL vs composite rolling FCI_ENDO_exCredit_AVG"),
  correlation = c(
    cor(d$FCI_OS_EXP, d$FCI_OS_ROLL, use = "complete.obs"),
    cor(d$FCI_OS_EXP, d$FCI_ENDO_exCredit_AVG, use = "complete.obs"),
    cor(d$FCI_OS_EXP, d$FCI_EXP_ENDO_exCredit, use = "complete.obs"),
    cor(d$FCI_OS_EXP, d$FCI_exCredit_ZSCORE, use = "complete.obs"),
    cor(d$FCI_OS_ROLL, d$FCI_ENDO_exCredit_AVG, use = "complete.obs")))
write_rev_csv(cons, "Rev_Aligned_Construction.csv")
print(cons)

# ---------------------------------------------------------------------------
# Part 1 — unit roots
# ---------------------------------------------------------------------------
cat("\n[1] Unit roots...\n")
ur_row <- function(x, label, sample) {
  x <- na.omit(x); if (length(x) < 40) return(NULL)
  a <- suppressWarnings(adf.test(x)); p <- suppressWarnings(pp.test(x))
  k <- suppressWarnings(kpss.test(x, null = "Level"))
  data.table(sample = sample, series = label, n = length(x),
             adf_p = round(a$p.value, 3), pp_p = round(p$p.value, 3),
             kpss_p = round(k$p.value, 3),
             verdict = ifelse(a$p.value > 0.10 & k$p.value <= 0.05, "I(1)-consistent",
                       ifelse(a$p.value <= 0.05 & k$p.value > 0.05, "I(0)-consistent",
                              "mixed/borderline")))
}
ur <- list()
for (s in names(samples)) {
  ds <- samples[[s]]
  ur[[paste(s, 1)]] <- ur_row(ds$FCI_OS_EXP,  "FCI_OS_EXP (level)", s)
  ur[[paste(s, 2)]] <- ur_row(c(NA, diff(ds$FCI_OS_EXP)), "FCI_OS_EXP (d1)", s)
  ur[[paste(s, 3)]] <- ur_row(ds$FCI_OS_ROLL, "FCI_OS_ROLL (level)", s)
  ur[[paste(s, 4)]] <- ur_row(ds$DXY, "DXY (level)", s)
}
ur <- rbindlist(ur)
write_rev_csv(ur, "Rev_Aligned_UnitRoots.csv")
print(ur)

# ---------------------------------------------------------------------------
# Part 2 — ARDL/PSS bounds on the SAME index (BIC + whitened, exact CVs)
# ---------------------------------------------------------------------------
cat("\n[2] ARDL bounds (same-index persistence assessment)...\n")
ardl_block <- function(ds, yv, xv, sample, case = 3) {
  cc <- as.data.frame(na.omit(ds[, c(yv, xv), with = FALSE]))
  if (nrow(cc) < 80) return(NULL)
  f <- as.formula(paste(yv, "~", paste(xv, collapse = "+")))
  fit <- tryCatch(ARDL::auto_ardl(f, data = cc, max_order = 8, selection = "BIC"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  bm <- fit$best_model
  bF  <- tryCatch(ARDL::bounds_f_test(bm, case = case), error = function(e) NULL)
  bFe <- tryCatch(ARDL::bounds_f_test(bm, case = case, exact = TRUE), error = function(e) NULL)
  bT  <- tryCatch(ARDL::bounds_t_test(bm, case = case), error = function(e) NULL)
  ue  <- tryCatch(ARDL::uecm(bm), error = function(e) NULL)
  bg12 <- tryCatch(lmtest::bgtest(ue, order = 12)$p.value, error = function(e) NA)
  mult <- tryCatch(as.data.frame(ARDL::multipliers(bm, type = "lr")), error = function(e) NULL)
  dxy <- if (!is.null(mult)) mult[grepl("DXY", mult[[1]]), ] else NULL
  data.table(sample = sample, y = yv, x = paste(xv, collapse = "+"),
             spec = "BIC-selected",
             order = paste(fit$best_order, collapse = ","), n = nrow(cc),
             F_stat = if (!is.null(bF)) round(as.numeric(bF$statistic), 3) else NA,
             F_p_asy = if (!is.null(bF)) round(as.numeric(bF$p.value), 4) else NA,
             F_p_exact = if (!is.null(bFe)) round(as.numeric(bFe$p.value), 4) else NA,
             t_stat = if (!is.null(bT)) round(as.numeric(bT$statistic), 3) else NA,
             t_p = if (!is.null(bT)) round(as.numeric(bT$p.value), 4) else NA,
             bg_p_order12 = round(bg12, 3),
             lr_dxy = if (!is.null(dxy) && nrow(dxy)) round(as.numeric(dxy[[2]]), 5) else NA,
             lr_dxy_p = if (!is.null(dxy) && nrow(dxy)) round(as.numeric(dxy[[5]]), 4) else NA)
}
whiten_fit <- function(ds, yv, xv, sample) {
  cc <- as.data.frame(na.omit(ds[, c(yv, xv), with = FALSE]))
  if (nrow(cc) < 80) return(NULL)
  f <- as.formula(paste(yv, "~", paste(xv, collapse = "+")))
  grid <- expand.grid(p = 2:12, q = 0:4)
  grid <- grid[order(grid$p + grid$q), ]
  best <- NULL; best_bg <- -Inf; best_ord <- c(NA, NA)
  for (i in seq_len(nrow(grid))) {
    m <- tryCatch(ARDL::ardl(f, data = cc, order = c(grid$p[i], grid$q[i])),
                  error = function(e) NULL)
    if (is.null(m)) next
    bg <- tryCatch(lmtest::bgtest(ARDL::uecm(m), order = 12)$p.value,
                   error = function(e) NA)
    if (is.na(bg)) next
    if (bg > best_bg) { best <- m; best_bg <- bg; best_ord <- c(grid$p[i], grid$q[i]) }
    if (bg > 0.10) break
  }
  if (is.null(best)) return(NULL)
  bF  <- ARDL::bounds_f_test(best, case = 3)
  bFe <- tryCatch(ARDL::bounds_f_test(best, case = 3, exact = TRUE), error = function(e) NULL)
  bT  <- ARDL::bounds_t_test(best, case = 3)
  mult <- as.data.frame(ARDL::multipliers(best, type = "lr"))
  dxy <- mult[grepl("DXY", mult[[1]]), ]
  data.table(sample = sample, y = yv, x = paste(xv, collapse = "+"),
             spec = "whitened (min lags passing BG12)",
             order = paste(best_ord, collapse = ","), n = nrow(cc),
             F_stat = round(as.numeric(bF$statistic), 3),
             F_p_asy = round(as.numeric(bF$p.value), 4),
             F_p_exact = if (!is.null(bFe)) round(as.numeric(bFe$p.value), 4) else NA,
             t_stat = round(as.numeric(bT$statistic), 3),
             t_p = round(as.numeric(bT$p.value), 4),
             bg_p_order12 = round(best_bg, 3),
             lr_dxy = round(as.numeric(dxy[[2]]), 5),
             lr_dxy_p = round(as.numeric(dxy[[5]]), 4))
}
ab <- list()
for (s in names(samples)) {
  ab[[paste(s, "bic")]]  <- ardl_block(samples[[s]], "FCI_OS_EXP", "DXY", s)
  ab[[paste(s, "wht")]]  <- whiten_fit(samples[[s]], "FCI_OS_EXP", "DXY", s)
  # rolling variant shown as a different estimand (expected: no levels relation)
  ab[[paste(s, "roll")]] <- ardl_block(samples[[s]], "FCI_OS_ROLL", "DXY", s)
}
ab <- rbindlist(ab, fill = TRUE)
write_rev_csv(ab, "Rev_Aligned_ARDL_Bounds.csv")
print(ab[, .(sample, y, spec, order, F_stat, F_p_asy, F_p_exact, t_p, bg_p_order12, lr_dxy, lr_dxy_p)])

# ---------------------------------------------------------------------------
# Part 3 — conditional DXY first stage on the aligned index
# ---------------------------------------------------------------------------
cat("\n[3] Conditional DXY first stage...\n")
fs_diag <- function(data, yvar, h, sample) {
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
  f_cls <- anova(fs_r, fs)$F[2]
  t_cls <- summary(fs)$coefficients["DXY", 3]
  pr2_cls <- t_cls^2 / (t_cls^2 + df.residual(fs))
  ct_hac <- lmtest::coeftest(fs, vcov = sandwich::NeweyWest(fs, lag = h + 1, prewhite = FALSE))
  data.table(sample = sample, index = yvar, horizon = h, n = nrow(rd),
             b_fs = round(coef(fs)["DXY"], 5),
             t_classical = round(t_cls, 2), F_classical = round(f_cls, 1),
             partialR2_classical = round(pr2_cls, 4),
             t_hac_nw = round(ct_hac["DXY", 3], 2),
             F_hac_nw = round(ct_hac["DXY", 3]^2, 1))
}
fs_out <- list()
for (s in names(samples)) for (h in c(6, 12, 18))
  for (v in c("FCI_OS_EXP", "FCI_OS_ROLL"))
    fs_out[[paste(s, v, h)]] <- fs_diag(samples[[s]], v, h, s)
fs_out <- rbindlist(fs_out)
write_rev_csv(fs_out, "Rev_Aligned_FirstStage.csv")
print(fs_out)

# ---------------------------------------------------------------------------
# Part 4 — reduced form, 2SLS, AR sets on the aligned index
# ---------------------------------------------------------------------------
cat("\n[4] IV-LP / AR sets (aligned index)...\n")
iv_out <- list()
for (s in names(samples)) for (h in c(6, 12, 18, 24)) {
  r <- iv_lp_h(samples[[s]], "Cred_Real_Total", "FCI_OS_EXP", "DXY", h,
               ar_grid = seq(-600, 300, by = 1))
  if (!is.null(r)) iv_out[[paste(s, h)]] <- data.table(sample = s, r)
}
iv_out <- rbindlist(iv_out)
# reduced form itself (per DXY point), with NW inference
rf_diag <- function(data, h, sample) {
  dh <- as.data.frame(data)
  dh$y_fwd  <- dplyr::lead(dh$Cred_Real_Total, h)
  dh$y_lag1 <- dplyr::lag(dh$Cred_Real_Total, 1)
  dh$y_lag2 <- dplyr::lag(dh$Cred_Real_Total, 2)
  dh$fci_lag1 <- dplyr::lag(dh$FCI_OS_EXP, 1)
  ex <- c("y_lag1", "y_lag2", "fci_lag1", EXOG_MACRO, intersect(ER_CTRLS, names(dh)))
  rd <- na.omit(dh[, c("y_fwd", ex, "DXY")])
  if (nrow(rd) < 50) return(NULL)
  rf <- lm(as.formula(paste("y_fwd ~ DXY +", paste(ex, collapse = "+"))), data = rd)
  ct <- lmtest::coeftest(rf, vcov = sandwich::NeweyWest(rf, lag = h + 1, prewhite = FALSE))
  data.table(sample = sample, horizon = h, n = nrow(rd),
             b_rf = round(ct["DXY", 1], 4), se_rf = round(ct["DXY", 2], 4),
             t_rf_nw = round(ct["DXY", 3], 2), p_rf_nw = round(ct["DXY", 4], 4))
}
rf_out <- list()
for (s in names(samples)) for (h in c(6, 12, 18, 24))
  rf_out[[paste(s, h)]] <- rf_diag(samples[[s]], h, s)
rf_out <- rbindlist(rf_out)
write_rev_csv(iv_out, "Rev_Aligned_IV_AR.csv")
write_rev_csv(rf_out, "Rev_Aligned_ReducedForm.csv")
print(iv_out[, .(sample, horizon, coef, p_value, first_stage_F, eff_F, ar_lo, ar_hi, n_obs)])
print(rf_out)

# ---------------------------------------------------------------------------
# Part 5 — Conley/CHR breakdown on the aligned index
# gamma_cal (trade channel a*b) is index-invariant: a = DXY->ToT leg,
# b = ToT->credit leg; read from script 49's calibration output.
# ---------------------------------------------------------------------------
cat("\n[5] Conley/CHR breakdown gamma* (aligned index)...\n")
cal <- fread(file.path(rev_paths$out_csv, "Rev_CHR_Calibration_Interval.csv"))
chr_out <- list()
for (h in c(6, 12, 18)) {
  base <- iv_out[sample == "full" & horizon == h]
  if (!nrow(base)) next
  # breakdown: smallest |gamma| at which the 90% NW CI for beta includes 0
  gseq <- seq(0, -1.5, by = -0.005)   # trade channel gamma_cal is small/positive;
                                       # breakdown searched on the adverse side
  gstar <- NA_real_
  for (g in gseq) {
    r <- iv_lp_h(samples$full, "Cred_Real_Total", "FCI_OS_EXP", "DXY", h,
                 gamma = g, ar_grid = c(0))   # AR grid unused here
    if (is.null(r)) break
    if (r$ci_lower <= 0 && r$ci_upper >= 0) { gstar <- g; break }
  }
  gcal <- cal[horizon == h, gamma_cal]
  chr_out[[as.character(h)]] <- data.table(
    horizon = h, beta_gamma0 = base$coef,
    ci_lo_gamma0 = base$ci_lower, ci_hi_gamma0 = base$ci_upper,
    gamma_breakdown = gstar,
    gamma_calibrated_trade = if (length(gcal)) gcal else NA_real_,
    ratio_cal_to_breakdown = if (length(gcal) && !is.na(gstar) && gstar != 0)
      abs(gcal / gstar) else NA_real_)
}
chr_out <- rbindlist(chr_out)
write_rev_csv(chr_out, "Rev_Aligned_CHR_Breakdown.csv")
print(chr_out)

# ---------------------------------------------------------------------------
# Part 6 — credit LP on the aligned index (and rolling variant)
# ---------------------------------------------------------------------------
cat("\n[6] Credit LP...\n")
lp_out <- list()
for (s in names(samples)) for (v in c("FCI_OS_EXP", "FCI_OS_ROLL")) {
  r <- run_lp(samples[[s]], "Cred_Real_Total", v, max_h = 18,
              controls = c("IMAEP_yoy", "IPC_yoy"))
  if (!nrow(r)) next
  cc <- na.omit(samples[[s]][, c("Cred_Real_Total", v), with = FALSE])
  sdx <- sd(cc[[v]])
  lp_out[[paste(s, v)]] <- data.table(sample = s, index = v, r,
                                      sd_x = sdx, coef_per_sd = r$coef * sdx)
}
lp_out <- rbindlist(lp_out)
write_rev_csv(lp_out, "Rev_Aligned_Credit_LP.csv")
print(lp_out[horizon %in% c(6, 12, 18),
             .(sample, index, horizon, coef, p_value, coef_per_sd, n_obs)])

# ---------------------------------------------------------------------------
# Part 7 — S1: null-imposed moving-block bootstrap for the DXY reduced form
# Persistence-preserving and JOINT (round-12 fix): the monthly DXY increment
# and the null-model regression disturbance dated the same month are
# resampled as a PAIR within each moving block, so any contemporaneous
# dependence between dollar innovations and credit disturbances is
# preserved under the null.  DXY* is rebuilt as a random walk from the
# resampled increments (preserves near-unit-root persistence), anchored at
# the sample's initial observed DXY level; y* = null-model fitted values +
# the paired resampled disturbances.  Blocks are overlapping with circular
# wrap-around, L = 18 months >= the maximum LP horizon (so the
# overlapping-horizon serial correlation of the h-step disturbance and the
# dollar's local persistence are both preserved within blocks).  H0: no DXY
# effect, imposed by construction (y* built from the null model that
# excludes DXY).  Two-sided p on the NW(h+1) t of DXY* in the re-estimated
# reduced form; controls held at observed values.
# ---------------------------------------------------------------------------
cat("\n[7] Null-imposed JOINT moving-block bootstrap (reduced form)...\n")
mbb_idx <- function(n, L) {
  k <- ceiling(n / L)
  starts <- sample.int(n, k, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) ((s:(s + L - 1) - 1) %% n) + 1))
  idx[1:n]
}
boot_rf <- function(data, h, sample, B = 999, L = 18) {
  dh <- as.data.frame(data)
  dh$y_fwd  <- dplyr::lead(dh$Cred_Real_Total, h)
  dh$y_lag1 <- dplyr::lag(dh$Cred_Real_Total, 1)
  dh$y_lag2 <- dplyr::lag(dh$Cred_Real_Total, 2)
  dh$fci_lag1 <- dplyr::lag(dh$FCI_OS_EXP, 1)
  ex <- c("y_lag1", "y_lag2", "fci_lag1", EXOG_MACRO, intersect(ER_CTRLS, names(dh)))
  rd <- na.omit(dh[, c("y_fwd", ex, "DXY")])
  n <- nrow(rd); if (n < 60) return(NULL)
  xs <- paste(ex, collapse = " + ")
  rf  <- lm(as.formula(paste("y_fwd ~ DXY +", xs)), data = rd)
  t_obs <- lmtest::coeftest(rf, vcov = sandwich::NeweyWest(rf, lag = h + 1,
                                                           prewhite = FALSE))["DXY", 3]
  m0 <- lm(as.formula(paste("y_fwd ~", xs)), data = rd)
  fit0 <- fitted(m0); e0 <- resid(m0)
  ddxy <- diff(rd$DXY)
  t_star <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- mbb_idx(n - 1, L)   # ONE index set: (dDXY_t, e_t) resampled jointly
    dxy_star <- rd$DXY[1] + c(0, cumsum(ddxy[idx]))
    y_star <- fit0 + c(e0[1], e0[idx + 1])
    rb <- rd; rb$DXY <- dxy_star; rb$y_fwd <- y_star
    mb <- lm(as.formula(paste("y_fwd ~ DXY +", xs)), data = rb)
    t_star[b] <- tryCatch(
      lmtest::coeftest(mb, vcov = sandwich::NeweyWest(mb, lag = h + 1,
                                                      prewhite = FALSE))["DXY", 3],
      error = function(e) NA_real_)
  }
  t_star <- t_star[!is.na(t_star)]
  data.table(sample = sample, horizon = h, n = n, B_effective = length(t_star),
             block_L = L, t_obs_nw = round(t_obs, 2),
             p_boot = round((1 + sum(abs(t_star) >= abs(t_obs))) /
                            (length(t_star) + 1), 4))
}
bb <- list()
for (h in c(6, 12, 18)) bb[[paste("full", h)]] <- boot_rf(samples$full, h, "full")
bb[["postIT 12"]] <- boot_rf(samples$postIT, 12, "postIT")
bb <- rbindlist(bb)
write_rev_csv(bb, "Rev_Aligned_RF_Bootstrap.csv")
print(bb)

# ---------------------------------------------------------------------------
# Part 8 — verdict (decision point for the framing rewrite)
# ---------------------------------------------------------------------------
cat("\n[8] Verdict...\n")
v_bounds  <- ab[y == "FCI_OS_EXP" & spec == "whitened (min lags passing BG12)"]
v_fs      <- fs_out[index == "FCI_OS_EXP" & horizon == 12]
v_iv      <- iv_out[horizon == 12]
v_rf      <- rf_out[horizon == 12]
v_bb      <- bb[horizon == 12]
v_lp      <- lp_out[index == "FCI_OS_EXP" & horizon == 12]
verdict <- data.table(
  leg = c("ARDL bounds full (whitened, exact p)",
          "ARDL bounds postIT (whitened, exact p)",
          "First stage full (classical F | HAC t)",
          "First stage postIT (classical F | HAC t)",
          "Reduced form h=12 full (NW p | bootstrap p)",
          "Reduced form h=12 postIT (NW p | bootstrap p)",
          "AR 90% set h=12 full",
          "Credit LP h=12 full (per-SD, p)",
          "Credit LP h=12 postIT (per-SD, p)"),
  value = c(
    sprintf("F=%.2f p=%.3f", v_bounds[sample == "full", F_stat],
            v_bounds[sample == "full", F_p_exact]),
    sprintf("F=%.2f p=%.3f", v_bounds[sample == "postIT", F_stat],
            v_bounds[sample == "postIT", F_p_exact]),
    sprintf("F=%.1f | t=%.2f", v_fs[sample == "full", F_classical],
            v_fs[sample == "full", t_hac_nw]),
    sprintf("F=%.1f | t=%.2f", v_fs[sample == "postIT", F_classical],
            v_fs[sample == "postIT", t_hac_nw]),
    sprintf("p=%.4f | p_boot=%.4f", v_rf[sample == "full", p_rf_nw],
            v_bb[sample == "full", p_boot]),
    sprintf("p=%.4f | p_boot=%.4f", v_rf[sample == "postIT", p_rf_nw],
            v_bb[sample == "postIT", p_boot]),
    sprintf("[%.0f, %.0f]", v_iv[sample == "full", ar_lo],
            v_iv[sample == "full", ar_hi]),
    sprintf("%.2f (p=%.4f)", v_lp[sample == "full", coef_per_sd],
            v_lp[sample == "full", p_value]),
    sprintf("%.2f (p=%.4f)", v_lp[sample == "postIT", coef_per_sd],
            v_lp[sample == "postIT", p_value])))
write_rev_csv(verdict, "Rev_Aligned_Verdict.csv")
print(verdict, row.names = FALSE)

cat(sprintf("\n=== 56 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
