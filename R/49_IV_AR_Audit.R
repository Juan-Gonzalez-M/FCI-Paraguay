# ============================================================================
# 49_IV_AR_Audit.R
# Full-sample IV/AR scaling audit (main-text Table 7 / Appendix K.1).
#
# Replicates the DXY-only IV-LP of script 23 (same data, controls, sample)
# for total real credit at h = 6, 12, 18 and decomposes the 2SLS estimate:
#
#   (1) identity check: 2SLS = reduced form / first stage (just-identified),
#       on the identical estimation sample;
#   (2) the conditional first-stage coefficient (DXY -> FCI within the full
#       control set), which is what scales per-FCI-unit IV magnitudes;
#   (3) the reduced form re-expressed per 1-SD DXY movement -- the
#       economically interpretable unit -- compared with the OLS-implied
#       effect of the same DXY movement (OLS coef x FS coef x SD(DXY));
#   (4) the AR p-value at the OLS point estimate, documenting that the AR
#       set excludes OLS because the reduced form is too large relative to
#       (FS x OLS), i.e., the LATE/exclusion tension, not a coding artifact.
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(sandwich); library(lmtest) })

d <- read.csv("../output/csv/New_External_Variables.csv")
d$fecha <- as.Date(d$fecha)
d <- d %>% arrange(fecha) %>%
  mutate(IMAEP_yoy_L1 = lag(IMAEP_yoy, 1), IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
         IPC_yoy_L1   = lag(IPC_yoy, 1),   IPC_yoy_L2   = lag(IPC_yoy, 2))

er  <- intersect(c("d_ToT", "US_10Y", "SP500"), names(d))
ex0 <- c("IMAEP_yoy","IMAEP_yoy_L1","IMAEP_yoy_L2","IPC_yoy","IPC_yoy_L1","IPC_yoy_L2", er)

cat("=== 49: full-sample IV/AR audit (Table 7 spec, total real credit) ===\n\n")
out <- data.frame()
for (h in c(6, 12, 18)) {
  dh <- d %>% mutate(y_fwd = lead(Cred_Real_Total, h),
                     y_lag1 = lag(Cred_Real_Total, 1),
                     y_lag2 = lag(Cred_Real_Total, 2),
                     fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1))
  ex <- c("y_lag1","y_lag2","fci_lag1", ex0)
  rd <- dh %>% dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG, all_of(ex), DXY) %>% na.omit()
  xs <- paste(ex, collapse = " + ")

  ols <- lm(as.formula(paste("y_fwd ~ FCI_ENDO_exCredit_AVG +", xs)), data = rd)
  rf  <- lm(as.formula(paste("y_fwd ~ DXY +", xs)), data = rd)
  fs  <- lm(as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY +", xs)), data = rd)

  b_ols <- coef(ols)["FCI_ENDO_exCredit_AVG"]
  b_rf  <- coef(rf)["DXY"]
  b_fs  <- coef(fs)["DXY"]
  b_iv_manual <- b_rf / b_fs

  iv <- ivreg::ivreg(as.formula(paste("y_fwd ~ FCI_ENDO_exCredit_AVG +", xs,
                                      "| DXY +", xs)), data = rd)
  b_iv <- coef(iv)["FCI_ENDO_exCredit_AVG"]

  # conditional first-stage F (classical, matching Table 7's Cond. F column)
  fs_r  <- lm(as.formula(paste("FCI_ENDO_exCredit_AVG ~", xs)), data = rd)
  f_cls <- anova(fs_r, fs)$F[2]

  # AR p-value at the OLS point estimate
  rd$y_t <- rd$y_fwd - as.numeric(b_ols) * rd$FCI_ENDO_exCredit_AVG
  ar_r <- lm(as.formula(paste("y_t ~", xs)), data = rd)
  ar_u <- lm(as.formula(paste("y_t ~ DXY +", xs)), data = rd)
  p_ar_at_ols <- anova(ar_r, ar_u)$`Pr(>F)`[2]

  sd_dxy <- sd(rd$DXY)
  rf_per_sd  <- b_rf * sd_dxy                       # credit pp per 1-SD DXY (IV view)
  ols_per_sd <- b_ols * b_fs * sd_dxy               # credit pp per 1-SD DXY (OLS view)

  cat(sprintf("h=%2d  N=%d  SD(DXY|sample)=%.2f\n", h, nrow(rd), sd_dxy))
  cat(sprintf("  OLS(FCI)          = %8.2f\n", b_ols))
  cat(sprintf("  2SLS (ivreg)      = %8.2f   manual RF/FS = %8.2f   (match: %s)\n",
              b_iv, b_iv_manual, ifelse(abs(b_iv - b_iv_manual) < 1e-6, "YES", "NO")))
  cat(sprintf("  Reduced form (DXY)= %8.4f pp/pt   First stage (DXY) = %+7.4f SD/pt   (cond. F = %5.1f)\n",
              b_rf, b_fs, f_cls))
  cat(sprintf("  Per 1-SD DXY: reduced form = %6.2f pp   OLS-implied = %6.2f pp   (wedge = %.1fx)\n",
              rf_per_sd, ols_per_sd, rf_per_sd / ols_per_sd))
  cat(sprintf("  AR p-value at OLS point (%.2f): %.4f\n\n", b_ols, p_ar_at_ols))

  out <- rbind(out, data.frame(horizon = h, n = nrow(rd), sd_dxy = sd_dxy,
    b_ols = b_ols, b_2sls = b_iv, b_2sls_manual = b_iv_manual,
    b_rf = b_rf, b_fs = b_fs, F_fs_classical = f_cls,
    rf_per_sd_dxy = rf_per_sd, ols_implied_per_sd_dxy = ols_per_sd,
    ar_p_at_ols = p_ar_at_ols))
}
write.csv(out, "../output/revision/csv/Rev_IV_AR_Audit.csv", row.names = FALSE)
cat("Saved: Rev_IV_AR_Audit.csv\n")

# ---------------------------------------------------------------------------
# Part 2:
#   (a) HAC-robust AR: the script-23 AR uses a classical anova F; recompute
#       the AR test and set with Newey-West(h+1) Wald inference -- an
#       independent implementation of the inversion;
#   (b) extended-grid check: is the AR set bounded below, or is the -250
#       edge a grid artifact (true set unbounded)?
#   (c) scaling: SD of the residualized DXY (projected on controls), and
#       per-10-point conditional reduced-form effects;
#   (d) persistence diagnostics: ADF/PP/KPSS for DXY level, 12m DXY change,
#       FCI_ENDO_exCredit; ADF on first-stage and reduced-form residuals;
#   (e) stationary instrument: 12-month DXY change, same spec -- conditional
#       F, reduced form, 2SLS, classical and HAC AR at h = 6, 12, 18.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(tseries))
cat("\n--- Part 2: HAC AR, boundedness, scaling, persistence, stationary IV ---\n\n")

d$DXY_d12 <- c(rep(NA, 12), diff(d$DXY, 12))

ar_p_hac <- function(rd, xs, d0, h) {
  rd$y_t <- rd$y_fwd - d0 * rd$FCI_ENDO_exCredit_AVG
  m <- lm(as.formula(paste("y_t ~ DXY +", xs)), data = rd)
  ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE))
  ct["DXY", 4]
}

aud2 <- data.frame()
for (h in c(6, 12, 18, 24)) {
  dh <- d %>% mutate(y_fwd = lead(Cred_Real_Total, h),
                     y_lag1 = lag(Cred_Real_Total, 1),
                     y_lag2 = lag(Cred_Real_Total, 2),
                     fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1))
  ex <- c("y_lag1","y_lag2","fci_lag1", ex0)
  rd <- dh %>% dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG, all_of(ex), DXY, DXY_d12) %>% na.omit()
  xs <- paste(ex, collapse = " + ")

  rf <- lm(as.formula(paste("y_fwd ~ DXY +", xs)), data = rd)
  ct_rf <- lmtest::coeftest(rf, vcov = sandwich::NeweyWest(rf, lag = h + 1, prewhite = FALSE))
  b_rf <- ct_rf["DXY", 1]; t_rf_nw <- ct_rf["DXY", 3]

  # residualized DXY
  dres <- resid(lm(as.formula(paste("DXY ~", xs)), data = rd))
  sd_res <- sd(dres)

  # HAC AR set on wide grid (step 1), plus deep-negative probes
  grid <- seq(-600, 100, by = 1)
  p_h <- sapply(grid, function(d0) ar_p_hac(rd, xs, d0, h))
  acc <- grid[p_h > 0.10]
  lo <- if (length(acc)) min(acc) else NA; up <- if (length(acc)) max(acc) else NA
  probe <- sapply(c(-1500, -5000), function(d0) ar_p_hac(rd, xs, d0, h))
  unbounded <- all(probe > 0.10) && !is.na(lo) && lo <= -599
  p_ar0_hac <- ar_p_hac(rd, xs, 0, h)

  cat(sprintf("h=%2d | RF=%.3f (NW t=%.2f) | SD(DXY_resid)=%.2f | RF x 10pt = %.2f pp | RF x 1SD-resid = %.2f pp\n",
              h, b_rf, t_rf_nw, sd_res, b_rf*10, b_rf*sd_res))
  cat(sprintf("      HAC-AR p(beta=0)=%.4f | HAC-AR 90%% set: [%s, %s]%s | deep probes p(-1500,-5000)=%.2f, %.2f\n",
              p_ar0_hac, ifelse(unbounded,"-Inf",as.character(lo)), as.character(up),
              ifelse(unbounded," (unbounded below)",""), probe[1], probe[2]))

  aud2 <- rbind(aud2, data.frame(horizon=h, rf=b_rf, rf_t_nw=t_rf_nw, sd_dxy_resid=sd_res,
    rf_per_10pt=b_rf*10, rf_per_1sd_resid=b_rf*sd_res, ar0_p_hac=p_ar0_hac,
    ar_hac_lo=ifelse(unbounded, -Inf, lo), ar_hac_up=up, unbounded_below=unbounded))
}
write.csv(aud2, "../output/revision/csv/Rev_IV_AR_Audit_HAC.csv", row.names = FALSE)

# (d) persistence diagnostics
cat("\nPersistence diagnostics (p-values):\n")
pers <- data.frame()
for (v in c("DXY", "DXY_d12", "FCI_ENDO_exCredit_AVG")) {
  x <- na.omit(d[[v]])
  suppressWarnings({
    p_adf <- adf.test(x)$p.value; p_pp <- pp.test(x)$p.value; p_kpss <- kpss.test(x)$p.value })
  cat(sprintf("  %-24s ADF p=%.3f  PP p=%.3f  KPSS p=%.3f\n", v, p_adf, p_pp, p_kpss))
  pers <- rbind(pers, data.frame(series=v, adf_p=p_adf, pp_p=p_pp, kpss_p=p_kpss))
}
# first-stage / reduced-form residual stationarity at h=12
dh <- d %>% mutate(y_fwd = lead(Cred_Real_Total, 12), y_lag1 = lag(Cred_Real_Total, 1),
                   y_lag2 = lag(Cred_Real_Total, 2), fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1))
ex <- c("y_lag1","y_lag2","fci_lag1", ex0)
rd <- dh %>% dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG, all_of(ex), DXY, DXY_d12) %>% na.omit()
xs <- paste(ex, collapse = " + ")
suppressWarnings({
  p_fs_res <- adf.test(resid(lm(as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY +", xs)), data = rd)))$p.value
  p_rf_res <- adf.test(resid(lm(as.formula(paste("y_fwd ~ DXY +", xs)), data = rd)))$p.value })
cat(sprintf("  first-stage residuals (h=12) ADF p=%.3f | reduced-form residuals ADF p=%.3f\n", p_fs_res, p_rf_res))
pers <- rbind(pers, data.frame(series=c("FS_residuals_h12","RF_residuals_h12"),
                               adf_p=c(p_fs_res,p_rf_res), pp_p=NA, kpss_p=NA))
write.csv(pers, "../output/revision/csv/Rev_IV_Persistence_Diagnostics.csv", row.names = FALSE)

# (e) stationary instrument: 12-month DXY change
cat("\nStationary instrument (12-month DXY change):\n")
st <- data.frame()
for (h in c(6, 12, 18)) {
  dh <- d %>% mutate(y_fwd = lead(Cred_Real_Total, h), y_lag1 = lag(Cred_Real_Total, 1),
                     y_lag2 = lag(Cred_Real_Total, 2), fci_lag1 = lag(FCI_ENDO_exCredit_AVG, 1))
  rd <- dh %>% dplyr::select(y_fwd, FCI_ENDO_exCredit_AVG, all_of(ex), DXY_d12) %>% na.omit()
  xs <- paste(ex, collapse = " + ")
  fs <- lm(as.formula(paste("FCI_ENDO_exCredit_AVG ~ DXY_d12 +", xs)), data = rd)
  ct_fs <- lmtest::coeftest(fs, vcov = sandwich::NeweyWest(fs, lag = h + 1, prewhite = FALSE))
  effF <- ct_fs["DXY_d12", 3]^2
  rf <- lm(as.formula(paste("y_fwd ~ DXY_d12 +", xs)), data = rd)
  ct_rf <- lmtest::coeftest(rf, vcov = sandwich::NeweyWest(rf, lag = h + 1, prewhite = FALSE))
  iv <- ivreg::ivreg(as.formula(paste("y_fwd ~ FCI_ENDO_exCredit_AVG +", xs, "| DXY_d12 +", xs)), data = rd)
  # HAC AR at 0 and HAC AR set
  arp <- function(d0) { rd$y_t <- rd$y_fwd - d0*rd$FCI_ENDO_exCredit_AVG
    m <- lm(as.formula(paste("y_t ~ DXY_d12 +", xs)), data = rd)
    lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE))["DXY_d12", 4] }
  grid <- seq(-600, 100, by = 1); pg <- sapply(grid, arp); acc <- grid[pg > 0.10]
  lo <- if (length(acc)) min(acc) else NA; up <- if (length(acc)) max(acc) else NA
  probe <- arp(-5000); unb <- !is.na(lo) && lo <= -599 && probe > 0.10
  cat(sprintf("  h=%2d: eff F=%5.1f | RF=%.3f (NW t=%.2f) | 2SLS=%.1f | AR p(0)=%.4f | AR set [%s, %s]%s\n",
              h, effF, ct_rf["DXY_d12",1], ct_rf["DXY_d12",3], coef(iv)["FCI_ENDO_exCredit_AVG"],
              arp(0), ifelse(unb,"-Inf",as.character(lo)), as.character(up), ifelse(unb," (unbounded)","")))
  st <- rbind(st, data.frame(horizon=h, effF_nw=effF, rf=ct_rf["DXY_d12",1], rf_t=ct_rf["DXY_d12",3],
    b2sls=coef(iv)["FCI_ENDO_exCredit_AVG"], ar0_p=arp(0),
    ar_lo=ifelse(unb,-Inf,lo), ar_up=up))
}
write.csv(st, "../output/revision/csv/Rev_IV_StationaryInstrument.csv", row.names = FALSE)
cat("\nSaved: Rev_IV_AR_Audit_HAC.csv, Rev_IV_Persistence_Diagnostics.csv, Rev_IV_StationaryInstrument.csv\n")

# ---------------------------------------------------------------------------
# Part 3: N.2 completion -- Z-Score individual-method LP, same convention as
# the PCA/VAR/DFM rows (replication check on PCA first).
# ---------------------------------------------------------------------------
cat("\n--- Part 3: individual-method LP (Z-Score completion) ---\n")
rv <- read.csv("../output/csv/FCI_Robustness_Versions.csv")
rv$fecha <- as.Date(rv$fecha)
dm <- d %>% left_join(rv %>% dplyr::select(fecha, FCI_exCredit_ZSCORE_norm,
  FCI_exCredit_PCA_norm, FCI_exCredit_VAR_norm, FCI_exCredit_DFM_norm), by = "fecha")

lp_h12 <- function(data, xvar, lagctrl = TRUE) {
  dh <- data %>% mutate(y_fwd = lead(Cred_Real_Total, 12),
                        y_lag1 = lag(Cred_Real_Total, 1), y_lag2 = lag(Cred_Real_Total, 2),
                        x_lag1 = lag(.data[[xvar]], 1))
  ctr <- if (lagctrl) c("IMAEP_yoy","IMAEP_yoy_L1","IMAEP_yoy_L2","IPC_yoy","IPC_yoy_L1","IPC_yoy_L2")
         else c("IMAEP_yoy","IPC_yoy")
  rhs <- c(xvar, "y_lag1", "y_lag2", "x_lag1", ctr)
  rd <- dh %>% dplyr::select(y_fwd, all_of(rhs)) %>% na.omit()
  m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
  ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = 13, prewhite = FALSE))
  c(b = ct[xvar,1], se = ct[xvar,2], p = ct[xvar,4], n = nrow(rd))
}
for (lc in c(TRUE, FALSE)) {
  cat(sprintf("  convention lagctrl=%s:\n", lc))
  for (v in c("FCI_exCredit_PCA_norm","FCI_exCredit_VAR_norm","FCI_exCredit_DFM_norm","FCI_exCredit_ZSCORE_norm")) {
    r <- lp_h12(dm, v, lc)
    cat(sprintf("    %-26s b=%7.2f  se=%.2f  p=%.4f  n=%d\n", v, r["b"], r["se"], r["p"], r["n"]))
  }
}

# ---------------------------------------------------------------------------
# Part 4: delta-method 90% interval for the CHR calibrated direct effect
# gamma_cal(h) = (DXY -> ToT growth) x (ToT growth -> credit at h),
# replicating script 38's components with their NW standard errors.
# Independence across the two regressions is assumed (stated in the note).
# ---------------------------------------------------------------------------
cat("\n--- Part 4: CHR calibration interval (delta method) ---\n")
m1 <- lm(d_ToT ~ DXY + IMAEP_yoy + IPC_yoy, data = d)
ct1 <- lmtest::coeftest(m1, vcov = sandwich::NeweyWest(m1, lag = 13, prewhite = FALSE))
a <- ct1["DXY",1]; se_a <- ct1["DXY",2]
cat(sprintf("  DXY -> ToT pass-through: a=%.5f (se %.5f)\n", a, se_a))
brk <- c(`6`=-0.135, `12`=-0.535, `18`=-0.554)
chr <- data.frame()
for (h in c(6,12,18)) {
  dh <- d %>% mutate(y_fwd = lead(Cred_Real_Total, h),
                     y_lag1 = lag(Cred_Real_Total,1), y_lag2 = lag(Cred_Real_Total,2),
                     x_lag1 = lag(d_ToT,1),
                     IM1=lag(IMAEP_yoy,1), IM2=lag(IMAEP_yoy,2),
                     IP1=lag(IPC_yoy,1), IP2=lag(IPC_yoy,2))
  rhs <- c("d_ToT","y_lag1","y_lag2","x_lag1","IMAEP_yoy","IM1","IM2","IPC_yoy","IP1","IP2")
  rd <- dh %>% dplyr::select(y_fwd, all_of(rhs)) %>% na.omit()
  m2 <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse="+"))), data = rd)
  ct2 <- lmtest::coeftest(m2, vcov = sandwich::NeweyWest(m2, lag = h+1, prewhite = FALSE))
  b <- ct2["d_ToT",1]; se_b <- ct2["d_ToT",2]
  g <- a*b; se_g <- sqrt(b^2*se_a^2 + a^2*se_b^2)
  lo <- g - 1.645*se_g; hi <- g + 1.645*se_g
  bmax <- max(abs(lo), abs(hi))
  cat(sprintf("  h=%2d: ToT->credit b=%.3f (se %.3f) | gamma_cal=%.4f, 90%% CI [%.4f, %.4f] | max|CI|/|breakdown| = %.1f%%\n",
              h, b, se_b, g, lo, hi, 100*bmax/abs(brk[as.character(h)])))
  chr <- rbind(chr, data.frame(horizon=h, a=a, se_a=se_a, b=b, se_b=se_b,
    gamma_cal=g, se_gamma=se_g, ci_lo=lo, ci_hi=hi, breakdown=brk[as.character(h)],
    max_ci_pct_of_breakdown=100*bmax/abs(brk[as.character(h)])))
}
write.csv(chr, "../output/revision/csv/Rev_CHR_Calibration_Interval.csv", row.names = FALSE)
cat("Saved: Rev_CHR_Calibration_Interval.csv\n")

# ---------------------------------------------------------------------------
# Part 5: JOINT moving-block bootstrap for the
# CHR calibrated direct effect. The delta method in Part 4 assumes the two
# calibration components are independent; here both regressions are
# re-estimated on the SAME resampled months (block length 18), so the
# bootstrap distribution of the product a*b embeds their joint dependence.
# ---------------------------------------------------------------------------
cat("\n--- Part 5: CHR calibration interval (joint block bootstrap) ---\n")
set.seed(20260712)
B_CHR <- 999; L_CHR <- 18
chr_bb <- data.frame()
for (h in c(6,12,18)) {
  dh <- d %>% mutate(y_fwd = lead(Cred_Real_Total, h),
                     y_lag1 = lag(Cred_Real_Total,1), y_lag2 = lag(Cred_Real_Total,2),
                     x_lag1 = lag(d_ToT,1),
                     IM1=lag(IMAEP_yoy,1), IM2=lag(IMAEP_yoy,2),
                     IP1=lag(IPC_yoy,1), IP2=lag(IPC_yoy,2))
  rhs2 <- c("d_ToT","y_lag1","y_lag2","x_lag1","IMAEP_yoy","IM1","IM2","IPC_yoy","IP1","IP2")
  need <- unique(c("y_fwd", rhs2, "DXY", "IMAEP_yoy", "IPC_yoy"))
  rd <- dh %>% dplyr::select(all_of(need)) %>% na.omit()
  n  <- nrow(rd); k <- ceiling(n / L_CHR); smax <- n - L_CHR + 1
  # point estimates on the joint sample (same rows for both regressions)
  a0 <- coef(lm(d_ToT ~ DXY + IMAEP_yoy + IPC_yoy, data = rd))["DXY"]
  b0 <- coef(lm(as.formula(paste("y_fwd ~", paste(rhs2, collapse="+"))), data = rd))["d_ToT"]
  gb <- rep(NA_real_, B_CHR)
  for (bb in seq_len(B_CHR)) {
    idx <- unlist(lapply(sample.int(smax, k, replace = TRUE),
                         function(s) s:(s + L_CHR - 1)))[1:n]
    rb <- rd[idx, ]
    ab <- tryCatch(coef(lm(d_ToT ~ DXY + IMAEP_yoy + IPC_yoy, data = rb))["DXY"],
                   error = function(e) NA_real_)
    bbb <- tryCatch(coef(lm(as.formula(paste("y_fwd ~", paste(rhs2, collapse="+"))),
                            data = rb))["d_ToT"], error = function(e) NA_real_)
    gb[bb] <- ab * bbb
  }
  ci <- quantile(gb, c(0.05, 0.95), na.rm = TRUE)
  bmax <- max(abs(ci))
  cat(sprintf("  h=%2d: gamma_cal(joint sample)=%.4f | joint-bootstrap 90%% CI [%.4f, %.4f] | max|CI|/|breakdown| = %.1f%% | includes 0: %s\n",
              h, a0*b0, ci[1], ci[2], 100*bmax/abs(brk[as.character(h)]),
              ifelse(ci[1] <= 0 & ci[2] >= 0, "yes", "no")))
  chr_bb <- rbind(chr_bb, data.frame(horizon=h, gamma_cal=a0*b0,
    ci_lo=ci[1], ci_hi=ci[2], breakdown=brk[as.character(h)],
    max_ci_pct_of_breakdown=100*bmax/abs(brk[as.character(h)]),
    includes_zero = (ci[1] <= 0 & ci[2] >= 0), B=B_CHR, block_len=L_CHR))
}
write.csv(chr_bb, "../output/revision/csv/Rev_CHR_Calibration_JointBootstrap.csv", row.names = FALSE)
cat("Saved: Rev_CHR_Calibration_JointBootstrap.csv\n")
