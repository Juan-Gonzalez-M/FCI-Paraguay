# ============================================================================
# 50_Persistence_Valid_Aggregate.R — Persistence-valid estimation of the
# aggregate dollar-credit relationship (colleague follow-up, Jul 2026).
# Motivation: the aggregate IV uses DXY/FCI LEVELS that Appendix K.1 finds
# nonstationary (DXY ADF p=0.63; FCI ADF p=0.68); stationary transformations
# have strong first stages but no reduced form (K.1c), so the causal signal
# lives in a shared persistent component. This script re-estimates the levels
# relationships in frameworks that are valid under the observed integration
# properties:
#   Block 1  unit-root battery (orders of integration; K.1 sanity gate)
#   Block 2  Engle-Granger residual-based cointegration (MacKinnon 2010 CVs)
#   Block 3  Johansen rank tests + VECM (cointegrating vector signs, alpha
#            loadings, weak-exogeneity of DXY)
#   Block 4  ARDL/PSS bounds tests (valid under mixed I(0)/I(1)) with
#            long-run multipliers and the PSS t-test on error correction
#   Block 5  verdict table + summary figure
#   Block 6  ex-TCN and rates-only expanding indices: does the DXY long-run
#            relation survive once the exchange rate is removed from the FCI?
#   Block 7  ARDL specification & stability audit for the headline models
#            (lags, serial correlation, RESET, CUSUM/supF, I(0)/I(1) bounds
#            incl. exact finite-sample, lag & deterministic sensitivity)
#   Block 8  pooled UECM with post-IT interactions: formal test of
#            H0: theta_preIT = theta_postIT for the long-run DXY multiplier
# PARALLEL/EXPLORATORY: not registered in RUN_REVISION.R; results feed a
# private decision memo only. ADDITIVE: writes only to output/revision/
# (Rev_PV_* CSVs, 292_*.png). Nothing else is read-modified.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
suppressMessages({library(tseries); library(urca); library(vars); library(ARDL)})
set.seed(20260713)
cat("=== 50: Persistence-valid aggregate estimation ===\n")

POSTIT <- as.Date("2011-05-01")

# ---------------------------------------------------------------------------
# Block 0 — data
# ---------------------------------------------------------------------------
d <- load_ext_data()

exp_fci <- fread(file.path(rev_paths$out_csv, "Rev_ExpandingFCI_Series.csv"))
exp_fci[, fecha := as.Date(fecha)]
d <- merge(d, exp_fci[, .(fecha, FCI_EXP_exCredit, FCI_EXP_ENDO_exCredit)],
           by = "fecha", all.x = TRUE)

macro <- as.data.table(readxl::read_excel(micro_paths$fci_xlsx, sheet = "Datos_macro"))
macro[, fecha := as.Date(Fecha)]
macro <- macro[, .(fecha, Creditos_deflactados)]
d <- merge(d, macro, by = "fecha", all.x = TRUE)
d[, lcred := log(Creditos_deflactados)]
setorder(d, fecha)

cat(sprintf("  Merged sample: %s to %s (n=%d); lcred non-NA: %d; FCI_EXP non-NA: %d\n",
            min(d$fecha), max(d$fecha), nrow(d),
            sum(!is.na(d$lcred)), sum(!is.na(d$FCI_EXP_ENDO_exCredit))))

samples <- list(full = d, postIT = d[fecha >= POSTIT])

# ---------------------------------------------------------------------------
# Block 1 — integration properties
# ---------------------------------------------------------------------------
cat("\n[1] Unit-root battery...\n")
ur_row <- function(x, label, sample) {
  x <- na.omit(x)
  if (length(x) < 40) return(NULL)
  a <- suppressWarnings(adf.test(x))
  p <- suppressWarnings(pp.test(x))
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
  specs <- list(
    DXY_level            = ds$DXY,
    DXY_d1               = c(NA, diff(ds$DXY)),
    lcred_level          = ds$lcred,
    lcred_d1             = c(NA, diff(ds$lcred)),
    FCI_roll_ENDOexC     = ds$FCI_ENDO_exCredit_AVG,
    FCI_exp_ENDOexC      = ds$FCI_EXP_ENDO_exCredit,
    FCI_exp_ENDOexC_d1   = c(NA, diff(ds$FCI_EXP_ENDO_exCredit)),
    CredGrowth_yoy       = ds$Cred_Real_Total,
    ToT_level            = ds$Paraguay_ToT,
    ToT_d1               = c(NA, diff(ds$Paraguay_ToT)))
  for (nm in names(specs)) ur[[paste(s, nm)]] <- ur_row(specs[[nm]], nm, s)
}
ur <- rbindlist(ur)
write_rev_csv(ur, "Rev_PV_UnitRoots.csv")

# Sanity gate vs published K.1 values (full sample: DXY ADF ~0.63, FCI ~0.68)
g1 <- ur[sample == "full" & series == "DXY_level", adf_p]
g2 <- ur[sample == "full" & series == "FCI_roll_ENDOexC", adf_p]
cat(sprintf("  K.1 gate: DXY ADF p=%.2f (published 0.63), rolling FCI ADF p=%.2f (published 0.68) -> %s\n",
            g1, g2, ifelse(abs(g1 - 0.63) < 0.1 & abs(g2 - 0.68) < 0.1, "PASS", "CHECK")))

# ---------------------------------------------------------------------------
# Block 2 — Engle-Granger residual-based cointegration
# MacKinnon (2010, QED WP 1227) asymptotic critical values, constant, no
# trend, N = variables in the cointegrating relation (standard ADF critical
# values are invalid for residual-based tests):
#   N=2: 1% -3.90, 5% -3.34, 10% -3.04
#   N=3: 1% -4.29, 5% -3.74, 10% -3.45
# ---------------------------------------------------------------------------
cat("\n[2] Engle-Granger residual cointegration...\n")
EG_CV <- list(`2` = c(p01 = -3.90, p05 = -3.34, p10 = -3.04),
              `3` = c(p01 = -4.29, p05 = -3.74, p10 = -3.45))
eg_row <- function(ds, yv, xv, sample) {
  cc <- na.omit(ds[, c(yv, xv), with = FALSE])
  if (nrow(cc) < 60) return(NULL)
  f <- as.formula(paste(yv, "~", paste(xv, collapse = "+")))
  m <- lm(f, cc)
  res <- residuals(m)
  # EG step 2: ADF on residuals, no deterministic terms, BIC lags
  adf <- ur.df(res, type = "none", selectlags = "BIC")
  stat <- adf@teststat[1]
  cv <- EG_CV[[as.character(length(xv) + 1)]]
  sig <- ifelse(stat < cv["p01"], "1%", ifelse(stat < cv["p05"], "5%",
                ifelse(stat < cv["p10"], "10%", "no")))
  data.table(sample = sample, y = yv, x = paste(xv, collapse = "+"),
             n = nrow(cc), beta_x1 = round(coef(m)[2], 5),
             eg_adf_stat = round(stat, 3),
             cv_5pct = cv["p05"], cv_10pct = cv["p10"],
             cointegrated_at = sig)
}
eg <- list()
for (s in names(samples)) {
  ds <- samples[[s]]
  eg[[paste(s, 1)]] <- eg_row(ds, "FCI_EXP_ENDO_exCredit", "DXY", s)
  eg[[paste(s, 2)]] <- eg_row(ds, "FCI_ENDO_exCredit_AVG", "DXY", s)
  eg[[paste(s, 3)]] <- eg_row(ds, "lcred", "DXY", s)
  eg[[paste(s, 4)]] <- eg_row(ds, "lcred", c("FCI_EXP_ENDO_exCredit", "DXY"), s)
}
eg <- rbindlist(eg)
write_rev_csv(eg, "Rev_PV_EngleGranger.csv")
print(eg[, .(sample, y, x, eg_adf_stat, cointegrated_at)])

# ---------------------------------------------------------------------------
# Block 3 — Johansen rank tests + VECM
# ---------------------------------------------------------------------------
cat("\n[3] Johansen + VECM...\n")
jo_block <- function(ds, vars_in, sys_label, sample, ecdet = "const") {
  cc <- na.omit(ds[, c("fecha", vars_in), with = FALSE])
  if (nrow(cc) < 80) return(NULL)
  X <- as.matrix(cc[, vars_in, with = FALSE])
  K <- tryCatch({
    vs <- vars::VARselect(X, lag.max = 8, type = "const")
    max(2, min(vs$selection["SC(n)"], vs$selection["HQ(n)"]))
  }, error = function(e) 2)
  cjo <- tryCatch(urca::ca.jo(X, type = "trace", ecdet = ecdet, K = K),
                  error = function(e) NULL)
  if (is.null(cjo)) return(NULL)
  ts_ <- cjo@teststat; cv <- cjo@cval   # rows: r<=p-1 ... r=0 (bottom = r=0)
  p <- length(vars_in)
  rank_rows <- data.table(
    sample = sample, system = sys_label, ecdet = ecdet, K = K, n = nrow(cc),
    hypothesis = rownames(cv),
    trace_stat = round(as.numeric(ts_), 2),
    cv_10 = cv[, "10pct"], cv_5 = cv[, "5pct"], cv_1 = cv[, "1pct"],
    reject_5 = as.numeric(ts_) > cv[, "5pct"],
    reject_10 = as.numeric(ts_) > cv[, "10pct"])
  # rank = number of rejected hypotheses starting from r=0 (last row upward)
  ord <- rev(seq_len(nrow(cv)))   # r=0 first
  r5 <- 0; for (i in ord) { if (as.numeric(ts_)[i] > cv[i, "5pct"]) r5 <- r5 + 1 else break }
  r10 <- 0; for (i in ord) { if (as.numeric(ts_)[i] > cv[i, "10pct"]) r10 <- r10 + 1 else break }
  vecm_rows <- NULL; wexo_rows <- NULL
  r_use <- max(r5, min(r10, 1))
  if (r_use >= 1 && r_use < p) {
    vecm <- urca::cajorls(cjo, r = r_use)
    beta <- vecm$beta            # normalized on first variable
    sm <- summary(vecm$rlm)
    alpha <- sapply(seq_len(p), function(i) coef(sm[[i]])["ect1", c("Estimate", "t value")])
    vecm_rows <- data.table(
      sample = sample, system = sys_label, ecdet = ecdet, rank_5 = r5, rank_10 = r10,
      r_used = r_use,
      beta_terms = paste(rownames(beta), collapse = "|"),
      beta_values = paste(round(beta[, 1], 4), collapse = "|"),
      alpha_eq = paste(vars_in, collapse = "|"),
      alpha_values = paste(round(alpha["Estimate", ], 4), collapse = "|"),
      alpha_tstats = paste(round(alpha["t value", ], 2), collapse = "|"))
    # weak exogeneity of DXY (last variable by construction below)
    idx_dxy <- which(vars_in == "DXY")
    A <- diag(p)[, -idx_dxy, drop = FALSE]
    wt <- tryCatch(urca::alrtest(z = cjo, A = A, r = r_use), error = function(e) NULL)
    if (!is.null(wt)) {
      lr <- wt@teststat; pv <- wt@pval[1]
      wexo_rows <- data.table(sample = sample, system = sys_label, ecdet = ecdet,
                              test = "DXY weakly exogenous (alpha_DXY = 0)",
                              LR_stat = round(as.numeric(lr), 3), p_value = round(pv, 3))
    }
  }
  list(rank = rank_rows, vecm = vecm_rows, wexo = wexo_rows)
}

systems <- list(
  `lcred_FCIexp_DXY` = c("lcred", "FCI_EXP_ENDO_exCredit", "DXY"),
  `lcred_DXY`        = c("lcred", "DXY"),
  `FCIexp_DXY`       = c("FCI_EXP_ENDO_exCredit", "DXY"),
  `FCIroll_DXY`      = c("FCI_ENDO_exCredit_AVG", "DXY"),
  `lcred_FCIroll_DXY`= c("lcred", "FCI_ENDO_exCredit_AVG", "DXY"))

jo_rank <- list(); jo_vecm <- list(); jo_wexo <- list()
for (s in names(samples)) for (sy in names(systems)) {
  for (ec in if (grepl("lcred", sy)) c("const", "trend") else "const") {
    res <- jo_block(samples[[s]], systems[[sy]], sy, s, ecdet = ec)
    if (is.null(res)) next
    key <- paste(s, sy, ec)
    jo_rank[[key]] <- res$rank; jo_vecm[[key]] <- res$vecm; jo_wexo[[key]] <- res$wexo
  }
}
jo_rank <- rbindlist(jo_rank); jo_vecm <- rbindlist(jo_vecm, fill = TRUE)
jo_wexo <- rbindlist(jo_wexo, fill = TRUE)
write_rev_csv(jo_rank, "Rev_PV_Johansen.csv")
write_rev_csv(jo_vecm, "Rev_PV_VECM.csv")
if (nrow(jo_wexo)) write_rev_csv(jo_wexo, "Rev_PV_WeakExogeneity.csv")
if (nrow(jo_vecm)) print(jo_vecm[, .(sample, system, ecdet, rank_5, rank_10, beta_values, alpha_tstats)])

# ---------------------------------------------------------------------------
# Block 4 — ARDL / Pesaran-Shin-Smith bounds
# Case 3 (unrestricted constant, no trend) throughout: the ARDL formula here
# carries no deterministic trend term, so case 5 is not defined for these
# models.  For the trending lcred this makes the bivariate DXY anchor a
# deliberately conservative specification; the PSS t-test on the lagged
# level guards against degenerate rejections driven by the trend (an F
# rejection WITHOUT a negative, significant lagged-level coefficient is the
# PSS degenerate case: no valid levels relationship, not evidence either way).
# ---------------------------------------------------------------------------
cat("\n[4] ARDL bounds tests...\n")
ardl_block <- function(ds, yv, xv, sample, case) {
  cc <- as.data.frame(na.omit(ds[, c(yv, xv), with = FALSE]))
  if (nrow(cc) < 80) return(NULL)
  f <- as.formula(paste(yv, "~", paste(xv, collapse = "+")))
  fit <- tryCatch(ARDL::auto_ardl(f, data = cc, max_order = 8, selection = "BIC"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  bm <- fit$best_model
  bF <- tryCatch(ARDL::bounds_f_test(bm, case = case), error = function(e) NULL)
  bT <- tryCatch(ARDL::bounds_t_test(bm, case = case), error = function(e) NULL)
  mult <- tryCatch(ARDL::multipliers(bm, type = "lr"), error = function(e) NULL)
  ue <- tryCatch(ARDL::uecm(bm), error = function(e) NULL)
  ect_coef <- NA_real_
  if (!is.null(ue)) {
    cn <- names(coef(ue)); lag_y <- grep(paste0("^L\\(", yv), cn, value = TRUE)
    if (length(lag_y)) ect_coef <- coef(ue)[lag_y[1]]
  }
  bounds <- data.table(
    sample = sample, y = yv, x = paste(xv, collapse = "+"), case = case,
    order = paste(fit$best_order, collapse = ","), n = nrow(cc),
    F_stat = if (!is.null(bF)) round(as.numeric(bF$statistic), 3) else NA,
    F_p = if (!is.null(bF)) round(as.numeric(bF$p.value), 4) else NA,
    t_stat = if (!is.null(bT)) round(as.numeric(bT$statistic), 3) else NA,
    t_p = if (!is.null(bT)) round(as.numeric(bT$p.value), 4) else NA,
    ect_coef = round(ect_coef, 4))
  lr <- NULL
  if (!is.null(mult)) {
    m <- as.data.frame(mult)   # cols: term, estimate, se, t, p (names vary by version)
    lr <- data.table(sample = sample, y = yv, x = paste(xv, collapse = "+"),
                     case = case, term = as.character(m[[1]]),
                     estimate = round(as.numeric(m[[2]]), 5),
                     se = round(as.numeric(m[[3]]), 5),
                     t = round(as.numeric(m[[4]]), 2),
                     p = round(as.numeric(m[[5]]), 4))
  }
  list(bounds = bounds, lr = lr)
}

ardl_specs <- list(
  list(y = "lcred", x = "DXY", cases = 3),
  list(y = "lcred", x = c("DXY", "Paraguay_ToT"), cases = 3),
  list(y = "lcred", x = "FCI_EXP_ENDO_exCredit", cases = 3),
  list(y = "FCI_EXP_ENDO_exCredit", x = "DXY", cases = 3),
  list(y = "FCI_ENDO_exCredit_AVG", x = "DXY", cases = 3))

ab <- list(); alr <- list()
for (s in names(samples)) for (sp in ardl_specs) for (cs in sp$cases) {
  res <- ardl_block(samples[[s]], sp$y, sp$x, s, cs)
  if (is.null(res)) next
  key <- paste(s, sp$y, paste(sp$x, collapse = "+"), cs)
  ab[[key]] <- res$bounds; alr[[key]] <- res$lr
}
ab <- rbindlist(ab); alr <- rbindlist(alr)
write_rev_csv(ab, "Rev_PV_ARDL_Bounds.csv")
write_rev_csv(alr, "Rev_PV_ARDL_LongRun.csv")
print(ab[, .(sample, y, x, case, F_stat, F_p, t_stat, t_p, ect_coef)])

# ---------------------------------------------------------------------------
# Block 5 — verdict + figure
# Expected supportive pattern: levels relationship exists (bounds/EG/Johansen);
# long-run DXY tightens FCI (+) and reduces credit (-); error correction
# negative and significant in the credit/FCI equation.
# ---------------------------------------------------------------------------
cat("\n[5] Verdict...\n")
vrows <- list()
for (i in seq_len(nrow(ab))) {
  r <- ab[i]
  lr_dxy <- alr[sample == r$sample & y == r$y & x == r$x & case == r$case & term == "DXY"]
  exp_sign <- if (r$y == "lcred") -1 else +1   # DXY should reduce credit, tighten FCI
  has_dxy <- nrow(lr_dxy) == 1
  f_rej   <- !is.na(r$F_p) & r$F_p < 0.10
  ect_ok  <- !is.na(r$t_p) & r$t_p < 0.10 & !is.na(r$ect_coef) & r$ect_coef < 0
  # PSS degenerate case: F rejects but the lagged-level coefficient is not
  # negative-significant -> no valid levels relationship (neither support
  # nor contradiction).
  degen   <- f_rej & !ect_ok
  vrows[[i]] <- data.table(
    design = "ARDL-PSS", sample = r$sample, y = r$y, x = r$x, case = r$case,
    levels_relationship = ifelse(f_rej & ect_ok,
                                 ifelse(r$F_p < 0.05, "yes (5%)", "yes (10%)"),
                          ifelse(degen, "degenerate (F rejects, no ECT)", "no")),
    lr_dxy_estimate = if (has_dxy) lr_dxy$estimate else NA,
    lr_dxy_p = if (has_dxy) lr_dxy$p else NA,
    lr_sign_as_expected = if (has_dxy) sign(lr_dxy$estimate) == exp_sign else NA,
    ect_neg_signif = ect_ok,
    supports_paper = {
      sgn <- if (has_dxy) (sign(lr_dxy$estimate) == exp_sign & lr_dxy$p < 0.10) else NA
      if (f_rej & ect_ok & isTRUE(sgn)) "YES"
      else if (degen) "NOT APPLICABLE (degenerate)"
      else if (f_rej & ect_ok) "PARTIAL"
      else "NO LEVELS RELATION"
    })
}
verdict <- rbindlist(vrows, fill = TRUE)
write_rev_csv(verdict, "Rev_PV_Verdict.csv")
print(verdict[, .(design, sample, y, x, case, levels_relationship,
                  lr_dxy_estimate, lr_dxy_p, supports_paper)])

# Figure: long-run DXY multipliers across ARDL specs
fig_dt <- alr[term == "DXY"]
if (nrow(fig_dt)) {
  fig_dt[, spec := paste0(y, " ~ ", x, " (case ", case, ", ", sample, ")")]
  fig_dt[, lo := estimate - Z90 * se]
  fig_dt[, hi := estimate + Z90 * se]
  fig_dt[, panel := ifelse(y == "lcred", "Long-run DXY effect on log real credit",
                           "Long-run DXY effect on FCI")]
  pF <- ggplot(fig_dt, aes(estimate, reorder(spec, estimate))) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, color = "steelblue4") +
    geom_point(size = 2.2, color = "steelblue4") +
    facet_wrap(~panel, scales = "free") +
    labs(title = "Persistence-valid long-run DXY effects (ARDL/PSS)",
         subtitle = "Long-run multipliers with 90% intervals; bounds-test details in Rev_PV_ARDL_Bounds.csv",
         x = "Long-run multiplier (per DXY point)", y = NULL) +
    theme_micro()
  save_rev_png(pF, "292_PV_LongRun_Summary.png", w = 12, h = 5.5)
}

cat("\nSUMMARY (ARDL designs):\n")
print(verdict[, .N, by = supports_paper])
cat("Johansen/EG results in Rev_PV_Johansen.csv / Rev_PV_EngleGranger.csv.\n")

# ---------------------------------------------------------------------------
# Block 6 — ex-TCN and rates-only expanding indices
# If the DXY long-run relation disappears once TCN is removed, the baseline
# cointegration primarily reflects a DXY-exchange-rate relation embedded in
# the FCI; if it survives, it demonstrates transmission into domestic
# financial prices beyond the exchange rate.  Indices are the transparent
# equal-weighted expanding z-scores, constructed exactly as in script 39.
# ---------------------------------------------------------------------------
cat("\n[6] Ex-TCN / rates-only expanding indices...\n")
inputs <- load_fci_inputs()
vars12 <- names(FCI_SIGNS)
Zexp <- sapply(vars12, function(v) expanding_z(inputs[[v]]) * FCI_SIGNS[v])
vars_extcn <- c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM",
                "Ratio_Cred_Depo", "Morosidad", "Rentabilidad", "Liquidez")
vars_rates <- c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM")
aux <- data.table(fecha = inputs$fecha,
                  FCI_EXP_ENDO_exC_exTCN = rowMeans(Zexp[, vars_extcn]),
                  FCI_EXP_RATES          = rowMeans(Zexp[, vars_rates]))
d <- merge(d, aux, by = "fecha", all.x = TRUE)
samples <- list(full = d, postIT = d[fecha >= POSTIT])

extcn_specs <- list(
  list(y = "FCI_EXP_ENDO_exC_exTCN", x = "DXY", cases = 3),
  list(y = "FCI_EXP_RATES",          x = "DXY", cases = 3))
xb <- list(); xlr <- list()
for (s in names(samples)) for (sp in extcn_specs) {
  res <- ardl_block(samples[[s]], sp$y, sp$x, s, 3)
  if (is.null(res)) next
  key <- paste(s, sp$y)
  xb[[key]] <- res$bounds; xlr[[key]] <- res$lr
}
xb <- rbindlist(xb); xlr <- rbindlist(xlr)
write_rev_csv(xb, "Rev_PV_ExTCN_ARDL_Bounds.csv")
write_rev_csv(xlr, "Rev_PV_ExTCN_ARDL_LongRun.csv")
cat("  Decision rule (colleague): survives ex-TCN -> more convincing;\n")
cat("  weakens-but-detectable -> TCN important-not-exclusive carrier;\n")
cat("  disappears -> exchange-rate-component evidence only.\n")
print(xb[, .(sample, y, F_stat, F_p, t_stat, t_p, ect_coef)])
print(xlr[term == "DXY", .(sample, y, estimate, se, t, p)])

# ---------------------------------------------------------------------------
# Block 7 — ARDL specification & stability audit (headline models)
# ---------------------------------------------------------------------------
cat("\n[7] ARDL specification & stability audit...\n")
audit_one <- function(ds, yv, xv, sample) {
  cc <- as.data.frame(na.omit(ds[, c(yv, xv), with = FALSE]))
  if (nrow(cc) < 80) return(NULL)
  f <- as.formula(paste(yv, "~", paste(xv, collapse = "+")))
  fit <- ARDL::auto_ardl(f, data = cc, max_order = 8, selection = "BIC")
  bm <- fit$best_model
  ue <- ARDL::uecm(bm)
  # plain design matrix for stability tests (avoid re-evaluating d()/L())
  X <- model.matrix(ue); y_ <- model.response(model.frame(ue))
  df <- data.frame(y_ = y_, X[, -1, drop = FALSE]); names(df) <- make.names(names(df))
  fs <- as.formula(paste("y_ ~", paste(setdiff(names(df), "y_"), collapse = "+")))
  # residual serial correlation and functional form
  bg1  <- lmtest::bgtest(ue, order = 1)$p.value
  bg12 <- lmtest::bgtest(ue, order = 12)$p.value
  rst  <- lmtest::resettest(ue, power = 2:3, type = "fitted")$p.value
  # parameter stability: Rec-CUSUM, OLS-CUSUM, supF (Andrews)
  cus  <- tryCatch(sctest(efp(fs, data = df, type = "Rec-CUSUM"))$p.value,
                   error = function(e) NA)
  ocus <- tryCatch(sctest(efp(fs, data = df, type = "OLS-CUSUM"))$p.value,
                   error = function(e) NA)
  supf <- tryCatch(sctest(Fstats(fs, data = df, from = 0.15))$p.value,
                   error = function(e) NA)
  # CUSUMSQ path (max deviation reported; formal tests above)
  csq <- tryCatch({
    w <- strucchange::recresid(fs, data = df)
    S <- cumsum(w^2) / sum(w^2)
    max(abs(S - seq_along(S) / length(S)))
  }, error = function(e) NA)
  # bounds with I(0)/I(1) critical values, asymptotic and exact finite-sample
  bFa <- ARDL::bounds_f_test(bm, case = 3, alpha = 0.05)
  bFe <- tryCatch(ARDL::bounds_f_test(bm, case = 3, alpha = 0.05, exact = TRUE),
                  error = function(e) NULL)
  bTa <- ARDL::bounds_t_test(bm, case = 3, alpha = 0.05)
  bTe <- tryCatch(ARDL::bounds_t_test(bm, case = 3, alpha = 0.05, exact = TRUE),
                  error = function(e) NULL)
  # sensitivity: second-best lag order and restricted-intercept case 2
  ord2 <- as.numeric(fit$top_orders[2, seq_len(1 + length(xv))])
  bm2 <- tryCatch(ARDL::ardl(f, data = cc, order = ord2), error = function(e) NULL)
  bF2 <- if (!is.null(bm2)) tryCatch(ARDL::bounds_f_test(bm2, case = 3),
                                     error = function(e) NULL) else NULL
  bFc2 <- tryCatch(ARDL::bounds_f_test(bm, case = 2), error = function(e) NULL)
  data.table(
    sample = sample, y = yv, x = paste(xv, collapse = "+"),
    max_order_considered = 8,
    order_selected = paste(fit$best_order, collapse = ","),
    n = nrow(cc),
    bg_p_order1 = round(bg1, 3), bg_p_order12 = round(bg12, 3),
    reset_p = round(rst, 3),
    rec_cusum_p = round(cus, 3), ols_cusum_p = round(ocus, 3),
    supF_p = round(supf, 3), cusumsq_maxdev = round(csq, 3),
    F_stat = round(as.numeric(bFa$statistic), 3),
    F_cv_I0_asy = round(bFa$tab[["Lower-bound I(0)"]], 3),
    F_cv_I1_asy = round(bFa$tab[["Upper-bound I(1)"]], 3),
    F_p_asy = round(as.numeric(bFa$p.value), 4),
    F_cv_I0_exact = if (!is.null(bFe)) round(bFe$tab[["Lower-bound I(0)"]], 3) else NA,
    F_cv_I1_exact = if (!is.null(bFe)) round(bFe$tab[["Upper-bound I(1)"]], 3) else NA,
    F_p_exact = if (!is.null(bFe)) round(as.numeric(bFe$p.value), 4) else NA,
    t_stat = round(as.numeric(bTa$statistic), 3),
    t_cv_I0_asy = round(bTa$tab[["Lower-bound I(0)"]], 3),
    t_cv_I1_asy = round(bTa$tab[["Upper-bound I(1)"]], 3),
    t_p_asy = round(as.numeric(bTa$p.value), 4),
    t_p_exact = if (!is.null(bTe)) round(as.numeric(bTe$p.value), 4) else NA,
    F_stat_order2 = if (!is.null(bF2)) round(as.numeric(bF2$statistic), 3) else NA,
    F_p_order2 = if (!is.null(bF2)) round(as.numeric(bF2$p.value), 4) else NA,
    order2 = paste(ord2, collapse = ","),
    F_stat_case2 = if (!is.null(bFc2)) round(as.numeric(bFc2$statistic), 3) else NA,
    F_p_case2 = if (!is.null(bFc2)) round(as.numeric(bFc2$p.value), 4) else NA)
}
headline <- list(c("FCI_EXP_ENDO_exCredit", "DXY"),
                 c("FCI_EXP_ENDO_exC_exTCN", "DXY"))
aud <- list()
for (s in names(samples)) for (hm in headline)
  aud[[paste(s, hm[1])]] <- audit_one(samples[[s]], hm[1], hm[2], s)
aud <- rbindlist(aud, fill = TRUE)
write_rev_csv(aud, "Rev_PV_ARDL_Audit.csv")
print(aud[, .(sample, y, order_selected, bg_p_order12, reset_p, rec_cusum_p,
              supF_p, F_stat, F_cv_I1_exact, F_p_exact, t_p_exact,
              F_p_order2, F_p_case2)])
cat("  Note: bootstrap bounds inference (bootCT/McNown et al.) not installed;\n")
cat("  exact finite-sample critical values reported instead.\n")

# --- 7b. Serial-correlation-robust (whitened) refit --------------------------
# The BIC-selected orders leave higher-order residual autocorrelation
# (BG(12) rejects), and the bounds F-test assumes white errors.  Refit with
# the minimal lag augmentation that passes BG(12) at 10% and re-report the
# bounds inference and long-run multiplier.  Also re-estimate the headline
# models on the pre-COVID post-IT sample (2011m5-2020m2), since supF flags
# instability within the post-IT window (COVID).
whiten_fit <- function(ds, yv, xv, sample) {
  cc <- as.data.frame(na.omit(ds[, c(yv, xv), with = FALSE]))
  if (nrow(cc) < 80) return(NULL)
  f <- as.formula(paste(yv, "~", paste(xv, collapse = "+")))
  grid <- expand.grid(p = 2:12, q = 0:4)
  grid <- grid[order(grid$p + grid$q), ]
  best <- NULL; best_bg <- -Inf
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
  bF <- ARDL::bounds_f_test(best, case = 3)
  bFe <- tryCatch(ARDL::bounds_f_test(best, case = 3, alpha = 0.05, exact = TRUE),
                  error = function(e) NULL)
  bT <- ARDL::bounds_t_test(best, case = 3)
  mult <- as.data.frame(ARDL::multipliers(best, type = "lr"))
  dxy <- mult[grepl("DXY", mult[[1]]), ]
  data.table(sample = sample, y = yv, order = paste(best_ord, collapse = ","),
             bg12_p = round(best_bg, 3),
             F_stat = round(as.numeric(bF$statistic), 3),
             F_p_asy = round(as.numeric(bF$p.value), 4),
             F_p_exact = if (!is.null(bFe)) round(as.numeric(bFe$p.value), 4) else NA,
             t_stat = round(as.numeric(bT$statistic), 3),
             t_p = round(as.numeric(bT$p.value), 4),
             lr_dxy = round(as.numeric(dxy[[2]]), 5),
             lr_dxy_p = round(as.numeric(dxy[[5]]), 4))
}
COVID_START <- as.Date("2020-03-01")
samples_w <- c(samples, list(postIT_preCOVID = d[fecha >= POSTIT & fecha < COVID_START]))
wh <- list()
for (s in names(samples_w)) for (hm in headline)
  wh[[paste(s, hm[1])]] <- whiten_fit(samples_w[[s]], hm[1], hm[2], s)
wh <- rbindlist(wh, fill = TRUE)
write_rev_csv(wh, "Rev_PV_ARDL_Whitened.csv")
print(wh)

# ---------------------------------------------------------------------------
# Block 8 — pooled UECM with post-IT break in the long-run multiplier
# Delta-method test of H0: theta_pre = theta_post, theta = -gamma/rho.
# ---------------------------------------------------------------------------
cat("\n[8] Pooled long-run break test (pre- vs post-IT)...\n")
break_one <- function(ds, yv, xv = "DXY", p = 2, q = 1, sr_interact = FALSE,
                      label = "") {
  cc <- na.omit(ds[, c("fecha", yv, xv), with = FALSE])
  setnames(cc, c("fecha", "y", "x"))
  cc[, D := as.numeric(fecha >= POSTIT)]
  cc[, dy := y - shift(y)][, dx := x - shift(x)]
  cc[, y_l1 := shift(y)][, x_l1 := shift(x)]
  cc[, Dy_l1 := D * y_l1][, Dx_l1 := D * x_l1]
  for (j in seq_len(p - 1)) cc[[paste0("dy_l", j)]] <- shift(cc$dy, j)
  for (j in seq_len(q) - 1) if (j > 0) cc[[paste0("dx_l", j)]] <- shift(cc$dx, j)
  rhs <- c("D", "y_l1", "x_l1", "Dy_l1", "Dx_l1", "dx",
           paste0("dy_l", seq_len(p - 1)),
           if (q > 1) paste0("dx_l", seq_len(q - 1)))
  if (sr_interact) {
    cc[, Ddx := D * dx]
    for (j in seq_len(p - 1)) cc[[paste0("Ddy_l", j)]] <- cc$D * cc[[paste0("dy_l", j)]]
    rhs <- c(rhs, "Ddx", paste0("Ddy_l", seq_len(p - 1)))
  }
  cc <- na.omit(cc[, c("dy", rhs), with = FALSE])
  m <- lm(as.formula(paste("dy ~", paste(rhs, collapse = "+"))), cc)
  V <- sandwich::NeweyWest(m, lag = 12, prewhite = FALSE)
  dm <- function(g) car::deltaMethod(m, g, vcov. = V)
  th_pre  <- dm("-x_l1/y_l1")
  th_post <- dm("-(x_l1 + Dx_l1)/(y_l1 + Dy_l1)")
  th_diff <- dm("-(x_l1 + Dx_l1)/(y_l1 + Dy_l1) + x_l1/y_l1")
  zdiff <- th_diff$Estimate / th_diff$SE
  # joint Wald: any change in the levels dynamics (gamma_D = rho_D = 0)
  R <- names(coef(m)) %in% c("Dy_l1", "Dx_l1")
  b <- coef(m)[R]; Vb <- V[R, R]
  wald <- as.numeric(t(b) %*% solve(Vb) %*% b)
  data.table(
    y = yv, spec = label, n = nrow(cc), p_lags = p, q_lags = q,
    rho_pre = round(coef(m)["y_l1"], 4),
    rho_post = round(coef(m)["y_l1"] + coef(m)["Dy_l1"], 4),
    theta_pre = round(th_pre$Estimate, 5), theta_pre_se = round(th_pre$SE, 5),
    theta_pre_lo90 = round(th_pre$Estimate - Z90 * th_pre$SE, 5),
    theta_pre_hi90 = round(th_pre$Estimate + Z90 * th_pre$SE, 5),
    theta_post = round(th_post$Estimate, 5), theta_post_se = round(th_post$SE, 5),
    theta_post_lo90 = round(th_post$Estimate - Z90 * th_post$SE, 5),
    theta_post_hi90 = round(th_post$Estimate + Z90 * th_post$SE, 5),
    equality_z = round(zdiff, 3),
    equality_p = round(2 * pnorm(-abs(zdiff)), 4),
    joint_wald_chi2 = round(wald, 3),
    joint_wald_p = round(pchisq(wald, df = 2, lower.tail = FALSE), 4))
}
brk <- rbindlist(list(
  break_one(d, "FCI_EXP_ENDO_exCredit", p = 2, q = 1, label = "baseline (2,1)"),
  break_one(d, "FCI_EXP_ENDO_exCredit", p = 3, q = 2, label = "adjacent lags (3,2)"),
  break_one(d, "FCI_EXP_ENDO_exCredit", p = 2, q = 1, sr_interact = TRUE,
            label = "short-run interactions"),
  break_one(d, "FCI_EXP_ENDO_exC_exTCN", p = 2, q = 1, label = "baseline (2,1)"),
  break_one(d, "FCI_EXP_ENDO_exC_exTCN", p = 3, q = 2, label = "adjacent lags (3,2)")),
  fill = TRUE)
write_rev_csv(brk, "Rev_PV_LongRun_Break.csv")
print(brk[, .(y, spec, theta_pre, theta_post, equality_p, joint_wald_p)])

# ---------------------------------------------------------------------------
# Block 9 — Gregory-Hansen (1996) break-robust cointegration (round-8 add)
# Residual-based cointegration allowing a one-time break: model C (level
# shift), C/T (level shift with trend), C/S (regime shift).  Reported both
# ways: (a) the standard GH scan (min ADF over 15%-trimmed candidate breaks,
# endogenous break location) against GH (1996, J. Econometrics, Table 1)
# asymptotic CVs for m = 1 regressor; (b) the ADF statistic at the IMPOSED
# 2011m5 IT break, evaluated conservatively against the same scan CVs.
# Interpreted narrowly: this concerns the DXY-FCI first-stage levels
# relation only, not the credit exclusion restriction.
# ---------------------------------------------------------------------------
cat("\n[9] Gregory-Hansen break-robust cointegration...\n")
GH_CV <- list(  # Gregory & Hansen (1996), Table 1, ADF*, m = 1
  C    = c(p01 = -5.13, p05 = -4.61, p10 = -4.34),
  CT   = c(p01 = -5.45, p05 = -4.99, p10 = -4.72),
  CS   = c(p01 = -5.47, p05 = -4.95, p10 = -4.68))

adf_resid_t <- function(res) {
  # ADF t on residuals, no deterministic terms, BIC lag selection (max 12)
  fit <- tryCatch(urca::ur.df(res, type = "none", lags = 12, selectlags = "BIC"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  as.numeric(fit@teststat[1])
}

gh_scan <- function(ds, yv, xv, sample) {
  cc <- as.data.frame(na.omit(ds[, c("fecha", yv, xv), with = FALSE]))
  n <- nrow(cc); if (n < 100) return(NULL)
  y <- cc[[yv]]; x <- cc[[xv]]; tr <- seq_len(n)
  lo <- floor(0.15 * n); hi <- ceiling(0.85 * n)
  out <- list()
  for (model in c("C", "CT", "CS")) {
    stats <- rep(NA_real_, n)
    for (b in lo:hi) {
      D <- as.numeric(tr > b)
      res <- switch(model,
        C  = resid(lm(y ~ D + x)),
        CT = resid(lm(y ~ D + tr + x)),
        CS = resid(lm(y ~ D + x + I(D * x))))
      stats[b] <- adf_resid_t(res)
    }
    i_min <- which.min(stats)
    # imposed IT break
    b_it <- max(which(cc$fecha < POSTIT))
    D <- as.numeric(tr > b_it)
    res_it <- switch(model,
      C  = resid(lm(y ~ D + x)),
      CT = resid(lm(y ~ D + tr + x)),
      CS = resid(lm(y ~ D + x + I(D * x))))
    t_it <- adf_resid_t(res_it)
    cv <- GH_CV[[model]]
    verdict <- function(t) ifelse(is.na(t), NA_character_,
      ifelse(t < cv["p01"], "reject at 1%", ifelse(t < cv["p05"], "reject at 5%",
      ifelse(t < cv["p10"], "reject at 10%", "no rejection"))))
    out[[model]] <- data.table(
      sample = sample, y = yv, x = xv, model = model, n = n,
      GH_stat = round(stats[i_min], 3),
      break_at_min = as.character(cc$fecha[i_min]),
      verdict_scan = verdict(stats[i_min]),
      t_at_2011m5 = round(t_it, 3),
      verdict_2011m5 = verdict(t_it),
      cv_1pct = cv["p01"], cv_5pct = cv["p05"], cv_10pct = cv["p10"])
  }
  rbindlist(out)
}

gh <- rbindlist(list(
  gh_scan(d, "FCI_EXP_ENDO_exCredit",  "DXY", "full"),
  gh_scan(d, "FCI_EXP_ENDO_exC_exTCN", "DXY", "full")), fill = TRUE)
write_rev_csv(gh, "Rev_PV_GregoryHansen.csv")
print(gh[, .(y, model, GH_stat, break_at_min, verdict_scan, t_at_2011m5, verdict_2011m5)])

# ---------------------------------------------------------------------------
# Block 10 — fixed-lag ARDL bounds grid (round-8 add)
# The headline bounds use BIC-selected orders (Block 4) and a whitened refit
# (Block 7b); referees' standard request is a fixed-lag grid.  Bounds F
# (case 3, asymptotic and exact finite-sample p) at imposed symmetric orders.
# ---------------------------------------------------------------------------
cat("\n[10] Fixed-lag ARDL bounds grid...\n")
grid_one <- function(ds, yv, xv, sample, p) {
  cc <- as.data.frame(na.omit(ds[, c(yv, xv), with = FALSE]))
  if (nrow(cc) < 80) return(NULL)
  f <- as.formula(paste(yv, "~", xv))
  bm <- tryCatch(ARDL::ardl(f, data = cc, order = c(p, p)), error = function(e) NULL)
  if (is.null(bm)) return(NULL)
  bFa <- tryCatch(ARDL::bounds_f_test(bm, case = 3), error = function(e) NULL)
  bFe <- tryCatch(ARDL::bounds_f_test(bm, case = 3, exact = TRUE), error = function(e) NULL)
  ue  <- tryCatch(ARDL::uecm(bm), error = function(e) NULL)
  bg12 <- tryCatch(lmtest::bgtest(ue, order = 12)$p.value, error = function(e) NA)
  ect <- NA_real_
  if (!is.null(ue)) {
    cn <- names(coef(ue)); ly <- grep(paste0("^L\\(", yv), cn, value = TRUE)
    if (length(ly)) ect <- coef(ue)[ly[1]]
  }
  data.table(sample = sample, y = yv, x = xv, order = sprintf("(%d,%d)", p, p),
             n = nrow(cc),
             F_stat = if (!is.null(bFa)) round(as.numeric(bFa$statistic), 3) else NA,
             F_p_asy = if (!is.null(bFa)) round(as.numeric(bFa$p.value), 4) else NA,
             F_p_exact = if (!is.null(bFe)) round(as.numeric(bFe$p.value), 4) else NA,
             ect_coef = round(ect, 4), bg_p_order12 = round(bg12, 3))
}
lag_grid <- list()
for (s in names(samples)) for (yv in c("FCI_EXP_ENDO_exCredit", "FCI_EXP_ENDO_exC_exTCN"))
  for (p in c(2, 4, 6, 8, 12))
    lag_grid[[paste(s, yv, p)]] <- grid_one(samples[[s]], yv, "DXY", s, p)
lag_grid <- rbindlist(lag_grid)
write_rev_csv(lag_grid, "Rev_PV_ARDL_LagGrid.csv")
print(lag_grid[, .(sample, y, order, F_stat, F_p_asy, F_p_exact, ect_coef, bg_p_order12)])

cat(sprintf("=== 50 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
