# ============================================================================
# 68_Symmetric_Gradient_Comparison.R
# SYMMETRIC dollar-vs-global-risk comparison at the FX-credit-exposure gradient.
#
# Motivation. Scripts 48 and 64 estimate the dollar and VIX exposure gradients
# in a matched estimating equation, but NOT with symmetrically constructed
# shocks: 30_Bank_Micro_Data_Construction.R:114-115 purges the DXY innovation
# against contemporaneous d_vix (among five other global factors), while
# 30:118 leaves the VIX innovation as a raw AR(3) residual. The manuscript
# discloses this and states that the fully symmetric comparison "is not
# estimated". Both referee-side audits identify it as the one required new
# estimation, because global financial conditions are Reviewer 1's leading
# alternative channel and the asymmetric comparator cannot adjudicate it.
#
# This script closes that gap. Two symmetric variants are estimated:
#
#   VARIANT "sym"    (primary): both innovations residualised on the SAME
#                    common factor set, with NEITHER purged against the other.
#                    Common set = dln_ipe, dln_tot, d_ffer, dln_sp500, d_us10y
#                    (i.e. d_vix is removed from the DXY purge, and the VIX
#                    innovation receives the five remaining factors it never
#                    had). This is the symmetric PLACEBO: the two carriers are
#                    treated identically, so a gradient difference cannot be an
#                    artefact of differential residualisation.
#
#   VARIANT "mutual" (secondary): each additionally purged against the other
#                    (DXY gets the common set + d_vix; VIX gets the common set
#                    + dln_dxy). This is the HORSE RACE: each carrier holds the
#                    other fixed. Reported because "symmetric" is ambiguous
#                    between these two readings and they answer different
#                    questions.
#
# Everything else is held at the published convention: identical AR lag order
# (3), identical standardisation rule and window (z on 2016m1-2025m12,
# 30:123-129), identical exposure construction and standardisation (64:56-61,
# observation-weighted after the merge), identical estimating equation, fixed
# effects, controls, per-horizon sample filter, inference battery and seeds
# (20260710 / 20260712).
#
# Parts
#   1. Innovation construction + correlation gate (auditability: how far the
#      symmetric carriers move from the published ones, and how correlated the
#      two symmetric carriers are with each other).
#   2. Standalone gradients per variant per shock, h = 3, 6, 12, 18.
#   3. Joint specification + gradient-difference Wald test, h = 12, 18, under
#      every vcov in the battery plus a reparameterised-carrier WCB.
#
# Nothing here overwrites an existing output: scripts 30/31/48/64 are untouched
# and micro_designA_panel.rds is read, never rewritten.
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
library(fixest)
set.seed(20260710)
cat("=== 68: symmetric DXY-vs-VIX gradient comparison ===\n")
cat(sprintf("WCB replications (MICRO_WCB_B): %d\n", WCB_B))

dA    <- read_rds_micro("micro_designA_panel.rds")
p2    <- read_rds_micro("micro_p2_sector.rds")
macro <- read_rds_micro("micro_macro.rds")

CTRLS <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
FE_A  <- "bst + bsc"

dA[, hedge_class := fifelse(sector %in% SECT_HEDGED,   "hedged",
                    fifelse(sector %in% SECT_UNHEDGED, "unhedged",
                    fifelse(sector %in% SECT_AMBIGUOUS, "ambiguous", "excluded")))]

# ============================================================================
# Part 1: symmetric innovations
# ar_innov is copied verbatim from 30:90-99 rather than refactored, because
# script 30 is a seeded canonical script whose outputs the appendix quotes.
# ============================================================================
cat("\n--- Part 1: symmetric innovation construction ---\n")

ar_innov <- function(x, extra = NULL, nlag = 3) {
  X <- sapply(1:nlag, function(l) shift(x, l))
  colnames(X) <- paste0("l", 1:nlag)
  df <- data.frame(y = x, X)
  if (!is.null(extra)) df <- cbind(df, extra)
  r <- rep(NA_real_, length(x))
  ok <- complete.cases(df)
  if (sum(ok) > 30) r[ok] <- resid(lm(y ~ ., data = df[ok, , drop = FALSE]))
  r
}

# Common factor set: the published DXY purge set MINUS the VIX, so that neither
# target variable appears in the other's residualisation.
COMMON <- c("dln_ipe", "dln_tot", "d_ffer", "dln_sp500", "d_us10y")

# Variant "sym": identical treatment, neither purged against the other.
macro[, shk_dxy_sym := ar_innov(dln_dxy, extra = macro[, ..COMMON])]
macro[, shk_vix_sym := ar_innov(d_vix,   extra = macro[, ..COMMON])]

# Variant "mutual": each additionally purged against the other.
macro[, shk_dxy_mut := ar_innov(dln_dxy, extra = macro[, c(COMMON, "d_vix"),   with = FALSE])]
macro[, shk_vix_mut := ar_innov(d_vix,   extra = macro[, c(COMMON, "dln_dxy"), with = FALSE])]

# Standardisation: identical rule and window to 30:123-129.
est_win <- macro$ym >= to_ym("2016-01-01") & macro$ym <= to_ym("2025-12-01")
NEWSHK  <- c("shk_dxy_sym", "shk_vix_sym", "shk_dxy_mut", "shk_vix_mut")
for (v in NEWSHK) {
  mu  <- mean(macro[[v]][est_win], na.rm = TRUE)
  sd_ <- sd(macro[[v]][est_win],   na.rm = TRUE)
  macro[, (v) := (get(v) - mu) / sd_]
}

# ---- Correlation gate ------------------------------------------------------
cg <- function(a, b) cor(macro[[a]][est_win], macro[[b]][est_win], use = "complete.obs")
corr_tab <- data.table(
  pair = c("published purged DXY vs symmetric DXY",
           "published raw VIX vs symmetric VIX",
           "symmetric DXY vs symmetric VIX",
           "mutual DXY vs mutual VIX",
           "published purged DXY vs published raw VIX",
           "published purged DXY vs mutual DXY",
           "published raw VIX vs mutual VIX"),
  corr = c(cg("shk_dxy_purged", "shk_dxy_sym"),
           cg("shk_vix",        "shk_vix_sym"),
           cg("shk_dxy_sym",    "shk_vix_sym"),
           cg("shk_dxy_mut",    "shk_vix_mut"),
           cg("shk_dxy_purged", "shk_vix"),
           cg("shk_dxy_purged", "shk_dxy_mut"),
           cg("shk_vix",        "shk_vix_mut")))
corr_tab[, corr := round(corr, 4)]
print(corr_tab)
write_micro_csv(corr_tab, "Micro_SymGradient_Correlations.csv")

# ---- Merge onto the existing panel (never rewritten) -----------------------
dA <- merge(dA, macro[, c("ym", NEWSHK), with = FALSE], by = "ym", all.x = TRUE)
for (v in NEWSHK) dA[, paste0(v, "_usd") := get(v) * usd]

# ---- Exposure: identical construction to 64:56-61 --------------------------
fx16 <- p2[ym <= to_ym("2016-12-01"), .(tot = sum(total)), by = .(sector, cur)]
fx16 <- dcast(fx16, sector ~ cur, value.var = "tot")
setnames(fx16, c("6200", "6900"), c("fc", "pyg"))
fx16[, fx_share_2016 := fc / (fc + pyg)]
dA <- merge(dA, fx16[, .(sector, fx_share_2016)], by = "sector", all.x = TRUE)
dA[, fxsh_z := zstd(fx_share_2016)]

for (v in NEWSHK) dA[, paste0(v, "_usd_fxsh") := get(paste0(v, "_usd")) * fxsh_z]

grad <- dA[hedge_class != "excluded"]
cat(sprintf("Gradient sample: %d rows, %d sectors, %d banks\n",
            nrow(grad), uniqueN(grad$sector), uniqueN(grad$bank)))

# ============================================================================
# Part 2: standalone gradients (mirrors 64 Part 1 exactly, per shock/variant)
# ============================================================================
cat("\n--- Part 2: standalone symmetric gradients ---\n")

SPECS <- list(
  list(variant = "sym",    shock = "DXY", stub = "shk_dxy_sym"),
  list(variant = "sym",    shock = "VIX", stub = "shk_vix_sym"),
  list(variant = "mutual", shock = "DXY", stub = "shk_dxy_mut"),
  list(variant = "mutual", shock = "VIX", stub = "shk_vix_mut"))

res <- list()
for (sp in SPECS) {
  lvl  <- paste0(sp$stub, "_usd")
  grd  <- paste0(sp$stub, "_usd_fxsh")
  for (h in c(3, 6, 12, 18)) {
    yv  <- paste0("g", h)
    ds  <- grad[get(paste0("ok", h)) == TRUE]
    rhs <- c(lvl, grd, CTRLS)
    m <- feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A)),
               data = ds, cluster = ~bank, notes = FALSE)
    ib <- inference_battery(m, grd, dk_lag = h + 1)
    p_wcb_month <- if (h %in% c(12, 18)) {
      wcb_pval(ds, yv, rhs, FE_A, grd, cluster = "ym", B = WCB_B)
    } else NA_real_
    res[[length(res) + 1]] <- data.table(
      variant = sp$variant, shock = sp$shock, h = h, b = ib$b,
      se_bank = ib$se_bank, p_bank = ib$p_bank,
      se_twoway = ib$se_twoway, p_twoway = ib$p_twoway,
      se_dk = ib$se_dk, p_dk = ib$p_dk,
      p_wcb_month = p_wcb_month, n = ib$n)
    cat(sprintf("%-6s %s h=%2d: b=%6.3f | p_bank=%.3f | p_twoway=%.3f | p_dk=%.3f | p_wcb_month=%s | n=%d\n",
                sp$variant, sp$shock, h, ib$b, ib$p_bank, ib$p_twoway, ib$p_dk,
                ifelse(is.na(p_wcb_month), "  -  ", sprintf("%.3f", p_wcb_month)), ib$n))
  }
}
out2 <- rbindlist(res)
write_micro_csv(out2, "Micro_SymGradient_TimeRobust.csv")

# ============================================================================
# Part 3: joint specification + gradient-difference test
# Structure and algebra copied from 64:205-269 (wald_diff at 64:225-232 and
# the reparameterised bootstrap carrier at 64:247).
# ============================================================================
cat("\n--- Part 3: symmetric DXY-vs-VIX gradient difference ---\n")
set.seed(20260712)

wald_diff <- function(m, V, nD, nV) {
  b  <- coef(m)
  d  <- unname(b[nD] - b[nV])
  v  <- V[nD, nD] + V[nV, nV] - 2 * V[nD, nV]
  se <- sqrt(v)
  c(d = d, se = se, p = 2 * pnorm(-abs(d / se)))
}

res3 <- list()
for (vr in c("sym", "mutual")) {
  sfx <- if (vr == "sym") c("shk_dxy_sym", "shk_vix_sym") else c("shk_dxy_mut", "shk_vix_mut")
  nD  <- paste0(sfx[1], "_usd_fxsh"); nV <- paste0(sfx[2], "_usd_fxsh")
  lD  <- paste0(sfx[1], "_usd");      lV <- paste0(sfx[2], "_usd")
  dA[, grad_sum := get(nD) + get(nV)]
  gr <- dA[hedge_class != "excluded"]
  for (h in c(12, 18)) {
    yv  <- paste0("g", h)
    ds  <- gr[get(paste0("ok", h)) == TRUE]
    rhs <- c(lD, lV, nD, nV, CTRLS)
    m   <- feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A)),
                 data = ds, cluster = ~bank, notes = FALSE)

    w_bank <- wald_diff(m, vcov(m, vcov = ~bank),      nD, nV)
    w_2way <- wald_diff(m, vcov(m, vcov = ~bank + ym), nD, nV)
    w_sec  <- wald_diff(m, vcov(m, vcov = ~sector),    nD, nV)
    w_dk   <- wald_diff(m, vcov(m, vcov = as.formula(sprintf("DK(%d) ~ ym", h + 1))),
                        nD, nV)

    rhs_rp <- c(lD, lV, "grad_sum", nV, CTRLS)
    p_wcb_month  <- wcb_pval(ds, yv, rhs_rp, FE_A, nV, cluster = "ym", B = WCB_B)
    p_wcb_sector <- wcb_pval(ds, yv, rhs_rp, FE_A, nV, cluster = "sector",
                             B = WCB_B, weights = "webb")

    res3[[length(res3) + 1]] <- data.table(
      variant = vr, h = h,
      b_dxy_gradient = unname(coef(m)[nD]),
      b_vix_gradient = unname(coef(m)[nV]),
      diff = unname(w_bank["d"]),
      se_bank = unname(w_bank["se"]), p_bank = unname(w_bank["p"]),
      se_twoway = unname(w_2way["se"]), p_twoway = unname(w_2way["p"]),
      se_sector = unname(w_sec["se"]), p_sector = unname(w_sec["p"]),
      se_dk = unname(w_dk["se"]), p_dk = unname(w_dk["p"]),
      p_wcb_month = p_wcb_month, p_wcb_webb_sector = p_wcb_sector,
      n = nobs(m))
    cat(sprintf("%-6s h=%2d: b_DXY=%6.3f b_VIX=%6.3f diff=%6.3f | p_bank=%.3f p_2way=%.3f p_sec=%.3f p_dk=%.3f | p_wcb_month=%.3f p_wcb_webb=%.3f\n",
                vr, h, coef(m)[nD], coef(m)[nV], w_bank["d"],
                w_bank["p"], w_2way["p"], w_sec["p"], w_dk["p"],
                p_wcb_month, p_wcb_sector))
  }
}
write_micro_csv(rbindlist(res3), "Micro_SymGradient_DXY_vs_VIX_Diff.csv")

cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
