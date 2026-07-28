# ============================================================================
# 69_Screening_Transformation_Matched.R
# TRANSFORMATION-MATCHED sensitivity for the candidate-instrument screening
# table (main-text Table 5; produced by 18_FCI_Improved_IV_LP.R SECTION A,
# output/csv/IV_First_Stage_Battery.csv).
#
# The objection. Table 5 ranks seven candidate external variables by their
# partial explanatory power for domestic financial conditions in a COMMON
# LEVELS specification. But the candidates do not receive a common
# transformation: six enter as levels (DXY, VIX, FFER, US_10Y, SP500,
# Selic_rate) while the terms-of-trade candidate enters as a 12-month growth
# rate (d_ToT, built at 17_FCI_New_External_Data.R:320). Partial-R^2 and t^2
# comparisons between near-I(1) and near-I(0) regressors quasi-mechanically
# favour the persistent candidate, so "the DXY has by far the largest in-sample
# low-frequency association" is vulnerable to the charge that it is an artefact
# of differential persistence rather than a fact about the dollar.
#
# What this script does. It re-runs the screening under three transformation
# grids applied IDENTICALLY to all seven candidates, and reports a persistence
# diagnostic beside every estimate so the comparison is self-documenting:
#
#   G1 "level"  — every candidate in levels. (This is the published
#                 specification except that the terms of trade enters as the
#                 LEVEL Paraguay_ToT rather than as d_ToT, which is the single
#                 change that makes the published grid internally consistent.)
#   G2 "d12"    — every candidate as a 12-month change: 100 x log-difference
#                 for price/index series (DXY, SP500, IPE, Paraguay_ToT),
#                 simple 12-month difference for rate/level series (FFER,
#                 US_10Y, Selic_rate, VIX). This imposes the d_ToT rule on
#                 every candidate and removes the persistence differential.
#   G3 "ma12"   — every candidate as a 12-month trailing moving average of the
#                 level: a common LOW-FREQUENCY filter that preserves
#                 persistence but equalises it across candidates.
#
# and two lag conventions:
#   "contemp"    — no own lags (the published Table 5 convention);
#   "lagmatched" — own L1 and L2 for every candidate, with the restriction
#                  tested jointly (anova F / HAC Wald) rather than as a t^2.
#
# Inference. Newey-West(3) is retained for exact comparability with Table 5,
# and Newey-West(13) is reported alongside because a fixed lag of 3 is far too
# short for near-unit-root levels (the bandwidth 44_Enhanced_Instrument_
# Robustness.R:85 uses). Classical partial R^2 and F follow 18:165-170.
#
# Sample. Held fixed WITHIN each grid by na.omit over the union of all
# transformed candidates and controls, so the 12-month filters' burn-in cannot
# differentially shorten one candidate's sample. N is reported per grid.
#
# This script writes only Rev_Screening_TransformationMatched.csv (+ a compact
# rank summary). It does NOT overwrite IV_First_Stage_Battery.csv: main-text
# Table 5 must remain reproducible from the script that produced it.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
suppressMessages(library(tseries))
cat("=== 69: transformation-matched candidate screening ===\n")

d <- load_ext_data()

DEP <- "FCI_ENDO_exCredit_AVG"   # 18:140-143

# Candidate -> (level column, transformation family for the d12 grid)
#   "logdiff" = price/index series, 100 * 12-month log difference
#   "diff"    = rate or index-level series, simple 12-month difference
CANDS <- list(
  DXY   = list(col = "DXY",          fam = "logdiff"),
  FFER  = list(col = "FFER",         fam = "diff"),
  US10Y = list(col = "US_10Y",       fam = "diff"),
  Selic = list(col = "Selic_rate",   fam = "diff"),
  ToT   = list(col = "Paraguay_ToT", fam = "logdiff"),
  SP500 = list(col = "SP500",        fam = "logdiff"),
  VIX   = list(col = "VIX",          fam = "diff"))

d12 <- function(x, fam) {
  if (fam == "logdiff") 100 * (log(x) - log(shift(x, 12))) else x - shift(x, 12)
}
ma12 <- function(x) {
  as.numeric(stats::filter(x, rep(1 / 12, 12), sides = 1))
}

GRIDS <- c("level", "d12", "ma12")

# ---- Build every transformed candidate column ------------------------------
for (nm in names(CANDS)) {
  cc <- CANDS[[nm]]$col; fam <- CANDS[[nm]]$fam
  d[, (paste0("X_level_", nm)) := get(cc)]
  d[, (paste0("X_d12_",   nm)) := d12(get(cc), fam)]
  d[, (paste0("X_ma12_",  nm)) := ma12(get(cc))]
}

# ---- Persistence diagnostics per transformed candidate ---------------------
pers_row <- function(x, grid, cand) {
  x <- na.omit(as.numeric(x))
  if (length(x) < 40) return(NULL)
  a  <- suppressWarnings(adf.test(x))
  pp <- suppressWarnings(pp.test(x))
  k  <- suppressWarnings(kpss.test(x, null = "Level"))
  ar1 <- suppressWarnings(tryCatch(coef(ar(x, order.max = 1, aic = FALSE,
                                           method = "ols"))[1],
                                   error = function(e) NA_real_))
  ar3 <- suppressWarnings(tryCatch(sum(coef(ar(x, order.max = 3, aic = FALSE,
                                               method = "ols"))[1:3]),
                                   error = function(e) NA_real_))
  data.table(grid = grid, candidate = cand, n_series = length(x),
             adf_p = round(a$p.value, 3), pp_p = round(pp$p.value, 3),
             kpss_p = round(k$p.value, 3),
             ar1 = round(as.numeric(ar1), 3),
             sum_ar3 = round(as.numeric(ar3), 3),
             verdict = ifelse(a$p.value > 0.10 & k$p.value <= 0.05, "I(1)-consistent",
                       ifelse(a$p.value <= 0.05 & k$p.value > 0.05, "I(0)-consistent",
                              "mixed/borderline")))
}

pers <- rbindlist(lapply(GRIDS, function(g)
  rbindlist(lapply(names(CANDS), function(nm)
    pers_row(d[[paste0("X_", g, "_", nm)]], g, nm)))), fill = TRUE)
write_rev_csv(pers, "Rev_Screening_Persistence.csv")
cat("\n--- Persistence of the transformed candidates ---\n")
print(pers)

# ---- Screening regressions -------------------------------------------------
cat("\n--- Screening regressions by grid and lag convention ---\n")
res <- list()

for (g in GRIDS) {
  xcols <- paste0("X_", g, "_", names(CANDS))
  for (lagconv in c("contemp", "lagmatched")) {
    # Lags of every candidate, so the fixed sample is common across candidates
    lagcols <- character(0)
    if (lagconv == "lagmatched") {
      for (xc in xcols) {
        d[, (paste0(xc, "_L1")) := shift(get(xc), 1)]
        d[, (paste0(xc, "_L2")) := shift(get(xc), 2)]
        lagcols <- c(lagcols, paste0(xc, c("_L1", "_L2")))
      }
    }
    keep <- c(DEP, EXOG_MACRO, xcols, lagcols)
    fsd  <- na.omit(d[, ..keep])
    Nfix <- nrow(fsd)

    ctrl <- paste(EXOG_MACRO, collapse = " + ")
    m_base <- lm(as.formula(paste(DEP, "~", ctrl)), data = fsd)
    ssr0   <- sum(resid(m_base)^2)

    for (nm in names(CANDS)) {
      xc   <- paste0("X_", g, "_", nm)
      regs <- if (lagconv == "lagmatched") paste0(xc, c("", "_L1", "_L2")) else xc
      m    <- lm(as.formula(paste(DEP, "~", paste(regs, collapse = " + "),
                                  "+", ctrl)), data = fsd)

      V3  <- sandwich::NeweyWest(m, lag = 3,  prewhite = FALSE)
      V13 <- sandwich::NeweyWest(m, lag = 13, prewhite = FALSE)

      # Own-coefficient HAC t on the contemporaneous term (comparable to Table 5)
      ct3  <- lmtest::coeftest(m, vcov. = V3)
      ct13 <- lmtest::coeftest(m, vcov. = V13)
      b    <- unname(coef(m)[xc])
      t3   <- unname(ct3[xc, 3]);  p3  <- unname(ct3[xc, 4])
      t13  <- unname(ct13[xc, 3]); p13 <- unname(ct13[xc, 4])

      # Classical partial R^2 on the full restriction (generalises 18:168)
      pr2 <- 1 - sum(resid(m)^2) / ssr0
      # Classical F: joint when lag-matched, t^2 when contemporaneous (18:170)
      f_cls <- if (lagconv == "lagmatched") {
        unname(anova(m_base, m)$F[2])
      } else {
        unname(summary(m)$coefficients[xc, 3]^2)
      }
      # HAC joint Wald on the same restriction
      f_hac3 <- suppressWarnings(tryCatch(
        unname(lmtest::waldtest(m_base, m, vcov = V3,  test = "F")$F[2]),
        error = function(e) NA_real_))
      f_hac13 <- suppressWarnings(tryCatch(
        unname(lmtest::waldtest(m_base, m, vcov = V13, test = "F")$F[2]),
        error = function(e) NA_real_))

      res[[length(res) + 1]] <- data.table(
        grid = g, lag_convention = lagconv, candidate = nm, N = Nfix,
        coef = b, t_nw3 = t3, p_nw3 = p3, t_nw13 = t13, p_nw13 = p13,
        partial_R2_classical = pr2, F_classical = f_cls,
        F_hac_nw3 = f_hac3, F_hac_nw13 = f_hac13)
    }
  }
}

out <- rbindlist(res)
# Rank within grid x lag convention by classical partial R^2 (Table 5's metric)
out[, rank_partialR2 := frank(-partial_R2_classical, ties.method = "min"),
    by = .(grid, lag_convention)]
out[, ratio_to_next := {
  s <- sort(partial_R2_classical, decreasing = TRUE)
  ifelse(rank_partialR2 == 1L, s[1] / s[2], NA_real_)
}, by = .(grid, lag_convention)]
write_rev_csv(out, "Rev_Screening_TransformationMatched.csv")

for (g in GRIDS) for (lc in c("contemp", "lagmatched")) {
  cat(sprintf("\n--- grid=%s, lags=%s (N=%d) ---\n", g, lc,
              out[grid == g & lag_convention == lc, unique(N)]))
  print(out[grid == g & lag_convention == lc][order(rank_partialR2),
    .(candidate, coef = round(coef, 4), t_nw3 = round(t_nw3, 2),
      t_nw13 = round(t_nw13, 2), pR2 = round(partial_R2_classical, 4),
      F_cls = round(F_classical, 1), rank = rank_partialR2)])
}

# ---- Compact summary: does the DXY still rank first? -----------------------
summ <- out[, .(top_candidate = candidate[which.min(rank_partialR2)],
                top_pR2 = round(max(partial_R2_classical), 4),
                dxy_pR2 = round(partial_R2_classical[candidate == "DXY"], 4),
                dxy_rank = rank_partialR2[candidate == "DXY"],
                dxy_t_nw3 = round(t_nw3[candidate == "DXY"], 2),
                dxy_t_nw13 = round(t_nw13[candidate == "DXY"], 2),
                dxy_ratio_to_next = round(ratio_to_next[candidate == "DXY"], 2),
                N = unique(N)),
            by = .(grid, lag_convention)]
write_rev_csv(summ, "Rev_Screening_RankSummary.csv")
cat("\n--- Rank summary ---\n")
print(summ)

cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
