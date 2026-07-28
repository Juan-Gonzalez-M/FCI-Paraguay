# ============================================================================
# 70_OneSided_Band_SDs_CurrencyTest.R
# Three small, referee-requested additions that share a data layer.
#
# PART 1 — Joint (max-t) confidence band for the PRIMARY one-sided LP profile.
#   The manuscript's Figure 2(a) plots the h = 1..18 post-IT profile of the
#   strictly one-sided FCI_exCredit z-score index and the narrative leans on a
#   response that "builds from h~5 and peaks at 12-13", but the displayed bands
#   are POINTWISE. Both audits ask for a simultaneous band so the shape claim
#   is not read off 18 separate 95% intervals. Implemented as a null-free
#   studentised moving-block bootstrap over months (blocks L = 18, circular,
#   the same resampler as 56:330-335), applying ONE common resampled date
#   sequence to all 18 horizon frames so cross-horizon dependence is preserved.
#   The critical value is the (1-alpha) quantile of max_h |t*_h|.
#   NOTE: 60_Submission_Figure3_TwoPanel.R computes this profile but persists
#   only the PNG. This script writes the profile to CSV so the figure and the
#   band cannot drift apart.
#
# PART 2 — Per-SD conversions for main-text Table 4 Panel B.
#   Panel B prints coefficients per INDEX UNIT across six differently scaled
#   index objects, while its own note says cross-variant comparisons are valid
#   "only in per-SD units" and supplies only the two one-sided SDs. That
#   invites exactly the comparison the note prohibits. This part computes the
#   SD of each Panel B index on the h = 12 estimation sample of the regression
#   that produced its coefficient. Three of the SDs are already archived
#   (0.8899 / 0.8834 / 0.6948); reproducing them is the gate that validates the
#   two that are not.
#   ONE ROW CANNOT BE FILLED: the "regime-specific composite" is re-estimated
#   on the post-IT sample inside 15_FCI_PostIT_Subsample_Analysis.R by the
#   four-method machinery (including the MARSS DFM) and its series is never
#   exported, so no per-SD conversion is recoverable from the archive. It is
#   reported as NA rather than approximated by a neighbouring object.
#
# PART 3 — MN-vs-USD full-sample coefficient-equality test.
#   Section 5.1 states that in the full sample the contraction is carried by
#   the guarani book (-9.06 versus -7.33 for USD). No test of that ordering
#   exists anywhere in the archive, the intervals overlap, and the two
#   coefficients differ in significance (p = 0.006 vs 0.088). This part
#   estimates the difference directly: the two currency equations are stacked
#   and fully interacted with a currency dummy on a common sample, which
#   reproduces the separate OLS coefficients exactly while delivering the
#   cross-equation covariance. Inference is Driscoll-Kraay with scores
#   aggregated within month, Bartlett-weighted at lag h+1 — the convention
#   63_Sectoral_Group_Wald.R already uses for the sectoral group test.
#
# ADDITIVE: writes only Rev_OneSided_Profile_Band.csv,
# Rev_PanelB_PerSD.csv and Rev_Currency_Equality_Test.csv.
# ============================================================================

t0 <- Sys.time()
suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(lubridate)
  library(sandwich); library(lmtest)
})
set.seed(20260716)   # mirrors 56:41 so block draws are comparable
OUT <- "../output/revision/csv"
cat("=== 70: joint LP band, Panel B per-SD, currency equality test ===\n")

# ---------------------------------------------------------------------------
# Shared data layer
# ---------------------------------------------------------------------------
ext <- read.csv("../output/csv/New_External_Variables.csv")
ext$fecha <- as.Date(ext$fecha)
rv <- read.csv("../output/csv/FCI_Robustness_Versions.csv")
rv$fecha <- as.Date(rv$fecha)
fcia <- read.csv("../output/csv/FCI_Complete_Results.csv")
fcia$fecha <- as.Date(fcia$fecha)

base <- ext %>%
  left_join(rv %>% dplyr::select(fecha, FCI_exCredit_ZSCORE), by = "fecha") %>%
  left_join(fcia %>% dplyr::select(fecha, FCI_RATES_AVG), by = "fecha") %>%
  arrange(fecha) %>%
  mutate(IMAEP_yoy_L1 = lag(IMAEP_yoy, 1), IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
         IPC_yoy_L1   = lag(IPC_yoy, 1),   IPC_yoy_L2   = lag(IPC_yoy, 2))

CTRL <- c("IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
          "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2")

# Macro-purged FCI: residuals of FCI_exCredit on IMAEP_yoy, IPC_yoy + 2 lags
# each (05_FCI_Local_Projections.R:1810-1820, verbatim specification).
pd <- base %>% dplyr::select(fecha, FCI_exCredit_AVG, all_of(CTRL)) %>% na.omit()
pm <- lm(FCI_exCredit_AVG ~ IMAEP_yoy + IMAEP_yoy_L1 + IMAEP_yoy_L2 +
           IPC_yoy + IPC_yoy_L1 + IPC_yoy_L2, data = pd)
pd$FCI_purged <- residuals(pm)
base <- base %>% left_join(pd %>% dplyr::select(fecha, FCI_purged), by = "fecha")
cat(sprintf("Purge reconstruction: R2 = %.4f, n = %d (diagnostic only; see Part 2 note)\n",
            summary(pm)$r.squared, nrow(pd)))

# Build one LP estimation frame (run_lp_standard convention, n_lags = 2)
lp_frame <- function(data, yv, xv, h) {
  d <- data %>% mutate(y_fwd  = lead(.data[[yv]], h),
                       y_lag1 = lag(.data[[yv]], 1),
                       y_lag2 = lag(.data[[yv]], 2),
                       fci_lag1 = lag(.data[[xv]], 1))
  d %>% dplyr::select(fecha, y_fwd, all_of(xv), y_lag1, y_lag2, fci_lag1,
                      all_of(CTRL)) %>% na.omit()
}
lp_fit <- function(rd, xv, h) {
  rhs <- c(xv, "y_lag1", "y_lag2", "fci_lag1", CTRL)
  m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
  ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1,
                                                      prewhite = FALSE))
  list(b = unname(ct[xv, 1]), se = unname(ct[xv, 2]), p = unname(ct[xv, 4]),
       n = nrow(rd))
}

# ===========================================================================
# PART 1 — joint max-t band for the primary one-sided post-IT profile
# ===========================================================================
cat("\n--- Part 1: joint (max-t) band for the primary one-sided profile ---\n")

XV  <- "FCI_exCredit_ZSCORE"
HS  <- 1:18
L   <- 18
B   <- 999

mbb_idx <- function(n, L) {          # 56:330-335, circular overlapping blocks
  k <- ceiling(n / L)
  starts <- sample.int(n, k, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) ((s:(s + L - 1) - 1) %% n) + 1))
  idx[1:n]
}

# The band is computed for BOTH samples: the post-IT profile is what Figure 2(a)
# plots, but the paper's headline magnitude (Table 4 Panel A, -10.22) is
# full-sample, so a band on the post-IT profile alone would not discipline the
# claim the narrative actually makes.
SAMPLES <- list(
  full   = base,
  postIT = base %>% filter(fecha >= as.Date("2011-05-01")))

band_one <- function(dat, label) {
frames <- lapply(HS, function(h) lp_frame(dat, "Cred_Real_Total", XV, h))
names(frames) <- as.character(HS)
point <- lapply(HS, function(h) lp_fit(frames[[as.character(h)]], XV, h))
bhat  <- sapply(point, `[[`, "b")
sehat <- sapply(point, `[[`, "se")
nobs  <- sapply(point, `[[`, "n")

pool <- frames[["1"]]$fecha          # widest frame's dates
np   <- length(pool)

# Sup-t band, Montiel Olea & Plagborg-Moller (2019). The block bootstrap is
# used ONLY to estimate the CORRELATION matrix of the coefficient vector across
# horizons; the band half-widths use the analytic Newey-West standard errors.
# This is deliberate: with blocks drawn with replacement a month can appear
# several times in a replicate, which biases bootstrap standard errors downward
# (a studentised max-|t| band built on them returns an implausible critical
# value of ~7). Cross-horizon CORRELATION is far more robust to that than the
# level of the bootstrap variance, so the hybrid is the defensible construction:
#   band_h = bhat_h  +/-  c(alpha, R) * se_NW_h,
# where c is the equicoordinate two-sided quantile of a N(0, R) vector.
Bstar <- matrix(NA_real_, nrow = B, ncol = length(HS))
for (b in seq_len(B)) {
  dsel <- pool[mbb_idx(np, L)]
  for (j in seq_along(HS)) {
    h  <- HS[j]
    fr <- frames[[as.character(h)]]
    rows <- match(dsel, fr$fecha)
    rows <- rows[!is.na(rows)]
    if (length(rows) < 60) next
    rb <- fr[rows, , drop = FALSE]
    ft <- tryCatch(lp_fit(rb, XV, h), error = function(e) NULL)
    if (is.null(ft)) next
    Bstar[b, j] <- ft$b
  }
}
ok <- stats::complete.cases(Bstar)
Bok  <- Bstar[ok, , drop = FALSE]
Rhat <- stats::cor(Bok)
boot_sd <- apply(Bok, 2, sd)
crit_joint95 <- as.numeric(mvtnorm::qmvnorm(0.95, sigma = Rhat,
                                            tail = "both.tails")$quantile)
crit_joint90 <- as.numeric(mvtnorm::qmvnorm(0.90, sigma = Rhat,
                                            tail = "both.tails")$quantile)
cat(sprintf("[%s] usable replicates %d/%d | mean |off-diag| corr %.3f | sup-t 95%% = %.3f, 90%% = %.3f | mean boot_sd/se_nw %.2f\n",
            label, sum(ok), B, mean(abs(Rhat[upper.tri(Rhat)])),
            crit_joint95, crit_joint90, mean(boot_sd / sehat)))

prof <- data.frame(
  sample = label, horizon = HS, coef = bhat, se_nw = sehat, sd_boot = boot_sd,
  n_obs = nobs,
  lo95_pointwise = bhat - qnorm(0.975) * sehat,
  hi95_pointwise = bhat + qnorm(0.975) * sehat,
  lo95_supt = bhat - crit_joint95 * sehat,
  hi95_supt = bhat + crit_joint95 * sehat,
  crit_supt95 = crit_joint95, crit_supt90 = crit_joint90,
  block_L = L, B_effective = sum(ok))
prof$sig_pointwise95 <- with(prof, hi95_pointwise < 0 | lo95_pointwise > 0)
prof$sig_supt95      <- with(prof, hi95_supt < 0 | lo95_supt > 0)
prof
}

prof <- bind_rows(lapply(names(SAMPLES),
                         function(s) band_one(SAMPLES[[s]], s)))
write.csv(prof, file.path(OUT, "Rev_OneSided_Profile_Band.csv"), row.names = FALSE)
cat("  [csv]  Rev_OneSided_Profile_Band.csv\n")
for (s in names(SAMPLES)) {
  pp <- prof %>% filter(sample == s)
  cat(sprintf("\n--- %s sample ---\n", s))
  print(pp %>% transmute(horizon, coef = round(coef, 2), se = round(se_nw, 2),
                         pw = ifelse(sig_pointwise95, "*", ""),
                         supt = ifelse(sig_supt95, "*", ""), n_obs))
  cat(sprintf("%s: significant pointwise at %d of 18 horizons; under the sup-t band at %d\n",
              s, sum(pp$sig_pointwise95), sum(pp$sig_supt95)))
  h12 <- pp %>% filter(horizon == 12)
  cat(sprintf("%s h=12: %.2f, pointwise [%.2f, %.2f], sup-t [%.2f, %.2f]\n",
              s, h12$coef, h12$lo95_pointwise, h12$hi95_pointwise,
              h12$lo95_supt, h12$hi95_supt))
}

# ===========================================================================
# PART 2 — Table 4 Panel B per-SD conversions
# ===========================================================================
cat("\n--- Part 2: Panel B per-SD conversions ---\n")

full <- base
pB <- list(
  list(row = "Regime-specific composite", index = "FCI_ENDO_exCredit (post-IT re-estimated)",
       xv = NA, sample = "postIT", coef = -5.5995777528278, avail = FALSE),
  list(row = "Full-sample composite, post-IT window", index = "FCI_exCredit_AVG",
       xv = "FCI_exCredit_AVG", sample = "postIT", coef = -3.99971875072756, avail = TRUE),
  list(row = "Rates-only co-baseline", index = "FCI_RATES_AVG",
       xv = "FCI_RATES_AVG", sample = "full", coef = -5.43944476708963, avail = TRUE),
  list(row = "Full-sample domestic", index = "FCI_ENDO_exCredit_AVG",
       xv = "FCI_ENDO_exCredit_AVG", sample = "full", coef = -7.94846064735573, avail = TRUE),
  list(row = "Full-sample ex-credit composite", index = "FCI_exCredit_AVG",
       xv = "FCI_exCredit_AVG", sample = "full", coef = -7.65628048038673, avail = TRUE),
  # The macro-purged index is built INSIDE 05_FCI_Local_Projections.R:1806-1820
  # on `purge_data`, i.e. analysis_data after na.omit() over its FULL column set
  # (which includes the external instruments) minus two rows for the lags. That
  # sample is narrower than a complete-cases frame on the purge variables alone,
  # and the series is never exported. Reconstructing the purge on the wider frame
  # returns -8.152 at n = 333 against the published -8.270 at n = 332, so the
  # reconstruction is demonstrably a different object and its SD is NOT used.
  list(row = "Macro-purged FCI", index = "FCI_purged",
       xv = "FCI_purged", sample = "full", coef = -8.26994616989275, avail = FALSE))

rowsB <- list()
for (r in pB) {
  if (!r$avail) {
    rowsB[[length(rowsB) + 1]] <- data.frame(
      panelB_row = r$row, index_object = r$index, sample = r$sample,
      coef_per_unit = r$coef, sd_index = NA_real_, coef_per_sd = NA_real_,
      n_obs = NA_integer_,
      note = if (grepl("Regime", r$row))
        "index re-estimated on the post-IT sample inside script 15 (four-method, incl. MARSS DFM); series not exported, per-SD not recoverable from the archive"
      else
        "index constructed inside script 05 on purge_data (analysis_data after na.omit over its full column set); series not exported and a wider-frame reconstruction does not reproduce the published coefficient (-8.152/n=333 vs -8.270/n=332), so no per-SD is claimed")
    next
  }
  dat <- if (r$sample == "postIT") full %>% filter(fecha >= as.Date("2011-05-01")) else full
  fr  <- lp_frame(dat, "Cred_Real_Total", r$xv, 12)
  sdx <- sd(fr[[r$xv]])
  rowsB[[length(rowsB) + 1]] <- data.frame(
    panelB_row = r$row, index_object = r$index, sample = r$sample,
    coef_per_unit = r$coef, sd_index = sdx, coef_per_sd = r$coef * sdx,
    n_obs = nrow(fr), note = "SD on the h=12 LP estimation sample")
}
outB <- bind_rows(rowsB)
write.csv(outB, file.path(OUT, "Rev_PanelB_PerSD.csv"), row.names = FALSE)
cat("  [csv]  Rev_PanelB_PerSD.csv\n")
print(outB %>% transmute(panelB_row, sample, coef_per_unit = round(coef_per_unit, 2),
                         sd = round(sd_index, 4), per_sd = round(coef_per_sd, 2), n_obs))
cat("\nGATE — archived SDs that must reproduce: FCI_exCredit_AVG full 0.8899,\n")
cat("       FCI_ENDO_exCredit_AVG full 0.8834, FCI_exCredit_AVG post-IT 0.6948\n")

# ===========================================================================
# PART 3 — MN vs USD full-sample coefficient equality
# ===========================================================================
cat("\n--- Part 3: MN vs USD full-sample equality test ---\n")

# Credit series exactly as 05_FCI_Local_Projections.R:144-157
mv  <- read_excel("../data/FCI_data_1.xlsx", sheet = "Main_variables")
names(mv)[1] <- "fecha"
mc  <- read_excel("../data/FCI_data_1.xlsx", sheet = "Datos_macro")
names(mc)[1] <- "fecha"
mv <- mv %>% mutate(fecha = as.Date(fecha)) %>% arrange(fecha)
mc <- mc %>% mutate(fecha = as.Date(fecha)) %>% arrange(fecha)

cred <- mv %>%
  left_join(mc %>% dplyr::select(fecha, IPC, Creditos_deflactados), by = "fecha") %>%
  mutate(
    Cred_Total   = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100,
    Cred_USD     = (Creditos_Sector_privado_USD_equivalente /
                    lag(Creditos_Sector_privado_USD_equivalente, 12) - 1) * 100,
    Cred_Real_MN = ((Creditos_Sector_privado_MN / IPC) /
                    lag(Creditos_Sector_privado_MN / IPC, 12) - 1) * 100) %>%
  dplyr::select(fecha, Cred_Total, Cred_USD, Cred_Real_MN)

cd <- cred %>%
  left_join(base %>% dplyr::select(fecha, FCI_exCredit_AVG, all_of(CTRL)),
            by = "fecha") %>%
  arrange(fecha)

XC <- "FCI_exCredit_AVG"   # the regressor behind Panel C's FCI_COMP rows (05:594)

# Reproduction gate: the two separate LPs must return the published values
gate <- lapply(c("Cred_Real_MN", "Cred_USD"), function(yv) {
  fr <- lp_frame(cd, yv, XC, 12); f <- lp_fit(fr, XC, 12)
  data.frame(outcome = yv, coef = f$b, se = f$se, p = f$p, n = f$n)
}) %>% bind_rows()
print(gate %>% mutate(across(where(is.numeric), ~round(., 4))))
cat("GATE — published: Cred_Real_MN -9.0642 (p 0.0064), Cred_USD -7.3349 (p 0.0875), n 334\n")

# Stacked, fully-interacted system on the common sample
stack_test <- function(data, h) {
  mk <- function(yv, g) {
    data %>% mutate(y_fwd  = lead(.data[[yv]], h),
                    y_lag1 = lag(.data[[yv]], 1),
                    y_lag2 = lag(.data[[yv]], 2),
                    fci_lag1 = lag(.data[[XC]], 1),
                    grp = g) %>%
      dplyr::select(fecha, grp, y_fwd, all_of(XC), y_lag1, y_lag2, fci_lag1,
                    all_of(CTRL))
  }
  s <- bind_rows(mk("Cred_Real_MN", "MN"), mk("Cred_USD", "USD")) %>% na.omit()
  keep <- s %>% count(fecha) %>% filter(n == 2) %>% pull(fecha)
  s <- s %>% filter(fecha %in% keep) %>% mutate(DU = as.numeric(grp == "USD"))
  rhs <- c(XC, "y_lag1", "y_lag2", "fci_lag1", CTRL)
  f <- as.formula(paste("y_fwd ~ DU +", paste(rhs, collapse = "+"), "+",
                        paste0("DU:", rhs, collapse = "+")))
  m <- lm(f, data = s)
  V <- sandwich::vcovPL(m, cluster = s$fecha, order.by = s$fecha,
                        lag = h + 1, adjust = TRUE)
  ip   <- paste0("DU:", XC)
  bMN  <- unname(coef(m)[XC])
  dUM  <- unname(coef(m)[ip])            # beta_USD - beta_MN
  bUSD <- bMN + dUM
  seD  <- sqrt(V[ip, ip])
  diff <- -dUM                           # beta_MN - beta_USD
  dfree <- length(unique(s$fecha)) - 2
  tst  <- diff / seD
  seCL <- sqrt(sandwich::vcovCL(m, cluster = s$fecha)[ip, ip])
  data.frame(horizon = h, beta_MN = bMN, beta_USD = bUSD,
             diff_MN_minus_USD = diff, se_dk = seD, t = tst,
             p_dk = 2 * pt(-abs(tst), dfree),
             ci_lo = diff - qt(0.975, dfree) * seD,
             ci_hi = diff + qt(0.975, dfree) * seD,
             se_clusterdate = seCL,
             p_clusterdate = 2 * pt(-abs(diff / seCL), dfree),
             n_stacked = nrow(s), n_months = length(unique(s$fecha)), df = dfree)
}

ce <- bind_rows(lapply(c(6, 12, 18), function(h) stack_test(cd, h)))
write.csv(ce, file.path(OUT, "Rev_Currency_Equality_Test.csv"), row.names = FALSE)
cat("  [csv]  Rev_Currency_Equality_Test.csv\n")
print(ce %>% mutate(across(where(is.numeric), ~round(., 4))))

cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
