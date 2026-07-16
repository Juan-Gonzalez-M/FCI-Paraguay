# ============================================================================
# 52_OneSided_FCI_Cobaseline.R
# Strictly one-sided (timing-clean) FCI co-baseline.
#
# The composite FCI averages four methods whose parameters (PCA loadings,
# VAR coefficients, smoothed DFM states) and normalizations use the full
# sample (Appendix A.5). The z-score method's INPUTS are strictly
# backward-looking (trailing 60-month standardization, script 01); its only
# full-sample step is the affine scale() rescaling, which cannot change LP
# t-statistics. This script promotes the raw (unscaled) z-score index as a
# declared timing-clean co-baseline and runs it through:
#
#   Part 1: construction checks -- correlation of raw vs. scale()d z-score
#           (must be 1 by affinity), and with the composite AVG.
#   Part 2: headline credit LP at h = 6, 12, 18 (script 49 conventions:
#           NW(h+1), lagged controls), full sample and post-IT, for
#           (a) one-sided z-score exCredit (raw),
#           (b) one-sided z-score comprehensive (raw),
#           (c) expanding-window recursively normalized index (script 39),
#           against the composite FCI_exCredit_AVG benchmark.
#   Part 3: DXY first stage for the one-sided index (levels), with
#           separately labeled diagnostics: classical partial R2, classical
#           conditional F, HAC(NW) t; plus 2SLS = RF/FS and AR p at 0.
#   Part 4: purging-window check -- script 23's globally-purged FCI is a
#           FULL-SAMPLE regression (23:528). Recompute the purge residuals
#           RECURSIVELY (expanding window, 60-obs burn-in) and compare the
#           h = 12 credit LP on both versions.
#
# Outputs: output/revision/csv/Rev_OneSided_LP.csv, Rev_OneSided_FirstStage.csv,
#          Rev_OneSided_PurgeCheck.csv, Rev_OneSided_Construction.csv
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(sandwich); library(lmtest) })

d <- read.csv("../output/csv/New_External_Variables.csv")
d$fecha <- as.Date(d$fecha)

comp <- read.csv("../output/csv/FCI_Complete_Results.csv")
comp$fecha <- as.Date(comp$fecha)
rv <- read.csv("../output/csv/FCI_Robustness_Versions.csv")
rv$fecha <- as.Date(rv$fecha)
exp_fci <- read.csv("../output/revision/csv/Rev_ExpandingFCI_Series.csv")
exp_fci$fecha <- as.Date(exp_fci$fecha)

d <- d %>%
  left_join(comp %>% dplyr::select(fecha, FCI_COMP_ZSCORE, FCI_COMP_ZSCORE_norm,
                                   FCI_ENDO_ZSCORE), by = "fecha") %>%
  left_join(rv %>% dplyr::select(fecha, FCI_exCredit_ZSCORE, FCI_exCredit_ZSCORE_norm),
            by = "fecha") %>%
  left_join(exp_fci %>% dplyr::select(fecha, FCI_EXP, FCI_EXP_exCredit), by = "fecha") %>%
  arrange(fecha) %>%
  mutate(IMAEP_yoy_L1 = lag(IMAEP_yoy, 1), IMAEP_yoy_L2 = lag(IMAEP_yoy, 2),
         IPC_yoy_L1   = lag(IPC_yoy, 1),   IPC_yoy_L2   = lag(IPC_yoy, 2))

POSTIT <- as.Date("2011-05-01")

cat("=== 52: strictly one-sided FCI co-baseline ===\n\n")

# ---------------------------------------------------------------------------
# Part 1: construction checks
# ---------------------------------------------------------------------------
cat("--- Part 1: construction checks ---\n")
cons <- data.frame(
  pair = c("ZSCORE raw vs scale()d (exCredit)",
           "ZSCORE raw vs scale()d (comprehensive)",
           "one-sided exCredit ZSCORE vs composite FCI_exCredit_AVG",
           "one-sided comprehensive ZSCORE vs composite FCI_COMP_AVG",
           "one-sided exCredit ZSCORE vs expanding-z exCredit (FCI_EXP_exCredit)"),
  correlation = c(
    cor(d$FCI_exCredit_ZSCORE, d$FCI_exCredit_ZSCORE_norm, use = "complete.obs"),
    cor(d$FCI_COMP_ZSCORE, d$FCI_COMP_ZSCORE_norm, use = "complete.obs"),
    cor(d$FCI_exCredit_ZSCORE, d$FCI_exCredit_AVG, use = "complete.obs"),
    cor(d$FCI_COMP_ZSCORE, d$FCI_COMP_AVG, use = "complete.obs"),
    cor(d$FCI_exCredit_ZSCORE, d$FCI_EXP_exCredit, use = "complete.obs")))
print(cons, row.names = FALSE)
write.csv(cons, "../output/revision/csv/Rev_OneSided_Construction.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Part 2: headline credit LP
# ---------------------------------------------------------------------------
lp_h <- function(data, xvar, h, lagctrl = TRUE) {
  dh <- data %>% mutate(y_fwd = lead(Cred_Real_Total, h),
                        y_lag1 = lag(Cred_Real_Total, 1),
                        y_lag2 = lag(Cred_Real_Total, 2),
                        x_lag1 = lag(.data[[xvar]], 1))
  ctr <- if (lagctrl) c("IMAEP_yoy","IMAEP_yoy_L1","IMAEP_yoy_L2",
                        "IPC_yoy","IPC_yoy_L1","IPC_yoy_L2")
         else c("IMAEP_yoy","IPC_yoy")
  rhs <- c(xvar, "y_lag1", "y_lag2", "x_lag1", ctr)
  rd <- dh %>% dplyr::select(y_fwd, all_of(rhs)) %>% na.omit()
  if (nrow(rd) < 50) return(rep(NA_real_, 5))
  m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
  ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE))
  c(b = unname(ct[xvar, 1]), se = unname(ct[xvar, 2]), p = unname(ct[xvar, 4]),
    n = nrow(rd), sd_x = sd(rd[[xvar]]))
}

cat("\n--- Part 2: headline credit LP (NW h+1, lagged controls) ---\n")
idx <- c(one_sided_exCredit      = "FCI_exCredit_ZSCORE",
         one_sided_comprehensive = "FCI_COMP_ZSCORE",
         expanding_z_exCredit    = "FCI_EXP_exCredit",
         composite_exCredit_AVG  = "FCI_exCredit_AVG")
lp <- data.frame()
for (smp in c("full", "postIT")) {
  ds <- if (smp == "postIT") d %>% filter(fecha >= POSTIT) else d
  for (v in names(idx)) {
    for (h in c(6, 12, 18)) {
      r <- lp_h(ds, unname(idx[v]), h)
      lp <- rbind(lp, data.frame(sample = smp, index = v, column = idx[v], horizon = h,
        b = r[1], se = r[2], p = r[3], n = r[4], sd_x = r[5],
        b_per_sd = r[1] * r[5]))
    }
  }
}
for (smp in c("full", "postIT")) {
  cat(sprintf("\n  [%s sample]\n", smp))
  sub <- lp[lp$sample == smp, ]
  for (v in names(idx)) {
    s <- sub[sub$index == v, ]
    cat(sprintf("  %-26s h6: %6.2f (p=%.3f)  h12: %6.2f (p=%.3f)  h18: %6.2f (p=%.3f)  n(h12)=%d\n",
        v, s$b[s$horizon==6], s$p[s$horizon==6], s$b[s$horizon==12], s$p[s$horizon==12],
        s$b[s$horizon==18], s$p[s$horizon==18], s$n[s$horizon==12]))
  }
}
write.csv(lp, "../output/revision/csv/Rev_OneSided_LP.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Part 3: DXY first stage for the one-sided index (script 49 control set)
# ---------------------------------------------------------------------------
cat("\n--- Part 3: DXY first stage, one-sided index ---\n")
er  <- intersect(c("d_ToT", "US_10Y", "SP500"), names(d))
ex0 <- c("IMAEP_yoy","IMAEP_yoy_L1","IMAEP_yoy_L2","IPC_yoy","IPC_yoy_L1","IPC_yoy_L2", er)

fs_diag <- function(data, yvar, h = 12) {
  dh <- data %>% mutate(y_fwd = lead(Cred_Real_Total, h),
                        y_lag1 = lag(Cred_Real_Total, 1),
                        y_lag2 = lag(Cred_Real_Total, 2),
                        fci_lag1 = lag(.data[[yvar]], 1))
  ex <- c("y_lag1","y_lag2","fci_lag1", ex0)
  rd <- dh %>% dplyr::select(y_fwd, all_of(yvar), all_of(ex), DXY) %>% na.omit()
  xs <- paste(ex, collapse = " + ")
  fs   <- lm(as.formula(paste(yvar, "~ DXY +", xs)), data = rd)
  fs_r <- lm(as.formula(paste(yvar, "~", xs)), data = rd)
  f_cls  <- anova(fs_r, fs)$F[2]
  t_cls  <- summary(fs)$coefficients["DXY", 3]
  # classical partial R2 = t^2 / (t^2 + df)
  df_res <- df.residual(fs)
  pr2_cls <- t_cls^2 / (t_cls^2 + df_res)
  ct_hac <- lmtest::coeftest(fs, vcov = sandwich::NeweyWest(fs, lag = h + 1, prewhite = FALSE))
  t_hac  <- ct_hac["DXY", 3]
  rf <- lm(as.formula(paste("y_fwd ~ DXY +", xs)), data = rd)
  b_iv <- coef(rf)["DXY"] / coef(fs)["DXY"]
  # AR p at 0 (classical)
  ar_r <- lm(as.formula(paste("y_fwd ~", xs)), data = rd)
  p_ar0 <- anova(ar_r, rf)$`Pr(>F)`[2]
  data.frame(index = yvar, horizon = h, n = nrow(rd),
             b_fs = coef(fs)["DXY"], t_classical = t_cls, F_classical = f_cls,
             partialR2_classical = pr2_cls, t_hac_nw = t_hac, F_hac_nw = t_hac^2,
             b_2sls = b_iv, ar_p_at_zero = p_ar0)
}
fs_out <- data.frame()
for (v in c("FCI_exCredit_ZSCORE", "FCI_COMP_ZSCORE", "FCI_exCredit_AVG")) {
  for (h in c(6, 12, 18)) fs_out <- rbind(fs_out, fs_diag(d, v, h))
}
print(fs_out %>% mutate(across(where(is.numeric), ~round(., 4))), row.names = FALSE)
write.csv(fs_out, "../output/revision/csv/Rev_OneSided_FirstStage.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Part 4: purging-window check (full-sample vs recursive purge, script 23 spec)
# ---------------------------------------------------------------------------
cat("\n--- Part 4: globally-purged FCI, full-sample vs recursive purge ---\n")
global_vars <- intersect(c("VIX", "DXY", "FFER", "US_10Y", "SP500", "Selic_rate"), names(d))
pd <- d
pvars <- c()
for (gv in global_vars) {
  pd[[paste0(gv, "_L1")]] <- lag(pd[[gv]], 1)
  pd[[paste0(gv, "_L2")]] <- lag(pd[[gv]], 2)
  pvars <- c(pvars, gv, paste0(gv, "_L1"), paste0(gv, "_L2"))
}
pf <- as.formula(paste("FCI_exCredit_AVG ~", paste(pvars, collapse = " + ")))
ok <- complete.cases(pd[, c("FCI_exCredit_AVG", pvars)])

# full-sample purge (replicates script 23)
pd$purged_full <- NA_real_
pd$purged_full[ok] <- residuals(lm(pf, data = pd[ok, ]))

# recursive purge: expanding window, min 60 complete obs, out-of-sample residual
idx_ok <- which(ok)
pd$purged_recursive <- NA_real_
for (i in seq_along(idx_ok)) {
  if (i <= 60) next
  train <- pd[idx_ok[1:(i - 1)], ]
  m <- lm(pf, data = train)
  row_i <- pd[idx_ok[i], , drop = FALSE]
  pd$purged_recursive[idx_ok[i]] <-
    row_i$FCI_exCredit_AVG - as.numeric(predict(m, newdata = row_i))
}
cat(sprintf("  corr(full-sample purge, recursive purge) = %.3f  (n=%d)\n",
    cor(pd$purged_full, pd$purged_recursive, use = "complete.obs"),
    sum(complete.cases(pd[, c("purged_full", "purged_recursive")]))))

pc <- data.frame()
for (v in c("purged_full", "purged_recursive")) {
  for (h in c(6, 12, 18)) {
    r <- lp_h(pd, v, h)
    pc <- rbind(pc, data.frame(purge = v, horizon = h,
      b = r[1], se = r[2], p = r[3], n = r[4], sd_x = r[5], b_per_sd = r[1] * r[5]))
  }
}
print(pc %>% mutate(across(where(is.numeric), ~round(., 4))), row.names = FALSE)
write.csv(pc, "../output/revision/csv/Rev_OneSided_PurgeCheck.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Gate verdict
# ---------------------------------------------------------------------------
g <- lp %>% filter(sample == "full", index == "one_sided_exCredit", horizon == 12)
cat(sprintf("\n=== GATE: one-sided exCredit z-score, full sample, h=12: b=%.2f (p=%.4f) ===\n",
    g$b, g$p))
cat(ifelse(!is.na(g$p) && g$p < 0.10 && g$b < 0,
    "GATE PASSED: timing-clean index confirms the headline result.\n",
    "GATE FAILED: timing-clean index does NOT confirm -- stop and reassess.\n"))

# ---------------------------------------------------------------------------
# Part 5 (round 9): recursive re-estimation of the BANK-LEVEL purged-DXY
# shock (script 30's shk_dxy_purged). The projection (AR(3) of the monthly
# log-DXY change on its own lags plus contemporaneous VIX, IPE, ToT, FFER,
# SP500, US10Y changes) is estimated on the full macro history; here it is
# re-estimated on expanding windows (minimum 60 complete observations,
# out-of-sample residual at each t) and the recursive innovation is
# correlated with the full-sample innovation over the 2016-2025 shock window.
# ---------------------------------------------------------------------------
cat("\n--- Part 5: recursive bank-shock purge check ---\n")
suppressPackageStartupMessages(library(readxl))
mv2  <- read_excel("../data/FCI_data_1.xlsx", sheet = "Main_variables")
gfc2 <- read_excel("../data/FCI_data_1.xlsx", sheet = "Global_Financial_Conditions")
names(gfc2)[1] <- "Fecha"; names(gfc2)[names(gfc2) == "S&P 500"] <- "SP500"
m5 <- merge(
  data.frame(fecha = as.Date(mv2$Fecha), VIX = mv2$VIX, FFER = mv2$FFER),
  data.frame(fecha = as.Date(gfc2$Fecha), DXY = gfc2$DXY, SP500 = gfc2$SP500,
             US10Y = gfc2$US_10Y, IPE = gfc2$IPE, ToT = gfc2$Paraguay_ToT),
  by = "fecha", all = TRUE)
m5 <- m5[order(m5$fecha), ]
dl <- function(x) 100 * c(NA, diff(log(x)))
df5 <- data.frame(fecha = m5$fecha,
  y = dl(m5$DXY), d_vix = c(NA, diff(m5$VIX)), dln_ipe = dl(m5$IPE),
  dln_tot = dl(m5$ToT), d_ffer = c(NA, diff(m5$FFER)),
  dln_sp500 = dl(m5$SP500), d_us10y = c(NA, diff(m5$US10Y)))
df5$l1 <- dplyr::lag(df5$y, 1); df5$l2 <- dplyr::lag(df5$y, 2); df5$l3 <- dplyr::lag(df5$y, 3)
ok5 <- complete.cases(df5[, -1])
f5 <- y ~ l1 + l2 + l3 + d_vix + dln_ipe + dln_tot + d_ffer + dln_sp500 + d_us10y
# full-sample innovation (replicates script 30's shk_dxy_purged up to scale)
df5$innov_full <- NA_real_
df5$innov_full[ok5] <- resid(lm(f5, data = df5[ok5, ]))
# recursive innovation: expanding window, min 60 obs, out-of-sample residual
idx5 <- which(ok5)
df5$innov_rec <- NA_real_
for (i in seq_along(idx5)) {
  if (i <= 60) next
  mtr <- lm(f5, data = df5[idx5[1:(i - 1)], ])
  df5$innov_rec[idx5[i]] <- df5$y[idx5[i]] -
    as.numeric(predict(mtr, newdata = df5[idx5[i], , drop = FALSE]))
}
win5 <- df5$fecha >= as.Date("2016-01-01") & df5$fecha <= as.Date("2025-12-01")
cc5 <- df5[win5 & complete.cases(df5[, c("innov_full", "innov_rec")]), ]
cat(sprintf("  corr(full-sample innovation, recursive innovation), 2016-2025: %.4f (n=%d)\n",
    cor(cc5$innov_full, cc5$innov_rec), nrow(cc5)))
write.csv(data.frame(n = nrow(cc5),
  correlation = cor(cc5$innov_full, cc5$innov_rec),
  sd_full = sd(cc5$innov_full), sd_rec = sd(cc5$innov_rec)),
  "../output/revision/csv/Rev_BankShock_RecursivePurge.csv", row.names = FALSE)
cat("Saved: Rev_BankShock_RecursivePurge.csv\n")
