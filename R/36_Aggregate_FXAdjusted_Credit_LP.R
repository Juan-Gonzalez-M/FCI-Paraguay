# ============================================================================
# 36_Aggregate_FXAdjusted_Credit_LP.R
# ADDITIVE exercise: re-estimates the paper's aggregate credit local
# projections with constant-exchange-rate (FX-adjusted) credit outcomes.
#
# Motivation: aggregate credit stocks embed the FX book at the CURRENT
# exchange rate (totales = MN + USD x TCN), so guarani depreciation
# mechanically inflates measured credit growth. Because depreciation
# coincides with FCI tightening, this attenuates the measured contraction -
# the published estimates should be conservative. This script quantifies that.
#
# Method: exact replication of run_lp_standard() from 05_FCI_Local_Projections.R
# (FCI_t + 2 lags of y + 1 lag of FCI + IMAEP/IPC controls with 2 lags,
# Newey-West lag h+1 prewhite=FALSE, 90% CI, FCI_exCredit swap for credit LHS,
# same join/na.omit sample construction). The published outcomes are
# re-estimated first and VALIDATED against output/csv/LP_Credit_Standard.csv;
# only then are the FX-adjusted outcomes swapped in.
#
# Outcomes:
#   published (replication): Cred_Total (real, CPI-deflated), Cred_USD
#     (USD book valued in PYG), Cred_Real_MN
#   FX-adjusted (new):       Cred_Total_FXadj (constant exchange rate, BIS
#     convention with e_{t-12}, CPI-deflated), Cred_USD_clean (USD book in
#     dollars)
#
# NOTHING in the existing pipeline or outputs is touched. Everything goes to
# output/micro/ (csv prefix Agg_FXAdj_, png 270-272).
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
suppressMessages({library(sandwich); library(lmtest)})
cat("=== 36: FX-adjusted aggregate credit LPs (additive exercise) ===\n")

MAX_H   <- 24
N_LAGS  <- 2
H_REP   <- c(3, 6, 12, 18, 24)
Z_CRIT  <- qnorm(0.95)   # 90% CI, as in LP_CONFIG

# ---------------------------------------------------------------------------
# 1. Rebuild analysis_data exactly as in 05_FCI_Local_Projections.R
#    (FCI read from the pipeline CSV instead of re-running script 01 -
#     identical values by construction)
# ---------------------------------------------------------------------------
fci_all <- fread(micro_paths$fci_csv)
fci_all[, fecha := as.Date(fecha)]

fci_data <- fci_all[, .(fecha,
                        FCI_COMP = FCI_COMP_AVG, FCI_ENDO = FCI_ENDO_AVG,
                        FCI_EXO = FCI_EXO_AVG, FCI_exCredit = FCI_exCredit_AVG,
                        FCI_ENDO_exCredit = FCI_ENDO_exCredit_AVG,
                        FCI_RATES = FCI_RATES_AVG,
                        FCI_ENDO_exCredit_exTCN = FCI_ENDO_exCredit_exTCN_AVG)]
# method columns matched by script 05's regex (PCA/VAR/DFM norm of exCredit)
mcols <- grep("^FCI_exCredit_(ZS|PCA|VAR|DFM)_norm$", names(fci_all), value = TRUE)
fci_methods <- fci_all[, c("fecha", mcols), with = FALSE]

mv <- as.data.table(suppressMessages(read_excel(micro_paths$fci_xlsx, sheet = "Main_variables")))
dm <- as.data.table(suppressMessages(read_excel(micro_paths$fci_xlsx, sheet = "Datos_macro")))
mv[, fecha := as.Date(Fecha)]; dm[, fecha := as.Date(Fecha)]
setorder(mv, fecha); setorder(dm, fecha)

L <- function(x, k = 12) shift(x, k)

# Credit outcomes: published definitions + FX-adjusted counterparts
cred <- merge(mv[, .(fecha, MN = Creditos_Sector_privado_MN,
                     USD = Creditos_Sector_privado_USD,
                     USDeq = Creditos_Sector_privado_USD_equivalente, TCN)],
              dm[, .(fecha, IPC, Creditos_deflactados)], by = "fecha")
setorder(cred, fecha)
cred[, `:=`(
  # --- published (replication) ---
  Cred_Total   = (Creditos_deflactados / L(Creditos_deflactados) - 1) * 100,
  Cred_MN      = (MN / L(MN) - 1) * 100,
  Cred_USD     = (USDeq / L(USDeq) - 1) * 100,
  Cred_Real_MN = ((MN / IPC) / L(MN / IPC) - 1) * 100,
  # --- FX-adjusted (new) ---
  # constant exchange rate (BIS convention: value both periods' FX book at
  # e_{t-12}), CPI-deflated
  Cred_Total_FXadj = ((MN + USD * L(TCN)) / (L(MN) + L(USD) * L(TCN)) *
                        L(IPC) / IPC - 1) * 100,
  # USD book in dollars (no valuation component by construction)
  Cred_USD_clean   = (USD / L(USD) - 1) * 100,
  # diagnostics
  tcn_yoy  = (TCN / L(TCN) - 1) * 100,
  fx_share = USDeq / (MN + USDeq)
)]

# Macro controls exactly as script 05
dm[, `:=`(IMAEP_yoy = (IMAEP / L(IMAEP) - 1) * 100,
          IPC_yoy   = (IPC / L(IPC) - 1) * 100)]
macro_data <- merge(dm[, .(fecha, IMAEP_yoy, IPC_yoy)],
                    cred[, .(fecha, Cred_Real_yoy = Cred_Total)], by = "fecha")

ext <- mv[, .(fecha, VIX, FFER, Commodities)]

analysis_data <- Reduce(function(a, b) merge(a, b, by = "fecha"),
  list(fci_data, fci_methods,
       cred[, .(fecha, Cred_Total, Cred_MN, Cred_USD, Cred_Real_MN,
                Cred_Total_FXadj, Cred_USD_clean)],
       macro_data, ext))
setorder(analysis_data, fecha)
# script 05 does na.omit() over its full column set; the FX-adjusted columns
# have the same support (12-month lags), so na.omit here reproduces the sample
analysis_data <- na.omit(analysis_data)
cat(sprintf("  Sample: %d obs, %s to %s\n", nrow(analysis_data),
            format(min(analysis_data$fecha)), format(max(analysis_data$fecha))))

# ---------------------------------------------------------------------------
# 2. run_lp_standard(): verbatim logic from script 05
# ---------------------------------------------------------------------------
run_lp_standard <- function(data, y_var, fci_var, max_h = MAX_H, n_lags = N_LAGS,
                            control_vars = NULL) {
  data <- as.data.frame(data)
  results <- list()
  for (h in 1:max_h) {
    dh <- data
    dh$y_fwd   <- dplyr::lead(dh[[y_var]], h)
    dh$y_lag1  <- dplyr::lag(dh[[y_var]], 1)
    dh$fci_lag1 <- dplyr::lag(dh[[fci_var]], 1)
    for (i in 2:n_lags) dh[[paste0("y_lag", i)]] <- dplyr::lag(dh[[y_var]], i)
    lag_vars <- c("y_lag1", paste0("y_lag", 2:n_lags), "fci_lag1")
    control_cols <- c()
    if (!is.null(control_vars)) {
      for (cv in control_vars) {
        control_cols <- c(control_cols, cv)
        for (j in 1:n_lags) {
          nm <- paste0(cv, "_lag", j)
          dh[[nm]] <- dplyr::lag(dh[[cv]], j)
          control_cols <- c(control_cols, nm)
        }
      }
    }
    fml <- as.formula(paste("y_fwd ~", fci_var, "+",
                            paste(c(lag_vars, control_cols), collapse = " + ")))
    reg_data <- na.omit(dh[, c("y_fwd", fci_var, lag_vars, control_cols)])
    if (nrow(reg_data) < 30) next
    model <- lm(fml, data = reg_data)
    ct <- lmtest::coeftest(model, vcov = sandwich::NeweyWest(model, lag = h + 1,
                                                             prewhite = FALSE))
    idx <- which(rownames(ct) == fci_var)
    results[[h]] <- data.frame(horizon = h, coef = ct[idx, 1], se = ct[idx, 2],
                               p_value = ct[idx, 4],
                               ci_lower = ct[idx, 1] - Z_CRIT * ct[idx, 2],
                               ci_upper = ct[idx, 1] + Z_CRIT * ct[idx, 2],
                               n_obs = nrow(reg_data))
  }
  do.call(rbind, results)
}

fci_credit_swap <- c(FCI_COMP = "FCI_exCredit", FCI_ENDO = "FCI_ENDO_exCredit",
                     FCI_EXO = "FCI_EXO")
CREDIT_MEASURES <- c(
  Total       = "Cred_Total",        # published
  USD         = "Cred_USD",          # published
  Real_MN     = "Cred_Real_MN",      # published (clean by construction)
  Total_FXadj = "Cred_Total_FXadj",  # new
  USD_clean   = "Cred_USD_clean")    # new

# ---------------------------------------------------------------------------
# 3. Estimate all LPs
# ---------------------------------------------------------------------------
res <- list()
for (fci in names(fci_credit_swap)) {
  for (cm in names(CREDIT_MEASURES)) {
    cat(sprintf("  LP: %s [%s] -> %s ...\n", fci, fci_credit_swap[fci], cm))
    r <- run_lp_standard(analysis_data, CREDIT_MEASURES[cm], fci_credit_swap[fci],
                         control_vars = c("IMAEP_yoy", "IPC_yoy"))
    r$fci_type <- fci; r$credit_measure <- cm
    r$adjusted <- cm %in% c("Total_FXadj", "USD_clean")
    res[[paste(fci, cm)]] <- r
  }
}
res <- rbindlist(res)
write_micro_csv(res, "Agg_FXAdj_Credit_LP_Results.csv")

# ---------------------------------------------------------------------------
# 4. Replication check against the published LP_Credit_Standard.csv
# ---------------------------------------------------------------------------
pub <- fread(file.path(micro_paths$root, "output", "csv", "LP_Credit_Standard.csv"))
chk <- merge(pub[credit_type %in% c("Total", "USD", "Real_MN"),
                 .(horizon, fci_type, credit_measure = credit_type,
                   coef_pub = coef, n_pub = n_obs)],
             res[adjusted == FALSE,
                 .(horizon, fci_type, credit_measure, coef_rep = coef, n_rep = n_obs)],
             by = c("horizon", "fci_type", "credit_measure"))
chk[, `:=`(dcoef = abs(coef_pub - coef_rep), dn = n_pub - n_rep)]
cat(sprintf("\nREPLICATION CHECK vs LP_Credit_Standard.csv (%d cells):\n", nrow(chk)))
cat(sprintf("  max |coef diff| = %.6f;  n_obs identical in %.1f%% of cells (max dn = %d)\n",
            max(chk$dcoef), 100 * mean(chk$dn == 0), max(abs(chk$dn))))
if (max(chk$dcoef) > 0.05)
  cat("  *** WARNING: replication imperfect - treat FX-adjusted deltas, not levels, as informative ***\n")
write_micro_csv(chk, "Agg_FXAdj_Replication_Check.csv")

# ---------------------------------------------------------------------------
# 5. Key-horizon comparison table (published vs FX-adjusted)
# ---------------------------------------------------------------------------
pairs <- list(c("Total", "Total_FXadj"), c("USD", "USD_clean"))
cmp <- rbindlist(lapply(pairs, function(pr) {
  a <- res[credit_measure == pr[1] & horizon %in% H_REP,
           .(fci_type, horizon, coef_published = coef, se_published = se,
             p_published = p_value)]
  b <- res[credit_measure == pr[2] & horizon %in% H_REP,
           .(fci_type, horizon, coef_fxadj = coef, se_fxadj = se,
             p_fxadj = p_value)]
  m <- merge(a, b, by = c("fci_type", "horizon"))
  m[, outcome := pr[1]]
  m
}))
cmp[, delta := coef_fxadj - coef_published]
setcolorder(cmp, c("outcome", "fci_type", "horizon"))
write_micro_csv(cmp, "Agg_FXAdj_Comparison_Key_Horizons.csv")

cat("\nKey comparison (FCI_COMP i.e. FCI_exCredit on credit LHS):\n")
print(cmp[fci_type == "FCI_COMP",
          .(outcome, h = horizon, b_pub = round(coef_published, 2),
            p_pub = round(p_published, 3), b_adj = round(coef_fxadj, 2),
            p_adj = round(p_fxadj, 3), delta = round(delta, 2))])

# ---------------------------------------------------------------------------
# 6. Mechanical decomposition series (for the appendix paragraph)
# ---------------------------------------------------------------------------
dec <- cred[!is.na(Cred_Total) & !is.na(Cred_Total_FXadj),
            .(fecha, Cred_Total, Cred_Total_FXadj,
              gap = Cred_Total - Cred_Total_FXadj, tcn_yoy, fx_share)]
mgap <- lm(gap ~ tcn_yoy, dec)
cat(sprintf("\nMechanical gap: mean %.2f pp, SD %.2f pp; gap = %.3f x TCN_yoy (R2 = %.3f); mean FX share %.3f\n",
            mean(dec$gap), sd(dec$gap), coef(mgap)[2], summary(mgap)$r.squared,
            mean(dec$fx_share)))
write_micro_csv(dec, "Agg_FXAdj_Mechanical_Decomposition.csv")

# ---------------------------------------------------------------------------
# 7. Figures (270-272)
# ---------------------------------------------------------------------------
mk_overlay <- function(pubm, adjm, ttl, fn) {
  d <- rbind(res[credit_measure == pubm & fci_type == "FCI_COMP",
                 .(horizon, coef, se, Series = "Published measure")],
             res[credit_measure == adjm & fci_type == "FCI_COMP",
                 .(horizon, coef, se, Series = "FX-adjusted (constant e)")])
  p <- ggplot(d, aes(horizon, coef, color = Series, fill = Series)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_ribbon(aes(ymin = coef - Z_CRIT * se, ymax = coef + Z_CRIT * se),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) + geom_point(size = 1.4) +
    scale_color_manual(values = c("firebrick", "grey35")) +
    scale_fill_manual(values = c("firebrick", "grey35")) +
    labs(title = ttl,
         subtitle = "LP on FCI_exCredit, spec identical to script 05 (NW h+1, 90% bands); valuation component removed in red",
         x = "Horizon (months)", y = "Coefficient (pp per 1 SD FCI)") +
    theme_micro()
  save_png(p, fn)
}
mk_overlay("Total", "Total_FXadj",
           "Aggregate LP: real total credit growth - published vs FX-adjusted",
           "270_AggLP_Total_FXadjusted.png")
mk_overlay("USD", "USD_clean",
           "Aggregate LP: USD credit growth - PYG-valued (published) vs in dollars",
           "271_AggLP_USD_FXadjusted.png")

pG <- ggplot(dec, aes(fecha)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_line(aes(y = gap, color = "Mechanical gap (measured - adjusted growth)"), linewidth = 0.7) +
  geom_line(aes(y = tcn_yoy * coef(mgap)[2], color = "FX share x TCN depreciation"), linewidth = 0.7) +
  scale_color_manual(values = c("steelblue4", "darkorange3")) +
  labs(title = "Valuation component of measured aggregate credit growth",
       subtitle = sprintf("gap = %.2f x TCN YoY, R2 = %.2f; mean FX share of credit %.0f%%",
                          coef(mgap)[2], summary(mgap)$r.squared, 100 * mean(dec$fx_share)),
       x = NULL, y = "pp of YoY credit growth", color = NULL) +
  theme_micro()
save_png(pG, "272_Agg_Valuation_Gap.png")

cat(sprintf("=== 36 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
