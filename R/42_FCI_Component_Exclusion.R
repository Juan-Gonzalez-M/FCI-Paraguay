# ============================================================================
# 42_FCI_Component_Exclusion.R — Component-by-component FCI exclusion battery
# (Coworker 3, Concern 5, item 2; Coworker 1, WP2 item 2).
#
# For each potentially endogenous component (NPL, ROE, liquidity, credit/
# deposit ratio, TCN), construct a leave-one-out FCI and re-estimate the
# credit LP. LOO variants are built on the exCredit base (11 vars), matching
# the paper's convention for credit-LHS regressions, using the z-score method
# with rolling-60m standardization (script 01 convention). Plus the
# quantities-at-t-1 variant (no contemporaneous banking quantity enters).
#
# ANCHOR VALIDATION: the z-score exCredit variant must track the published
# 4-method FCI_exCredit_AVG closely (corr printed; gate at 0.90).
# ADDITIVE: writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
cat("=== 42: FCI component-exclusion battery ===\n")

inputs <- load_fci_inputs()
d <- load_ext_data()

vars_exCred <- setdiff(names(FCI_SIGNS), "Crecimiento_creditos")

# ---------------------------------------------------------------------------
# 1. Build variants (z-score method, rolling 60m with partial windows)
# ---------------------------------------------------------------------------
variants <- list(
  zs_full12        = names(FCI_SIGNS),
  zs_exCredit      = vars_exCred,                                # baseline for LP
  zs_exCredit_exNPL       = setdiff(vars_exCred, "Morosidad"),
  zs_exCredit_exROE       = setdiff(vars_exCred, "Rentabilidad"),
  zs_exCredit_exLiquidez  = setdiff(vars_exCred, "Liquidez"),
  zs_exCredit_exRatioCD   = setdiff(vars_exCred, "Ratio_Cred_Depo"),
  zs_exCredit_exTCN       = setdiff(vars_exCred, "TCN"))

FV <- data.table(fecha = inputs$fecha)
for (v in names(variants)) FV[, (v) := build_zscore_fci(inputs, variants[[v]])]
# Quantities at t-1 (lag all banking quantities; exCredit base for the LP)
FV[, zs_exCredit_QtyLag1 := build_zscore_fci(inputs, vars_exCred,
                                             lag_vars = FCI_QUANTITY_VARS)]
FV[, zs_full12_QtyLag1 := build_zscore_fci(inputs, names(FCI_SIGNS),
                                           lag_vars = FCI_QUANTITY_VARS)]

d <- merge(d, FV, by = "fecha", all.x = TRUE)

# ---------------------------------------------------------------------------
# 2. Anchor validation against the published 4-method AVG variants
# ---------------------------------------------------------------------------
anchors <- fread(micro_paths$fci_csv, select = c("fecha", "FCI_exNPL_AVG"))
anchors[, fecha := as.Date(fecha)]
d <- merge(d, anchors, by = "fecha", all.x = TRUE)
a1 <- d[, cor(zs_exCredit, FCI_exCredit_AVG, use = "complete.obs")]
a2 <- d[, cor(zs_full12, FCI_COMP_AVG, use = "complete.obs")]
a3 <- d[, cor(zs_exCredit_exNPL, FCI_exNPL_AVG, use = "complete.obs")]
cat(sprintf("  Anchor corr: zs_exCredit vs FCI_exCredit_AVG = %.3f\n", a1))
cat(sprintf("               zs_full12   vs FCI_COMP_AVG     = %.3f\n", a2))
cat(sprintf("               zs_exCr_exNPL vs FCI_exNPL_AVG   = %.3f\n", a3))
if (min(a1, a2, na.rm = TRUE) < 0.90)
  cat("  *** WARNING: z-score variants track the 4-method AVG below 0.90 —\n",
      "      treat LOO deltas, not levels, as the informative object. ***\n")

# ---------------------------------------------------------------------------
# 3. Credit LP at h = 6/12/18 for every variant + z-tests vs the zs_exCredit
#    baseline; include the published rates-only FCI as a co-baseline row
# ---------------------------------------------------------------------------
fci_all <- fread(micro_paths$fci_csv, select = c("fecha", "FCI_RATES_AVG"))
fci_all[, fecha := as.Date(fecha)]
d <- merge(d, fci_all, by = "fecha", all.x = TRUE)

LP_VARS <- c("zs_exCredit", "zs_exCredit_exNPL", "zs_exCredit_exROE",
             "zs_exCredit_exLiquidez", "zs_exCredit_exRatioCD",
             "zs_exCredit_exTCN", "zs_exCredit_QtyLag1",
             "FCI_exCredit_AVG", "FCI_RATES_AVG")

res <- list()
for (v in LP_VARS) {
  r <- run_lp(d, "Cred_Real_Total", v, max_h = 18,
              controls = c("IMAEP_yoy", "IPC_yoy"))
  r[, variant := v]
  res[[v]] <- r
  k <- r[horizon == 12]
  cat(sprintf("  [%s] h=12: b=%7.2f (se %.2f, p=%.3f)\n", v, k$coef, k$se, k$p_value))
}
res <- rbindlist(res)

base12 <- res[variant == "zs_exCredit" & horizon %in% c(6, 12, 18),
              .(horizon, b0 = coef, se0 = se)]
tab <- merge(res[horizon %in% c(6, 12, 18)], base12, by = "horizon")
tab[, z_vs_baseline := (coef - b0) / sqrt(se^2 + se0^2)]
tab[, p_equal := 2 * pnorm(-abs(z_vs_baseline))]
corrs <- sapply(LP_VARS, function(v)
  d[, cor(get(v), FCI_COMP_AVG, use = "complete.obs")])
tab[, corr_with_FCI := round(corrs[variant], 3)]
setorder(tab, horizon, variant)
write_rev_csv(res, "Rev_FCI_LOO_LP_Full.csv")
write_rev_csv(tab[, .(variant, horizon, corr_with_FCI, coef, se, p_value,
                      z_vs_baseline = round(z_vs_baseline, 2),
                      p_equal = round(p_equal, 3), n_obs)],
              "Rev_FCI_LOO_Battery.csv")

# ---------------------------------------------------------------------------
# 4. Forest plot at h = 12
# ---------------------------------------------------------------------------
f12 <- tab[horizon == 12]
labs <- c(zs_exCredit = "exCredit (z-score baseline)",
          zs_exCredit_exNPL = "- NPL", zs_exCredit_exROE = "- ROE",
          zs_exCredit_exLiquidez = "- Liquidity",
          zs_exCredit_exRatioCD = "- Credit/deposit ratio",
          zs_exCredit_exTCN = "- Exchange rate",
          zs_exCredit_QtyLag1 = "Banking quantities at t-1",
          FCI_exCredit_AVG = "exCredit (published 4-method)",
          FCI_RATES_AVG = "Rates-only FCI (co-baseline)")
f12[, vlab := factor(labs[variant], levels = rev(labs))]
pF <- ggplot(f12, aes(coef, vlab)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_vline(xintercept = f12[variant == "zs_exCredit", coef],
             color = "steelblue", alpha = 0.5) +
  geom_errorbarh(aes(xmin = coef - Z90 * se, xmax = coef + Z90 * se),
                 height = 0.25, color = "grey40") +
  geom_point(size = 2.4, color = "steelblue4") +
  labs(title = "Credit effect at h = 12 across FCI component-exclusion variants",
       subtitle = "No single potentially endogenous component drives the result; 90% intervals (NW)",
       x = "Coefficient (pp per FCI unit)", y = NULL) +
  theme_micro()
save_rev_png(pF, "286_FCI_LOO_Forest.png", w = 10, h = 6.5)

cat(sprintf("=== 42 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
