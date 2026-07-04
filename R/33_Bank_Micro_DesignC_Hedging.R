# ============================================================================
# 33_Bank_Micro_DesignC_Hedging.R
# Design C (plan §4): hedging heterogeneity - the theory test.
#
#   Δʰ ln L_{b,s,c,t+h} = β₁ₕ (Shockₜ × USD_c) + β₂ₕ (Shockₜ × USD_c × Unhedgedₛ)
#                          + α_{b×s×t} + γ_{b×s×c} + ε
#
# Distinguishing prediction of the dual-engine framework: the USD contraction
# concentrates in UNHEDGED (guarani-earning) sectors; naturally hedged dollar
# earners (agriculture, cattle) are flat or mildly positive. Pure GFC or pure
# trade-channel stories do not predict this sign pattern (response to R1 #4).
#
# Sharp test: hedged vs unhedged sectors only (ambiguous excluded).
# Gradient test: continuous treatment = sector FX share in 2016 (proxy for
# dollar-earning intensity; national-accounts export intensity not available
# at this sector breakdown - documented deviation from plan §4.1).
# Caveats built in (plan §4.3): Consumo/Vivienda tiny FX books -> mid-
# dollarization contrast as primary; agriculture seasonality -> month-of-year
# x sector x currency FE robustness.
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
library(fixest)
cat("=== 33: Design C - hedging heterogeneity ===\n")

dA <- read_rds_micro("micro_designA_panel.rds")   # from script 31
p2 <- read_rds_micro("micro_p2_sector.rds")

CTRLS <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
FE_A  <- "bst + bsc"

# Sector classification (plan §4.1)
dA[, hedge_class := fifelse(sector %in% SECT_HEDGED,   "hedged",
                    fifelse(sector %in% SECT_UNHEDGED, "unhedged",
                    fifelse(sector %in% SECT_AMBIGUOUS, "ambiguous", "excluded")))]
dA[, unhedged := as.integer(hedge_class == "unhedged")]
dA[, shk_usd_unh := shk_usd * unhedged]

# Mid-dollarization unhedged subset (primary contrast, plan §4.3 caveat)
SECT_UNH_MID <- c("COMERCIO AL POR MENOR", "SERVICIOS", "CONSTRUCCION")
# Tiny-FX sectors admitted only above the 5 bn threshold
SECT_TINY_FX <- c("CONSUMO", "VIVIENDA")

# Continuous gradient: sector FX share in 2016 (higher = more dollar-earning)
fx16 <- p2[ym <= to_ym("2016-12-01"), .(tot = sum(total)), by = .(sector, cur)]
fx16 <- dcast(fx16, sector ~ cur, value.var = "tot")
setnames(fx16, c("6200", "6900"), c("fc", "pyg"))
fx16[, fx_share_2016 := fc / (fc + pyg)]
dA <- merge(dA, fx16[, .(sector, fx_share_2016)], by = "sector", all.x = TRUE)
dA[, fxsh_z := zstd(fx_share_2016)]
dA[, shk_usd_fxsh := shk_usd * fxsh_z]

# Month-of-year FE for seasonality robustness
dA[, moy := (ym - 1L) %% 12L + 1L]
dA[, scm := paste(sector, cur, moy)]

# ---------------------------------------------------------------------------
# 1. Estimation variants
# ---------------------------------------------------------------------------
run_triple <- function(d, tripvar, label, fe = FE_A, wcb_h = H_KEY, B = WCB_B) {
  res <- list()
  for (h in HORIZONS) {
    yv  <- paste0("g", h)
    ds  <- d[get(paste0("ok", h)) == TRUE]
    rhs <- c("shk_usd", tripvar, CTRLS)
    m <- tryCatch(
      feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", fe)),
            data = ds, cluster = ~bank, notes = FALSE),
      error = function(e) NULL)
    if (is.null(m)) next
    for (pv in c("shk_usd", tripvar)) {
      ib <- inference_battery(m, pv, dk_lag = h + 1)
      pw <- if (h %in% wcb_h && pv == tripvar) wcb_pval(ds, yv, rhs, fe, pv, B = B) else NA_real_
      res[[length(res) + 1]] <- data.table(
        variant = label, param = pv, h = h, b = ib$b, se_bank = ib$se_bank,
        p_bank = ib$p_bank, p_wcb = pw, se_dk = ib$se_dk, p_dk = ib$p_dk, n = ib$n)
    }
    cat(sprintf("  [%s] h=%2d: b1=%7.3f  b2(%s)=%7.3f (p=%.3f)\n", label, h,
                coef(m)["shk_usd"], tripvar, coef(m)[tripvar],
                fixest::coeftable(m)[tripvar, 4]))
  }
  rbindlist(res)
}

# (a) Sharp test: hedged vs unhedged, tiny-FX sectors only above 5 bn
sharp <- dA[hedge_class %in% c("hedged", "unhedged") &
            (!sector %in% SECT_TINY_FX | total > 5000)]
cat("\n(a) Sharp test: hedged vs unhedged\n")
resC_sharp <- run_triple(sharp, "shk_usd_unh", "sharp")

# (b) Primary contrast: hedged vs mid-dollarization unhedged (plan §4.3)
mid <- dA[sector %in% c(SECT_HEDGED, SECT_UNH_MID)]
cat("\n(b) Primary contrast: hedged vs mid-dollarization unhedged\n")
resC_mid <- run_triple(mid, "shk_usd_unh", "mid_dollarization")

# (c) Gradient: continuous 2016 sector FX share (all sectors incl. ambiguous)
grad <- dA[hedge_class != "excluded"]
cat("\n(c) Gradient: continuous sector FX share\n")
resC_grad <- run_triple(grad, "shk_usd_fxsh", "gradient_fxshare")

# (d) Seasonality robustness: add sector x currency x month-of-year FE
cat("\n(d) Sharp test + month-of-year x sector x currency FE\n")
resC_moy <- run_triple(sharp, "shk_usd_unh", "sharp_moyFE",
                       fe = paste(FE_A, "+ scm"), wcb_h = 12)

resC <- rbind(resC_sharp, resC_mid, resC_grad, resC_moy)
write_micro_csv(resC, "Micro_DesignC_Hedging_Results.csv")

# ---------------------------------------------------------------------------
# 2. Split-sample view: beta_h estimated separately on hedged / unhedged
# ---------------------------------------------------------------------------
cat("\nSplit-sample estimates...\n")
split_res <- list()
for (grp in c("hedged", "unhedged")) {
  dg <- dA[hedge_class == grp & (!sector %in% SECT_TINY_FX | total > 5000)]
  for (h in HORIZONS) {
    ds <- dg[get(paste0("ok", h)) == TRUE]
    m <- tryCatch(
      feols(as.formula(paste0("g", h, " ~ shk_usd + ", paste(CTRLS, collapse = "+"),
                              " | ", FE_A)), data = ds, cluster = ~bank, notes = FALSE),
      error = function(e) NULL)
    if (is.null(m)) next
    cr <- coef_row(m, "shk_usd")
    split_res[[length(split_res) + 1]] <- data.table(
      group = grp, h = h, b = cr$b, se_bank = cr$se, p_bank = cr$p, n = cr$n)
  }
}
split_res <- rbindlist(split_res)
write_micro_csv(split_res, "Micro_DesignC_Split_Sample.csv")

# ---------------------------------------------------------------------------
# 3. Figures
# ---------------------------------------------------------------------------
d2 <- resC[variant == "sharp" & param == "shk_usd_unh"]
plot_dynamic_beta(d2,
  "Design C: triple interaction - Shock x USD x Unhedged sector",
  "Contraction concentrated in unhedged USD credit is the dual-engine's distinguishing prediction; 90% bands",
  "265_DesignC_Triple_Interaction.png")

pS <- ggplot(split_res, aes(h, b, color = group, fill = group)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = b - 1.645 * se_bank, ymax = b + 1.645 * se_bank),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  scale_color_manual(values = c(hedged = "forestgreen", unhedged = "firebrick"),
                     labels = c("Hedged (dollar earners)", "Unhedged (guarani earners)")) +
  scale_fill_manual(values = c(hedged = "forestgreen", unhedged = "firebrick"),
                    labels = c("Hedged (dollar earners)", "Unhedged (guarani earners)")) +
  labs(title = "Design C: within-cell USD-vs-PYG differential by borrower hedging status",
       subtitle = "Split-sample beta_h; hedged = agriculture + cattle; unhedged = guarani-earning sectors; 90% bands",
       x = "Horizon (months)", y = expression(beta[h] ~ "(pp)"),
       color = NULL, fill = NULL) +
  theme_micro()
save_png(pS, "266_DesignC_Hedged_vs_Unhedged.png")

cat(sprintf("=== 33 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
