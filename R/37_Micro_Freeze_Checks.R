# ============================================================================
# 37_Micro_Freeze_Checks.R — Tier-1 pre-freeze checks on the micro results
# (Coworker 1's micro assessment):
#  1. COVID-excluded and tightening-only versions of the Design C GRADIENT and
#     hedged/unhedged split — the tests that will carry the micro section.
#  2. h = 1-3 TCN-timing audit (is the fast h=1 effect a valuation-rate
#     timing artifact?).
#  3. Reconciliation wedge decomposition (is the 0.926 micro-vs-aggregate
#     correlation fully attributable to financieras coverage?).
# ADDITIVE: reads micro rds caches; writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
library(fixest)
cat("=== 37: Micro freeze checks ===\n")

dA    <- read_rds_micro("micro_designA_panel.rds")
p2    <- read_rds_micro("micro_p2_sector.rds")
p1    <- read_rds_micro("micro_p1_carteras.rds")
p4    <- read_rds_micro("micro_p4_bankchars.rds")
macro <- read_rds_micro("micro_macro.rds")

CTRLS <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
FE_A  <- "bst + bsc"

# ---------------------------------------------------------------------------
# 1. Design C gradient & split: COVID-excluded and tightening-only
# ---------------------------------------------------------------------------
cat("\n[1] Design C gradient/split robustness...\n")

# Reconstruct script 33's classification and gradient treatment
dA[, hedge_class := fifelse(sector %in% SECT_HEDGED,   "hedged",
                    fifelse(sector %in% SECT_UNHEDGED, "unhedged",
                    fifelse(sector %in% SECT_AMBIGUOUS, "ambiguous", "excluded")))]
fx16 <- p2[ym <= to_ym("2016-12-01"), .(tot = sum(total)), by = .(sector, cur)]
fx16 <- dcast(fx16, sector ~ cur, value.var = "tot")
setnames(fx16, c("6200", "6900"), c("fc", "pyg"))
fx16[, fx_share_2016 := fc / (fc + pyg)]
dA <- merge(dA, fx16[, .(sector, fx_share_2016)], by = "sector", all.x = TRUE)
dA[, fxsh_z := zstd(fx_share_2016)]
dA[, shk_pos := pmax(shk_dxy_purged, 0)]
dA[, `:=`(shk_usd_fxsh = shk_usd * fxsh_z,
          shkpos_usd = shk_pos * usd)]
dA[, `:=`(shkpos_usd_fxsh = shkpos_usd * fxsh_z)]
dA[, covid := ym >= COVID_START & ym <= COVID_END]

SECT_TINY_FX <- c("CONSUMO", "VIVIENDA")
sharp <- quote(hedge_class %in% c("hedged", "unhedged") &
               (!sector %in% SECT_TINY_FX | total > 5000))
grad  <- quote(hedge_class != "excluded")

run_var <- function(d, label, shockvar, tripvar, sub, drop_covid = FALSE,
                    wcb_h = c(6, 12, 18)) {
  res <- list()
  for (h in HORIZONS) {
    yv <- paste0("g", h)
    ds <- d[get(paste0("ok", h)) == TRUE][eval(sub)]
    if (drop_covid) ds <- ds[covid == FALSE]
    rhs <- c(shockvar, tripvar, CTRLS)
    m <- tryCatch(feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"),
                                         "|", FE_A)),
                        data = ds, cluster = ~bank, notes = FALSE),
                  error = function(e) NULL)
    if (is.null(m) || !tripvar %in% names(coef(m))) next
    cr <- coef_row(m, tripvar)
    pw <- if (h %in% wcb_h) wcb_pval(ds, yv, rhs, FE_A, tripvar) else NA_real_
    res[[length(res) + 1]] <- data.table(variant = label, param = tripvar, h = h,
                                         b = cr$b, se_bank = cr$se, p_bank = cr$p,
                                         p_wcb = pw, n = cr$n)
  }
  rbindlist(res)
}

split_var <- function(d, label, shockvar, drop_covid = FALSE) {
  res <- list()
  for (grp in c("hedged", "unhedged")) {
    dg <- d[hedge_class == grp & (!sector %in% SECT_TINY_FX | total > 5000)]
    if (drop_covid) dg <- dg[covid == FALSE]
    for (h in HORIZONS) {
      ds <- dg[get(paste0("ok", h)) == TRUE]
      m <- tryCatch(feols(as.formula(paste0("g", h, " ~ ", shockvar, " + ",
                                            paste(CTRLS, collapse = "+"), " | ", FE_A)),
                          data = ds, cluster = ~bank, notes = FALSE),
                    error = function(e) NULL)
      if (is.null(m) || !shockvar %in% names(coef(m))) next
      cr <- coef_row(m, shockvar)
      res[[length(res) + 1]] <- data.table(variant = label, group = grp, h = h,
                                           b = cr$b, se_bank = cr$se,
                                           p_bank = cr$p, n = cr$n)
    }
  }
  rbindlist(res)
}

gradient_res <- rbind(
  run_var(dA, "gradient_baseline",   "shk_usd",    "shk_usd_fxsh",    grad),
  run_var(dA, "gradient_exCOVID",    "shk_usd",    "shk_usd_fxsh",    grad, drop_covid = TRUE),
  run_var(dA, "gradient_tightening", "shkpos_usd", "shkpos_usd_fxsh", grad),
  run_var(dA, "gradient_tight_exCOVID", "shkpos_usd", "shkpos_usd_fxsh", grad,
          drop_covid = TRUE))
dA[, shk_usd_unh := shk_usd * as.integer(hedge_class == "unhedged")]
dA[, shkpos_usd_unh := shkpos_usd * as.integer(hedge_class == "unhedged")]
triple_res <- rbind(
  run_var(dA, "sharp_baseline",   "shk_usd",    "shk_usd_unh",    sharp),
  run_var(dA, "sharp_exCOVID",    "shk_usd",    "shk_usd_unh",    sharp, drop_covid = TRUE),
  run_var(dA, "sharp_tightening", "shkpos_usd", "shkpos_usd_unh", sharp))
split_res <- rbind(
  split_var(dA, "split_baseline", "shk_usd"),
  split_var(dA, "split_exCOVID",  "shk_usd", drop_covid = TRUE),
  split_var(dA, "split_tightening", "shkpos_usd"))

out1 <- rbind(gradient_res, triple_res, fill = TRUE)
write_rev_csv(out1, "Rev_DesignC_Gradient_Robustness.csv")
write_rev_csv(split_res, "Rev_DesignC_Split_Robustness.csv")

cat("\nGradient at h=12 across variants:\n")
print(gradient_res[h == 12, .(variant, b = round(b, 3), p_bank = round(p_bank, 3),
                              p_wcb = round(p_wcb, 3), n)])
cat("\nHedged split at h=12:\n")
print(split_res[h == 12, .(variant, group, b = round(b, 3), p = round(p_bank, 3))])

gplot <- gradient_res[variant %in% c("gradient_baseline", "gradient_exCOVID",
                                     "gradient_tightening")]
gplot[, Variant := fcase(variant == "gradient_baseline", "Baseline",
                         variant == "gradient_exCOVID", "COVID excluded",
                         variant == "gradient_tightening", "Tightening only")]
pG <- ggplot(gplot, aes(h, b, color = Variant, fill = Variant)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = b - 1.645 * se_bank, ymax = b + 1.645 * se_bank),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
  scale_color_manual(values = c("steelblue4", "firebrick", "darkorange3")) +
  scale_fill_manual(values = c("steelblue4", "firebrick", "darkorange3")) +
  labs(title = "Design C hedging gradient: robustness of the load-bearing test",
       subtitle = "Shock x USD x sector FX-share (per SD); 90% bands, cluster by bank",
       x = "Horizon (months)", y = "Gradient coefficient (pp)") +
  theme_micro()
save_rev_png(pG, "280_DesignC_Gradient_Robustness.png")

# ---------------------------------------------------------------------------
# 2. TCN-timing audit for h = 1-3
# ---------------------------------------------------------------------------
cat("\n[2] TCN-timing audit (h = 1-3)...\n")
fxv <- macro[, .(ym, TCN, IPC)]
setorder(fxv, ym)
fxv[, `:=`(TCN_l1 = shift(TCN, 1))]
fxv[, TCN_avg := (TCN + TCN_l1) / 2]

audit <- list()
for (conv in c("baseline_t", "lag1", "avg_t_t1")) {
  d <- merge(p2[, .(bank, sector, cur, ym, total, drop_event)],
             fxv, by = "ym", all.x = TRUE)
  e_use <- switch(conv, baseline_t = quote(TCN), lag1 = quote(TCN_l1),
                  avg_t_t1 = quote(TCN_avg))
  d[, lval := fifelse(cur == CUR_FC, log(total / eval(e_use)),
                      log(total / IPC * 100))]
  d[total <= 0, lval := NA_real_]
  add_lead_growth(d, "lval", by = c("bank", "sector", "cur"), horizons = 1:3)
  add_lag_growth(d, "lval", by = c("bank", "sector", "cur"))
  setorderv(d, c("bank", "sector", "cur", "ym"))
  for (h in 1:3) {
    d[, paste0("ok", h) := {
      idx <- match(ym + h, ym)
      total > CELL_MIN_BASE & !is.na(total[idx]) & total[idx] > CELL_MIN_BASE
    }, by = .(bank, sector, cur)]
    d[, paste0("g", h) := winsorize(get(paste0("g", h))), by = cur]
  }
  d <- merge(d, macro[, .(ym, shk_dxy_purged)], by = "ym", all.x = TRUE)
  d[, usd := as.integer(cur == CUR_FC)]
  d <- merge(d, p4[, .(bank, ym, tier1_l12_z, size_l12_z, fc_dep_share_l12_z)],
             by = c("bank", "ym"), all.x = TRUE)
  d[, `:=`(shk_usd = shk_dxy_purged * usd, glag_usd = g_lag12 * usd,
           tier1_usd = tier1_l12_z * usd, size_usd = size_l12_z * usd,
           fcdep_usd = fc_dep_share_l12_z * usd)]
  d <- d[ym >= to_ym("2016-01-01") & ym <= to_ym("2025-12-01") &
         drop_event == FALSE & !is.na(shk_dxy_purged)]
  d[, n_cur := uniqueN(cur[total > CELL_MIN_BASE]), by = .(bank, sector, ym)]
  d <- d[n_cur == 2]
  d[, `:=`(bst = paste(bank, sector, ym), bsc = paste(bank, sector, cur))]
  for (h in 1:3) {
    ds <- d[get(paste0("ok", h)) == TRUE]
    m <- feols(as.formula(paste0("g", h, " ~ shk_usd + ",
                                 paste(CTRLS, collapse = "+"), " | ", FE_A)),
               data = ds, cluster = ~bank, notes = FALSE)
    cr <- coef_row(m, "shk_usd")
    audit[[length(audit) + 1]] <- data.table(convention = conv, h = h,
                                             b = cr$b, se = cr$se, p = cr$p, n = cr$n)
    cat(sprintf("  [%s] h=%d: b=%.3f (p=%.3f)\n", conv, h, cr$b, cr$p))
  }
}
audit <- rbindlist(audit)
write_rev_csv(audit, "Rev_TCN_Timing_Audit.csv")

# ---------------------------------------------------------------------------
# 3. Reconciliation wedge decomposition
# ---------------------------------------------------------------------------
cat("\n[3] Reconciliation wedge decomposition...\n")
dm <- as.data.table(suppressMessages(read_excel(micro_paths$fci_xlsx,
                                                sheet = "Datos_macro")))
dm[, ym := to_ym(Fecha)]
sys <- p1[, .(micro = sum(vigente + vencida, na.rm = TRUE)), by = ym]
rec <- merge(sys, dm[, .(ym, agg = Creditos)], by = "ym")
setorder(rec, ym)
rec[, `:=`(ratio = micro / agg,
           micro_yoy = 100 * (micro / shift(micro, 12) - 1),
           agg_yoy   = 100 * (agg / shift(agg, 12) - 1))]
rec <- merge(rec, macro[, .(ym, TCN)], by = "ym", all.x = TRUE)
rec[, `:=`(gap = micro_yoy - agg_yoy, tcn_yoy = 100 * (TCN / shift(TCN, 12) - 1),
           trend = seq_len(.N))]

# (a) coverage ratio path (financieras share of the aggregate = 1 - ratio)
# (b) growth correlation by subperiod
per <- rbind(
  data.table(period = "full",       corr = rec[, cor(micro_yoy, agg_yoy, use = "complete.obs")]),
  data.table(period = "2017-2019",  corr = rec[ym <= to_ym("2019-12-01"),
                                               cor(micro_yoy, agg_yoy, use = "complete.obs")]),
  data.table(period = "2020-2022",  corr = rec[ym >= to_ym("2020-01-01") & ym <= to_ym("2022-12-01"),
                                               cor(micro_yoy, agg_yoy, use = "complete.obs")]),
  data.table(period = "2023-2025",  corr = rec[ym >= to_ym("2023-01-01"),
                                               cor(micro_yoy, agg_yoy, use = "complete.obs")]))
# (c) what explains the gap: coverage trend vs FX valuation (both books are
# nominal PYG so TCN should NOT drive the gap unless coverage differs by currency)
mg <- lm(gap ~ trend + tcn_yoy, rec)
ct <- coeftest(mg, vcov. = NeweyWest(mg, lag = 13, prewhite = FALSE))
# (d) correlation after removing the smooth coverage drift
rec[, gap_detr := resid(lm(gap ~ poly(trend, 3), rec, na.action = na.exclude))]
corr_detr <- rec[, cor(micro_yoy - gap + gap_detr, agg_yoy, use = "complete.obs")]

cat(sprintf("  Coverage ratio: %.3f (2017) -> %.3f (2025); financieras share of aggregate: %.1f%% -> %.1f%%\n",
            rec[!is.na(micro_yoy)][1, ratio], rec[.N, ratio],
            100 * (1 - rec[!is.na(micro_yoy)][1, ratio]), 100 * (1 - rec[.N, ratio])))
cat("  YoY correlation by period:\n"); print(per)
cat(sprintf("  gap ~ trend + TCN_yoy: b_trend=%.4f (p=%.3f), b_tcn=%.3f (p=%.3f)\n",
            ct["trend", 1], ct["trend", 4], ct["tcn_yoy", 1], ct["tcn_yoy", 4]))
cat(sprintf("  Correlation after removing smooth coverage drift: %.3f (raw %.3f)\n",
            corr_detr, per[period == "full", corr]))

write_rev_csv(rec[, .(date = ym_date(ym), micro, agg, ratio, micro_yoy, agg_yoy,
                      gap, tcn_yoy)], "Rev_Reconciliation_Decomposition.csv")
write_rev_csv(per, "Rev_Reconciliation_Correlations.csv")

pR <- ggplot(rec[!is.na(gap)], aes(ym_date(ym))) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_line(aes(y = gap, color = "Micro - aggregate YoY growth gap"), linewidth = 0.7) +
  geom_line(aes(y = 100 * (1 - ratio), color = "Financieras share of aggregate (%)"),
            linewidth = 0.7) +
  scale_color_manual(values = c("darkorange3", "steelblue4")) +
  labs(title = "Reconciliation wedge: coverage drift, not construction error",
       subtitle = "The growth gap tracks the declining financieras share; TCN plays no role (both books nominal PYG)",
       x = NULL, y = "pp / %", color = NULL) +
  theme_micro()
save_rev_png(pR, "281_Reconciliation_Wedge.png")

cat(sprintf("=== 37 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
