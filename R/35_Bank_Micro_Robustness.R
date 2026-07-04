# ============================================================================
# 35_Bank_Micro_Robustness.R
# Phase 4 (plan §5): robustness battery for Design A (appendix).
#
# Multiple-testing discipline (plan §5): the single pre-specified headline is
# beta_12 from Design A with the purged DXY shock and baseline filters
# (script 31). Everything here is supporting evidence.
#
# Variants: COVID exclusion/interaction; alternative shocks; merger
# treatments (no-drop / balanced 13-bank / pro-forma); cell-size thresholds;
# deflation variants (incl. the deliberately UNADJUSTED book to show the sign
# of the valuation bias); asymmetry; extensive margin; cell-specific trends;
# winsorization sensitivity. Wild-bootstrap p at h = 12 for every variant.
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
library(fixest)
cat("=== 35: Robustness battery (Design A) ===\n")

macro <- read_rds_micro("micro_macro.rds")
p2    <- read_rds_micro("micro_p2_sector.rds")
p4    <- read_rds_micro("micro_p4_bankchars.rds")
meta  <- read_rds_micro("micro_meta.rds")

H_ROB  <- c(6, 12, 18)
CTRLS  <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
FE_A   <- "bst + bsc"

# ---------------------------------------------------------------------------
# 0. Panel builder (parameterized version of script 31's)
# ---------------------------------------------------------------------------
build_panel <- function(p2x = p2, min_cell = CELL_MIN_BASE, winsor = WINSOR_P,
                        lval_var = "lval", drop_mode = "baseline") {
  d <- copy(p2x)
  if (lval_var == "lval_nomPYG") {  # FC in USD, PYG nominal (no CPI deflation)
    d[, lval_nomPYG := fifelse(cur == CUR_FC, log(total / TCN), log(total))]
    d[total <= 0, lval_nomPYG := NA_real_]
  }
  add_lead_growth(d, lval_var, by = c("bank", "sector", "cur"))
  add_lag_growth(d, lval_var, by = c("bank", "sector", "cur"))
  setorderv(d, c("bank", "sector", "cur", "ym"))
  for (h in H_ROB) {
    d[, paste0("ok", h) := {
      idx <- match(ym + h, ym)
      total > min_cell & !is.na(total[idx]) & total[idx] > min_cell
    }, by = .(bank, sector, cur)]
  }
  if (!is.null(winsor)) {
    for (h in H_ROB) d[, paste0("g", h) := winsorize(get(paste0("g", h)), winsor), by = cur]
  }
  d <- merge(d, macro[, .(ym, shk_dxy_purged, shk_dxy_raw, shk_vix, shk_fci)],
             by = "ym", all.x = TRUE)
  d[, usd := as.integer(cur == CUR_FC)]
  d <- merge(d, p4[, .(bank, ym, tier1_l12_z, size_l12_z, fc_dep_share_l12_z)],
             by = c("bank", "ym"), all.x = TRUE)
  d[, `:=`(shk_usd = shk_dxy_purged * usd, shkraw_usd = shk_dxy_raw * usd,
           vix_usd = shk_vix * usd, fci_usd = shk_fci * usd,
           glag_usd = g_lag12 * usd, tier1_usd = tier1_l12_z * usd,
           size_usd = size_l12_z * usd, fcdep_usd = fc_dep_share_l12_z * usd)]
  d <- d[ym >= to_ym("2016-01-01") & ym <= to_ym("2025-12-01") & !is.na(shk_dxy_purged)]
  if (drop_mode == "baseline")   d <- d[drop_event == FALSE]
  if (drop_mode == "balanced13") d <- d[bank %in% meta$continuous_banks]
  d[, n_cur := uniqueN(cur[total > min_cell]), by = .(bank, sector, ym)]
  d <- d[n_cur == 2]
  d[, `:=`(bst = paste(bank, sector, ym), bsc = paste(bank, sector, cur))]
  d
}

# Pro-forma merged panel: fold each exiting bank into its absorber over the
# entire sample, then rebuild (plan §1.3 robustness (i))
p2_proforma <- local({
  px <- copy(p2)
  ev <- meta$events[classified == "merger" & !is.na(absorber)]
  for (i in seq_len(nrow(ev))) px[bank == ev$exiting_bank[i], bank := ev$absorber[i]]
  px[, .(total = sum(total, na.rm = TRUE), TCN = TCN[1], IPC = IPC[1],
         drop_event = FALSE),
     by = .(bank, sector, cur, ym)][
    , lval := fifelse(cur == CUR_FC, log(total / TCN), log(total / IPC * 100))][
      total <= 0, lval := NA_real_][]
})

# ---------------------------------------------------------------------------
# 1. Variant runner
# ---------------------------------------------------------------------------
run_variant <- function(d, variant, shockvar = "shk_usd", extra_rhs = NULL,
                        fe = FE_A, subset_expr = NULL, B = WCB_B) {
  res <- list()
  for (h in H_ROB) {
    yv <- paste0("g", h)
    ds <- d[get(paste0("ok", h)) == TRUE]
    if (!is.null(subset_expr)) ds <- ds[eval(subset_expr)]
    rhs <- c(shockvar, extra_rhs, CTRLS)
    m <- tryCatch(
      feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", fe)),
            data = ds, cluster = ~bank, notes = FALSE),
      error = function(e) NULL)
    if (is.null(m)) next
    params <- c(shockvar, extra_rhs)
    for (pv in params) {
      cr <- coef_row(m, pv)
      pw <- if (h == 12) wcb_pval(ds, yv, rhs, fe, pv, B = B) else NA_real_
      res[[length(res) + 1]] <- data.table(variant = variant, param = pv, h = h,
                                           b = cr$b, se_bank = cr$se,
                                           p_bank = cr$p, p_wcb = pw, n = cr$n)
    }
    cat(sprintf("  [%s] h=%2d: b=%7.3f (p=%.3f)\n", variant, h,
                coef(m)[shockvar], fixest::coeftable(m)[shockvar, 4]))
  }
  rbindlist(res)
}

RES <- list()
base <- build_panel()

# 1. Baseline (reference row, matches script 31)
RES$base <- run_variant(base, "baseline_purgedDXY")

# 2. COVID (plan §5.1)
RES$cvd  <- run_variant(base, "covid_excluded",
                        subset_expr = quote(ym < COVID_START | ym > COVID_END))
base[, covid := as.integer(ym >= COVID_START & ym <= COVID_END)]
base[, shk_usd_cov := shk_usd * covid]
RES$cvdi <- run_variant(base, "covid_interaction", extra_rhs = "shk_usd_cov")

# 3. Alternative shocks (plan §5.2; Fed broad dollar index not in database -
#    documented as unavailable)
RES$raw <- run_variant(base, "raw_DXY", shockvar = "shkraw_usd")
RES$fci <- run_variant(base, "FCI_shock", shockvar = "fci_usd")
RES$vix <- run_variant(base, "VIX", shockvar = "vix_usd")

# 4. Merger treatments (plan §5.3)
RES$nodrop <- run_variant(build_panel(drop_mode = "none"), "no_event_drop")
RES$bal13  <- run_variant(build_panel(drop_mode = "balanced13"), "balanced_13banks")
RES$profo  <- run_variant(build_panel(p2x = p2_proforma), "proforma_mergers")

# 5. Cell-size thresholds (plan §5.4)
RES$th5  <- run_variant(build_panel(min_cell = 5000),  "threshold_5bn")
RES$th10 <- run_variant(build_panel(min_cell = 10000), "threshold_10bn")

# 6. Deflation / valuation variants (plan §5.5; US CPI not in database, so the
#    nominal-USD FC book is the baseline and PYG-nominal is the variant)
RES$nomP <- run_variant(build_panel(lval_var = "lval_nomPYG"), "FCusd_PYGnominal")
RES$unadj <- run_variant(build_panel(lval_var = "lval_nom"), "UNADJUSTED_valuation_bias_demo")

# 7. Asymmetry (plan §5.6)
base[, `:=`(shk_pos_usd = pmax(shk_dxy_purged, 0) * usd,
            shk_neg_usd = pmin(shk_dxy_purged, 0) * usd)]
RES$asym <- run_variant(base, "asymmetric_pos", shockvar = "shk_pos_usd",
                        extra_rhs = "shk_neg_usd")

# 8. Cell-specific linear trends (plan §5.8): de-dollarization must not drive beta
RES$trend <- run_variant(base, "cell_linear_trends", fe = "bst + bsc[[ym]]")

# 9. Winsorization sensitivity (plan §1.4)
RES$w25 <- run_variant(build_panel(winsor = c(0.025, 0.975)), "winsor_2.5pct")
RES$w0  <- run_variant(build_panel(winsor = NULL), "no_winsorizing")

rob <- rbindlist(RES)
write_micro_csv(rob, "Micro_Robustness_Battery.csv")

# ---------------------------------------------------------------------------
# 2. Extensive margin (plan §5.7): FX-lending activity indicator (LPM)
# ---------------------------------------------------------------------------
cat("\nExtensive margin (LPM on FX-cell activity)...\n")
fcgrid <- p2[cur == CUR_FC]
fcgrid[, ever := any(total > CELL_MIN_BASE), by = .(bank, sector)]
fcgrid <- fcgrid[ever == TRUE]
setorderv(fcgrid, c("bank", "sector", "ym"))
for (h in H_ROB) {
  fcgrid[, paste0("act", h) := {
    idx <- match(ym + h, ym)
    as.integer(!is.na(total[idx]) & total[idx] > CELL_MIN_BASE)
  }, by = .(bank, sector)]
}
fcgrid <- merge(fcgrid, macro[, .(ym, shk = shk_dxy_purged)], by = "ym", all.x = TRUE)
fcgrid[, `:=`(bsec = paste(bank, sector), moy = (ym - 1L) %% 12L + 1L)]
fcgrid <- fcgrid[ym >= to_ym("2016-01-01") & ym <= to_ym("2025-12-01") &
                 drop_event == FALSE & !is.na(shk)]
ext <- list()
for (h in H_ROB) {
  m <- feols(as.formula(paste0("act", h, " ~ shk | bsec + moy")),
             data = fcgrid, cluster = ~bank, notes = FALSE)
  cr <- coef_row(m, "shk")
  ext[[length(ext) + 1]] <- data.table(h = h, b = cr$b, se_bank = cr$se,
                                       p_bank = cr$p, n = cr$n)
  cat(sprintf("  h=%2d: b=%.4f (p=%.3f)\n", h, cr$b, cr$p))
}
write_micro_csv(rbindlist(ext), "Micro_Extensive_Margin_LPM.csv")

# ---------------------------------------------------------------------------
# 3. Forest plot of beta_12 across variants
# ---------------------------------------------------------------------------
f12 <- rob[h == 12 & param %in% c("shk_usd", "shkraw_usd", "fci_usd", "vix_usd",
                                  "shk_pos_usd")]
f12[, variant := factor(variant, levels = rev(unique(variant)))]
pF <- ggplot(f12, aes(x = b, y = variant)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_vline(xintercept = rob[variant == "baseline_purgedDXY" & h == 12 &
                              param == "shk_usd", b],
             color = "steelblue", alpha = 0.5) +
  geom_errorbarh(aes(xmin = b - 1.645 * se_bank, xmax = b + 1.645 * se_bank),
                 height = 0.25, color = "grey40") +
  geom_point(size = 2.4, color = "steelblue4") +
  labs(title = expression("Robustness battery: Design A " * beta[12] * " across variants"),
       subtitle = "90% intervals, cluster by bank; blue line = baseline estimate",
       x = expression(beta[12] ~ "(pp)"), y = NULL) +
  theme_micro()
save_png(pF, "269_Robustness_Forest_Beta12.png", w = 10, h = 7)

cat(sprintf("=== 35 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
