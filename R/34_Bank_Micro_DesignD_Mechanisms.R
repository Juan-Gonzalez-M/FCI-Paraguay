# ============================================================================
# 34_Bank_Micro_DesignD_Mechanisms.R
# Design D (plan §3.3): mechanism corroboration - deposits and NPLs by currency.
#
# (i)  Deposit flows: Δʰ ln Deposits_{b,c,t+h} on Shock x USD with bank x time
#      FE - does the public fly INTO dollar deposits (funding squeeze on the
#      PYG side) or OUT (dollar liquidity squeeze)? Either signs the funding
#      mechanism.
# (ii) NPL response: Vencida/(Vigente+Vencida) by bank x sector x currency, in
#      pp changes (bounded ratio - no logs). A dollar shock should raise
#      delinquency in USD books of unhedged sectors with a 6-12 month lag.
#      This answers R1 #5: the FCI's banking components are OUTCOMES of the
#      identified shock, which is exactly why the leave-one-out and rates-only
#      co-baseline FCIs are used.
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
library(fixest)
cat("=== 34: Design D - mechanisms (deposits and NPLs by currency) ===\n")

macro <- read_rds_micro("micro_macro.rds")
dep   <- read_rds_micro("micro_deposits.rds")
dA    <- read_rds_micro("micro_designA_panel.rds")

# ---------------------------------------------------------------------------
# 1. Deposit flows by currency
# ---------------------------------------------------------------------------
add_lead_growth(dep, "lval", by = c("bank", "cur"))
for (h in HORIZONS) dep[, paste0("g", h) := winsorize(get(paste0("g", h))), by = cur]
dep <- merge(dep, macro[, .(ym, shk = shk_dxy_purged)], by = "ym", all.x = TRUE)
dep[, usd := as.integer(cur == CUR_FC)]
dep[, shk_usd := shk * usd]
dep[, `:=`(bt = paste(bank, ym), bc = paste(bank, cur))]
depE <- dep[ym >= to_ym("2016-01-01") & ym <= to_ym("2025-12-01") &
            drop_event == FALSE & !is.na(shk) & deposits > CELL_MIN_BASE]
cat(sprintf("  Deposit panel: %d rows\n", nrow(depE)))

res_dep <- list()
for (h in HORIZONS) {
  yv <- paste0("g", h)
  m <- tryCatch(feols(as.formula(paste(yv, "~ shk_usd | bt + bc")),
                      data = depE, cluster = ~bank, notes = FALSE),
                error = function(e) NULL)
  if (is.null(m)) next
  ib <- inference_battery(m, "shk_usd", dk_lag = h + 1)
  pw <- if (h %in% H_KEY) wcb_pval(depE, yv, "shk_usd", "bt + bc", "shk_usd") else NA_real_
  res_dep[[h]] <- data.table(outcome = "deposits", h = h, b = ib$b,
                             se_bank = ib$se_bank, p_bank = ib$p_bank, p_wcb = pw,
                             se_dk = ib$se_dk, p_dk = ib$p_dk, n = ib$n)
  cat(sprintf("  [deposits] h=%2d: b=%7.3f  p_bank=%.3f\n", h, ib$b, ib$p_bank))
}
res_dep <- rbindlist(res_dep)

# ---------------------------------------------------------------------------
# 2. NPL response (pp changes) - overall USD differential and unhedged triple
# ---------------------------------------------------------------------------
# dA already carries npl, hedging class comes from script 33's classification
dA[, hedge_class := fifelse(sector %in% SECT_HEDGED,   "hedged",
                    fifelse(sector %in% SECT_UNHEDGED, "unhedged",
                    fifelse(sector %in% SECT_AMBIGUOUS, "ambiguous", "excluded")))]
dA[, unhedged := as.integer(hedge_class == "unhedged")]

setorderv(dA, c("bank", "sector", "cur", "ym"))
for (h in HORIZONS) {
  dA[, paste0("dnpl", h) := {
    idx <- match(ym + h, ym)
    npl[idx] - npl
  }, by = .(bank, sector, cur)]
  dA[, paste0("dnpl", h) := winsorize(get(paste0("dnpl", h))), by = cur]
}
dA[, shk_usd_unh := shk_usd * unhedged]

res_npl <- list()
for (h in HORIZONS) {
  yv <- paste0("dnpl", h)
  ds <- dA[get(paste0("ok", h)) == TRUE]
  # (a) average USD differential
  m1 <- tryCatch(feols(as.formula(paste(yv, "~ shk_usd | bst + bsc")),
                       data = ds, cluster = ~bank, notes = FALSE),
                 error = function(e) NULL)
  # (b) triple with unhedged (sharp sectors only)
  ds2 <- ds[hedge_class %in% c("hedged", "unhedged")]
  m2 <- tryCatch(feols(as.formula(paste(yv, "~ shk_usd + shk_usd_unh | bst + bsc")),
                       data = ds2, cluster = ~bank, notes = FALSE),
                 error = function(e) NULL)
  if (!is.null(m1)) {
    ib <- inference_battery(m1, "shk_usd", dk_lag = h + 1)
    pw <- if (h %in% H_KEY) wcb_pval(ds, yv, "shk_usd", "bst + bsc", "shk_usd") else NA_real_
    res_npl[[length(res_npl) + 1]] <- data.table(
      outcome = "npl_pp", param = "shk_usd", h = h, b = ib$b,
      se_bank = ib$se_bank, p_bank = ib$p_bank, p_wcb = pw,
      se_dk = ib$se_dk, p_dk = ib$p_dk, n = ib$n)
  }
  if (!is.null(m2)) {
    for (pv in c("shk_usd", "shk_usd_unh")) {
      cr <- coef_row(m2, pv)
      res_npl[[length(res_npl) + 1]] <- data.table(
        outcome = "npl_pp_triple", param = pv, h = h, b = cr$b,
        se_bank = cr$se, p_bank = cr$p, p_wcb = NA_real_,
        se_dk = NA_real_, p_dk = NA_real_, n = cr$n)
    }
  }
  if (!is.null(m1))
    cat(sprintf("  [NPL] h=%2d: dUSD=%6.3f pp (p=%.3f)%s\n", h,
                coef(m1)["shk_usd"], fixest::coeftable(m1)["shk_usd", 4],
                if (!is.null(m2)) sprintf("  x unhedged=%6.3f (p=%.3f)",
                    coef(m2)["shk_usd_unh"], fixest::coeftable(m2)["shk_usd_unh", 4]) else ""))
}
res_npl <- rbindlist(res_npl)

write_micro_csv(rbind(res_dep[, .(outcome, param = "shk_usd", h, b, se_bank,
                                  p_bank, p_wcb, se_dk, p_dk, n)], res_npl),
                "Micro_DesignD_Mechanisms.csv")

# ---------------------------------------------------------------------------
# 3. Figures
# ---------------------------------------------------------------------------
plot_dynamic_beta(res_dep,
  "Design D: differential USD-vs-PYG deposit response to a 1 SD dollar shock",
  "Bank x time FE; positive = flight INTO dollar deposits, negative = dollar funding squeeze; 90% bands",
  "267_DesignD_Deposits_by_Currency.png")

npl_cmp <- rbind(
  res_npl[outcome == "npl_pp" & param == "shk_usd",
          .(h, b, se = se_bank, Series = "USD differential (all sectors)")],
  res_npl[outcome == "npl_pp_triple" & param == "shk_usd_unh",
          .(h, b, se = se_bank, Series = "x Unhedged (additional)")])
pN <- ggplot(npl_cmp, aes(h, b, color = Series, fill = Series)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = b - 1.645 * se, ymax = b + 1.645 * se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  scale_color_manual(values = c("firebrick", "darkorange3")) +
  scale_fill_manual(values = c("firebrick", "darkorange3")) +
  labs(title = "Design D: NPL response by currency (pp) after a 1 SD dollar shock",
       subtitle = "Rising USD delinquency in unhedged sectors signs the borrower balance-sheet channel (R1 #5)",
       x = "Horizon (months)", y = "Delta NPL ratio (pp)") +
  theme_micro()
save_png(pN, "268_DesignD_NPL_by_Currency.png")

cat(sprintf("=== 34 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
