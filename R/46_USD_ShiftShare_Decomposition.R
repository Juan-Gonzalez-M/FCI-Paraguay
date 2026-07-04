# ============================================================================
# 46_USD_ShiftShare_Decomposition.R — Reconciling the aggregate USD credit
# contraction with the positive within-cell currency differential
# (Colleague 1, feedback on changes & response, point 2 — "mandatory").
#
# Apparent tension: aggregate USD credit (in dollars) contracts MORE than
# total credit after tightening (script 36), while WITHIN bank x sector cells
# USD books outperform PYG books after a dollar shock (micro Design A).
# Reconciliation: aggregation. Decompose the aggregate USD-vs-PYG growth
# differential at horizon h into
#   WITHIN  = sum_c w_c * (g^USD_c - g^PYG_c)   (common cell weights)
#   BETWEEN = residual composition term (USD credit concentrated in cells
#             that contract more overall)
# and show which component responds to the dollar shock. Expected: within
# responds positively (the micro result), between negatively and dominantly.
# Also estimated on the FCI (the aggregate LP treatment) to document the
# treatment difference explicitly.
# ADDITIVE: writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
cat("=== 46: USD shift-share decomposition (within vs between cells) ===\n")

macro <- read_rds_micro("micro_macro.rds")
p2    <- read_rds_micro("micro_p2_sector.rds")

H_SET <- c(6, 12, 18)

# Cell-level valuation-adjusted log levels already in p2 (lval): FC in USD,
# PYG CPI-deflated. Compute h-horizon growth per cell and currency.
d <- p2[drop_event == FALSE & total > CELL_MIN_BASE,
        .(bank, sector, cur, ym, total, lval)]
setorderv(d, c("bank", "sector", "cur", "ym"))
for (h in H_SET) {
  d[, paste0("g", h) := {
    idx <- match(ym + h, ym)
    100 * (lval[idx] - lval)
  }, by = .(bank, sector, cur)]
  d[, paste0("g", h) := winsorize(get(paste0("g", h))), by = cur]
}

# Reshape to cell (bank x sector) x month with both currency books
wide <- dcast(d, bank + sector + ym ~ cur,
              value.var = c("total", paste0("g", H_SET)))
setnames(wide, gsub("_6200", "_FC", gsub("_6900", "_PYG", names(wide))))

series <- list()
for (h in H_SET) {
  gF <- paste0("g", h, "_FC"); gP <- paste0("g", h, "_PYG")
  w <- wide[!is.na(get(gF)) & !is.na(get(gP))]        # dual-growth cells
  s <- w[, {
    wF <- total_FC / sum(total_FC)                     # USD-book weights
    wP <- total_PYG / sum(total_PYG)                   # PYG-book weights
    wC <- (total_FC + total_PYG) / sum(total_FC + total_PYG)  # common weights
    agg_usd <- sum(wF * get(gF)); agg_pyg <- sum(wP * get(gP))
    within  <- sum(wC * (get(gF) - get(gP)))
    .(h = h, diff_agg = agg_usd - agg_pyg, within = within,
      between = (agg_usd - agg_pyg) - within, n_cells = .N)
  }, by = ym]
  series[[as.character(h)]] <- s
}
series <- rbindlist(series)
series <- merge(series, macro[, .(ym, shk = shk_dxy_purged)], by = "ym")
fci <- fread(micro_paths$fci_csv, select = c("fecha", "FCI_exCredit_AVG"))
fci[, ym := to_ym(fecha)]
series <- merge(series, fci[, .(ym, FCI = FCI_exCredit_AVG)], by = "ym", all.x = TRUE)
series <- series[ym >= to_ym("2016-01-01") & ym <= to_ym("2025-12-01")]

cat("Average decomposition (2016-2025):\n")
print(series[, .(diff_agg = round(mean(diff_agg), 2), within = round(mean(within), 2),
                 between = round(mean(between), 2)), by = h])

# Time-series LP of each component on the dollar shock and on the FCI
lp_comp <- function(yvar, xvar) {
  out <- list()
  for (hh in H_SET) {
    s <- series[h == hh]
    setorder(s, ym)
    s[, `:=`(y_l1 = shift(get(yvar), 1), x_l1 = shift(get(xvar), 1))]
    reg <- na.omit(s[, .(y = get(yvar), x = get(xvar), y_l1, x_l1)])
    m <- lm(y ~ x + y_l1 + x_l1, reg)
    ct <- coeftest(m, vcov. = NeweyWest(m, lag = hh + 1, prewhite = FALSE))
    out[[length(out) + 1]] <- data.table(component = yvar, treatment = xvar,
                                         h = hh, b = ct["x", 1], se = ct["x", 2],
                                         p = ct["x", 4], n = nrow(reg))
  }
  rbindlist(out)
}
res <- rbindlist(lapply(c("diff_agg", "within", "between"), function(y)
  rbind(lp_comp(y, "shk"), lp_comp(y, "FCI"))))
cat("\nComponent responses:\n")
print(dcast(res[, .(component, treatment, h, b = round(b, 3), p = round(p, 3))],
            component + treatment ~ h, value.var = c("b", "p")))

write_rev_csv(series[, .(date = ym_date(ym), h, diff_agg, within, between, n_cells)],
              "Rev_ShiftShare_Series.csv")
write_rev_csv(res, "Rev_ShiftShare_USD_Decomposition.csv")

# Figure: component responses to the dollar shock
rp <- res[treatment == "shk"]
rp[, clab := factor(component, levels = c("diff_agg", "within", "between"),
                    labels = c("Aggregate USD - PYG differential",
                               "Within-cell component", "Between-cell composition"))]
pS <- ggplot(rp, aes(factor(h), b, fill = clab)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = b - 1.645 * se, ymax = b + 1.645 * se),
                position = position_dodge(width = 0.8), width = 0.2, color = "grey30") +
  scale_fill_manual(values = c("grey40", "steelblue4", "firebrick")) +
  labs(title = "Reconciling aggregate and micro: shift-share decomposition of the USD-PYG differential",
       subtitle = "Response to a 1 SD dollar shock; within-cell reallocation is positive, between-cell composition negative; 90% bars (NW)",
       x = "Horizon (months)", y = "Response (pp)", fill = NULL) +
  theme_micro()
save_rev_png(pS, "291_ShiftShare_Decomposition.png", w = 10.5, h = 6)

cat(sprintf("=== 46 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
