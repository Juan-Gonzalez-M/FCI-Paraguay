# ============================================================================
# 45_COVID_Reclassification_Check.R — The COVID-reclassification variant
# (Colleague 1's "one remaining mandatory exercise").
#
# Concern: during 2020-22 the bulletin moved reprogrammed loans into special
# categories ("Medida Excepcional COVID 19 - Vigente/Vencida", "Medidas
# transitorias") that sit OUTSIDE the vigente+vencida stocks the micro panel
# uses — and the aggregate-vs-micro gap loads on TCN, i.e., the missing chunk
# is currency-asymmetric. If reclassification removed balances asymmetrically
# by currency, it could mechanically generate a within-cell currency
# differential concentrated in the COVID window.
#
# Data constraint: the COVID categories exist in the Carteras sheet at
# bank x currency level ONLY (Credito Sector has no COVID categories), so an
# exact bank x sector x currency rebuild is impossible. Three-part design:
#   (a) SIGN THE BIAS: currency composition and size of COVID balances by
#       month — which currency book is missing more?
#   (b) BANK x CURRENCY TEST: Design A-style LP of the within-bank currency
#       differential (bank x time + bank x currency FE) on the baseline book
#       vs the reclassification-INCLUSIVE book (COVID categories added back).
#   (c) GRADIENT SENSITIVITY: apportion each bank x currency COVID balance to
#       sectors proportionally to that bank-currency's sector composition,
#       rebuild cell stocks, re-run the hedging gradient and sharp triple.
# ADDITIVE: writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
library(fixest)
cat("=== 45: COVID-reclassification check ===\n")

macro <- read_rds_micro("micro_macro.rds")
p1    <- read_rds_micro("micro_p1_carteras.rds")     # carries drop_event
p2    <- read_rds_micro("micro_p2_sector.rds")
p4    <- read_rds_micro("micro_p4_bankchars.rds")

# ---------------------------------------------------------------------------
# 0. COVID-category balances from the Carteras sheet (bank x currency x month)
# ---------------------------------------------------------------------------
cart <- as.data.table(suppressMessages(read_excel(micro_paths$bank_xlsx,
                                                  sheet = "Carteras")))
setnames(cart, c("fecha", "bank", "cuenta", "cur", "importe"))
cart[, ym := to_ym(fecha)]
COVID_CTAS <- c("Medida Excepcional COVID 19 - Vigente",
                "Medida Excepcional COVID 19 - Vencida",
                "Medidas transitorias")
cov <- cart[cuenta %in% COVID_CTAS,
            .(covid_bal = sum(importe, na.rm = TRUE)), by = .(bank, cur, ym)]
cat(sprintf("  COVID-category balances: %d bank x currency x month cells, %s to %s\n",
            nrow(cov), format(ym_date(min(cov$ym))), format(ym_date(max(cov$ym)))))

# ---------------------------------------------------------------------------
# (a) Sign the bias: size and currency composition of the missing balances
# ---------------------------------------------------------------------------
base_bc <- p1[, .(bank, cur, ym, book = total, drop_event)]
bc <- merge(base_bc, cov, by = c("bank", "cur", "ym"), all.x = TRUE)
bc[is.na(covid_bal), covid_bal := 0]
bc[, book_incl := book + covid_bal]

comp <- bc[, .(book = sum(book), covid = sum(covid_bal)), by = .(cur, ym)]
comp[, share_pct := 100 * covid / (book + covid)]
comp_w <- dcast(comp, ym ~ cur, value.var = "share_pct")
setnames(comp_w, c("6200", "6900"), c("share_FC", "share_PYG"))
peak <- comp_w[which.max(pmax(share_FC, share_PYG, na.rm = TRUE))]
cat(sprintf("  Peak COVID-category share of book: FC %.1f%% / PYG %.1f%% (at %s)\n",
            peak$share_FC, peak$share_PYG, format(ym_date(peak$ym))))
tot_cov <- bc[covid_bal > 0, .(covid = sum(covid_bal)), by = cur]
tot_cov[, pct := 100 * covid / sum(covid)]
cat(sprintf("  Cumulative COVID balances by currency: FC %.1f%% vs PYG %.1f%%\n",
            tot_cov[cur == CUR_FC, pct], tot_cov[cur == CUR_PYG, pct]))
write_rev_csv(comp_w[, .(date = ym_date(ym), share_FC, share_PYG)],
              "Rev_COVID_Reclass_Composition.csv")

# ---------------------------------------------------------------------------
# (b) Bank x currency Design A: baseline book vs reclassification-inclusive
# ---------------------------------------------------------------------------
cat("\n(b) Bank x currency LP: baseline vs add-back...\n")
fxv <- macro[, .(ym, TCN, IPC, shk_dxy_purged)]

run_bc <- function(bookvar, label) {
  d <- merge(bc, fxv, by = "ym", all.x = TRUE)
  d[, lval := fifelse(cur == CUR_FC, log(get(bookvar) / TCN),
                      log(get(bookvar) / IPC * 100))]
  d[get(bookvar) <= 0, lval := NA_real_]
  add_lead_growth(d, "lval", by = c("bank", "cur"))
  for (h in HORIZONS) d[, paste0("g", h) := winsorize(get(paste0("g", h))), by = cur]
  d[, usd := as.integer(cur == CUR_FC)]
  d[, shk_usd := shk_dxy_purged * usd]
  d[, `:=`(bt = paste(bank, ym), bcur = paste(bank, cur))]
  d <- d[ym >= to_ym("2016-01-01") & ym <= to_ym("2025-12-01") &
         drop_event == FALSE & !is.na(shk_dxy_purged)]
  res <- list()
  for (h in HORIZONS) {
    yv <- paste0("g", h)
    m <- tryCatch(feols(as.formula(paste(yv, "~ shk_usd | bt + bcur")),
                        data = d, cluster = ~bank, notes = FALSE),
                  error = function(e) NULL)
    if (is.null(m)) next
    cr <- coef_row(m, "shk_usd")
    pw <- if (h %in% c(6, 12, 18)) wcb_pval(d, yv, "shk_usd", "bt + bcur", "shk_usd")
          else NA_real_
    res[[h]] <- data.table(book = label, h = h, b = cr$b, se_bank = cr$se,
                           p_bank = cr$p, p_wcb = pw, n = cr$n)
  }
  rbindlist(res)
}
res_bc <- rbind(run_bc("book", "baseline_excl_COVID_cats"),
                run_bc("book_incl", "inclusive_addback"))
cmp_bc <- dcast(res_bc, h ~ book, value.var = c("b", "p_bank", "p_wcb"))
cat("  Bank x currency differential at h = 6/12/18 (baseline vs add-back):\n")
print(res_bc[h %in% c(6, 12, 18),
             .(book, h, b = round(b, 3), p_bank = round(p_bank, 3),
               p_wcb = round(p_wcb, 3))])
write_rev_csv(res_bc, "Rev_COVID_Reclass_BankCur_LP.csv")

# ---------------------------------------------------------------------------
# (c) Gradient sensitivity under proportional sector apportionment
# ---------------------------------------------------------------------------
cat("\n(c) Gradient with COVID balances apportioned to sectors...\n")
p2x <- copy(p2)[, .(bank, sector, cur, ym, total, TCN, IPC, drop_event)]
p2x[, sec_share := total / sum(total), by = .(bank, cur, ym)]
p2x <- merge(p2x, cov, by = c("bank", "cur", "ym"), all.x = TRUE)
p2x[is.na(covid_bal), covid_bal := 0]
p2x[, total_incl := total + covid_bal * sec_share]

build_grad_panel <- function(totvar) {
  d <- copy(p2x)
  d[, lval := fifelse(cur == CUR_FC, log(get(totvar) / TCN),
                      log(get(totvar) / IPC * 100))]
  d[get(totvar) <= 0, lval := NA_real_]
  add_lead_growth(d, "lval", by = c("bank", "sector", "cur"))
  add_lag_growth(d, "lval", by = c("bank", "sector", "cur"))
  setorderv(d, c("bank", "sector", "cur", "ym"))
  for (h in c(6, 12, 18)) {
    d[, paste0("ok", h) := {
      idx <- match(ym + h, ym)
      get(totvar) > CELL_MIN_BASE & !is.na(get(totvar)[idx]) &
        get(totvar)[idx] > CELL_MIN_BASE
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
  d[, n_cur := uniqueN(cur[get(totvar) > CELL_MIN_BASE]), by = .(bank, sector, ym)]
  d <- d[n_cur == 2]
  d[, `:=`(bst = paste(bank, sector, ym), bsc = paste(bank, sector, cur))]
  # hedging treatments (as in scripts 33/37)
  fx16 <- p2[ym <= to_ym("2016-12-01"), .(tot = sum(total)), by = .(sector, cur)]
  fx16 <- dcast(fx16, sector ~ cur, value.var = "tot")
  setnames(fx16, c("6200", "6900"), c("fc", "pyg"))
  fx16[, fx_share_2016 := fc / (fc + pyg)]
  d <- merge(d, fx16[, .(sector, fx_share_2016)], by = "sector", all.x = TRUE)
  d[, fxsh_z := zstd(fx_share_2016)]
  d[, shk_usd_fxsh := shk_usd * fxsh_z]
  d[, hedge_class := fifelse(sector %in% SECT_HEDGED, "hedged",
                     fifelse(sector %in% SECT_UNHEDGED, "unhedged",
                     fifelse(sector %in% SECT_AMBIGUOUS, "ambiguous", "excluded")))]
  d[, shk_usd_unh := shk_usd * as.integer(hedge_class == "unhedged")]
  d
}

CTRLS <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
FE_A  <- "bst + bsc"
grad_res <- list()
for (tv in c("total", "total_incl")) {
  lab <- ifelse(tv == "total", "baseline_excl_COVID_cats", "inclusive_apportioned")
  dg <- build_grad_panel(tv)
  for (h in c(6, 12, 18)) {
    yv <- paste0("g", h)
    # gradient
    ds <- dg[get(paste0("ok", h)) == TRUE & hedge_class != "excluded"]
    rhs <- c("shk_usd", "shk_usd_fxsh", CTRLS)
    m <- feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A)),
               data = ds, cluster = ~bank, notes = FALSE)
    cr <- coef_row(m, "shk_usd_fxsh")
    pw <- wcb_pval(ds, yv, rhs, FE_A, "shk_usd_fxsh")
    grad_res[[length(grad_res) + 1]] <- data.table(
      book = lab, test = "gradient", h = h, b = cr$b, se_bank = cr$se,
      p_bank = cr$p, p_wcb = pw, n = cr$n)
    # sharp triple
    ds2 <- dg[get(paste0("ok", h)) == TRUE &
              hedge_class %in% c("hedged", "unhedged") &
              (!sector %in% c("CONSUMO", "VIVIENDA") | get(tv) > 5000)]
    rhs2 <- c("shk_usd", "shk_usd_unh", CTRLS)
    m2 <- feols(as.formula(paste(yv, "~", paste(rhs2, collapse = "+"), "|", FE_A)),
                data = ds2, cluster = ~bank, notes = FALSE)
    cr2 <- coef_row(m2, "shk_usd_unh")
    grad_res[[length(grad_res) + 1]] <- data.table(
      book = lab, test = "sharp_triple", h = h, b = cr2$b, se_bank = cr2$se,
      p_bank = cr2$p, p_wcb = NA_real_, n = cr2$n)
    cat(sprintf("  [%s] h=%2d: gradient=%.3f (p_wcb=%.3f)  triple=%.3f (p=%.3f)\n",
                lab, h, cr$b, pw, cr2$b, cr2$p))
  }
}
grad_res <- rbindlist(grad_res)
write_rev_csv(grad_res, "Rev_COVID_Reclass_Gradient.csv")

# ---------------------------------------------------------------------------
# Figures and verdict
# ---------------------------------------------------------------------------
pC <- ggplot(melt(comp_w, id.vars = "ym", variable.name = "book",
                  value.name = "share")[!is.na(share)],
             aes(ym_date(ym), share, color = book)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(share_FC = "firebrick", share_PYG = "steelblue4"),
                     labels = c("Foreign currency book", "Guarani book")) +
  labs(title = "COVID-reclassified balances as a share of each currency book",
       subtitle = "Medida Excepcional COVID-19 + Medidas transitorias, share of (book + reclassified)",
       x = NULL, y = "% of book", color = NULL) +
  theme_micro()
save_rev_png(pC, "289_COVID_Reclass_Composition.png")

g12 <- grad_res[test == "gradient"]
pG <- ggplot(g12, aes(factor(h), b, color = book)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_pointrange(aes(ymin = b - 1.645 * se_bank, ymax = b + 1.645 * se_bank),
                  position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("steelblue4", "darkorange3"),
                     labels = c("Baseline (excl. COVID categories)",
                                "Inclusive (apportioned add-back)")) +
  labs(title = "Hedging gradient: baseline vs COVID-reclassification-inclusive book",
       subtitle = "Shock x USD x sector FX share; 90% intervals, cluster by bank",
       x = "Horizon (months)", y = "Gradient coefficient (pp)", color = NULL) +
  theme_micro()
save_rev_png(pG, "290_COVID_Reclass_Gradient.png")

b0 <- grad_res[test == "gradient" & book == "baseline_excl_COVID_cats" & h == 12, b]
b1 <- grad_res[test == "gradient" & book == "inclusive_apportioned" & h == 12, b]
cat(sprintf("\nVERDICT: gradient at h=12 baseline %.3f vs inclusive %.3f (change %.0f%%)\n",
            b0, b1, 100 * (b1 - b0) / abs(b0)))

cat(sprintf("=== 45 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
