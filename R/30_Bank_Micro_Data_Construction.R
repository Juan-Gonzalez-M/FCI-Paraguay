# ============================================================================
# 30_Bank_Micro_Data_Construction.R
# Phase 0 of the bank-level micro plan (Bank_Level_Micro_Plan_IREF_Revision.md)
#
# Builds, from the BCP Boletin de Bancos workbook (Jan 2016 - May 2026):
#   P1  credit by bank x currency x month                (sheet Carteras)
#   P2  credit by bank x sector x currency x month        (sheet Credito Sector)
#   P3  crop detail by bank x activity x currency x month (sheet Credito Actividad)
#   P4  bank characteristics / exposures                  (sheets EEFF + Ratios)
#   DEP deposits by bank x currency x month               (sheet EEFF)
# plus the shock series (purged DXY innovation etc.), the valuation adjustment
# (plan §1.2), entry/exit & merger protocol (§1.3), and the aggregate
# reconciliation check (§1.2 verification).
#
# NEW PHASE: reads existing project data but modifies nothing outside
# output/micro/.
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
cat("=== 30: Bank-level micro data construction ===\n")

# ============================================================================
# 0. USER-EXTENSION SLOT: TCN and CPI for Jan-May 2026
# ----------------------------------------------------------------------------
# The paper database (FCI_data_1.xlsx) ends Dec 2025 while the bank panel runs
# to May 2026. The valuation adjustment needs end-of-month TCN (PYG/USD) and
# PY CPI for those months. Paste them here when available; until then, FC-book
# outcomes dated after Dec 2025 are dropped (NA), which only shortens the
# usable leads for the last few shock months.
# Example: TCN_2026 <- c(`2026-01` = 7850, `2026-02` = 7900, ...)
TCN_2026 <- c()   # end-of-period PYG/USD, Jan..May 2026
IPC_2026 <- c()   # CPI index (same base as Datos_macro$IPC), Jan..May 2026
# ============================================================================

# ---------------------------------------------------------------------------
# 1. Read bank workbook
# ---------------------------------------------------------------------------
fx  <- micro_paths$bank_xlsx
cat("Reading workbook:", basename(fx), "\n")

cart <- as.data.table(suppressMessages(read_excel(fx, sheet = "Carteras")))
csec <- as.data.table(suppressMessages(read_excel(fx, sheet = "Credito Sector")))
cact <- as.data.table(suppressMessages(read_excel(fx, sheet = "Credito Actividad")))
eeff <- as.data.table(suppressMessages(read_excel(fx, sheet = "EEFF")))
rati <- as.data.table(suppressMessages(read_excel(fx, sheet = "Ratios")))

setnames(cart, c("fecha", "bank", "cuenta", "cur", "importe"))
setnames(csec, c("fecha", "bank", "cur", "sector", "vencida", "vigente"))
setnames(cact, c("fecha", "bank", "actividad", "cur", "vencida", "vigente"))
setnames(eeff, c("fecha", "bank", "rubro", "cur", "importe", "reporte"))
setnames(rati, c("fecha", "bank", "rubro", "total"))

for (d in list(cart, csec, cact, eeff, rati)) d[, ym := to_ym(fecha)]

# ---------------------------------------------------------------------------
# 2. Macro series and shock construction
# ---------------------------------------------------------------------------
cat("Building macro/shock series...\n")
fdat <- micro_paths$fci_xlsx

mv  <- as.data.table(suppressMessages(read_excel(fdat, sheet = "Main_variables")))
dm  <- as.data.table(suppressMessages(read_excel(fdat, sheet = "Datos_macro")))
gfc <- as.data.table(suppressMessages(read_excel(fdat, sheet = "Global_Financial_Conditions")))
setnames(gfc, 1, "Fecha")
setnames(gfc, "S&P 500", "SP500")

macro <- merge(
  mv[, .(ym = to_ym(Fecha), TCN, VIX, FFER)],
  dm[, .(ym = to_ym(Fecha), IPC)], by = "ym", all = TRUE)
macro <- merge(macro,
  gfc[, .(ym = to_ym(Fecha), DXY, SP500, US10Y = US_10Y,
          Selic = Selic_rate, IPE, ToT = Paraguay_ToT)], by = "ym", all = TRUE)

# Append user-supplied 2026 months if provided
if (length(TCN_2026)) {
  add <- data.table(ym = to_ym(paste0(names(TCN_2026), "-01")), TCN = as.numeric(TCN_2026))
  if (length(IPC_2026)) add[, IPC := as.numeric(IPC_2026)[match(ym, to_ym(paste0(names(IPC_2026), "-01")))]]
  macro <- rbind(macro, add, fill = TRUE)
  setorder(macro, ym)
}

# FCI series from the existing pipeline output
fci <- fread(micro_paths$fci_csv, select = c("fecha", "FCI_COMP_AVG", "FCI_ENDO_AVG"))
fci[, ym := to_ym(fecha)]
macro <- merge(macro, fci[, .(ym, FCI = FCI_COMP_AVG)], by = "ym", all.x = TRUE)
setorder(macro, ym)

# Innovations: residuals of monthly changes on 3 own lags (AR pre-whitening)
ar_innov <- function(x, extra = NULL, nlag = 3) {
  X <- sapply(1:nlag, function(l) shift(x, l))
  colnames(X) <- paste0("l", 1:nlag)
  df <- data.frame(y = x, X)
  if (!is.null(extra)) df <- cbind(df, extra)
  r <- rep(NA_real_, length(x))
  ok <- complete.cases(df)
  if (sum(ok) > 30) r[ok] <- resid(lm(y ~ ., data = df[ok, , drop = FALSE]))
  r
}

macro[, dln_dxy   := 100 * c(NA, diff(log(DXY)))]
macro[, dln_sp500 := 100 * c(NA, diff(log(SP500)))]
macro[, dln_ipe   := 100 * c(NA, diff(log(IPE)))]
macro[, dln_tot   := 100 * c(NA, diff(log(ToT)))]
macro[, d_vix     := c(NA, diff(VIX))]
macro[, d_ffer    := c(NA, diff(FFER))]
macro[, d_us10y   := c(NA, diff(US10Y))]

# (a) Raw DXY innovation: AR(3) residual of the monthly log change
macro[, shk_dxy_raw := ar_innov(dln_dxy)]

# (b) Purged DXY innovation (primary shock, plan §2.1): residual after removing
#     contemporaneous trade/commodity/global-financial channels named by R1
macro[, shk_dxy_purged := ar_innov(
  dln_dxy, extra = macro[, .(d_vix, dln_ipe, dln_tot, d_ffer, dln_sp500, d_us10y)])]

# (c) VIX innovation (placebo, plan §2.4)
macro[, shk_vix := ar_innov(d_vix)]

# (d) Aggregate FCI shock ("reduced-form treatment" variant)
macro[, shk_fci := ar_innov(FCI)]

# Standardize all shocks to SD 1 over the estimation window 2016m1-2025m12
est_win <- macro$ym >= to_ym("2016-01-01") & macro$ym <= to_ym("2025-12-01")
for (v in c("shk_dxy_raw", "shk_dxy_purged", "shk_vix", "shk_fci")) {
  mu <- mean(macro[[v]][est_win], na.rm = TRUE)
  sd_ <- sd(macro[[v]][est_win], na.rm = TRUE)
  macro[, (v) := (get(v) - mu) / sd_]
}

shocks <- macro[ym >= to_ym("2015-01-01"),
                .(date = ym_date(ym), ym, TCN, IPC,
                  shk_dxy_raw, shk_dxy_purged, shk_vix, shk_fci)]
write_micro_csv(shocks, "Micro_Shock_Series.csv")

cat(sprintf("  Purged vs raw DXY innovation corr (2016-25): %.3f\n",
            cor(macro[est_win, shk_dxy_raw], macro[est_win, shk_dxy_purged],
                use = "complete.obs")))

# ---------------------------------------------------------------------------
# 3. Panels P1-P3 (credit) and DEP (deposits)
# ---------------------------------------------------------------------------
cat("Assembling panels...\n")

# P1: bank x currency x month, wide over the portfolio accounts of interest
p1_items <- c("Cartera Vigente", "Cartera Vencida", "Renovados",
              "Refinanciados", "Reestructurados")
p1 <- dcast(cart[cuenta %in% p1_items],
            bank + cur + ym ~ cuenta, value.var = "importe", fun.aggregate = sum)
setnames(p1, p1_items, c("vigente", "vencida", "renovados", "refinanciados", "reestructurados"))

# P2: bank x sector x currency x month (main analysis panel)
p2 <- csec[, .(vigente = sum(vigente, na.rm = TRUE),
               vencida = sum(vencida, na.rm = TRUE)), by = .(bank, sector, cur, ym)]

# P3: crop detail
p3 <- cact[, .(vigente = sum(vigente, na.rm = TRUE),
               vencida = sum(vencida, na.rm = TRUE)), by = .(bank, actividad, cur, ym)]

# DEP: deposits by currency from EEFF (BG items)
dep_items <- c("Depósitos a la Vista", "Depósitos a Plazo Fijo",
               "Depósitos CDA", "Depósitos Cta. Cte.")
dep <- eeff[rubro %in% dep_items,
            .(deposits = sum(importe, na.rm = TRUE)), by = .(bank, cur, ym)]

# ---------------------------------------------------------------------------
# 4. Aggregate reconciliation (plan §1.2 verification) + unit check
# ---------------------------------------------------------------------------
cat("Reconciling against the aggregate credit series...\n")
sys_credit <- p1[, .(micro_total = sum(vigente + vencida, na.rm = TRUE)), by = ym]
agg <- merge(sys_credit, dm[, .(ym = to_ym(Fecha), Creditos)], by = "ym")
agg[, ratio := micro_total / Creditos]
agg[, micro_yoy := 100 * (micro_total / shift(micro_total, 12) - 1)]
agg[, agg_yoy   := 100 * (Creditos    / shift(Creditos, 12) - 1)]
rec_cor <- cor(agg$micro_yoy, agg$agg_yoy, use = "complete.obs")
cat(sprintf("  Level ratio micro/aggregate: median %.3f (range %.3f-%.3f)\n",
            median(agg$ratio), min(agg$ratio), max(agg$ratio)))
cat(sprintf("  YoY growth correlation: %.3f (plan expects > 0.98; wedge = financieras not in boletin)\n",
            rec_cor))
write_micro_csv(agg[, .(date = ym_date(ym), micro_total, Creditos, ratio, micro_yoy, agg_yoy)],
                "Micro_Aggregate_Reconciliation.csv")

# ---------------------------------------------------------------------------
# 5. Entry/exit and merger protocol (plan §1.3)
# ---------------------------------------------------------------------------
cat("Entry/exit and merger protocol...\n")
span <- p1[, .(first = min(ym), last = max(ym),
               n_months = uniqueN(ym)), by = bank]
ym_min <- min(p1$ym); ym_max <- max(p1$ym)
span[, `:=`(entrant = first > ym_min, exiter = last < ym_max)]
continuous_banks <- span[first == ym_min & last == ym_max, bank]
cat(sprintf("  %d banks, %d continuous over the full 125 months\n",
            nrow(span), length(continuous_banks)))

# For each exit, identify the absorbing bank (largest total-book jump in the
# exit month or the following one) and classify merger vs wind-down
bank_tot <- p1[, .(tot = sum(vigente + vencida, na.rm = TRUE)), by = .(bank, ym)]
setorder(bank_tot, bank, ym)
bank_tot[, dln := c(NA, diff(log(tot))), by = bank]

events <- list()
for (b in span[exiter == TRUE, bank]) {
  ex_ym <- span[bank == b, last]
  ex_size <- bank_tot[bank == b & ym == ex_ym, tot]
  cand <- bank_tot[ym %in% (ex_ym + 1:2) & bank != b & !is.na(dln)]
  cand <- cand[, .(max_jump = max(dln, na.rm = TRUE)), by = bank][order(-max_jump)]
  absorber <- cand$bank[1]; jump <- cand$max_jump[1]
  # Merger if the absorber's jump is commensurate with the exiting book
  is_merger <- is.finite(jump) && jump > 0.10
  events[[length(events) + 1]] <- data.table(
    exiting_bank = b, exit_ym = ex_ym, exit_date = ym_date(ex_ym),
    exit_book = ex_size, absorber = ifelse(is_merger, absorber, NA_integer_),
    absorber_jump_pct = 100 * jump, classified = ifelse(is_merger, "merger", "wind-down"))
}
events <- rbindlist(events)
entries <- span[entrant == TRUE, .(bank, entry_ym = first, entry_date = ym_date(first))]
print(events); print(entries)
write_micro_csv(events,  "Micro_Exit_Events.csv")
write_micro_csv(entries, "Micro_Entry_Events.csv")

# Baseline drop flags: exiting + absorbing entities in a ±3-month window
# around each exit; entrants excluded for their first 12 months
drop_flags <- rbindlist(list(
  events[, .(bank = exiting_bank, ym_from = exit_ym - 3, ym_to = exit_ym + 3)],
  events[!is.na(absorber), .(bank = absorber, ym_from = exit_ym - 3, ym_to = exit_ym + 3)],
  entries[, .(bank, ym_from = entry_ym, ym_to = entry_ym + 11)]
))

flag_drop <- function(dt) {
  dt[, drop_event := FALSE]
  for (i in seq_len(nrow(drop_flags))) {
    dt[bank == drop_flags$bank[i] & ym >= drop_flags$ym_from[i] &
       ym <= drop_flags$ym_to[i], drop_event := TRUE]
  }
  invisible(dt)
}

# ---------------------------------------------------------------------------
# 6. Valuation adjustment (plan §1.2) — adjusted log levels
#    FC (6200): reported in PYG at current e_t  =>  lval = ln(L / e_t)  [USD]
#    PYG (6900): lval = ln(L / IPC * 100)                              [real]
# ---------------------------------------------------------------------------
cat("Applying valuation adjustment...\n")
fxci <- macro[, .(ym, TCN, IPC)]

add_lval <- function(dt, levelvar) {
  dt <- merge(dt, fxci, by = "ym", all.x = TRUE)
  dt[, lval := fifelse(
        cur == CUR_FC,  log(get(levelvar) / TCN),
        log(get(levelvar) / IPC * 100))]
  dt[get(levelvar) <= 0, lval := NA_real_]
  # Unadjusted variant (nominal PYG) kept for diagnostics/audit
  dt[, lval_nom := fifelse(get(levelvar) > 0, log(get(levelvar)), NA_real_)]
  dt
}

p2[, total := vigente + vencida]
p2 <- add_lval(p2, "total")
p1[, total := vigente + vencida]
p1 <- add_lval(p1, "total")
dep <- add_lval(dep, "deposits")

# NPL ratio (bounded, pp — no logs; plan §3.3)
p2[, npl := 100 * vencida / total]
p2[total <= 0, npl := NA_real_]

# ---------------------------------------------------------------------------
# 7. Bank characteristics / exposures (plan §1.5), P4
# ---------------------------------------------------------------------------
cat("Building exposure variables...\n")

# FC deposit share and total deposits
dep_w <- dcast(dep, bank + ym ~ cur, value.var = "deposits")
setnames(dep_w, c("6200", "6900"), c("dep_fc", "dep_pyg"))
dep_w[, dep_total := rowSums(cbind(dep_fc, dep_pyg), na.rm = TRUE)]
dep_w[, fc_dep_share := dep_fc / dep_total]

# External funding (EEFF BG 'Externo'); share of (Externo + deposits)
extf <- eeff[rubro == "Externo" & reporte == "BG",
             .(externo = sum(importe, na.rm = TRUE)), by = .(bank, ym)]
p4 <- merge(dep_w, extf, by = c("bank", "ym"), all.x = TRUE)
p4[is.na(externo), externo := 0]
p4[, ext_fund_share := externo / (externo + dep_total)]

# FC loan share from P1
loan_w <- dcast(p1, bank + ym ~ cur, value.var = "total")
setnames(loan_w, c("6200", "6900"), c("loan_fc", "loan_pyg"))
loan_w[, fc_loan_share := loan_fc / rowSums(cbind(loan_fc, loan_pyg), na.rm = TRUE)]
p4 <- merge(p4, loan_w[, .(bank, ym, fc_loan_share, loan_fc, loan_pyg)],
            by = c("bank", "ym"), all.x = TRUE)

# Ratios sheet: Tier 1, liquidity, risk-weighted assets (size)
get_ratio <- function(name, newname) {
  r <- rati[rubro == name, .(bank, ym, v = total)]
  setnames(r, "v", newname); r
}
p4 <- Reduce(function(a, b) merge(a, b, by = c("bank", "ym"), all.x = TRUE), list(
  p4,
  get_ratio("Relación entre TIER 1/ACPR", "tier1"),
  get_ratio("Disponible + Inversiones Temporales/Depósitos", "liquidity"),
  get_ratio("ACTIVOS Y CONTINGENTES PONDERADOS", "rwa")))
p4[, size := log(rwa)]   # ln risk-weighted assets (total assets not in bulletin)

# 12-month lags (baseline timing) and 2016 fixed averages (robustness)
setorder(p4, bank, ym)
expo_vars <- c("fc_dep_share", "ext_fund_share", "fc_loan_share", "tier1", "liquidity", "size")
for (v in expo_vars) {
  p4[, paste0(v, "_l12") := {idx <- match(ym - 12, ym); get(v)[idx]}, by = bank]
  p4[, paste0(v, "_2016") := mean(get(v)[ym <= to_ym("2016-12-01")], na.rm = TRUE), by = bank]
}

# Standardized versions (mean 0, SD 1 across bank-months in estimation window)
for (v in paste0(rep(expo_vars, 2), rep(c("_l12", "_2016"), each = length(expo_vars)))) {
  p4[, paste0(v, "_z") := zstd(get(v))]
}

# ---------------------------------------------------------------------------
# 8. Filters, flags, and sample descriptives
# ---------------------------------------------------------------------------
flag_drop(p1); flag_drop(p2); flag_drop(dep)

# Dual-currency cell coverage (plan verified fact: ~85%)
cell_cov <- p2[total > CELL_MIN_BASE,
               .(n_cur = uniqueN(cur)), by = .(bank, sector, ym)]
cat(sprintf("  Cells > 1 bn PYG: %d; share dual-currency: %.1f%%\n",
            nrow(cell_cov), 100 * mean(cell_cov$n_cur == 2)))

# Sector FX shares (latest month) — should match plan §0 (Agri 0.89, ...)
fxsh <- p2[ym == max(ym), .(tot = sum(total)), by = .(sector, cur)]
fxsh <- dcast(fxsh, sector ~ cur, value.var = "tot")
setnames(fxsh, c("6200", "6900"), c("fc", "pyg"))
fxsh[, fx_share := fc / (fc + pyg)]
setorder(fxsh, -fx_share)
print(fxsh[, .(sector, fx_share = round(fx_share, 2))])
write_micro_csv(fxsh, "Micro_Sector_FX_Shares.csv")

# System FX share over time (52% -> 44% check)
sys_fx <- p1[, .(tot = sum(total)), by = .(cur, ym)]
sys_fx <- dcast(sys_fx, ym ~ cur, value.var = "tot")
setnames(sys_fx, c("6200", "6900"), c("fc", "pyg"))
sys_fx[, fx_share := fc / (fc + pyg)]
setorder(sys_fx, ym)
fx_first <- first(sys_fx$fx_share); fx_last <- last(sys_fx$fx_share)
cat(sprintf("  System FX share: %.3f (2016m1) -> %.3f (last)\n", fx_first, fx_last))

# ---------------------------------------------------------------------------
# 9. Save
# ---------------------------------------------------------------------------
meta <- list(events = events, entries = entries, drop_flags = drop_flags,
             continuous_banks = continuous_banks, span = span,
             rec_cor = rec_cor, built = Sys.time())

write_rds_micro(macro, "micro_macro.rds")
write_rds_micro(p1,    "micro_p1_carteras.rds")
write_rds_micro(p2,    "micro_p2_sector.rds")
write_rds_micro(p3,    "micro_p3_actividad.rds")
write_rds_micro(p4,    "micro_p4_bankchars.rds")
write_rds_micro(dep,   "micro_deposits.rds")
write_rds_micro(meta,  "micro_meta.rds")

samp <- data.table(
  item = c("banks", "continuous_banks", "months", "p2_rows", "dual_currency_share",
           "reconciliation_yoy_corr", "fx_share_first", "fx_share_last",
           "tcn_2026_supplied"),
  value = c(nrow(span), length(continuous_banks), uniqueN(p2$ym), nrow(p2),
            round(mean(cell_cov$n_cur == 2), 3), round(rec_cor, 3),
            round(fx_first, 3), round(fx_last, 3), length(TCN_2026)))
write_micro_csv(samp, "Micro_Sample_Info.csv")

cat(sprintf("=== 30 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
