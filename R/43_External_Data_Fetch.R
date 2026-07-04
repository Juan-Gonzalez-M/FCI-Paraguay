# ============================================================================
# 43_External_Data_Fetch.R — Fetch and cache external series (FRED + Atlanta
# Fed), per user-approved scope. Downloads once into output/revision/external/
# (data/ is never touched); builds the monthly merged file consumed by 44.
#
# Series:
#   DTWEXB   Fed nominal broad dollar index, daily, 1995-2019 (discontinued)
#   DTWEXBGS Fed nominal broad dollar index, daily, 2006-     (current)
#            -> spliced at the 2006-2019 overlap (ratio of overlap means)
#   DEXBZUS  BRL per USD, daily (Brazil channel for the enhanced purge)
#   DEXUSEU  USD per EUR, daily, 1999- (placebo instrument)
#   TEDRATE  TED spread, daily, ends Jan 2022 (dollar-funding stress)
#   Wu-Xia shadow federal funds rate (Atlanta Fed xlsx; best-effort, optional)
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
cat("=== 43: External data fetch (FRED, cached) ===\n")

options(timeout = 300)

bd_old <- fred_monthly(fred_get("DTWEXB"),   "BD_old")
bd_new <- fred_monthly(fred_get("DTWEXBGS"), "BD_new")
brl    <- fred_monthly(fred_get("DEXBZUS"),  "BRL")
eur    <- fred_monthly(fred_get("DEXUSEU"),  "EUR")
gbp    <- fred_monthly(fred_get("DEXUSUK"),  "GBP")   # for the dollar-free EUR/GBP cross
ted    <- fred_monthly(fred_get("TEDRATE"),  "TED")

# ---- Splice the broad dollar index ------------------------------------------
broad <- NULL
if (!is.null(bd_old) && !is.null(bd_new)) {
  ov <- merge(bd_old, bd_new, by = "fecha")
  k  <- mean(ov$BD_old / ov$BD_new)
  cat(sprintf("  Broad dollar splice: %d overlap months, ratio %.4f\n", nrow(ov), k))
  broad <- rbind(bd_old[fecha < min(ov$fecha)],
                 bd_new[, .(fecha, BD_old = BD_new * k)])
  setnames(broad, "BD_old", "BroadDollar")
  setorder(broad, fecha)
} else if (!is.null(bd_new)) {
  broad <- setnames(copy(bd_new), "BD_new", "BroadDollar")
  cat("  NOTE: DTWEXB unavailable; broad dollar starts 2006 (DTWEXBGS only)\n")
}

# ---- Wu-Xia shadow rate (best-effort) ----------------------------------------
wuxia <- NULL
wx_cache <- file.path(rev_paths$external, "WuXiaShadowRate.xls")
if (!file.exists(wx_cache)) {
  for (u in c("https://www.atlantafed.org/-/media/documents/datafiles/cqer/research/wu-xia-shadow-federal-funds-rate/WuXiaShadowRate.xls",
              "https://www.atlantafed.org/-/media/documents/datafiles/cqer/research/wu-xia-shadow-federal-funds-rate/WuXiaShadowRate.xlsx")) {
    ok <- tryCatch({download.file(u, wx_cache, quiet = TRUE, mode = "wb"); TRUE},
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && file.exists(wx_cache) && file.size(wx_cache) > 5000) break
    if (file.exists(wx_cache)) file.remove(wx_cache)
  }
}
if (file.exists(wx_cache)) {
  wuxia <- tryCatch({
    w <- as.data.table(suppressMessages(read_excel(wx_cache)))
    setnames(w, 1:2, c("period", "WuXia"))
    # period is typically YYYYMM numeric or a date
    if (is.numeric(w$period)) {
      w[, fecha := as.Date(sprintf("%d-%02d-01", period %/% 100, period %% 100))]
    } else w[, fecha := as.Date(paste0(substr(as.character(period), 1, 7), "-01"))]
    w[!is.na(WuXia) & !is.na(fecha), .(fecha, WuXia = as.numeric(WuXia))]
  }, error = function(e) NULL)
}
cat(if (is.null(wuxia)) "  Wu-Xia shadow rate: unavailable (FFER used as-is downstream)\n"
    else sprintf("  Wu-Xia shadow rate: %d months\n", nrow(wuxia)))

# ---- Merge to a monthly grid -------------------------------------------------
grid <- data.table(fecha = seq(as.Date("1995-01-01"), as.Date("2025-12-01"), by = "month"))
ext <- Reduce(function(a, b) if (is.null(b)) a else merge(a, b, by = "fecha", all.x = TRUE),
              list(grid, broad, brl, eur, gbp, ted, wuxia))
avail <- sapply(setdiff(names(ext), "fecha"),
                function(v) sum(!is.na(ext[[v]])))
cat("  Monthly observations per series:\n"); print(avail)
write_rev_csv(ext, "Rev_External_Series.csv")

cat(sprintf("=== 43 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
