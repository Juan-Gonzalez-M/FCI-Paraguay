# ============================================================================
# 47_FXAdj_PostIT_Asymmetric_Check.R — Does the valuation adjustment hurt the
# paper ANYWHERE? Completes the FX-adjustment audit beyond the full-sample
# credit LPs of script 36:
#   (a) POST-IT subsample (2011m5+): measured vs constant-exchange-rate
#       outcomes (Total real, USD) on FCI_exCredit and FCI_ENDO_exCredit.
#   (b) ASYMMETRIC LP (tightening vs easing, script-05 convention):
#       measured vs adjusted, full sample and post-IT.
# Expected direction: depreciation coincides with tightening, so the
# valuation component inflates measured credit growth in tightening episodes,
# attenuating the measured tightening coefficient. Adjustment should
# strengthen (or leave unchanged) the published pattern. Verify, not assume.
# ADDITIVE: writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
cat("=== 47: FX-adjustment audit — post-IT and asymmetric ===\n")

d <- load_ext_data()
fxadj <- fread(file.path(micro_paths$out_csv, "Agg_FXAdj_Mechanical_Decomposition.csv"))
fxadj[, fecha := as.Date(fecha)]
mv <- as.data.table(suppressMessages(read_excel(micro_paths$fci_xlsx,
                                                sheet = "Main_variables")))
mv[, fecha := as.Date(Fecha)]
setorder(mv, fecha)
mv[, `:=`(Cred_USD_meas = 100 * (Creditos_Sector_privado_USD_equivalente /
                              shift(Creditos_Sector_privado_USD_equivalente, 12) - 1),
          Cred_USD_clean = 100 * (Creditos_Sector_privado_USD /
                              shift(Creditos_Sector_privado_USD, 12) - 1))]
d <- merge(d, fxadj[, .(fecha, Cred_Total_FXadj)], by = "fecha", all.x = TRUE)
d <- merge(d, mv[, .(fecha, Cred_USD_meas, Cred_USD_clean)], by = "fecha", all.x = TRUE)
POSTIT <- as.Date("2011-05-01")

OUTCOMES <- list(
  c(y = "Cred_Real_Total",  lab = "Total (measured, published)"),
  c(y = "Cred_Total_FXadj", lab = "Total (FX-adjusted)"),
  c(y = "Cred_USD_meas",    lab = "USD (measured PYG-valued, published)"),
  c(y = "Cred_USD_clean",   lab = "USD (in dollars)"))
FCIS <- c("FCI_exCredit_AVG", "FCI_ENDO_exCredit_AVG")

# ---------------------------------------------------------------------------
# (a) Post-IT subsample, standard LP (script-05/36 spec)
# ---------------------------------------------------------------------------
cat("\n(a) Post-IT subsample LPs...\n")
resA <- list()
for (f in FCIS) for (oc in OUTCOMES) {
  r <- run_lp(d[fecha >= POSTIT], oc["y"], f, max_h = 18,
              controls = c("IMAEP_yoy", "IPC_yoy"))
  if (!nrow(r)) next
  r[, `:=`(sample = "post-IT", fci = f, outcome = oc["lab"])]
  resA[[paste(f, oc["y"])]] <- r
  k <- r[horizon == 12]
  cat(sprintf("  [%s | %s] h=12: b=%7.2f (se %.2f, p=%.3f, n=%d)\n",
              sub("_AVG", "", f), oc["lab"], k$coef, k$se, k$p_value, k$n_obs))
}
resA <- rbindlist(resA)

# ---------------------------------------------------------------------------
# (b) Asymmetric LP (FCI_pos / FCI_neg, script-05 convention), h = 12
# ---------------------------------------------------------------------------
cat("\n(b) Asymmetric LPs (h = 12)...\n")
run_asym <- function(data, yvar, fvar, h = 12) {
  dd <- as.data.frame(data)
  dd$pos <- pmax(dd[[fvar]], 0); dd$neg <- abs(pmin(dd[[fvar]], 0))
  dd$y_fwd <- dplyr::lead(dd[[yvar]], h)
  dd$y_l1 <- dplyr::lag(dd[[yvar]], 1); dd$y_l2 <- dplyr::lag(dd[[yvar]], 2)
  dd$f_l1 <- dplyr::lag(dd[[fvar]], 1)
  ctrl <- c("IMAEP_yoy", "IPC_yoy")
  for (cv in ctrl) for (j in 1:2) dd[[paste0(cv, "_l", j)]] <- dplyr::lag(dd[[cv]], j)
  cols <- c("y_fwd", "pos", "neg", "y_l1", "y_l2", "f_l1", ctrl,
            paste0(rep(ctrl, each = 2), "_l", 1:2))
  reg <- na.omit(dd[, cols])
  if (nrow(reg) < 40) return(NULL)
  m <- lm(y_fwd ~ ., reg)
  ct <- coeftest(m, vcov. = NeweyWest(m, lag = h + 1, prewhite = FALSE))
  wald <- tryCatch({   # H0: pos = -(-neg) i.e. symmetric effect
    L <- rep(0, length(coef(m))); names(L) <- names(coef(m))
    L["pos"] <- 1; L["neg"] <- 1
    v <- as.numeric(t(L) %*% NeweyWest(m, lag = h + 1, prewhite = FALSE) %*% L)
    z <- as.numeric(L %*% coef(m)) / sqrt(v)
    2 * pnorm(-abs(z))
  }, error = function(e) NA_real_)
  data.table(b_tight = ct["pos", 1], p_tight = ct["pos", 4],
             b_ease = -ct["neg", 1], p_ease = ct["neg", 4],
             wald_p_sym = wald, n = nrow(reg))
}
resB <- list()
for (smp in c("full", "post-IT")) for (oc in OUTCOMES[c(1, 2)]) {
  ds <- if (smp == "full") d else d[fecha >= POSTIT]
  r <- run_asym(ds, oc["y"], "FCI_exCredit_AVG")
  if (is.null(r)) next
  r[, `:=`(sample = smp, outcome = oc["lab"])]
  resB[[paste(smp, oc["y"])]] <- r
  cat(sprintf("  [%s | %s] tightening %7.2f (p=%.3f)  easing %7.2f (p=%.3f)  Wald sym p=%.3f\n",
              smp, oc["lab"], r$b_tight, r$p_tight, r$b_ease, r$p_ease, r$wald_p_sym))
}
resB <- rbindlist(resB)

write_rev_csv(resA, "Rev_FXAdj_PostIT_LP.csv")
write_rev_csv(resB, "Rev_FXAdj_Asymmetric_LP.csv")

# Verdict lines
cmp <- dcast(resA[horizon %in% c(6, 12, 18) & fci == "FCI_ENDO_exCredit_AVG"],
             horizon ~ outcome, value.var = "coef")
cat("\nPost-IT comparison (FCI_ENDO_exCredit):\n"); print(cmp)
cat(sprintf("=== 47 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
