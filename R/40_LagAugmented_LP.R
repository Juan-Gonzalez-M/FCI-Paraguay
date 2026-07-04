# ============================================================================
# 40_LagAugmented_LP.R — Lag-augmented local projections
# (Montiel Olea & Plagborg-Moller 2021; Coworker 1 WP6).
#
# Re-estimates the headline credit LPs with lag augmentation (one extra lag of
# the outcome and of the FCI, plus a third control lag) and Eicker-Huber-White
# robust SEs, side by side with the Newey-West baseline. Outcomes: published
# real total credit, FX-adjusted total (constant exchange rate), and USD book
# in dollars (both from the script-36 exercise). Regressor: FCI_exCredit_AVG.
# ADDITIVE: writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
cat("=== 40: Lag-augmented LP (MOP-PM 2021) ===\n")

d <- load_ext_data()

# FX-adjusted outcomes: Total_FXadj from the script-36 output; USD_clean rebuilt
fxadj <- fread(file.path(micro_paths$out_csv, "Agg_FXAdj_Mechanical_Decomposition.csv"))
fxadj[, fecha := as.Date(fecha)]
mv <- as.data.table(suppressMessages(read_excel(micro_paths$fci_xlsx,
                                                sheet = "Main_variables")))
mv[, fecha := as.Date(Fecha)]
setorder(mv, fecha)
mv[, Cred_USD_clean := 100 * (Creditos_Sector_privado_USD /
                                shift(Creditos_Sector_privado_USD, 12) - 1)]
d <- merge(d, fxadj[, .(fecha, Cred_Total_FXadj)], by = "fecha", all.x = TRUE)
d <- merge(d, mv[, .(fecha, Cred_USD_clean)], by = "fecha", all.x = TRUE)

OUTCOMES <- c(Total = "Cred_Real_Total", Total_FXadj = "Cred_Total_FXadj",
              USD_clean = "Cred_USD_clean")
XVAR <- "FCI_exCredit_AVG"

res <- list()
for (oc in names(OUTCOMES)) {
  for (la in c(FALSE, TRUE)) {
    r <- run_lp(d, OUTCOMES[oc], XVAR, max_h = 24,
                controls = c("IMAEP_yoy", "IPC_yoy"), lag_augment = la)
    r[, `:=`(outcome = oc, inference = ifelse(la, "lag-augmented + EHW",
                                              "baseline + Newey-West"))]
    res[[paste(oc, la)]] <- r
    cat(sprintf("  [%s | %s] h=12: b=%.2f (se %.2f, p=%.3f)\n",
                oc, ifelse(la, "lag-aug/EHW", "NW"), r[horizon == 12, coef],
                r[horizon == 12, se], r[horizon == 12, p_value]))
  }
}
res <- rbindlist(res)
write_rev_csv(res, "Rev_LagAugmented_LP.csv")

key <- dcast(res[horizon %in% c(3, 6, 12, 18, 24)],
             outcome + horizon ~ inference, value.var = c("coef", "se", "p_value"))
write_rev_csv(key, "Rev_LagAugmented_LP_KeyHorizons.csv")

pl <- res[outcome == "Total"]
pL <- ggplot(pl, aes(horizon, coef, color = inference, fill = inference)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
  scale_color_manual(values = c("steelblue4", "darkorange3")) +
  scale_fill_manual(values = c("steelblue4", "darkorange3")) +
  labs(title = "Headline credit LP: Newey-West baseline vs lag-augmented (MOP-PM 2021)",
       subtitle = "Real total credit growth on FCI_exCredit; 90% bands",
       x = "Horizon (months)", y = "Coefficient (pp)") +
  theme_micro()
save_rev_png(pL, "284_LagAugmented_LP.png")

cat(sprintf("=== 40 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
