# ============================================================================
# 41_Falsification_Placebos.R — Structurally insulated placebo outcomes
# (Coworker 3, Concern 1, item 2).
#
# If the balance-sheet/credit channel is the operative mechanism, FCI and
# dollar shocks should NOT predict structurally insulated activities:
#   - Electricidad y agua (binational hydro, contracts in USD, insulated)
#   - Consumo Publico (government consumption)
#   - Agricultura (self-financed offshore; the paper's existing falsification)
# while the "should respond" contrasts are credit and (weakly) PIB.
# Quarterly LPs, h = 1-8 quarters, NW(h+1); regressors: quarterly-averaged
# FCI_COMP_AVG and, separately, the purged-DXY shock (from the micro phase).
# ADDITIVE: writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
cat("=== 41: Falsification placebos (quarterly) ===\n")

qsa <- as.data.table(suppressMessages(read_excel(micro_paths$fci_xlsx,
                                                 sheet = "Quarterly_SA")))
qsa[, qtr := as.Date(cut(as.Date(Date), "quarter"))]
setorder(qsa, qtr)
setnames(qsa, c("Electricidad y agua", "Consumo Público", "Consumo Privado"),
         c("Electricidad", "Consumo_Publico", "Consumo_Privado"))

PLACEBOS <- c(Electricidad = "Electricidad", Consumo_Publico = "Consumo_Publico",
              Agricultura = "Agricultura")
CONTRASTS <- c(PIB = "PIB")
for (v in c(PLACEBOS, CONTRASTS)) qsa[, paste0(v, "_yoy") :=
                                        100 * (get(v) / shift(get(v), 4) - 1)]

# Quarterly regressors: FCI, purged-DXY shock, credit (contrast outcome)
d   <- load_ext_data()
shk <- fread(file.path(micro_paths$out_csv, "Micro_Shock_Series.csv"))
shk[, fecha := as.Date(date)]
dm  <- merge(d[, .(fecha, FCI_COMP_AVG, Cred_Real_Total, IMAEP_yoy, IPC_yoy)],
             shk[, .(fecha, shk_dxy_purged)], by = "fecha", all.x = TRUE)
dq <- dm[, .(FCI = mean(FCI_COMP_AVG, na.rm = TRUE),
             shk_dxy = mean(shk_dxy_purged, na.rm = TRUE),
             Cred = mean(Cred_Real_Total, na.rm = TRUE),
             IPC_yoy = mean(IPC_yoy, na.rm = TRUE)),
         by = .(qtr = as.Date(cut(fecha, "quarter")))]
setorder(dq, qtr)

q <- merge(qsa[, c("qtr", paste0(c(PLACEBOS, CONTRASTS), "_yoy")), with = FALSE],
           dq, by = "qtr")
q[, PIB_ctrl := PIB_yoy]

# Quarterly LP: y_{t+h} ~ x + 2 lags y + lag x + IPC control (NW h+1)
lp_q <- function(y, x, max_h = 8) {
  out <- list()
  for (h in 1:max_h) {
    dh <- as.data.frame(q)
    dh$y_fwd <- dplyr::lead(dh[[y]], h)
    dh$y_l1 <- dplyr::lag(dh[[y]], 1); dh$y_l2 <- dplyr::lag(dh[[y]], 2)
    dh$x_l1 <- dplyr::lag(dh[[x]], 1)
    reg <- na.omit(dh[, c("y_fwd", x, "y_l1", "y_l2", "x_l1", "IPC_yoy")])
    if (nrow(reg) < 25) next
    m <- lm(as.formula(paste("y_fwd ~", x, "+ y_l1 + y_l2 + x_l1 + IPC_yoy")), reg)
    ct <- coeftest(m, vcov. = NeweyWest(m, lag = h + 1, prewhite = FALSE))
    out[[h]] <- data.table(horizon = h, coef = ct[x, 1], se = ct[x, 2],
                           p_value = ct[x, 4], n_obs = nrow(reg))
  }
  rbindlist(out)
}

res <- list()
for (y in c(names(PLACEBOS), names(CONTRASTS), "Cred")) {
  yv <- if (y == "Cred") "Cred" else paste0(y, "_yoy")
  for (x in c("FCI", "shk_dxy")) {
    r <- lp_q(yv, x)
    r[, `:=`(outcome = y, shock = x,
             role = fifelse(y %in% names(PLACEBOS), "placebo", "contrast"))]
    res[[paste(y, x)]] <- r
    k <- r[horizon == 4]
    if (nrow(k)) cat(sprintf("  [%s ~ %s] h=4: b=%7.3f (p=%.3f)\n", y, x, k$coef, k$p_value))
  }
}
res <- rbindlist(res)
write_rev_csv(res, "Rev_Falsification_Placebos.csv")

# Summary at h = 4 quarters (~ 12 months)
summ <- res[horizon == 4, .(outcome, shock, role, coef = round(coef, 3),
                            p = round(p_value, 3))]
write_rev_csv(summ, "Rev_Falsification_Summary_h4.csv")

res[, olab := factor(outcome, levels = c("Electricidad", "Consumo_Publico",
                                         "Agricultura", "PIB", "Cred"),
                     labels = c("Electricity/water (binational)", "Public consumption",
                                "Agriculture", "GDP", "Real credit"))]
res[, slab := ifelse(shock == "FCI", "FCI (comprehensive)", "Purged-DXY shock")]
pF <- ggplot(res, aes(horizon, coef)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = coef - Z90 * se, ymax = coef + Z90 * se,
                  fill = role), alpha = 0.25) +
  geom_line(aes(color = role), linewidth = 0.8) +
  facet_grid(slab ~ olab, scales = "free_y") +
  scale_color_manual(values = c(contrast = "firebrick", placebo = "grey40")) +
  scale_fill_manual(values = c(contrast = "firebrick", placebo = "grey40")) +
  labs(title = "Falsification battery: insulated placebos vs responding outcomes",
       subtitle = "Quarterly LPs, 90% bands (NW). Placebos (grey) should be flat; credit (red) should respond",
       x = "Horizon (quarters)", y = "Coefficient") +
  theme_micro() + theme(legend.position = "none")
save_rev_png(pF, "285_Falsification_Placebos.png", w = 13, h = 6.5)

cat(sprintf("=== 41 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
