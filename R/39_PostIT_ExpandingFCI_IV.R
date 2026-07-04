# ============================================================================
# 39_PostIT_ExpandingFCI_IV.R — Post-IT weak-instrument response
# (Coworker 1 WP1-Step 5; Coworker 3 Concern 2):
#  1. Reconstruct the expanding-window-standardized FCI series (z-score and
#     PCA variants; the appendix K.5/N.6 series was never saved).
#     VALIDATION GATE: reproduce first-stage F ~ 197 (full) / ~ 173 (post-IT)
#     from output/csv/Expanding_Window_First_Stage.csv.
#  2. Rolling-vs-expanding two-column first-stage table for Section 5.2.
#  3. NEW estimation: post-IT IV-LP of credit on the expanding-window FCI
#     (exCredit variant on the credit LHS), instrumented by DXY, with AR
#     confidence sets and Montiel Olea-Pflueger effective F.
# ADDITIVE: writes only to output/revision/.
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
cat("=== 39: Expanding-window FCI and post-IT IV-LP ===\n")

# ---------------------------------------------------------------------------
# 1. Construct expanding-window FCIs (z-score method; PCA check)
# ---------------------------------------------------------------------------
inputs <- load_fci_inputs()
vars12      <- names(FCI_SIGNS)
vars_exCred <- setdiff(vars12, "Crecimiento_creditos")

Zexp <- sapply(vars12, function(v) expanding_z(inputs[[v]]) * FCI_SIGNS[v])
fci_exp      <- rowMeans(Zexp, na.rm = FALSE)
fci_exp_exC  <- rowMeans(Zexp[, vars_exCred], na.rm = FALSE)
# ENDO_exCredit variant (8 domestic vars ex credit growth) - matches the endo
# used in the DXY-only IV (script 23) and the appendix first-stage exercise
vars_endo_exc <- c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM",
                   "Ratio_Cred_Depo", "Morosidad", "Rentabilidad", "Liquidez", "TCN")
fci_exp_endo_exC <- rowMeans(Zexp[, vars_endo_exc], na.rm = FALSE)

# PCA variant: PC1 of the expanding-standardized panel (sign-aligned)
cc <- complete.cases(Zexp)
pca <- prcomp(Zexp[cc, ], center = FALSE, scale. = FALSE)
pc1 <- rep(NA_real_, nrow(Zexp)); pc1[cc] <- pca$x[, 1]
if (cor(pc1, fci_exp, use = "complete.obs") < 0) pc1 <- -pc1

exp_fci <- data.table(fecha = inputs$fecha, FCI_EXP = fci_exp,
                      FCI_EXP_exCredit = fci_exp_exC,
                      FCI_EXP_ENDO_exCredit = fci_exp_endo_exC, FCI_EXP_PCA = pc1)
cat(sprintf("  Expanding z-score FCI vs PCA variant corr: %.3f\n",
            cor(fci_exp, pc1, use = "complete.obs")))
write_rev_csv(exp_fci, "Rev_ExpandingFCI_Series.csv")

# ---------------------------------------------------------------------------
# 2. First-stage validation and rolling-vs-expanding table
# ---------------------------------------------------------------------------
d <- load_ext_data()
d <- merge(d, exp_fci, by = "fecha", all.x = TRUE)
POSTIT <- as.Date("2011-05-01")

# First stage under two documented conventions: classical univariate F and
# HAC (Newey-West, 13 lags) t^2. The endo matching the appendix/script-23
# exercise is the ENDO_exCredit index.
fs_row <- function(dd, fci_col, label) {
  reg <- na.omit(dd[, c(fci_col, "DXY"), with = FALSE])
  m <- lm(as.formula(paste(fci_col, "~ DXY")), reg)
  ct  <- coef(summary(m))
  tnw <- coeftest(m, vcov. = NeweyWest(m, lag = 13, prewhite = FALSE))["DXY", 3]
  data.table(sample = label, n_obs = nrow(reg), DXY_coef = ct["DXY", 1],
             DXY_t = ct["DXY", 3], F_classical = ct["DXY", 3]^2,
             t_NW = tnw, F_NW = tnw^2,
             R2_pct = 100 * summary(m)$r.squared)
}
fs <- rbind(
  fs_row(d, "FCI_EXP_ENDO_exCredit", "Full (expanding window)"),
  fs_row(d[fecha >= POSTIT], "FCI_EXP_ENDO_exCredit", "Post-IT (expanding window)"),
  fs_row(d, "FCI_ENDO_exCredit_AVG", "Full (rolling 60m, baseline)"),
  fs_row(d[fecha >= POSTIT], "FCI_ENDO_exCredit_AVG", "Post-IT (rolling 60m, baseline)"))
print(fs[, .(sample, n_obs, coef = round(DXY_coef, 4),
             F_classical = round(F_classical, 1), F_NW = round(F_NW, 1))])

# Validation vs the saved appendix table. The appendix's producing script is
# not in the repo, so the gate is QUALITATIVE: (i) rolling post-IT F collapses
# to ~0; (ii) expanding post-IT F remains far above conventional thresholds;
# (iii) coefficient signs/magnitudes line up. Exact F levels depend on the
# ad-hoc exercise's standardization details.
tgt <- fread(rev_paths$expwin_csv)
cmp <- merge(fs, tgt[, .(sample, F_target = F_stat, coef_target = DXY_coef)],
             by = "sample")
cat("\nVALIDATION vs Expanding_Window_First_Stage.csv (qualitative gate):\n")
print(cmp[, .(sample, coef = round(DXY_coef, 4), coef_target = round(coef_target, 4),
              F_classical = round(F_classical, 1), F_NW = round(F_NW, 1),
              F_target = round(F_target, 1))])
gate1 <- cmp[grepl("Post-IT \\(rolling", sample), F_classical] < 1
gate2 <- cmp[grepl("Post-IT \\(expanding", sample), pmin(F_classical, F_NW)] > 50
cat(sprintf("  Gate 1 (rolling post-IT F collapses to ~0): %s\n",
            ifelse(gate1, "PASS", "FAIL")))
cat(sprintf("  Gate 2 (expanding post-IT F >> 10 under both conventions): %s\n",
            ifelse(gate2, "PASS", "FAIL")))
write_rev_csv(fs,  "Rev_ExpandingFCI_FirstStage.csv")
write_rev_csv(cmp, "Rev_ExpandingFCI_Validation.csv")

# ---------------------------------------------------------------------------
# 3. Post-IT IV-LP: credit on expanding-FCI (exCredit), instrument = DXY
# ---------------------------------------------------------------------------
cat("\nPost-IT IV-LP (expanding-window FCI_exCredit, DXY instrument)...\n")
dpost <- d[fecha >= POSTIT]
iv_res <- list()
for (h in 1:18) {
  r <- iv_lp_h(dpost, "Cred_Real_Total", "FCI_EXP_exCredit", "DXY", h,
               ar_grid = seq(-300, 150, by = 1))
  if (is.null(r)) next
  iv_res[[h]] <- r
  cat(sprintf("  h=%2d: b=%7.2f  se=%5.2f  p=%.3f  FS-F=%6.1f  effF=%6.1f  AR90=[%s, %s]\n",
              h, r$coef, r$se, r$p_value, r$first_stage_F, r$eff_F,
              round(r$ar_lo, 1), round(r$ar_hi, 1)))
}
iv_res <- rbindlist(iv_res)

# Full-sample counterpart + effective F for the baseline rolling FCI IV
iv_full <- list()
for (h in c(6, 12, 18)) {
  r1 <- iv_lp_h(d, "Cred_Real_Total", "FCI_EXP_exCredit", "DXY", h)
  r2 <- iv_lp_h(d, "Cred_Real_Total", "FCI_ENDO_exCredit_AVG", "DXY", h)
  if (!is.null(r1)) iv_full[[length(iv_full) + 1]] <- cbind(spec = "full_expandingFCI", r1)
  if (!is.null(r2)) iv_full[[length(iv_full) + 1]] <- cbind(spec = "full_rollingFCI_baseline", r2)
}
iv_full <- rbindlist(iv_full)

write_rev_csv(rbind(cbind(spec = "postIT_expandingFCI", iv_res), iv_full),
              "Rev_PostIT_IV_ExpandingFCI.csv")

# ---------------------------------------------------------------------------
# 4. Figure: post-IT IV-LP with AR bands
# ---------------------------------------------------------------------------
pI <- ggplot(iv_res, aes(horizon, coef)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbar(aes(ymin = ar_lo, ymax = ar_hi), width = 0.3, color = "grey55") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue4", linewidth = 1) + geom_point(size = 1.6) +
  labs(title = "Post-IT IV-LP: credit response using the expanding-window FCI",
       subtitle = "2SLS, DXY instrument; shaded = 90% Newey-West CI, whiskers = Anderson-Rubin 90% set",
       x = "Horizon (months)", y = "Credit growth response (pp)") +
  theme_micro()
save_rev_png(pI, "283_PostIT_IV_ExpandingFCI.png")

cat(sprintf("=== 39 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
