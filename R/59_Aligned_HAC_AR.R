# ============================================================================
# 59_Aligned_HAC_AR.R — HAC Anderson-Rubin sets for the aligned one-sided
# ex-credit ex-TCN domestic index (script 56 object), matching main-text
# Table 7's HAC-inversion convention (script 49 Part 2).  Full sample,
# h = 6, 12, 18, 24; deep-negative probes verify boundedness.
# ADDITIVE: writes output/revision/csv/Rev_Aligned_HAC_AR.csv.
# ============================================================================
source("revision_helpers.R")
cat("=== 59: HAC AR sets, aligned index ===\n")

d <- load_ext_data()
inputs <- load_fci_inputs()
VARS_OS <- c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM",
             "Ratio_Cred_Depo", "Morosidad", "Rentabilidad", "Liquidez")
aux <- data.table(fecha = inputs$fecha,
                  FCI_OS_EXP = build_zscore_fci(inputs, VARS_OS, std_fn = expanding_z))
d <- merge(d, aux, by = "fecha", all.x = TRUE)
setorder(d, fecha)

ar_p_hac <- function(rd, xs, d0, h) {
  rd$y_t <- rd$y_fwd - d0 * rd$FCI_OS_EXP
  m <- lm(as.formula(paste("y_t ~ DXY +", xs)), data = rd)
  lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1,
                                                 prewhite = FALSE))["DXY", 4]
}

out <- data.frame()
for (h in c(6, 12, 18, 24)) {
  dh <- as.data.frame(d)
  dh$y_fwd  <- dplyr::lead(dh$Cred_Real_Total, h)
  dh$y_lag1 <- dplyr::lag(dh$Cred_Real_Total, 1)
  dh$y_lag2 <- dplyr::lag(dh$Cred_Real_Total, 2)
  dh$fci_lag1 <- dplyr::lag(dh$FCI_OS_EXP, 1)
  ex <- c("y_lag1", "y_lag2", "fci_lag1", EXOG_MACRO, intersect(ER_CTRLS, names(dh)))
  rd <- na.omit(dh[, c("y_fwd", "FCI_OS_EXP", ex, "DXY")])
  xs <- paste(ex, collapse = " + ")
  grid <- seq(-800, 200, by = 1)
  p_h <- vapply(grid, function(d0) ar_p_hac(rd, xs, d0, h), numeric(1))
  acc <- grid[p_h > 0.10]
  lo <- if (length(acc)) min(acc) else NA; up <- if (length(acc)) max(acc) else NA
  probe <- vapply(c(-1500, -5000), function(d0) ar_p_hac(rd, xs, d0, h), numeric(1))
  unbounded <- all(probe > 0.10) && !is.na(lo) && lo <= -799
  out <- rbind(out, data.frame(horizon = h, n = nrow(rd),
    ar0_p_hac = ar_p_hac(rd, xs, 0, h),
    ar_hac_lo = ifelse(unbounded, -Inf, lo), ar_hac_up = up,
    unbounded_below = unbounded))
  cat(sprintf("h=%2d: AR p(0)=%.4f  90%% set [%s, %s]%s\n", h,
              out$ar0_p_hac[nrow(out)], ifelse(unbounded, "-Inf", as.character(lo)),
              as.character(up), ifelse(unbounded, " UNBOUNDED", "")))
}
write_rev_csv(out, "Rev_Aligned_HAC_AR.csv")
