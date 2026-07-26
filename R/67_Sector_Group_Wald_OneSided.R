# ============================================================================
# 67_Sector_Group_Wald_OneSided.R — Round-17 colleague item: apply the paper's
# own timing standard to the sectoral group-equality test.
#
# PROBLEM.  Table 10 Panel B and Appendix F.8 estimate the stacked group-Wald
# test on FCI_COMP_AVG, the four-method composite whose PCA/VAR/DFM parameters
# are estimated on the FULL sample.  That is exactly the look-ahead the paper's
# promotion of a strictly one-sided index is designed to avoid.  The asymmetry
# is visible to any referee who has just been taught the timing standard.
#
# FIX.  Re-estimate the identical stacked system on timing-clean indices:
#   FCI_COMP_AVG        four-method composite, full-sample parameters (GATE:
#                       must reproduce Rev_Sector_Group_Wald.csv exactly)
#   FCI_COMP_ZSCORE     one-sided trailing 60m z-score, SAME 12 variables
#                       (the strict apples-to-apples timing-clean counterpart)
#   FCI_exCredit_ZSCORE one-sided trailing 60m z-score, ex-credit (the paper's
#                       headline predictive index)
#
# Estimator, controls, horizons and inference are copied verbatim from
# 63_Sectoral_Group_Wald.R; only the index column varies.  Index scales differ
# (quarter-end SD ~0.93 vs ~0.46), so raw coefficients are NOT comparable
# across rows; sd_x and diff_per_sd are reported for that reason.  t-statistics
# and p-values are scale-invariant and carry the comparison.
#
# ADDITIVE: writes only output/revision/csv/Rev_Sector_Group_Wald_OneSided.csv.
# ============================================================================

t0 <- Sys.time()
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(lubridate)
  library(sandwich); library(lmtest)
})
cat("=== 67: Sectoral group equality under the one-sided timing standard ===\n")

OUT <- "../output/revision/csv"

INDICES <- c(
  "FCI_COMP_AVG"        = "composite, full-sample parameters (published baseline)",
  "FCI_COMP_ZSCORE"     = "one-sided trailing z-score, same 12 variables",
  "FCI_exCredit_ZSCORE" = "one-sided trailing z-score, ex-credit (headline index)"
)

# --- Data (script 63 conventions, verbatim) --------------------------------
fci <- read.csv("../output/csv/FCI_Complete_Results.csv")
fci$fecha <- as.Date(fci$fecha)
stopifnot(all(names(INDICES) %in% names(fci)))

macro <- read_excel("../data/FCI_data_1.xlsx", sheet = "Datos_macro")
names(macro)[1] <- "fecha"
macro <- macro %>% mutate(fecha = as.Date(fecha)) %>% arrange(fecha) %>%
  mutate(IPC_yoy = (IPC / lag(IPC, 12) - 1) * 100)
ctrl_q <- macro %>% filter(month(fecha) %% 3 == 0) %>%
  transmute(q_date = fecha, IPC_yoy)

qd <- read_excel("../data/FCI_data_1.xlsx", sheet = "Quarterly_SA")
names(qd)[1] <- "fecha"
qd <- qd %>% mutate(fecha = as.Date(fecha)) %>% arrange(fecha)
col_mapping <- c("Agricultura" = "Agricultura",
                 "Ganadería forestal,  pesca y minería" = "Ganaderia_fp",
                 "Manufactura" = "Manufactura",
                 "Electricidad y agua" = "Electricidad_y_agua",
                 "Construcción" = "Construccion",
                 "Servicios" = "Servicios",
                 "PIB" = "PIB")
for (old in names(col_mapping))
  if (old %in% names(qd)) names(qd)[names(qd) == old] <- col_mapping[old]

SECT6 <- c("Agricultura", "Ganaderia_fp", "Manufactura",
           "Electricidad_y_agua", "Construccion", "Servicios")
stopifnot(all(SECT6 %in% names(qd)), "PIB" %in% names(qd))

qd <- qd %>% mutate(Impuestos = PIB - rowSums(across(all_of(SECT6))))
GRP_FIN <- c("Ganaderia_fp", "Construccion", "Servicios")
GRP_INS <- c("Agricultura", "Manufactura", "Electricidad_y_agua", "Impuestos")
qd <- qd %>% mutate(FinDep    = rowSums(across(all_of(GRP_FIN))),
                    Insulated = rowSums(across(all_of(GRP_INS))))
for (v in c("FinDep", "Insulated"))
  qd[[paste0("AG_", v)]] <- (qd[[v]] / dplyr::lag(qd[[v]], 4) - 1) * 100
qd <- qd %>% mutate(q_date = fecha)

# --- Stacked, fully-interacted system (verbatim from 63) -------------------
build <- function(data, h) {
  d <- data %>% mutate(
    IPC_yoy_lag1 = lag(IPC_yoy, 1), IPC_yoy_lag2 = lag(IPC_yoy, 2),
    fci_lag1     = lag(FCI_COMP, 1))
  mk <- function(yv, g) {
    d %>% mutate(y_fwd  = lead(.data[[yv]], h),
                 y_lag1 = lag(.data[[yv]], 1),
                 y_lag2 = lag(.data[[yv]], 2),
                 grp    = g) %>%
      dplyr::select(q_date, grp, y_fwd, FCI_COMP, y_lag1, y_lag2,
                    fci_lag1, IPC_yoy, IPC_yoy_lag1, IPC_yoy_lag2)
  }
  s <- bind_rows(mk("AG_FinDep", "F"), mk("AG_Insulated", "I")) %>% na.omit()
  keep <- s %>% count(q_date) %>% filter(n == 2) %>% pull(q_date)
  s %>% filter(q_date %in% keep) %>% mutate(DI = as.numeric(grp == "I"))
}

rhs <- c("FCI_COMP", "y_lag1", "y_lag2", "fci_lag1",
         "IPC_yoy", "IPC_yoy_lag1", "IPC_yoy_lag2")

group_wald <- function(index_col) {
  fci_q <- fci %>% filter(month(fecha) %% 3 == 0) %>%
    transmute(q_date = fecha, FCI_COMP = .data[[index_col]])
  az <- fci_q %>% inner_join(qd, by = "q_date") %>%
    left_join(ctrl_q, by = "q_date") %>% arrange(q_date)

  out <- data.frame()
  for (h in 1:8) {
    s <- build(az, h)
    if (nrow(s) < 60) next
    f <- as.formula(paste("y_fwd ~ DI +", paste(rhs, collapse = "+"),
                          "+", paste0("DI:", rhs, collapse = "+")))
    m <- lm(f, data = s)
    V <- sandwich::vcovPL(m, cluster = s$q_date, order.by = s$q_date,
                          lag = h + 1, adjust = TRUE)
    bF   <- unname(coef(m)["FCI_COMP"])
    dIF  <- unname(coef(m)["DI:FCI_COMP"])
    bI   <- bF + dIF
    seD  <- sqrt(V["DI:FCI_COMP", "DI:FCI_COMP"])
    diff <- -dIF
    dfree <- length(unique(s$q_date)) - 2
    tstat <- diff / seD
    pval  <- 2 * pt(-abs(tstat), dfree)
    tc    <- qt(0.975, dfree)
    seCL <- sqrt(sandwich::vcovCL(m, cluster = s$q_date)["DI:FCI_COMP",
                                                         "DI:FCI_COMP"])
    # SD of the index on the estimation sample actually used (one obs per date)
    sdx <- sd(s$FCI_COMP[!duplicated(s$q_date)])
    out <- rbind(out, data.frame(
      index = index_col, horizon = h,
      beta_FinDep = bF, beta_Insulated = bI,
      diff_F_minus_I = diff, se = seD, t = tstat, p_value = pval,
      ci_lo = diff - tc * seD, ci_hi = diff + tc * seD,
      se_clusterdate = seCL,
      p_clusterdate = 2 * pt(-abs(diff / seCL), dfree),
      sd_x = sdx, diff_per_sd = diff * sdx, se_per_sd = seD * sdx,
      ci_lo_per_sd = (diff - tc * seD) * sdx,
      ci_hi_per_sd = (diff + tc * seD) * sdx,
      n_stacked = nrow(s), n_quarters = length(unique(s$q_date)), df = dfree))
  }
  out
}

res <- bind_rows(lapply(names(INDICES), group_wald))

# --- GATE: the composite block must reproduce the published table ----------
pub_path <- file.path(OUT, "Rev_Sector_Group_Wald.csv")
if (file.exists(pub_path)) {
  pub <- read.csv(pub_path)
  new <- res %>% filter(index == "FCI_COMP_AVG") %>% arrange(horizon)
  pub <- pub %>% arrange(horizon)
  cmp <- c("beta_FinDep", "beta_Insulated", "diff_F_minus_I", "se",
           "p_value", "n_quarters")
  maxdev <- max(abs(as.matrix(new[, cmp]) - as.matrix(pub[, cmp])))
  cat(sprintf("\nGATE: max |deviation| from Rev_Sector_Group_Wald.csv = %.3e\n",
              maxdev))
  if (maxdev > 1e-8)
    stop("GATE FAILED: composite block does not reproduce the published table.")
  cat("GATE PASSED: composite block reproduces the published table exactly.\n")
} else {
  cat("\nGATE SKIPPED: Rev_Sector_Group_Wald.csv not found.\n")
}

write.csv(res, file.path(OUT, "Rev_Sector_Group_Wald_OneSided.csv"),
          row.names = FALSE)

cat("\n--- Group-equality test (beta_F - beta_I) by index, h = 1..8 ---\n")
for (ix in names(INDICES)) {
  cat(sprintf("\n%s  [%s]\n", ix, INDICES[[ix]]))
  r <- res %>% filter(index == ix)
  cat(sprintf("  SD of index on the quarterly estimation sample: %.4f\n",
              r$sd_x[1]))
  for (i in seq_len(nrow(r)))
    cat(sprintf("  h=%dQ  diff=%7.3f  [%7.3f, %7.3f]  DK p=%.4f  qtr p=%.4f  per-SD=%7.3f  N_q=%d\n",
                r$horizon[i], r$diff_F_minus_I[i], r$ci_lo[i], r$ci_hi[i],
                r$p_value[i], r$p_clusterdate[i], r$diff_per_sd[i],
                r$n_quarters[i]))
}

cat(sprintf("\nDone in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
