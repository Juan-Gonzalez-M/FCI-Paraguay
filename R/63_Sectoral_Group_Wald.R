# ============================================================================
# 63_Sectoral_Group_Wald.R — Round-14 colleague item: genuine coefficient-
# equality test between the two pre-classified sector groups.
#
# PROBLEM.  Script 57 reports an "equality test" that is in fact a separate
# local projection on the growth DIFFERENTIAL (AG_FinDep - AG_Insulated),
# whose controls are lags of the differential.  That is a different estimand
# from beta_F - beta_I: each group LP conditions on lags of its OWN dependent
# variable, so the differential-LP coefficient does not equal the difference
# of the two group coefficients (-0.56 vs -1.91-(-1.16) = -0.75 at h=2Q;
# -2.28 vs -2.52 at h=4Q).  It is therefore not a Wald test of equal responses.
#
# FIX.  Estimate the two group equations jointly as a fully-interacted stacked
# system on a COMMON sample.  Full interaction of every regressor with the
# group dummy makes the stacked point estimates numerically identical to the
# two separate OLS regressions, while delivering the cross-equation covariance
# needed to test H0: beta_F = beta_I.  Inference is Driscoll-Kraay
# (cross-sectionally correlated + serially correlated errors), lag h+1,
# matching the NW(h+1) convention used elsewhere.
#
# Specification otherwise identical to script 57 / main-text Table 10.
# ADDITIVE: writes only output/revision/csv/Rev_Sector_Group_Wald.csv.
# ============================================================================

t0 <- Sys.time()
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(lubridate)
  library(sandwich); library(lmtest)
})
cat("=== 61: Sectoral group coefficient-equality (stacked Wald) ===\n")

OUT <- "../output/revision/csv"

# --- Data (script 57 conventions, verbatim) --------------------------------
fci <- read.csv("../output/csv/FCI_Complete_Results.csv")
fci$fecha <- as.Date(fci$fecha)
fci_q <- fci %>% filter(month(fecha) %% 3 == 0) %>%
  transmute(q_date = fecha, FCI_COMP = FCI_COMP_AVG)

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

az <- fci_q %>% inner_join(qd, by = "q_date") %>%
  left_join(ctrl_q, by = "q_date") %>% arrange(q_date)

# --- Stacked, fully-interacted system --------------------------------------
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
  # common sample: keep only dates present for BOTH groups
  keep <- s %>% count(q_date) %>% filter(n == 2) %>% pull(q_date)
  s %>% filter(q_date %in% keep) %>% mutate(DI = as.numeric(grp == "I"))
}

rhs <- c("FCI_COMP", "y_lag1", "y_lag2", "fci_lag1",
         "IPC_yoy", "IPC_yoy_lag1", "IPC_yoy_lag2")
res <- data.frame()
for (h in 1:8) {
  s <- build(az, h)
  if (nrow(s) < 60) next
  f <- as.formula(paste("y_fwd ~ DI +", paste(rhs, collapse = "+"),
                        "+", paste0("DI:", rhs, collapse = "+")))
  m <- lm(f, data = s)
  # Driscoll-Kraay: scores aggregated within each date (cross-equation
  # dependence), Bartlett-weighted across dates (serial correlation).
  # NOTE: order.by must be supplied explicitly; without it vcovPL degenerates.
  V <- sandwich::vcovPL(m, cluster = s$q_date, order.by = s$q_date,
                        lag = h + 1, adjust = TRUE)
  bF   <- unname(coef(m)["FCI_COMP"])                 # finance-dependent
  dIF  <- unname(coef(m)["DI:FCI_COMP"])              # beta_I - beta_F
  bI   <- bF + dIF
  seD  <- sqrt(V["DI:FCI_COMP", "DI:FCI_COMP"])
  diff <- -dIF                                        # beta_F - beta_I
  dfree <- length(unique(s$q_date)) - 2               # DK dof (n clusters)
  tstat <- diff / seD
  pval  <- 2 * pt(-abs(tstat), dfree)
  tc    <- qt(0.975, dfree)
  # cross-check: cluster-robust (by date) without the serial-correlation kernel
  seCL <- sqrt(sandwich::vcovCL(m, cluster = s$q_date)["DI:FCI_COMP",
                                                       "DI:FCI_COMP"])
  res <- rbind(res, data.frame(
    horizon = h, beta_FinDep = bF, beta_Insulated = bI,
    diff_F_minus_I = diff, se = seD, t = tstat, p_value = pval,
    ci_lo = diff - tc * seD, ci_hi = diff + tc * seD,
    se_clusterdate = seCL,
    p_clusterdate = 2 * pt(-abs(diff / seCL), dfree),
    n_stacked = nrow(s), n_quarters = length(unique(s$q_date)), df = dfree))
}

print(format(res, digits = 4))
write.csv(res, file.path(OUT, "Rev_Sector_Group_Wald.csv"), row.names = FALSE)

cat("\n--- Headline horizons ---\n")
for (h in c(2, 4)) {
  r <- res[res$horizon == h, ]
  cat(sprintf("h=%dQ: beta_F=%.2f  beta_I=%.2f  diff=%.2f  [%.2f, %.2f]  p=%.4f  (N_q=%d)\n",
              h, r$beta_FinDep, r$beta_Insulated, r$diff_F_minus_I,
              r$ci_lo, r$ci_hi, r$p_value, r$n_quarters))
}
cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
