# ============================================================================
# revision_helpers.R — Shared utilities for the revision-extras phase
# (scripts 37-44). Implements coworker-advised exercises. ADDITIVE: reuses
# micro_helpers.R; writes only to output/revision/.
# ============================================================================

source("micro_helpers.R")   # data.table, readxl, ggplot2, paths, to_ym, etc.
suppressMessages({library(sandwich); library(lmtest); library(ivreg)})

rev_paths <- list(
  out      = file.path(micro_paths$root, "output", "revision"),
  out_csv  = file.path(micro_paths$root, "output", "revision", "csv"),
  out_png  = file.path(micro_paths$root, "output", "revision", "png"),
  external = file.path(micro_paths$root, "output", "revision", "external"),
  ext_csv  = file.path(micro_paths$root, "output", "csv", "New_External_Variables.csv"),
  lp_credit_csv = file.path(micro_paths$root, "output", "csv", "LP_Credit_Standard.csv"),
  expwin_csv = file.path(micro_paths$root, "output", "csv", "Expanding_Window_First_Stage.csv")
)
invisible(lapply(rev_paths[c("out", "out_csv", "out_png", "external")],
                 dir.create, recursive = TRUE, showWarnings = FALSE))

write_rev_csv <- function(dt, name) {
  fwrite(as.data.table(dt), file.path(rev_paths$out_csv, name))
  cat("  [csv] ", name, "\n")
}
save_rev_png <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(rev_paths$out_png, name), p, width = w, height = h, dpi = 150)
  cat("  [png] ", name, "\n")
}

Z90 <- qnorm(0.95)

# ---- Data loaders -----------------------------------------------------------
# ext_data: script 17's consolidated monthly dataset (FCIs, DXY, controls,
# Cred_Real_Total) — the same input scripts 18-23 use.
load_ext_data <- function() {
  d <- fread(rev_paths$ext_csv)
  d[, fecha := as.Date(fecha)]
  setorder(d, fecha)
  d[, `:=`(IMAEP_yoy_L1 = shift(IMAEP_yoy, 1), IMAEP_yoy_L2 = shift(IMAEP_yoy, 2),
           IPC_yoy_L1 = shift(IPC_yoy, 1), IPC_yoy_L2 = shift(IPC_yoy, 2))]
  d
}

# Script 23's exogenous-control convention (macro controls + ER controls)
EXOG_MACRO <- c("IMAEP_yoy", "IMAEP_yoy_L1", "IMAEP_yoy_L2",
                "IPC_yoy", "IPC_yoy_L1", "IPC_yoy_L2")
ER_CTRLS   <- c("d_ToT", "US_10Y", "SP500")

# ---- FRED fetch with cache --------------------------------------------------
fred_get <- function(id) {
  cache <- file.path(rev_paths$external, paste0(id, ".csv"))
  if (!file.exists(cache)) {
    url <- paste0("https://fred.stlouisfed.org/graph/fredgraph.csv?id=", id)
    ok <- tryCatch({
      download.file(url, cache, quiet = TRUE, method = "libcurl"); TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok || !file.exists(cache) || file.size(cache) < 100) {
      if (file.exists(cache)) file.remove(cache)
      cat("  [fred]  FAILED:", id, "\n"); return(NULL)
    }
    cat("  [fred]  downloaded:", id, "\n")
  } else cat("  [fred]  cached:", id, "\n")
  d <- fread(cache)
  setnames(d, c("date", "value"))
  d[, date := as.Date(date)]
  d[, value := suppressWarnings(as.numeric(value))]
  d[!is.na(value)]
}

# Monthly average of a daily/weekly FRED series
fred_monthly <- function(d, name) {
  if (is.null(d)) return(NULL)
  m <- d[, .(v = mean(value, na.rm = TRUE)),
         by = .(fecha = as.Date(format(date, "%Y-%m-01")))]
  setnames(m, "v", name)
  setorder(m, fecha)
  m
}

# ---- LP engine (script 05/23 convention, verified in script 36) --------------
# y_{t+h} ~ x_t + y_lag1 + y_lag2 + x_lag1 + controls(+2 lags each); NW(h+1).
# lag_augment: adds y_lag3 + x_lag2 + 3rd control lag and uses HC1 SEs
# (Montiel Olea & Plagborg-Moller 2021).
run_lp <- function(data, y_var, x_var, max_h = 24, controls = c("IMAEP_yoy", "IPC_yoy"),
                   lag_augment = FALSE, min_n = 30, subset_expr = NULL) {
  data <- as.data.frame(data)
  n_lags <- if (lag_augment) 3 else 2
  out <- list()
  for (h in 1:max_h) {
    dh <- data
    dh$y_fwd  <- dplyr::lead(dh[[y_var]], h)
    for (i in 1:n_lags) dh[[paste0("y_lag", i)]] <- dplyr::lag(dh[[y_var]], i)
    dh$x_lag1 <- dplyr::lag(dh[[x_var]], 1)
    if (lag_augment) dh$x_lag2 <- dplyr::lag(dh[[x_var]], 2)
    rhs <- c(x_var, paste0("y_lag", 1:n_lags), "x_lag1", if (lag_augment) "x_lag2")
    for (cv in controls) {
      if (!cv %in% names(dh)) next
      rhs <- c(rhs, cv)
      for (j in 1:n_lags) {
        nm <- paste0(cv, "_lag", j)
        dh[[nm]] <- dplyr::lag(dh[[cv]], j)
        rhs <- c(rhs, nm)
      }
    }
    if (!is.null(subset_expr)) dh <- dh[eval(subset_expr, dh), ]
    reg <- na.omit(dh[, c("y_fwd", rhs)])
    if (nrow(reg) < min_n) next
    m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = " + "))), data = reg)
    V <- if (lag_augment) sandwich::vcovHC(m, type = "HC1")
         else sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE)
    ct <- lmtest::coeftest(m, vcov = V)
    out[[h]] <- data.frame(horizon = h, coef = ct[x_var, 1], se = ct[x_var, 2],
                           p_value = ct[x_var, 4],
                           ci_lower = ct[x_var, 1] - Z90 * ct[x_var, 2],
                           ci_upper = ct[x_var, 1] + Z90 * ct[x_var, 2],
                           n_obs = nrow(reg))
  }
  rbindlist(out)
}

# ---- IV-LP engine (script 23 convention) -------------------------------------
# 2SLS of y_{t+h} on endo, instrumented by `instr`; exog = y lags, endo lag,
# macro + ER controls. Returns coef, NW inference, conventional first-stage F,
# MOP effective F (single instrument: HAC-robust squared t of the instrument
# in the first stage), and the AR 90% confidence set by grid inversion.
# gamma: assumed direct effect of the instrument on the outcome (CHR 2012);
# outcome is replaced by y - gamma*Z before estimation.
iv_lp_h <- function(data, y_var, endo, instr, h,
                    er_controls = ER_CTRLS, gamma = 0,
                    ar_grid = seq(-250, 100, by = 0.5), min_n = 50) {
  dh <- as.data.frame(data)
  dh$y_fwd  <- dplyr::lead(dh[[y_var]], h)
  dh$y_lag1 <- dplyr::lag(dh[[y_var]], 1)
  dh$y_lag2 <- dplyr::lag(dh[[y_var]], 2)
  dh$fci_lag1 <- dplyr::lag(dh[[endo]], 1)
  exog <- c("y_lag1", "y_lag2", "fci_lag1", EXOG_MACRO,
            intersect(er_controls, names(dh)))
  reg <- na.omit(dh[, c("y_fwd", endo, exog, instr)])
  if (nrow(reg) < min_n) return(NULL)
  reg$y_adj <- reg$y_fwd - gamma * reg[[instr]]
  exog_str <- paste(exog, collapse = " + ")

  res <- tryCatch({
    ivm <- ivreg::ivreg(as.formula(paste("y_adj ~", endo, "+", exog_str,
                                         "|", instr, "+", exog_str)), data = reg)
    V   <- sandwich::NeweyWest(ivm, lag = h + 1, prewhite = FALSE)
    ct  <- lmtest::coeftest(ivm, vcov = V)
    dg  <- summary(ivm, diagnostics = TRUE)$diagnostics
    fsF <- tryCatch(dg["Weak instruments", "statistic"], error = function(e) NA)

    # First stage + effective F (robust)
    fs  <- lm(as.formula(paste(endo, "~", instr, "+", exog_str)), data = reg)
    tfs <- lmtest::coeftest(fs, vcov = sandwich::NeweyWest(fs, lag = h + 1,
                                                           prewhite = FALSE))
    effF <- tfs[instr, 3]^2
    r2f  <- summary(fs)$r.squared
    fs0  <- lm(as.formula(paste(endo, "~", exog_str)), data = reg)
    pR2  <- 1 - sum(resid(fs)^2) / sum(resid(fs0)^2)

    # AR confidence set (just-identified; grid inversion, homoskedastic F as
    # in script 23)
    ar_p <- vapply(ar_grid, function(b0) {
      yt <- reg$y_adj - b0 * reg[[endo]]
      fr <- lm(as.formula(paste("yt ~", exog_str)), data = reg)
      fu <- lm(as.formula(paste("yt ~", instr, "+", exog_str)), data = reg)
      anova(fr, fu)$`Pr(>F)`[2]
    }, numeric(1))
    acc <- ar_grid[ar_p > 0.10]
    data.frame(horizon = h, gamma = gamma,
               coef = ct[endo, 1], se = ct[endo, 2], p_value = ct[endo, 4],
               ci_lower = ct[endo, 1] - Z90 * ct[endo, 2],
               ci_upper = ct[endo, 1] + Z90 * ct[endo, 2],
               first_stage_F = fsF, eff_F = effF, partial_R2 = pR2,
               ar_lo = if (length(acc)) min(acc) else NA_real_,
               ar_hi = if (length(acc)) max(acc) else NA_real_,
               n_obs = nrow(reg))
  }, error = function(e) NULL)
  res
}

# ---- Rolling-window standardization (script 01 convention: 60m window,
# partial windows allowed for initial observations) ---------------------------
roll_z <- function(x, window = 60, min_obs = 12) {
  n <- length(x); z <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    w <- x[max(1, t - window + 1):t]
    w <- w[!is.na(w)]
    if (length(w) >= min_obs && sd(w) > 0 && !is.na(x[t]))
      z[t] <- (x[t] - mean(w)) / sd(w)
  }
  z
}

expanding_z <- function(x, min_obs = 36) {
  n <- length(x); z <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    w <- x[1:t]; w <- w[!is.na(w)]
    if (length(w) >= min_obs && sd(w) > 0 && !is.na(x[t]))
      z[t] <- (x[t] - mean(w)) / sd(w)
  }
  z
}

# FCI sign conventions (R/CLAUDE.md): positive = tighter
FCI_SIGNS <- c(TPM = 1, Spread_activas_pasivas = 1, Spread_mercado_TPM = 1,
               Crecimiento_creditos = -1, Ratio_Cred_Depo = -1, Morosidad = 1,
               Rentabilidad = -1, Liquidez = -1, TCN = 1, Commodities = -1,
               FFER = 1, VIX = 1)
FCI_QUANTITY_VARS <- c("Crecimiento_creditos", "Ratio_Cred_Depo", "Morosidad",
                       "Rentabilidad", "Liquidez")

# Load the 12 FCI input variables (with credit growth computed as in script 01)
load_fci_inputs <- function() {
  mv <- as.data.table(suppressMessages(read_excel(micro_paths$fci_xlsx,
                                                  sheet = "Main_variables")))
  mv[, fecha := as.Date(Fecha)]
  setorder(mv, fecha)
  mv[, Crecimiento_creditos :=
       (Creditos_Sector_privado_totales / shift(Creditos_Sector_privado_totales, 12) - 1) * 100]
  mv[, c("fecha", names(FCI_SIGNS)), with = FALSE]
}

# Build a z-score-method FCI from a variable subset under a standardization fn
build_zscore_fci <- function(inputs, vars, std_fn = roll_z, lag_vars = NULL) {
  X <- copy(inputs)
  if (!is.null(lag_vars)) {
    for (v in intersect(lag_vars, vars)) X[, (v) := shift(get(v), 1)]
  }
  Z <- sapply(vars, function(v) std_fn(X[[v]]) * FCI_SIGNS[v])
  rowMeans(Z, na.rm = FALSE)
}
