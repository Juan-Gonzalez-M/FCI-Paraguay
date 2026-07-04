# ============================================================================
# micro_helpers.R — Shared utilities for the bank-level micro phase (scripts 30-35)
# Implements Bank_Level_Micro_Plan_IREF_Revision.md. NEW PHASE: does not touch
# or depend on modifications to the existing aggregate pipeline (scripts 01-23).
# ============================================================================

suppressMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
})

# ---- Paths -----------------------------------------------------------------
# Resolve project root robustly: helpers live in <root>/R/ and scripts are run
# from R/ (same convention as the aggregate pipeline).
micro_paths <- local({
  root <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  if (!file.exists(file.path(root, "data", "FCI_data_1.xlsx"))) {
    root <- normalizePath(getwd(), mustWork = FALSE)  # fallback if run from root
  }
  list(
    root     = root,
    bank_xlsx = file.path(root, "data", "Boletin_Bancos_May2026.xlsx"),
    fci_xlsx  = file.path(root, "data", "FCI_data_1.xlsx"),
    fci_csv   = file.path(root, "output", "csv", "FCI_Complete_Results.csv"),
    out       = file.path(root, "output", "micro"),
    out_csv   = file.path(root, "output", "micro", "csv"),
    out_png   = file.path(root, "output", "micro", "png"),
    out_rds   = file.path(root, "output", "micro", "rds")
  )
})
invisible(lapply(micro_paths[c("out", "out_csv", "out_png", "out_rds")],
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# ---- Month index -----------------------------------------------------------
# ym: integer month index (year*12 + month) so leads/lags are exact integers
to_ym   <- function(d) year(as.IDate(d)) * 12L + month(as.IDate(d))
ym_date <- function(ym) as.Date(sprintf("%d-%02d-01", (ym - 1L) %/% 12L, (ym - 1L) %% 12L + 1L))

# ---- Constants (plan §1.3, §4.1) -------------------------------------------
CUR_FC  <- 6200  # foreign currency (verified: FX share 52% in 2016 -> 44% in 2026)
CUR_PYG <- 6900  # guaranies

SECT_HEDGED    <- c("AGRICULTURA", "GANADERIA")
SECT_UNHEDGED  <- c("COMERCIO AL POR MENOR", "SERVICIOS", "CONSTRUCCION",
                    "CONSUMO", "VIVIENDA", "ACTIVIDADES INMOBILIARIAS")
SECT_AMBIGUOUS <- c("INDUSTRIA", "COMERCIO AL POR MAYOR", "SECTOR FINANCIERO")
SECT_EXCLUDED  <- c("ADMINISTRACION PUBLICA", "OTROS")

COVID_START <- to_ym("2020-03-01"); COVID_END <- to_ym("2021-06-01")

HORIZONS <- 1:18
H_KEY    <- c(3, 6, 12)   # main-text horizons per plan §6.1

# Cell-size filter (plan §1.4): balances are in millions of PYG (verified in
# script 30 reconciliation), so 1 bn PYG = 1000.
CELL_MIN_BASE <- 1000     # 1 bn PYG
CELL_MIN_ALT  <- c(5000, 10000)

WINSOR_P <- c(0.01, 0.99)

# ---- Small utilities --------------------------------------------------------
winsorize <- function(x, p = WINSOR_P) {
  q <- quantile(x, p, na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

zstd <- function(x) as.numeric(scale(x))

# Long-difference growth over horizon h on an adjusted log level (plan §1.4):
# y_h(t) = 100 * (lval_{t+h} - lval_t), computed within groups on a complete
# ym grid so gaps never produce silent mis-leads.
add_lead_growth <- function(dt, lvar, by, horizons = HORIZONS, prefix = "g") {
  setorderv(dt, c(by, "ym"))
  for (h in horizons) {
    dt[, paste0(prefix, h) := {
      idx <- match(ym + h, ym)
      100 * (get(lvar)[idx] - get(lvar))
    }, by = by]
  }
  invisible(dt)
}

add_lag_growth <- function(dt, lvar, by, h = 12, name = "g_lag12") {
  setorderv(dt, c(by, "ym"))
  dt[, (name) := {
    idx <- match(ym - h, ym)
    100 * (get(lvar) - get(lvar)[idx])
  }, by = by]
  invisible(dt)
}

# ---- Inference -----------------------------------------------------------
# ~17-18 bank clusters => wild cluster bootstrap is the mandated baseline
# (plan §2.2). Implemented natively (restricted WCB-t, Rademacher weights,
# Cameron-Gelbach-Miller 2008 / MacKinnon-Webb): does not depend on
# fwildclusterboot compiling, and handles multi-way high-dimensional FEs.

WCB_B <- as.integer(Sys.getenv("MICRO_WCB_B", "999"))   # override for quick runs

# Restricted wild cluster bootstrap-t p-value for H0: coefficient on `param`= 0.
#   dt      : estimation data (rows with NAs in used vars are fine - dropped
#             identically in full and restricted fits via obs alignment)
#   yvar    : outcome column name
#   rhs     : character vector of regressors (must contain `param`)
#   fe      : fixed-effects part as string, e.g. "bank^sector^ym + bank^sector^cur"
#   cluster : cluster variable name
wcb_pval <- function(dt, yvar, rhs, fe, param, cluster = "bank", B = WCB_B) {
  out <- tryCatch({
    f_full <- as.formula(paste(yvar, "~", paste(rhs, collapse = " + "), "|", fe))
    m1 <- fixest::feols(f_full, data = dt, cluster = as.formula(paste0("~", cluster)),
                        notes = FALSE)
    t_obs <- fixest::coeftable(m1)[param, 3]
    keep  <- fixest::obs(m1)
    ds    <- dt[keep]
    rhs_r <- setdiff(rhs, param)
    f_res <- as.formula(paste(yvar, "~",
                              if (length(rhs_r)) paste(rhs_r, collapse = " + ") else "1",
                              "|", fe))
    m0    <- fixest::feols(f_res, data = ds, notes = FALSE)
    fit0  <- fitted(m0); e0 <- resid(m0)
    cl    <- ds[[cluster]]; ucl <- unique(cl); G <- length(ucl)
    cl_i  <- match(cl, ucl)
    fcl   <- as.formula(paste0("~", cluster))
    tb    <- numeric(B)
    for (b in seq_len(B)) {
      w  <- sample(c(-1, 1), G, replace = TRUE)[cl_i]
      ds[, .wcb_y := fit0 + e0 * w]
      mb <- fixest::feols(as.formula(paste(".wcb_y ~", paste(rhs, collapse = " + "), "|", fe)),
                          data = ds, cluster = fcl, notes = FALSE)
      tb[b] <- fixest::coeftable(mb)[param, 3]
    }
    mean(abs(tb) >= abs(t_obs))
  }, error = function(e) NA_real_)
  out
}

# Extract coefficient row (estimate, analytic-cluster se/p) from a fixest model
coef_row <- function(model, param) {
  ct <- tryCatch(fixest::coeftable(model), error = function(e) NULL)
  if (is.null(ct) || !param %in% rownames(ct)) {
    return(list(b = NA_real_, se = NA_real_, p = NA_real_, n = NA_integer_))
  }
  list(b = ct[param, 1], se = ct[param, 2], p = ct[param, 4],
       n = model$nobs)
}

# Analytic inference battery for one coefficient (plan §2.2, §5.9):
# cluster(bank) + two-way(bank,time) + Driscoll-Kraay(bandwidth h+1)
inference_battery <- function(model, param, dk_lag = 13) {
  base <- coef_row(model, param)
  se2  <- tryCatch(sqrt(diag(vcov(model, vcov = ~ bank + ym)))[param],
                   error = function(e) NA_real_)
  seDK <- tryCatch(sqrt(diag(vcov(model,
                    vcov = as.formula(sprintf("DK(%d) ~ ym", dk_lag)))))[param],
                   error = function(e) NA_real_)
  p2   <- if (is.na(se2))  NA_real_ else 2 * pnorm(-abs(base$b / se2))
  pDK  <- if (is.na(seDK)) NA_real_ else 2 * pnorm(-abs(base$b / seDK))
  list(b = base$b, se_bank = base$se, p_bank = base$p,
       se_twoway = se2, p_twoway = p2, se_dk = seDK, p_dk = pDK, n = base$n)
}

# ---- Plot theme -------------------------------------------------------------
theme_micro <- function() {
  theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(color = "grey35", size = 10),
          legend.position = "bottom")
}

save_png <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(micro_paths$out_png, name), p, width = w, height = h, dpi = 150)
  cat("  [png] ", name, "\n")
}

write_micro_csv <- function(dt, name) {
  fwrite(as.data.table(dt), file.path(micro_paths$out_csv, name))
  cat("  [csv] ", name, "\n")
}

read_rds_micro  <- function(name) {
  obj <- readRDS(file.path(micro_paths$out_rds, name))
  # restore data.table over-allocation so := works by reference after readRDS
  if (data.table::is.data.table(obj)) obj <- data.table::setDT(data.table::copy(obj))
  obj
}
write_rds_micro <- function(obj, name) {
  saveRDS(obj, file.path(micro_paths$out_rds, name))
  cat("  [rds] ", name, "\n")
}

# Dynamic beta_h figure with 90% bands
plot_dynamic_beta <- function(res, title, subtitle, fname,
                              bcol = "b", secol = "se_bank") {
  res <- as.data.table(res)
  res[, lo := get(bcol) - 1.645 * get(secol)]
  res[, hi := get(bcol) + 1.645 * get(secol)]
  p <- ggplot(res, aes(x = h, y = get(bcol))) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
    geom_line(color = "steelblue4", linewidth = 1) +
    geom_point(color = "steelblue4", size = 1.6) +
    labs(title = title, subtitle = subtitle,
         x = "Horizon (months)", y = expression(beta[h] ~ "(pp)")) +
    theme_micro()
  save_png(p, fname)
}
