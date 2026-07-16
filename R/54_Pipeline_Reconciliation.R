# ============================================================================
# 54_Pipeline_Reconciliation.R — Canonical-pipeline reconciliation (round 9).
#
# Colleague item: closely related post-IT headline estimates differ
# (-5.60 / -5.31 / -4.13 / -4.00 total credit; -12.84 / -7.24 / -6.72 USD).
# This script establishes that the differences come from the INDEX object
# (post-IT re-estimated vs full-sample-constructed) and the index VARIANT
# (ENDO_exCredit vs exCredit), NOT from different outcome data or controls:
#
#   Part 1  outcome-series identity: the credit outcome used by scripts 05/15
#           (YoY of Creditos_deflactados from FCI_data_1.xlsx) vs the outcome
#           used by scripts 17/47/52 (Cred_Real_Total in
#           New_External_Variables.csv) -- correlation, mean absolute
#           difference, max discrepancy on the common sample.
#   Part 2  machinery validation + label audit: the canonical LP spec
#           (y lags 1-2, FCI lag 1, IMAEP_yoy and IPC_yoy each with 2 lags,
#           NW(h+1)) is run with candidate index columns to identify which
#           index reproduces the published full-sample h=12 estimate
#           (-7.65628, LP_Credit_Standard.csv, labeled "FCI_COMP") -- this
#           settles the Appendix E.0c index-label question.
#   Part 3  matched panel (the promoted one-sided primary + composite
#           comparator under ONE pipeline): full sample and post-IT,
#           h = 6/12/18, per-unit, per-SD, 95% CI, N.
#   Part 4  attribution table: every post-IT total-credit headline variant
#           with its index object, variant, and pipeline stated -- the
#           reconciliation exhibit for Appendix E.3/O.8.
#
# Canonical declaration implemented here (round-9 decision): the strictly
# one-sided trailing z-score index (full-sample-constructed, hence complete
# 60-month windows at every post-IT date; no estimated parameters, no future
# information) is the PRIMARY causal-timing specification; the four-method
# composite is the descriptive/historical index and a measurement-combination
# robustness check; the post-IT re-estimated composite of script 15 is the
# regime-specific measurement variant.
# Outputs: Rev_Canonical_Pipeline.csv, Rev_OneSided_Matched.csv,
#          Rev_Outcome_Reconciliation.csv (output/revision/csv/)
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(sandwich); library(lmtest); library(readxl) })

POSTIT <- as.Date("2011-05-01")

d <- read.csv("../output/csv/New_External_Variables.csv")
d$fecha <- as.Date(d$fecha)
comp <- read.csv("../output/csv/FCI_Complete_Results.csv"); comp$fecha <- as.Date(comp$fecha)
rv <- read.csv("../output/csv/FCI_Robustness_Versions.csv"); rv$fecha <- as.Date(rv$fecha)
d <- d %>%
  left_join(comp %>% dplyr::select(fecha, FCI_COMP_ZSCORE), by = "fecha") %>%
  left_join(rv %>% dplyr::select(fecha, FCI_exCredit_ZSCORE), by = "fecha") %>%
  arrange(fecha)

cat("=== 54: canonical-pipeline reconciliation ===\n\n")

# ---------------------------------------------------------------------------
# Part 1 — outcome-series identity
# ---------------------------------------------------------------------------
mac <- read_excel("../data/FCI_data_1.xlsx", sheet = "Datos_macro")
fcol <- names(mac)[1]
mac <- mac %>% rename(fecha = !!sym(fcol)) %>% mutate(fecha = as.Date(fecha)) %>%
  arrange(fecha) %>%
  mutate(Cred_xlsx = (Creditos_deflactados / lag(Creditos_deflactados, 12) - 1) * 100)
cmp <- d %>% dplyr::select(fecha, Cred_Real_Total) %>%
  inner_join(mac %>% dplyr::select(fecha, Cred_xlsx), by = "fecha") %>% na.omit()
rec <- data.frame(
  n = nrow(cmp),
  correlation = cor(cmp$Cred_Real_Total, cmp$Cred_xlsx),
  mean_abs_diff = mean(abs(cmp$Cred_Real_Total - cmp$Cred_xlsx)),
  max_abs_diff = max(abs(cmp$Cred_Real_Total - cmp$Cred_xlsx)))
cat("--- Part 1: outcome-series identity (scripts 05/15 vs 17/47/52) ---\n")
print(rec, row.names = FALSE)
write.csv(rec, "../output/revision/csv/Rev_Outcome_Reconciliation.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Canonical LP machinery (script 05/15 convention: y lags 1-2, FCI lag 1,
# IMAEP_yoy + 2 lags, IPC_yoy + 2 lags, NW(h+1))
# ---------------------------------------------------------------------------
lp_can <- function(data, yvar, xvar, h) {
  dh <- data %>% mutate(
    y_fwd = lead(.data[[yvar]], h),
    y_lag1 = lag(.data[[yvar]], 1), y_lag2 = lag(.data[[yvar]], 2),
    fci_lag1 = lag(.data[[xvar]], 1),
    IMAEP_lag1 = lag(IMAEP_yoy, 1), IMAEP_lag2 = lag(IMAEP_yoy, 2),
    IPC_lag1 = lag(IPC_yoy, 1), IPC_lag2 = lag(IPC_yoy, 2))
  rhs <- c(xvar, "y_lag1", "y_lag2", "fci_lag1",
           "IMAEP_yoy", "IMAEP_lag1", "IMAEP_lag2",
           "IPC_yoy", "IPC_lag1", "IPC_lag2")
  rd <- dh %>% dplyr::select(y_fwd, all_of(rhs)) %>% na.omit()
  if (nrow(rd) < 30) return(NULL)
  m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
  ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE))
  list(b = unname(ct[xvar, 1]), se = unname(ct[xvar, 2]), p = unname(ct[xvar, 4]),
       n = nrow(rd), sd_x = sd(rd[[xvar]]))
}

# ---------------------------------------------------------------------------
# Part 2 — machinery validation + index-label audit
# ---------------------------------------------------------------------------
cat("\n--- Part 2: which index reproduces the published -7.65628? ---\n")
pub <- -7.65628048038672
for (v in c("FCI_exCredit_AVG", "FCI_COMP_AVG")) {
  r <- lp_can(d, "Cred_Real_Total", v, 12)
  cat(sprintf("  %-20s h=12: b=%.5f  (published -7.65628: %s)\n", v, r$b,
              ifelse(abs(r$b - pub) < 1e-3, "MATCH", "no")))
}

# ---------------------------------------------------------------------------
# Part 3 — matched panel under the canonical pipeline
# ---------------------------------------------------------------------------
cat("\n--- Part 3: matched panel (one-sided primary vs composite comparator) ---\n")
z95 <- qnorm(0.975)
panel <- data.frame()
specs <- c(one_sided_exCredit = "FCI_exCredit_ZSCORE",
           one_sided_comprehensive = "FCI_COMP_ZSCORE",
           composite_exCredit = "FCI_exCredit_AVG",
           composite_ENDO_exCredit = "FCI_ENDO_exCredit_AVG")
for (smp in c("full", "postIT")) {
  ds <- if (smp == "postIT") d %>% filter(fecha >= POSTIT) else d
  for (v in names(specs)) for (h in c(6, 12, 18)) {
    r <- lp_can(ds, "Cred_Real_Total", unname(specs[v]), h)
    if (is.null(r)) next
    panel <- rbind(panel, data.frame(sample = smp, index = v, column = unname(specs[v]),
      horizon = h, b = r$b, se = r$se, p = r$p, n = r$n, sd_x = r$sd_x,
      b_per_sd = r$b * r$sd_x,
      ci95_lo = r$b - z95 * r$se, ci95_hi = r$b + z95 * r$se))
  }
}
print(panel %>% mutate(across(where(is.numeric), ~round(., 3))), row.names = FALSE)
write.csv(panel, "../output/revision/csv/Rev_OneSided_Matched.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Part 4 — attribution table for the post-IT total-credit variants
# ---------------------------------------------------------------------------
attribution <- data.frame(
  estimate = c(-5.60, -5.31,
               round(panel$b[panel$sample=="postIT" & panel$index=="composite_ENDO_exCredit" & panel$horizon==12], 2),
               round(panel$b[panel$sample=="postIT" & panel$index=="composite_exCredit" & panel$horizon==12], 2),
               round(panel$b[panel$sample=="postIT" & panel$index=="one_sided_exCredit" & panel$horizon==12], 2)),
  index_object = c("post-IT re-estimated composite", "post-IT re-estimated composite",
                   "full-sample-constructed composite", "full-sample-constructed composite",
                   "one-sided trailing z-score (primary)"),
  variant = c("ENDO_exCredit (8 vars)", "exCredit (11 vars)",
              "ENDO_exCredit (9->8 vars)", "exCredit (11 vars)",
              "exCredit (11 vars)"),
  source = c("script 15 / Table 4", "script 15 / App E.3",
             "script 54 (canonical spec)", "script 54 (canonical spec)",
             "script 54 (canonical spec, PRIMARY)"),
  outcome_and_controls = "identical (Part 1; canonical spec)")
cat("\n--- Part 4: post-IT total-credit attribution ---\n")
print(attribution, row.names = FALSE)
write.csv(attribution, "../output/revision/csv/Rev_Canonical_Pipeline.csv", row.names = FALSE)
cat("\nSaved: Rev_Outcome_Reconciliation.csv, Rev_OneSided_Matched.csv, Rev_Canonical_Pipeline.csv\n")
