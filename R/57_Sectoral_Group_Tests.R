# ============================================================================
# 57_Sectoral_Group_Tests.R — Round-11 colleague item E3/E4.
#
# Both round-11 audits object that the "segmented economy" claim rests on
# sector-by-sector significance counting (non-rejection is not insulation)
# with no group-level or multiplicity-adjusted inference.  Required:
#   (1) PRE-DEFINED group aggregates built from published value-added
#       levels (additive at constant 2014 prices):
#         finance-dependent = Livestock+ + Construction + Services (~56% GDP)
#         insulated         = Agriculture + Manufacturing + Electricity
#                             + Taxes residual (~44% GDP)
#       Group definitions follow the paper's Section 2.3 classification
#       (colleague 2's suggested composition), fixed before estimation.
#   (2) Equality test between the two groups: LP on the growth DIFFERENCE
#       (finance-dependent minus insulated), same specification; the
#       coefficient is the response differential and its NW p-value tests
#       H0: equal responses directly.
#   (3) Family-level multiplicity: Benjamini-Hochberg FDR q-values over the
#       sector x horizon family (6 production sectors x h = 1..8Q).
#   (4) Aggregation check: GDP-share-weighted sum of sector responses vs
#       the directly estimated aggregate GDP response.
#
# Specification: identical to script 27 / main-text Table 10 (total
# effects): FCI_COMP (end-of-quarter), controls IPC_yoy (+2 lags), 2 lags
# of dep. var, 1 FCI lag, NW(h+1).  Full sample.
#
# ADDITIVE: writes only output/revision/csv/Rev_Sector_Group*.csv.
# ============================================================================

t0 <- Sys.time()
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(lubridate)
  library(sandwich); library(lmtest); library(data.table)
})
cat("=== 57: Sectoral group-level tests ===\n")

OUT <- "../output/revision/csv"

# ---------------------------------------------------------------------------
# Data (script 27 conventions)
# ---------------------------------------------------------------------------
fci <- read.csv("../output/csv/FCI_Complete_Results.csv")
fci$fecha <- as.Date(fci$fecha)
fci_q <- fci %>%
  filter(month(fecha) %% 3 == 0) %>%
  transmute(q_date = fecha, FCI_COMP = FCI_COMP_AVG)

macro <- read_excel("../data/FCI_data_1.xlsx", sheet = "Datos_macro")
names(macro)[1] <- "fecha"
macro <- macro %>%
  mutate(fecha = as.Date(fecha)) %>% arrange(fecha) %>%
  mutate(IPC_yoy = (IPC / lag(IPC, 12) - 1) * 100)
ctrl_q <- macro %>%
  filter(month(fecha) %% 3 == 0) %>%
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

# Taxes residual (constant-price levels are additive)
qd <- qd %>% mutate(Impuestos = PIB - rowSums(across(all_of(SECT6))))

# Pre-defined groups (Section 2.3 classification)
GRP_FIN <- c("Ganaderia_fp", "Construccion", "Servicios")
GRP_INS <- c("Agricultura", "Manufactura", "Electricidad_y_agua", "Impuestos")
qd <- qd %>% mutate(FinDep    = rowSums(across(all_of(GRP_FIN))),
                    Insulated = rowSums(across(all_of(GRP_INS))))

# Average GDP shares (full sample, constant 2014 prices)
shares <- sapply(c(SECT6, "Impuestos", "FinDep", "Insulated"),
                 function(v) mean(qd[[v]] / qd$PIB, na.rm = TRUE))
cat("\nAverage GDP shares:\n"); print(round(100 * shares, 1))

# YoY growth
for (v in c(SECT6, "Impuestos", "FinDep", "Insulated", "PIB"))
  qd[[paste0("AG_", v)]] <- (qd[[v]] / dplyr::lag(qd[[v]], 4) - 1) * 100
qd <- qd %>% mutate(AG_GroupDiff = AG_FinDep - AG_Insulated, q_date = fecha)

az <- fci_q %>% inner_join(qd, by = "q_date") %>%
  left_join(ctrl_q, by = "q_date") %>% arrange(q_date)
cat(sprintf("\nMerged quarterly sample: %d obs, %s to %s\n",
            nrow(az), min(az$q_date), max(az$q_date)))

# ---------------------------------------------------------------------------
# LP engine (script 27 convention)
# ---------------------------------------------------------------------------
lp_q <- function(data, y_var, max_h = 8) {
  out <- data.frame()
  for (h in 1:max_h) {
    dh <- data %>%
      mutate(y_fwd = lead(.data[[y_var]], h),
             y_lag1 = lag(.data[[y_var]], 1), y_lag2 = lag(.data[[y_var]], 2),
             fci_lag1 = lag(FCI_COMP, 1),
             IPC_yoy_lag1 = lag(IPC_yoy, 1), IPC_yoy_lag2 = lag(IPC_yoy, 2))
    rhs <- c("FCI_COMP", "y_lag1", "y_lag2", "fci_lag1",
             "IPC_yoy", "IPC_yoy_lag1", "IPC_yoy_lag2")
    rd <- dh %>% dplyr::select(y_fwd, all_of(rhs)) %>% na.omit()
    if (nrow(rd) < 30) next
    m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
    ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1,
                                                         prewhite = FALSE))
    out <- rbind(out, data.frame(
      y = y_var, horizon = h, coef = ct["FCI_COMP", 1], se = ct["FCI_COMP", 2],
      p_value = ct["FCI_COMP", 4], n_obs = nrow(rd)))
  }
  out
}

# ---------------------------------------------------------------------------
# 1. Group aggregates + equality test
# ---------------------------------------------------------------------------
cat("\n[1] Group aggregates and equality test...\n")
grp <- rbind(lp_q(az, "AG_FinDep"), lp_q(az, "AG_Insulated"),
             lp_q(az, "AG_GroupDiff"), lp_q(az, "AG_PIB"))
# NOTE: shares[...] carries its own names, so naming the elements of c() again
# produced "FinDep.FinDep"-style names and the lookup below silently returned NA
# for both group rows. unname() is required for the match to resolve.
grp$share_pct <- round(100 * c(FinDep    = unname(shares["FinDep"]),
                               Insulated = unname(shares["Insulated"]),
                               GroupDiff = NA, PIB = 1)[
  sub("AG_", "", grp$y)], 1)
fwrite(grp, file.path(OUT, "Rev_Sector_GroupTests.csv"))
print(grp %>% filter(horizon %in% c(2, 4)) %>%
        mutate(across(where(is.numeric), ~round(., 3))), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Sector x horizon family with BH-FDR q-values
# ---------------------------------------------------------------------------
cat("\n[2] Sector family FDR...\n")
fam <- do.call(rbind, lapply(paste0("AG_", SECT6), function(v) lp_q(az, v)))
fam$q_bh <- p.adjust(fam$p_value, method = "BH")
fwrite(fam, file.path(OUT, "Rev_Sector_Family_FDR.csv"))
cat("  Rejections at 10%: nominal p:", sum(fam$p_value < 0.10),
    "of", nrow(fam), "| BH q:", sum(fam$q_bh < 0.10), "\n")
print(fam %>% filter(q_bh < 0.20 | (horizon %in% c(2, 4) & p_value < 0.10)) %>%
        mutate(across(where(is.numeric), ~round(., 3))), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Weighted-aggregation check
# ---------------------------------------------------------------------------
cat("\n[3] Weighted sector responses vs aggregate GDP...\n")
w7 <- shares[c(SECT6, "Impuestos")]
fam7 <- do.call(rbind, lapply(paste0("AG_", c(SECT6, "Impuestos")),
                              function(v) lp_q(az, v)))
agg_check <- do.call(rbind, lapply(c(2, 4), function(h) {
  fh <- fam7[fam7$horizon == h, ]
  fh$w <- w7[sub("AG_", "", fh$y)]
  data.frame(horizon = h,
             weighted_sector_sum = sum(fh$w * fh$coef),
             gdp_direct = grp$coef[grp$y == "AG_PIB" & grp$horizon == h],
             gdp_direct_p = grp$p_value[grp$y == "AG_PIB" & grp$horizon == h])
}))
fwrite(agg_check, file.path(OUT, "Rev_Sector_Aggregation_Check.csv"))
print(round(agg_check, 3), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Verdict
# ---------------------------------------------------------------------------
cat("\n[4] Verdict...\n")
gd <- grp %>% filter(y == "AG_GroupDiff", horizon %in% c(2, 4))
verdict <- data.frame(
  test = c("Group equality h=2Q (diff, p)", "Group equality h=4Q (diff, p)",
           "FinDep h=2Q (coef, p)", "FinDep h=4Q (coef, p)",
           "Insulated h=2Q (coef, p)", "Insulated h=4Q (coef, p)",
           "BH-FDR rejections at q<0.10 (of 48)"),
  value = c(sprintf("%.2f (p=%.3f)", gd$coef[gd$horizon == 2], gd$p_value[gd$horizon == 2]),
            sprintf("%.2f (p=%.3f)", gd$coef[gd$horizon == 4], gd$p_value[gd$horizon == 4]),
            sprintf("%.2f (p=%.3f)", grp$coef[grp$y == "AG_FinDep" & grp$horizon == 2],
                    grp$p_value[grp$y == "AG_FinDep" & grp$horizon == 2]),
            sprintf("%.2f (p=%.3f)", grp$coef[grp$y == "AG_FinDep" & grp$horizon == 4],
                    grp$p_value[grp$y == "AG_FinDep" & grp$horizon == 4]),
            sprintf("%.2f (p=%.3f)", grp$coef[grp$y == "AG_Insulated" & grp$horizon == 2],
                    grp$p_value[grp$y == "AG_Insulated" & grp$horizon == 2]),
            sprintf("%.2f (p=%.3f)", grp$coef[grp$y == "AG_Insulated" & grp$horizon == 4],
                    grp$p_value[grp$y == "AG_Insulated" & grp$horizon == 4]),
            sprintf("%d", sum(fam$q_bh < 0.10))))
fwrite(verdict, file.path(OUT, "Rev_Sector_Group_Verdict.csv"))
print(verdict, row.names = FALSE)

cat(sprintf("\n=== 57 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
