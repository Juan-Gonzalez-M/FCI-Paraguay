# ============================================================================
# 65_Sector_Group_Stability.R — Round-15 colleague item: is the sector-group
# contrast driven by one large component or by a convenient classification?
#
# CONCERN.  The finance-dependent block (livestock+, construction, services)
# is ~56% of GDP, of which services alone is ~47% of GDP; manufacturing sits
# in the insulated block despite being partly agro-processing; and the taxes
# residual has no behavioural interpretation.  The blocks were pre-classified
# but not pre-registered, and were constructed during the revision.  The
# equality rejection at h=4Q is also convention-dependent (Driscoll-Kraay
# p = 0.008 vs quarter-clustered p = 0.090).
#
# THIS SCRIPT re-runs the stacked coefficient-equality test of script 63 under:
#   (a) leave-one-component-out: each of the 7 components dropped in turn
#       (including the "excluding services" row the referees would ask for);
#   (b) alternative classifications: manufacturing moved to the finance-
#       dependent block; the taxes residual excluded entirely;
#   (c) the post-IT subsample (May 2011+), which the manuscript's conclusion
#       currently asserts without a displayed group test.
#
# The estimator, controls, FE-free stacked design, Driscoll-Kraay lag h+1 and
# quarter-clustered cross-check are IDENTICAL to script 63; only the group
# membership and the sample vary.
#
# CONSTRUCTION RULE.  Impuestos is always PIB minus the full six-sector set,
# so dropping a component from a GROUP never redefines the residual.  A
# dropped component leaves the estimation entirely (it is in neither block).
#
# ADDITIVE: writes only output/revision/csv/Rev_Sector_Group_Stability.csv.
# ============================================================================

t0 <- Sys.time()
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(lubridate)
  library(sandwich); library(lmtest)
})
cat("=== 65: Sectoral group stability (LOO / alt-classification / post-IT) ===\n")

OUT    <- "../output/revision/csv"
POSTIT <- as.Date("2011-05-01")

# --- Data (script 63 conventions, verbatim) ---------------------------------
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

# Residual is ALWAYS defined against the full six-sector set (see header).
qd <- qd %>% mutate(Impuestos = PIB - rowSums(across(all_of(SECT6))))
qd <- qd %>% mutate(q_date = fecha)

base_panel <- fci_q %>% inner_join(qd, by = "q_date") %>%
  left_join(ctrl_q, by = "q_date") %>% arrange(q_date)

# Average GDP shares (full sample), reported so each variant's block weight is
# auditable rather than asserted.
COMPONENTS <- c(SECT6, "Impuestos")
shares <- sapply(COMPONENTS, function(v)
  mean(base_panel[[v]] / base_panel$PIB, na.rm = TRUE))
cat("\nAverage GDP shares (full sample):\n")
print(round(shares, 4))

# --- Stacked, fully-interacted system (script 63 estimator) -----------------
RHS <- c("FCI_COMP", "y_lag1", "y_lag2", "fci_lag1",
         "IPC_yoy", "IPC_yoy_lag1", "IPC_yoy_lag2")

group_wald <- function(data, grp_fin, grp_ins, label, horizons = 1:8) {
  d0 <- data %>%
    mutate(FinDep    = rowSums(across(all_of(grp_fin))),
           Insulated = rowSums(across(all_of(grp_ins))))
  for (v in c("FinDep", "Insulated"))
    d0[[paste0("AG_", v)]] <- (d0[[v]] / dplyr::lag(d0[[v]], 4) - 1) * 100

  build <- function(dd, h) {
    d <- dd %>% mutate(
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

  out <- data.frame()
  for (h in horizons) {
    s <- build(d0, h)
    if (nrow(s) < 60) next
    f <- as.formula(paste("y_fwd ~ DI +", paste(RHS, collapse = "+"),
                          "+", paste0("DI:", RHS, collapse = "+")))
    m <- lm(f, data = s)
    # order.by must be explicit or vcovPL degenerates (round-14 gotcha).
    V <- sandwich::vcovPL(m, cluster = s$q_date, order.by = s$q_date,
                          lag = h + 1, adjust = TRUE)
    bF   <- unname(coef(m)["FCI_COMP"])
    dIF  <- unname(coef(m)["DI:FCI_COMP"])
    bI   <- bF + dIF
    seD  <- sqrt(V["DI:FCI_COMP", "DI:FCI_COMP"])
    diff <- -dIF
    dfree <- length(unique(s$q_date)) - 2
    tc    <- qt(0.975, dfree)
    seCL <- sqrt(sandwich::vcovCL(m, cluster = s$q_date)["DI:FCI_COMP",
                                                         "DI:FCI_COMP"])
    out <- rbind(out, data.frame(
      variant = label, horizon = h,
      beta_FinDep = bF, beta_Insulated = bI, diff_F_minus_I = diff,
      se = seD, t = diff / seD, p_value = 2 * pt(-abs(diff / seD), dfree),
      ci_lo = diff - tc * seD, ci_hi = diff + tc * seD,
      se_clusterdate = seCL,
      p_clusterdate = 2 * pt(-abs(diff / seCL), dfree),
      share_fin = sum(shares[grp_fin]), share_ins = sum(shares[grp_ins]),
      n_stacked = nrow(s), n_quarters = length(unique(s$q_date)), df = dfree))
  }
  out
}

GRP_FIN <- c("Ganaderia_fp", "Construccion", "Servicios")
GRP_INS <- c("Agricultura", "Manufactura", "Electricidad_y_agua", "Impuestos")

variants <- list()
variants[["baseline"]] <- list(fin = GRP_FIN, ins = GRP_INS, data = base_panel)

# (a) leave-one-component-out
for (v in COMPONENTS) {
  variants[[paste0("drop_", v)]] <- list(
    fin  = setdiff(GRP_FIN, v),
    ins  = setdiff(GRP_INS, v),
    data = base_panel)
}

# (b) alternative classifications
variants[["alt_manufactura_to_findep"]] <- list(
  fin  = c(GRP_FIN, "Manufactura"),
  ins  = setdiff(GRP_INS, "Manufactura"),
  data = base_panel)
variants[["alt_impuestos_excluded"]] <- list(
  fin  = GRP_FIN,
  ins  = setdiff(GRP_INS, "Impuestos"),
  data = base_panel)

# (c) post-IT subsample, baseline classification
variants[["postIT_baseline"]] <- list(
  fin  = GRP_FIN, ins = GRP_INS,
  data = base_panel %>% filter(q_date >= POSTIT))

res <- do.call(rbind, lapply(names(variants), function(nm) {
  v <- variants[[nm]]
  cat(sprintf("  [%-28s] fin=%d comp (%.0f%% GDP)  ins=%d comp (%.0f%% GDP)\n",
              nm, length(v$fin), 100 * sum(shares[v$fin]),
              length(v$ins), 100 * sum(shares[v$ins])))
  group_wald(v$data, v$fin, v$ins, nm)
}))

write.csv(res, file.path(OUT, "Rev_Sector_Group_Stability.csv"), row.names = FALSE)
cat(sprintf("  [csv]  Rev_Sector_Group_Stability.csv (%d rows)\n", nrow(res)))

cat("\n--- h = 4Q across variants (the headline horizon) ---\n")
r4 <- res[res$horizon == 4, ]
print(format(r4[, c("variant", "beta_FinDep", "beta_Insulated", "diff_F_minus_I",
                    "ci_lo", "ci_hi", "p_value", "p_clusterdate",
                    "share_fin", "n_quarters")], digits = 3), row.names = FALSE)

cat("\n--- h = 2Q across variants ---\n")
r2 <- res[res$horizon == 2, ]
print(format(r2[, c("variant", "diff_F_minus_I", "p_value", "p_clusterdate",
                    "n_quarters")], digits = 3), row.names = FALSE)

cat("\n--- Post-IT full horizon profile (the currently undisplayed claim) ---\n")
rp <- res[res$variant == "postIT_baseline", ]
if (nrow(rp)) {
  print(format(rp[, c("horizon", "beta_FinDep", "beta_Insulated",
                      "diff_F_minus_I", "ci_lo", "ci_hi", "p_value",
                      "p_clusterdate", "n_quarters")], digits = 3),
        row.names = FALSE)
} else {
  cat("  post-IT sample did not clear the n_stacked >= 60 guard at any horizon\n")
}

cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
