# ============================================================================
# 48_Gradient_TimeRobust_Inference.R
# Time-dependence-robust inference for the hedging gradient (Design C).
#
# Replicates the gradient specification of script 33 exactly (same panel,
# controls, fixed effects, sample filters) and reports, per horizon:
#   - bank-clustered analytic p (baseline, replication check vs script 33)
#   - two-way bank x month clustered p
#   - Driscoll-Kraay p (lag h+1)
#   - wild cluster bootstrap p clustered by MONTH (Rademacher, B=999),
#     pricing the common time-series shock directly.
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
library(fixest)
set.seed(20260710)   # fix bootstrap draws so reported p-values are exactly reproducible
cat("=== 48: gradient time-robust inference ===\n")

dA <- read_rds_micro("micro_designA_panel.rds")
p2 <- read_rds_micro("micro_p2_sector.rds")

CTRLS <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
FE_A  <- "bst + bsc"

dA[, hedge_class := fifelse(sector %in% SECT_HEDGED,   "hedged",
                    fifelse(sector %in% SECT_UNHEDGED, "unhedged",
                    fifelse(sector %in% SECT_AMBIGUOUS, "ambiguous", "excluded")))]

fx16 <- p2[ym <= to_ym("2016-12-01"), .(tot = sum(total)), by = .(sector, cur)]
fx16 <- dcast(fx16, sector ~ cur, value.var = "tot")
setnames(fx16, c("6200", "6900"), c("fc", "pyg"))
fx16[, fx_share_2016 := fc / (fc + pyg)]
dA <- merge(dA, fx16[, .(sector, fx_share_2016)], by = "sector", all.x = TRUE)
dA[, fxsh_z := zstd(fx_share_2016)]
dA[, shk_usd_fxsh := shk_usd * fxsh_z]

grad <- dA[hedge_class != "excluded"]

res <- list()
for (h in c(3, 6, 12, 18)) {
  yv  <- paste0("g", h)
  ds  <- grad[get(paste0("ok", h)) == TRUE]
  rhs <- c("shk_usd", "shk_usd_fxsh", CTRLS)
  m <- feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A)),
             data = ds, cluster = ~bank, notes = FALSE)
  ib <- inference_battery(m, "shk_usd_fxsh", dk_lag = h + 1)
  p_wcb_month <- if (h %in% c(12, 18)) {
    wcb_pval(ds, yv, rhs, FE_A, "shk_usd_fxsh", cluster = "ym", B = WCB_B)
  } else NA_real_
  res[[length(res) + 1]] <- data.table(
    h = h, b = ib$b,
    se_bank = ib$se_bank, p_bank = ib$p_bank,
    se_twoway = ib$se_twoway, p_twoway = ib$p_twoway,
    se_dk = ib$se_dk, p_dk = ib$p_dk,
    p_wcb_month = p_wcb_month, n = ib$n)
  cat(sprintf("h=%2d: b=%6.3f | p_bank=%.3f | p_twoway=%.3f | p_dk=%.3f | p_wcb_month=%s | n=%d\n",
              h, ib$b, ib$p_bank, ib$p_twoway, ib$p_dk,
              ifelse(is.na(p_wcb_month), "  -  ", sprintf("%.3f", p_wcb_month)), ib$n))
}
out <- rbindlist(res)
write_micro_csv(out, "Micro_Gradient_TimeRobust.csv")
cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---------------------------------------------------------------------------
# Part 2: SECTOR-robust and time-block inference.
# The identifying exposure (2016 sector FX share) varies at the sector level
# (13 sectors), so sector-level dependence is the first-order clustering
# concern for the gradient. Adds: sector-clustered SE, two-way sector x month,
# wild cluster bootstrap by sector, and a moving-block bootstrap over months
# (block length 18 >= horizon, pseudo-time relabeling so bank x sector x month
# FE do not merge across resampled blocks).
# ---------------------------------------------------------------------------
cat("\n--- Part 2: sector-robust and time-block inference ---\n")

block_boot_p <- function(ds, yv, rhs, param, L = 18, B = 499) {
  f <- as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "| bst + bsc"))
  m0 <- feols(f, data = ds, cluster = ~bank, notes = FALSE)
  b0 <- coef(m0)[param]; t0 <- fixest::coeftable(m0)[param, 3]
  ms <- sort(unique(ds$ym)); nm <- length(ms); k <- ceiling(nm / L)
  starts_max <- nm - L + 1
  tb <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    st <- sample.int(starts_max, k, replace = TRUE)
    sel <- unlist(lapply(st, function(s) ms[s:(s + L - 1)]))[1:nm]
    parts <- vector("list", nm)
    for (j in seq_len(nm)) {
      pj <- ds[ym == sel[j]]
      pj[, bst := paste(bank, sector, j)]
      parts[[j]] <- pj
    }
    db <- rbindlist(parts)
    mb <- tryCatch(feols(f, data = db, cluster = ~bank, notes = FALSE),
                   error = function(e) NULL)
    if (is.null(mb) || !(param %in% rownames(fixest::coeftable(mb)))) next
    ct <- fixest::coeftable(mb)
    tb[b] <- (ct[param, 1] - b0) / ct[param, 2]
  }
  mean(abs(tb) >= abs(t0), na.rm = TRUE)
}

res2 <- list()
for (h in c(12, 18)) {
  yv  <- paste0("g", h)
  ds  <- grad[get(paste0("ok", h)) == TRUE]
  rhs <- c("shk_usd", "shk_usd_fxsh", CTRLS)
  f   <- as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "| bst + bsc"))
  m   <- feols(f, data = ds, cluster = ~bank, notes = FALSE)
  b   <- coef(m)["shk_usd_fxsh"]
  se_sec  <- sqrt(diag(vcov(m, vcov = ~sector)))["shk_usd_fxsh"]
  se_2wsm <- sqrt(diag(vcov(m, vcov = ~sector + ym)))["shk_usd_fxsh"]
  p_sec   <- 2 * pnorm(-abs(b / se_sec))
  p_2wsm  <- 2 * pnorm(-abs(b / se_2wsm))
  p_wcb_sector <- wcb_pval(ds, yv, rhs, FE_A, "shk_usd_fxsh",
                           cluster = "sector", B = WCB_B)
  p_block <- block_boot_p(ds, yv, rhs, "shk_usd_fxsh", L = 18, B = 499)
  res2[[length(res2) + 1]] <- data.table(
    h = h, b = b, se_sector = se_sec, p_sector = p_sec,
    p_twoway_sector_month = p_2wsm, p_wcb_sector = p_wcb_sector,
    p_blockboot_18m = p_block, n = nobs(m))
  cat(sprintf("h=%2d: b=%6.3f | p_sector=%.3f | p_2way(sec,ym)=%.3f | p_wcb_sector=%.3f | p_block18=%.3f\n",
              h, b, p_sec, p_2wsm, p_wcb_sector, p_block))
}
out2 <- rbindlist(res2)
write_micro_csv(out2, "Micro_Gradient_SectorRobust.csv")
cat(sprintf("\nDone (all parts) in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---------------------------------------------------------------------------
# Part 3: wild cluster bootstrap for the hedged/unhedged SPLIT-SAMPLE rows of
# main-text Table 8 (previously cluster-robust only; header said "WCB p").
# Replicates script 33's split spec exactly.
# ---------------------------------------------------------------------------
cat("\n--- Part 3: split-sample WCB (Table 8 rows) ---\n")
SECT_TINY_FX <- c("CONSUMO", "VIVIENDA")
res3 <- list()
for (grp in c("hedged", "unhedged")) {
  dg <- dA[hedge_class == grp & (!sector %in% SECT_TINY_FX | total > 5000)]
  for (h in c(6, 12)) {
    yv <- paste0("g", h); ds <- dg[get(paste0("ok", h)) == TRUE]
    rhs <- c("shk_usd", CTRLS)
    m <- feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A)),
               data = ds, cluster = ~bank, notes = FALSE)
    cr <- coef_row(m, "shk_usd")
    pw <- wcb_pval(ds, yv, rhs, FE_A, "shk_usd", B = WCB_B)
    res3[[length(res3) + 1]] <- data.table(group = grp, h = h, b = cr$b,
      p_bank = cr$p, p_wcb = pw, n = cr$n)
    cat(sprintf("  %-9s h=%2d: b=%6.2f  p_bank=%.3f  p_wcb=%.3f\n", grp, h, cr$b, cr$p, pw))
  }
}
write_micro_csv(rbindlist(res3), "Micro_Split_WCB.csv")
cat("Part 3 done\n")

# ---------------------------------------------------------------------------
# Part 4: design-based few-cluster inference.
# With 13 sector clusters, Rademacher weights can be unreliable; adds:
#   (a) Webb 6-point wild cluster bootstrap by sector (MacKinnon-Webb 2018)
#   (b) CR2 (bias-reduced) cluster-robust SE with Satterthwaite df, computed
#       on the FWL-partialled regression (variables demeaned by bst + bsc)
#   (c) leave-one-sector-out stability: gradient re-estimated dropping each
#       sector in turn (13 runs per horizon).
# NO permutation of sector FX shares (exposure
# is not randomly assigned; exchangeability indefensible).
# ---------------------------------------------------------------------------
cat("\n--- Part 4: few-cluster (Webb / CR2 / leave-one-sector-out) ---\n")
set.seed(20260712)

res4 <- list(); loso <- list()
for (h in c(12, 18)) {
  yv  <- paste0("g", h)
  ds  <- grad[get(paste0("ok", h)) == TRUE]
  rhs <- c("shk_usd", "shk_usd_fxsh", CTRLS)

  # (a) Webb-weight WCB by sector
  p_webb <- wcb_pval(ds, yv, rhs, FE_A, "shk_usd_fxsh",
                     cluster = "sector", B = WCB_B, weights = "webb")

  # (b) CR2 + Satterthwaite on the FWL-partialled regression
  f   <- as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A))
  m   <- feols(f, data = ds, cluster = ~bank, notes = FALSE)
  du  <- ds[fixest::obs(m)]
  dm  <- fixest::demean(du[, c(yv, rhs), with = FALSE],
                        f = du[, .(bst, bsc)])
  dm  <- as.data.frame(dm); dm$sector <- du$sector
  mlm <- lm(as.formula(paste(yv, "~ 0 +", paste(rhs, collapse = "+"))), data = dm)
  ct2 <- clubSandwich::coef_test(mlm,
           vcov = clubSandwich::vcovCR(mlm, cluster = dm$sector, type = "CR2"),
           coefs = "shk_usd_fxsh")
  b_cr2 <- ct2$beta; se_cr2 <- ct2$SE; df_cr2 <- ct2$df_Satt; p_cr2 <- ct2$p_Satt

  res4[[length(res4) + 1]] <- data.table(
    h = h, b = coef(m)["shk_usd_fxsh"],
    p_wcb_webb_sector = p_webb,
    se_cr2 = se_cr2, df_satt = df_cr2, p_cr2_satt = p_cr2)
  cat(sprintf("h=%2d: b=%6.3f | p_webb(sector)=%.3f | CR2 se=%.3f df=%.1f p=%.3f\n",
              h, coef(m)["shk_usd_fxsh"], p_webb, se_cr2, df_cr2, p_cr2))

  # (c) leave-one-sector-out
  for (sec in sort(unique(as.character(ds$sector)))) {
    dl <- ds[sector != sec]
    ml <- tryCatch(feols(f, data = dl, cluster = ~bank, notes = FALSE),
                   error = function(e) NULL)
    if (is.null(ml) || !"shk_usd_fxsh" %in% rownames(fixest::coeftable(ml))) next
    cl <- fixest::coeftable(ml)["shk_usd_fxsh", ]
    loso[[length(loso) + 1]] <- data.table(
      h = h, dropped_sector = sec, b = cl[1], se_bank = cl[2], p_bank = cl[4],
      n = nobs(ml))
  }
}
out4 <- rbindlist(res4); write_micro_csv(out4, "Micro_Gradient_FewCluster.csv")
outL <- rbindlist(loso); write_micro_csv(outL, "Micro_Gradient_LOSO.csv")
for (hh in c(12, 18)) {
  ol <- outL[h == hh]
  cat(sprintf("LOSO h=%d: b range [%.3f, %.3f]; p_bank range [%.3f, %.3f]; sign flips: %d/%d\n",
              hh, min(ol$b), max(ol$b), min(ol$p_bank), max(ol$p_bank),
              sum(sign(ol$b) != sign(median(ol$b))), nrow(ol)))
}
cat("Part 4 done\n")

# ---------------------------------------------------------------------------
# Part 5: seeded ex-COVID gradient WCB (main-text Table 8 row).
# COVID window = COVID_START..COVID_END (micro_helpers.R: 2020-03..2021-06).
# ---------------------------------------------------------------------------
cat("\n--- Part 5: ex-COVID gradient (seeded WCB) ---\n")
set.seed(20260712)
ds5 <- grad[ok12 == TRUE & (ym < COVID_START | ym > COVID_END)]
rhs5 <- c("shk_usd", "shk_usd_fxsh", CTRLS)
m5 <- feols(as.formula(paste("g12 ~", paste(rhs5, collapse = "+"), "|", FE_A)),
            data = ds5, cluster = ~bank, notes = FALSE)
ct5 <- fixest::coeftable(m5)["shk_usd_fxsh", ]
pw5 <- wcb_pval(ds5, "g12", rhs5, FE_A, "shk_usd_fxsh", B = WCB_B)
out5 <- data.table(h = 12, b = ct5[1], se_bank = ct5[2], p_bank = ct5[4],
                   p_wcb = pw5, n = nobs(m5))
write_micro_csv(out5, "Micro_Gradient_ExCovid.csv")
cat(sprintf("ex-COVID gradient h=12: b=%.3f se=%.3f p_bank=%.3f p_wcb=%.3f n=%d\n",
            ct5[1], ct5[2], ct5[4], pw5, nobs(m5)))
