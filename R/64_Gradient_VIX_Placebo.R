# ============================================================================
# 64_Gradient_VIX_Placebo.R
# Matched placebo for the FX-credit-exposure gradient: replaces the dollar
# shock with VIX innovations in the SAME triple-interaction specification,
# sample, fixed effects, controls, horizons and inference battery as
# 48_Gradient_TimeRobust_Inference.R.
#
# Motivation. The manuscript's dollar-specificity claim is stated at the
# GRADIENT estimand (Shock x USD x 2016 sector FX-credit share) but the
# displayed VIX placebo (script 31, Micro_DesignA_Main.csv) is estimated on
# the AVERAGE USD-PYG differential. A null average VIX differential does not
# establish a null VIX gradient. This script estimates the matched object.
#
# Parts
#   1. VIX gradient, h = 3, 6, 12, 18: bank-cluster, two-way bank x month,
#      Driscoll-Kraay(h+1), month-clustered wild cluster bootstrap.
#   2. Sector-level and time-block inference (h = 12, 18): sector cluster,
#      two-way sector x month, sector WCB, moving-block bootstrap (L = 18).
#   3. Few-cluster checks (h = 12, 18): Webb-weight sector WCB, CR2 +
#      Satterthwaite on the FWL-partialled regression, leave-one-sector-out.
#   4. Joint DXY-vs-VIX difference test: both gradients in one regression,
#      H0: beta_DXY_gradient = beta_VIX_gradient, under every vcov in the
#      battery plus a wild cluster bootstrap on a reparameterised carrier.
#   5. Export of the 2016 baseline sector FX-credit shares (the exposure used
#      by the continuous gradient), with ranking and standardisation.
#
# Two facts to disclose wherever the placebo is reported:
#   (a) shk_dxy_purged is orthogonalised against contemporaneous d_vix by
#       construction (30_Bank_Micro_Data_Construction.R:114-115), so this is a
#       clean placebo but NOT a horse race between correlated shocks.
#   (b) shk_vix is the raw AR(3) innovation of d_VIX (30:118) - the same
#       convention already used for the average-differential placebo, so the
#       two VIX exercises are internally consistent.
#
# Seeds mirror script 48 (20260710 / 20260712) so the placebo and the dollar
# gradient are compared under identical bootstrap draws.
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
library(fixest)
set.seed(20260710)
cat("=== 64: matched VIX gradient placebo ===\n")

dA <- read_rds_micro("micro_designA_panel.rds")
p2 <- read_rds_micro("micro_p2_sector.rds")

CTRLS <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
FE_A  <- "bst + bsc"

dA[, hedge_class := fifelse(sector %in% SECT_HEDGED,   "hedged",
                    fifelse(sector %in% SECT_UNHEDGED, "unhedged",
                    fifelse(sector %in% SECT_AMBIGUOUS, "ambiguous", "excluded")))]

# ---- Exposure: identical construction to 48:30-35 ---------------------------
fx16 <- p2[ym <= to_ym("2016-12-01"), .(tot = sum(total)), by = .(sector, cur)]
fx16 <- dcast(fx16, sector ~ cur, value.var = "tot")
setnames(fx16, c("6200", "6900"), c("fc", "pyg"))
fx16[, fx_share_2016 := fc / (fc + pyg)]
dA <- merge(dA, fx16[, .(sector, fx_share_2016)], by = "sector", all.x = TRUE)
dA[, fxsh_z := zstd(fx_share_2016)]

dA[, shk_usd_fxsh := shk_usd * fxsh_z]   # dollar gradient (script 48 object)
dA[, vix_usd_fxsh := vix_usd * fxsh_z]   # matched VIX gradient (this script)

grad <- dA[hedge_class != "excluded"]

# ============================================================================
# Part 1: VIX gradient, analytic + month-clustered WCB
# ============================================================================
cat("\n--- Part 1: VIX gradient, time-robust inference ---\n")
res <- list()
for (h in c(3, 6, 12, 18)) {
  yv  <- paste0("g", h)
  ds  <- grad[get(paste0("ok", h)) == TRUE]
  rhs <- c("vix_usd", "vix_usd_fxsh", CTRLS)
  m <- feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A)),
             data = ds, cluster = ~bank, notes = FALSE)
  ib <- inference_battery(m, "vix_usd_fxsh", dk_lag = h + 1)
  p_wcb_month <- if (h %in% c(12, 18)) {
    wcb_pval(ds, yv, rhs, FE_A, "vix_usd_fxsh", cluster = "ym", B = WCB_B)
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
write_micro_csv(rbindlist(res), "Micro_VIXGradient_TimeRobust.csv")

# ============================================================================
# Part 2: sector-robust and time-block inference
# block_boot_p is copied verbatim from 48:76-100. Script 48 is a seeded
# canonical script quoted in the appendix; it is deliberately NOT refactored.
# ============================================================================
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
  rhs <- c("vix_usd", "vix_usd_fxsh", CTRLS)
  f   <- as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "| bst + bsc"))
  m   <- feols(f, data = ds, cluster = ~bank, notes = FALSE)
  b   <- coef(m)["vix_usd_fxsh"]
  se_sec  <- sqrt(diag(vcov(m, vcov = ~sector)))["vix_usd_fxsh"]
  se_2wsm <- sqrt(diag(vcov(m, vcov = ~sector + ym)))["vix_usd_fxsh"]
  p_sec   <- 2 * pnorm(-abs(b / se_sec))
  p_2wsm  <- 2 * pnorm(-abs(b / se_2wsm))
  p_wcb_sector <- wcb_pval(ds, yv, rhs, FE_A, "vix_usd_fxsh",
                           cluster = "sector", B = WCB_B)
  p_block <- block_boot_p(ds, yv, rhs, "vix_usd_fxsh", L = 18, B = 499)
  res2[[length(res2) + 1]] <- data.table(
    h = h, b = b, se_sector = se_sec, p_sector = p_sec,
    p_twoway_sector_month = p_2wsm, p_wcb_sector = p_wcb_sector,
    p_blockboot_18m = p_block, n = nobs(m))
  cat(sprintf("h=%2d: b=%6.3f | p_sector=%.3f | p_2way(sec,ym)=%.3f | p_wcb_sector=%.3f | p_block18=%.3f\n",
              h, b, p_sec, p_2wsm, p_wcb_sector, p_block))
}
write_micro_csv(rbindlist(res2), "Micro_VIXGradient_SectorRobust.csv")

# ============================================================================
# Part 3: few-cluster design-based checks (Webb WCB / CR2 / leave-one-sector-out)
# ============================================================================
cat("\n--- Part 3: few-cluster (Webb / CR2 / leave-one-sector-out) ---\n")
set.seed(20260712)

res3 <- list(); loso <- list()
for (h in c(12, 18)) {
  yv  <- paste0("g", h)
  ds  <- grad[get(paste0("ok", h)) == TRUE]
  rhs <- c("vix_usd", "vix_usd_fxsh", CTRLS)

  p_webb <- wcb_pval(ds, yv, rhs, FE_A, "vix_usd_fxsh",
                     cluster = "sector", B = WCB_B, weights = "webb")

  f   <- as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A))
  m   <- feols(f, data = ds, cluster = ~bank, notes = FALSE)
  du  <- ds[fixest::obs(m)]
  dm  <- fixest::demean(du[, c(yv, rhs), with = FALSE], f = du[, .(bst, bsc)])
  dm  <- as.data.frame(dm); dm$sector <- du$sector
  mlm <- lm(as.formula(paste(yv, "~ 0 +", paste(rhs, collapse = "+"))), data = dm)
  ct2 <- clubSandwich::coef_test(mlm,
           vcov = clubSandwich::vcovCR(mlm, cluster = dm$sector, type = "CR2"),
           coefs = "vix_usd_fxsh")

  res3[[length(res3) + 1]] <- data.table(
    h = h, b = coef(m)["vix_usd_fxsh"],
    p_wcb_webb_sector = p_webb,
    se_cr2 = ct2$SE, df_satt = ct2$df_Satt, p_cr2_satt = ct2$p_Satt)
  cat(sprintf("h=%2d: b=%6.3f | p_webb(sector)=%.3f | CR2 se=%.3f df=%.1f p=%.3f\n",
              h, coef(m)["vix_usd_fxsh"], p_webb, ct2$SE, ct2$df_Satt, ct2$p_Satt))

  for (sec in sort(unique(as.character(ds$sector)))) {
    dl <- ds[sector != sec]
    ml <- tryCatch(feols(f, data = dl, cluster = ~bank, notes = FALSE),
                   error = function(e) NULL)
    if (is.null(ml) || !"vix_usd_fxsh" %in% rownames(fixest::coeftable(ml))) next
    cl <- fixest::coeftable(ml)["vix_usd_fxsh", ]
    loso[[length(loso) + 1]] <- data.table(
      h = h, dropped_sector = sec, b = cl[1], se_bank = cl[2], p_bank = cl[4],
      n = nobs(ml))
  }
}
write_micro_csv(rbindlist(res3), "Micro_VIXGradient_FewCluster.csv")
outL <- rbindlist(loso); write_micro_csv(outL, "Micro_VIXGradient_LOSO.csv")
for (hh in c(12, 18)) {
  ol <- outL[h == hh]
  cat(sprintf("LOSO h=%d: b range [%.3f, %.3f]; p_bank range [%.3f, %.3f]; sign flips: %d/%d\n",
              hh, min(ol$b), max(ol$b), min(ol$p_bank), max(ol$p_bank),
              sum(sign(ol$b) != sign(median(ol$b))), nrow(ol)))
}

# ============================================================================
# Part 4: joint DXY-vs-VIX gradient difference test
#
# Both gradients enter one regression:
#   g_h ~ shk_usd + vix_usd + shk_usd_fxsh + vix_usd_fxsh + CTRLS | bst + bsc
# and H0: beta(shk_usd_fxsh) = beta(vix_usd_fxsh) is tested as
#   d = b_D - b_V,  Var(d) = V_DD + V_VV - 2 V_DV
# under each vcov in the battery.
#
# For a bootstrap version, the model is reparameterised: with
# grad_sum = shk_usd_fxsh + vix_usd_fxsh, the regressors (grad_sum,
# vix_usd_fxsh) span the same space and the coefficient on vix_usd_fxsh
# equals beta_V - beta_D. A restricted wild cluster bootstrap on that carrier
# is therefore a bootstrap test of equality.
# ============================================================================
cat("\n--- Part 4: DXY vs VIX gradient difference ---\n")
set.seed(20260712)

dA[, grad_sum := shk_usd_fxsh + vix_usd_fxsh]
grad <- dA[hedge_class != "excluded"]

wald_diff <- function(m, V) {
  b <- coef(m)
  d <- unname(b["shk_usd_fxsh"] - b["vix_usd_fxsh"])
  v <- V["shk_usd_fxsh", "shk_usd_fxsh"] + V["vix_usd_fxsh", "vix_usd_fxsh"] -
       2 * V["shk_usd_fxsh", "vix_usd_fxsh"]
  se <- sqrt(v)
  c(d = d, se = se, p = 2 * pnorm(-abs(d / se)))
}

res4 <- list()
for (h in c(12, 18)) {
  yv  <- paste0("g", h)
  ds  <- grad[get(paste0("ok", h)) == TRUE]
  rhs <- c("shk_usd", "vix_usd", "shk_usd_fxsh", "vix_usd_fxsh", CTRLS)
  m   <- feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A)),
               data = ds, cluster = ~bank, notes = FALSE)

  w_bank <- wald_diff(m, vcov(m, vcov = ~bank))
  w_2way <- wald_diff(m, vcov(m, vcov = ~bank + ym))
  w_sec  <- wald_diff(m, vcov(m, vcov = ~sector))
  w_dk   <- wald_diff(m, vcov(m, vcov = as.formula(sprintf("DK(%d) ~ ym", h + 1))))

  rhs_rp <- c("shk_usd", "vix_usd", "grad_sum", "vix_usd_fxsh", CTRLS)
  p_wcb_month  <- wcb_pval(ds, yv, rhs_rp, FE_A, "vix_usd_fxsh",
                           cluster = "ym", B = WCB_B)
  p_wcb_sector <- wcb_pval(ds, yv, rhs_rp, FE_A, "vix_usd_fxsh",
                           cluster = "sector", B = WCB_B, weights = "webb")

  res4[[length(res4) + 1]] <- data.table(
    h = h,
    b_dxy_gradient = unname(coef(m)["shk_usd_fxsh"]),
    b_vix_gradient = unname(coef(m)["vix_usd_fxsh"]),
    diff = unname(w_bank["d"]),
    se_bank = unname(w_bank["se"]), p_bank = unname(w_bank["p"]),
    se_twoway = unname(w_2way["se"]), p_twoway = unname(w_2way["p"]),
    se_sector = unname(w_sec["se"]), p_sector = unname(w_sec["p"]),
    se_dk = unname(w_dk["se"]), p_dk = unname(w_dk["p"]),
    p_wcb_month = p_wcb_month, p_wcb_webb_sector = p_wcb_sector,
    n = nobs(m))
  cat(sprintf("h=%2d: b_DXY=%6.3f b_VIX=%6.3f diff=%6.3f | p_bank=%.3f p_2way=%.3f p_sec=%.3f p_dk=%.3f | p_wcb_month=%.3f p_wcb_webb=%.3f\n",
              h, coef(m)["shk_usd_fxsh"], coef(m)["vix_usd_fxsh"], w_bank["d"],
              w_bank["p"], w_2way["p"], w_sec["p"], w_dk["p"],
              p_wcb_month, p_wcb_sector))
}
write_micro_csv(rbindlist(res4), "Micro_Gradient_DXY_vs_VIX_Diff.csv")

# ============================================================================
# Part 5: 2016 baseline sector FX-credit shares (the exposure variable)
# Printed so the exposure ranking used by the continuous gradient is auditable
# (the appendix previously printed May-2026 shares only).
# ============================================================================
cat("\n--- Part 5: 2016 baseline sector FX-credit shares ---\n")
sh16 <- copy(fx16)[, .(sector, fc_2016 = fc, pyg_2016 = pyg, fx_share_2016)]
# Two standardisations are reported because they are not the same object:
#   fxsh_z_est    - the one the gradient actually uses. Script 48 (and this
#                   script) standardise AFTER merging onto the panel, so the
#                   mean and SD are observation-weighted across bank x sector x
#                   month x currency cells, not equal-weighted across sectors.
#   fxsh_z_sector - equal-weighted across the 13 sectors, reported so the
#                   exposure ranking can be read without the panel.
mu_obs <- mean(dA$fx_share_2016, na.rm = TRUE)
sd_obs <- sd(dA$fx_share_2016, na.rm = TRUE)
sh16[, fxsh_z_est    := (fx_share_2016 - mu_obs) / sd_obs]
sh16[, fxsh_z_sector := zstd(fx_share_2016)]
setorder(sh16, -fx_share_2016)
sh16[, rank_2016 := seq_len(.N)]
sh16[, in_gradient_sample := sector %in% unique(as.character(grad$sector))]
sh16[, rank_in_sample := NA_integer_]
sh16[in_gradient_sample == TRUE, rank_in_sample := seq_len(.N)]
write_micro_csv(sh16, "Micro_Sector_FX_Shares_2016.csv")
cat(sprintf("Observation-weighted exposure moments: mean=%.4f sd=%.4f\n", mu_obs, sd_obs))
print(sh16[, .(rank_2016, sector, fx_share_2016 = round(fx_share_2016, 3),
               fxsh_z_est = round(fxsh_z_est, 2), in_gradient_sample, rank_in_sample)])

cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
