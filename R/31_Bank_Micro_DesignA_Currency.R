# ============================================================================
# 31_Bank_Micro_DesignA_Currency.R
# Design A (plan §2): within-bank, within-sector currency-denomination test.
#
#   Δʰ ln L_{b,s,c,t+h} = βₕ (Shockₜ × USD_c) + α_{b×s×t} + γ_{b×s×c}
#                          + λ (X_{b,t-12} × USD_c) + ε
#
# bank×sector×time FE absorb ALL aggregate shocks and bank-sector demand;
# identification is the within-cell currency contrast. Primary shock: purged
# DXY innovation. Inference: restricted wild cluster bootstrap by bank
# (baseline), two-way and Driscoll-Kraay as secondary.
# Falsifications (plan §2.4): VIX placebo, random-timing permutation placebo.
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
library(fixest)
cat("=== 31: Design A - within-bank currency-denomination test ===\n")

macro <- read_rds_micro("micro_macro.rds")
p2    <- read_rds_micro("micro_p2_sector.rds")
p4    <- read_rds_micro("micro_p4_bankchars.rds")

# ---------------------------------------------------------------------------
# 1. Build the Design A estimation panel
# ---------------------------------------------------------------------------
build_designA_panel <- function(min_cell = CELL_MIN_BASE, winsor = WINSOR_P,
                                lval_var = "lval") {
  d <- copy(p2)

  # Lead growths on the valuation-adjusted log level, and 12m lagged growth
  add_lead_growth(d, lval_var, by = c("bank", "sector", "cur"))
  add_lag_growth(d, lval_var, by = c("bank", "sector", "cur"))

  # Cell-size filter at t; the lead side is enforced by requiring the cell to
  # pass the filter at t+h as well (plan §1.4: stock > threshold at both ends)
  setorderv(d, c("bank", "sector", "cur", "ym"))
  for (h in HORIZONS) {
    d[, paste0("ok", h) := {
      idx <- match(ym + h, ym)
      total > min_cell & !is.na(total[idx]) & total[idx] > min_cell
    }, by = .(bank, sector, cur)]
  }

  # Winsorize growth rates within currency (plan §1.4)
  for (h in HORIZONS) {
    d[, paste0("g", h) := winsorize(get(paste0("g", h)), winsor), by = cur]
  }

  # Shock, treatment, controls
  d <- merge(d, macro[, .(ym, shk_dxy_purged, shk_dxy_raw, shk_vix, shk_fci)],
             by = "ym", all.x = TRUE)
  d[, usd := as.integer(cur == CUR_FC)]
  d <- merge(d, p4[, .(bank, ym, tier1_l12_z, size_l12_z, fc_dep_share_l12_z)],
             by = c("bank", "ym"), all.x = TRUE)

  # Interactions (levels absorbed by the FE structure)
  d[, `:=`(shk_usd   = shk_dxy_purged * usd,
           shkraw_usd = shk_dxy_raw   * usd,
           vix_usd   = shk_vix        * usd,
           fci_usd   = shk_fci        * usd,
           glag_usd  = g_lag12        * usd,
           tier1_usd = tier1_l12_z    * usd,
           size_usd  = size_l12_z     * usd,
           fcdep_usd = fc_dep_share_l12_z * usd)]

  # Sample: shock months 2016m1-2025m12, event windows dropped,
  # dual-currency cells only (the contrast needs both books)
  d <- d[ym >= to_ym("2016-01-01") & ym <= to_ym("2025-12-01") &
          drop_event == FALSE & !is.na(shk_dxy_purged)]
  d[, n_cur := uniqueN(cur[total > min_cell]), by = .(bank, sector, ym)]
  d <- d[n_cur == 2]
  d[, `:=`(bst = paste(bank, sector, ym), bsc = paste(bank, sector, cur))]
  d
}

dA <- build_designA_panel()
cat(sprintf("  Design A panel: %s cell-month rows, %d banks, %d sectors\n",
            format(nrow(dA), big.mark = ","), uniqueN(dA$bank), uniqueN(dA$sector)))

CTRLS <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
FE_A  <- "bst + bsc"

# ---------------------------------------------------------------------------
# 2. Main LP loop, h = 1..18: purged DXY (primary) + VIX placebo
# ---------------------------------------------------------------------------
run_designA <- function(d, shockvar, label, wcb_h = HORIZONS, B = WCB_B) {
  res <- list()
  for (h in HORIZONS) {
    yv  <- paste0("g", h)
    ds  <- d[get(paste0("ok", h)) == TRUE]
    rhs <- c(shockvar, CTRLS)
    m <- tryCatch(
      feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", FE_A)),
            data = ds, cluster = ~bank, notes = FALSE),
      error = function(e) NULL)
    if (is.null(m)) next
    ib <- inference_battery(m, shockvar, dk_lag = h + 1)
    pw <- if (h %in% wcb_h) wcb_pval(ds, yv, rhs, FE_A, shockvar, B = B) else NA_real_
    res[[h]] <- data.table(shock = label, h = h, b = ib$b,
                           se_bank = ib$se_bank, p_bank = ib$p_bank,
                           p_wcb = pw, se_twoway = ib$se_twoway,
                           p_twoway = ib$p_twoway, se_dk = ib$se_dk,
                           p_dk = ib$p_dk, n = ib$n,
                           n_banks = uniqueN(ds[fixest::obs(m), bank]))
    cat(sprintf("  [%s] h=%2d: b=%7.3f  se=%.3f  p_bank=%.3f  p_wcb=%s  n=%d\n",
                label, h, ib$b, ib$se_bank, ib$p_bank,
                ifelse(is.na(pw), "-", sprintf("%.3f", pw)), ib$n))
  }
  rbindlist(res)
}

cat("\nPrimary: purged DXY innovation x USD\n")
resA <- run_designA(dA, "shk_usd", "purged_DXY")

cat("\nFalsification: VIX placebo x USD (plan §2.4)\n")
resV <- run_designA(dA, "vix_usd", "VIX_placebo", wcb_h = H_KEY)

write_micro_csv(rbind(resA, resV), "Micro_DesignA_Main.csv")

# ---------------------------------------------------------------------------
# 3. Random-timing permutation placebo at h = 12 (500 draws, plan §2.4)
# ---------------------------------------------------------------------------
cat("\nPermutation placebo (500 draws, h = 12)...\n")
set.seed(20260703)
h12 <- dA[ok12 == TRUE & !is.na(g12)]
shk_by_ym <- unique(h12[, .(ym, shk_dxy_purged)])
b_obs <- resA[h == 12, b]

perm_b <- replicate(500, {
  perm <- copy(shk_by_ym)
  perm[, shk_perm := sample(shk_dxy_purged)]
  dp <- merge(h12, perm[, .(ym, shk_perm)], by = "ym")
  dp[, perm_usd := shk_perm * usd]
  m <- feols(as.formula(paste("g12 ~ perm_usd +", paste(CTRLS, collapse = "+"), "|", FE_A)),
             data = dp, notes = FALSE)
  coef(m)["perm_usd"]
})
p_perm <- mean(abs(perm_b) >= abs(b_obs))
cat(sprintf("  beta12 observed = %.3f; permutation p = %.3f\n", b_obs, p_perm))
write_micro_csv(data.table(draw = seq_along(perm_b), beta12_perm = perm_b,
                           beta12_obs = b_obs, p_perm = p_perm),
                "Micro_DesignA_Placebo_Permutation.csv")

# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------
plot_dynamic_beta(resA,
  "Design A: differential USD-vs-PYG credit response to a 1 SD dollar shock",
  "Within bank x sector x month; purged DXY innovation; 90% bands (cluster by bank); valuation-adjusted growth",
  "260_DesignA_Dynamic_Beta.png")

cmp <- rbind(resA[, .(h, b, se = se_bank, Shock = "Purged DXY")],
             resV[, .(h, b, se = se_bank, Shock = "VIX placebo")])
pV <- ggplot(cmp, aes(h, b, color = Shock, fill = Shock)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = b - 1.645 * se, ymax = b + 1.645 * se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  scale_color_manual(values = c("steelblue4", "grey40")) +
  scale_fill_manual(values = c("steelblue4", "grey40")) +
  labs(title = "Design A falsification: dollar shock vs VIX placebo",
       subtitle = "A pure global-risk shock should not open a within-cell currency differential",
       x = "Horizon (months)", y = expression(beta[h] ~ "(pp)")) +
  theme_micro()
save_png(pV, "261_DesignA_VIX_Placebo.png")

pP <- ggplot(data.table(b = perm_b), aes(b)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white") +
  geom_vline(xintercept = b_obs, color = "firebrick", linewidth = 1) +
  labs(title = "Design A: random-timing placebo distribution (500 permutations, h = 12)",
       subtitle = sprintf("Observed beta12 = %.3f (red); permutation p = %.3f", b_obs, p_perm),
       x = expression(beta[12] ~ "under permuted shock dates"), y = "Count") +
  theme_micro()
save_png(pP, "262_DesignA_Permutation_Placebo.png")

write_rds_micro(dA, "micro_designA_panel.rds")
cat(sprintf("=== 31 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
