# ============================================================================
# 32_Bank_Micro_DesignB_Exposure.R
# Design B (plan §3): bank lending channel via exposure heterogeneity.
#
# (i)  Bank level:        Δʰ ln L_{b,t+h} = βₕ (Shockₜ × Exposure_{b,t-12})
#                                            + α_b + τ_t + φ X_{b,t-12} + ε
# (ii) Bank×sector level: Δʰ ln L_{b,s,t+h} = βₕ (Shockₜ × Exposure_b)
#                                            + α_{s×t} + γ_{b×s} + ε   (preferred)
#      sector×time FE absorb sector credit demand (Degryse et al. 2019 JFI
#      industry-time FE as the Khwaja-Mian substitute).
#
# Time FE absorb the aggregate shock => no instrument needed; the post-IT
# weak-IV concern (R1 #2) is irrelevant by construction. Outcome: constant-
# exchange-rate, CPI-deflated total credit (FX-adjusted aggregation).
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
library(fixest)
cat("=== 32: Design B - bank lending channel via exposure heterogeneity ===\n")

macro <- read_rds_micro("micro_macro.rds")
p1    <- read_rds_micro("micro_p1_carteras.rds")
p2    <- read_rds_micro("micro_p2_sector.rds")
p4    <- read_rds_micro("micro_p4_bankchars.rds")

EXPOSURES <- c(fc_dep_share = "FC deposit share",
               ext_fund_share = "External funding share",
               tier1 = "Tier 1 ratio",
               size = "Size (ln RWA)",
               liquidity = "Liquidity ratio")

e_ref  <- macro[ym == to_ym("2016-01-01"), TCN]

# ---------------------------------------------------------------------------
# 1. Outcomes: constant-FX real credit (bank level and bank x sector level)
# ---------------------------------------------------------------------------
mk_constfx <- function(dt, by_cols) {
  w <- dcast(dt, as.formula(paste(paste(c(by_cols, "ym"), collapse = "+"), "~ cur")),
             value.var = "total", fun.aggregate = sum)
  setnames(w, c("6200", "6900"), c("fc", "pyg"))
  w <- merge(w, macro[, .(ym, TCN, IPC)], by = "ym", all.x = TRUE)
  w[, constfx := pyg + fc * (e_ref / TCN)]           # FC book at constant e
  w[, lval := log(constfx / IPC * 100)]              # CPI-deflated
  w[constfx <= 0 | is.na(TCN), lval := NA_real_]
  w
}

bk <- mk_constfx(p1[, .(bank, cur, ym, total, drop_event)], "bank")
bk <- merge(bk, unique(p1[, .(bank, ym, drop_event)]), by = c("bank", "ym"))
add_lead_growth(bk, "lval", by = "bank")
add_lag_growth(bk, "lval", by = "bank")
for (h in HORIZONS) bk[, paste0("g", h) := winsorize(get(paste0("g", h)))]

bs <- mk_constfx(p2[, .(bank, sector, cur, ym, total)], c("bank", "sector"))
bs <- merge(bs, unique(p2[, .(bank, ym, drop_event)]), by = c("bank", "ym"))
add_lead_growth(bs, "lval", by = c("bank", "sector"))
for (h in HORIZONS) bs[, paste0("g", h) := winsorize(get(paste0("g", h)))]
bs <- bs[fc + pyg > CELL_MIN_BASE]

# Attach shock and exposures (12-month lag, standardized; plan §1.5)
expo_z <- p4[, c("bank", "ym", paste0(names(EXPOSURES), "_l12_z"),
                 paste0(names(EXPOSURES), "_2016_z")), with = FALSE]
prep <- function(d) {
  d <- merge(d, macro[, .(ym, shk = shk_dxy_purged)], by = "ym", all.x = TRUE)
  d <- merge(d, expo_z, by = c("bank", "ym"), all.x = TRUE)
  for (v in names(EXPOSURES)) {
    d[, paste0("shk_", v) := shk * get(paste0(v, "_l12_z"))]
    d[, paste0("shk_", v, "_16") := shk * get(paste0(v, "_2016_z"))]
  }
  d[ym >= to_ym("2016-01-01") & ym <= to_ym("2025-12-01") &
    drop_event == FALSE & !is.na(shk)]
}
bk <- prep(bk); bs <- prep(bs)
bs[, `:=`(st = paste(sector, ym), bsec = paste(bank, sector))]
cat(sprintf("  Bank-level: %d rows; bank x sector: %d rows\n", nrow(bk), nrow(bs)))

# ---------------------------------------------------------------------------
# 2. Estimation: one exposure at a time, then jointly (h = 1..18)
# ---------------------------------------------------------------------------
BANK_CTRL <- c("tier1_l12_z", "size_l12_z", "liquidity_l12_z")

run_designB <- function(level = c("bank", "bank_sector"), expos, joint = FALSE,
                        wcb_h = H_KEY, B = WCB_B) {
  level <- match.arg(level)
  d  <- if (level == "bank") bk else bs
  fe <- if (level == "bank") "bank + ym" else "st + bsec"
  res <- list()
  for (h in HORIZONS) {
    yv <- paste0("g", h)
    for (v in expos) {
      param <- paste0("shk_", v)
      rhs <- if (joint) paste0("shk_", expos) else param
      if (level == "bank") rhs <- c(rhs, setdiff(BANK_CTRL, paste0(v, "_l12_z")))
      m <- tryCatch(
        feols(as.formula(paste(yv, "~", paste(rhs, collapse = "+"), "|", fe)),
              data = d, cluster = ~bank, notes = FALSE),
        error = function(e) NULL)
      if (is.null(m)) next
      ib <- inference_battery(m, param, dk_lag = h + 1)
      pw <- if (h %in% wcb_h) wcb_pval(d, yv, rhs, fe, param, B = B) else NA_real_
      res[[length(res) + 1]] <- data.table(
        level = level, spec = ifelse(joint, "joint", "single"),
        exposure = v, h = h, b = ib$b, se_bank = ib$se_bank,
        p_bank = ib$p_bank, p_wcb = pw, se_twoway = ib$se_twoway,
        p_twoway = ib$p_twoway, se_dk = ib$se_dk, p_dk = ib$p_dk, n = ib$n)
      if (joint) break   # joint spec: one model per h covers all exposures
    }
    if (joint) {
      # collect remaining coefficients from the joint model
      for (v in expos[-1]) {
        param <- paste0("shk_", v)
        ib <- inference_battery(m, param, dk_lag = h + 1)
        res[[length(res) + 1]] <- data.table(
          level = level, spec = "joint", exposure = v, h = h, b = ib$b,
          se_bank = ib$se_bank, p_bank = ib$p_bank, p_wcb = NA_real_,
          se_twoway = ib$se_twoway, p_twoway = ib$p_twoway,
          se_dk = ib$se_dk, p_dk = ib$p_dk, n = ib$n)
      }
    }
  }
  rbindlist(res)
}

cat("\n(ii) Bank x sector level with sector-time FE (preferred)...\n")
resB2 <- run_designB("bank_sector", names(EXPOSURES))
cat("\n(i) Bank level...\n")
resB1 <- run_designB("bank", names(EXPOSURES))
cat("\nJoint specifications...\n")
resBj <- rbind(run_designB("bank_sector", names(EXPOSURES), joint = TRUE, wcb_h = 12),
               run_designB("bank", names(EXPOSURES), joint = TRUE, wcb_h = 12))

resB <- rbind(resB1, resB2, resBj)
write_micro_csv(resB, "Micro_DesignB_Exposure_Results.csv")

key <- resB[h == 12 & spec == "single" & level == "bank_sector"]
cat("\nKey results (bank x sector, h = 12, single exposures):\n")
print(key[, .(exposure, b = round(b, 3), se = round(se_bank, 3),
              p_bank = round(p_bank, 3), p_wcb = round(p_wcb, 3))])

# 2016-fixed exposure robustness (primary exposure only)
cat("\nRobustness: exposures fixed at 2016 average (FC deposit share)...\n")
res16 <- list()
for (h in H_KEY) {
  m <- feols(as.formula(paste0("g", h, " ~ shk_fc_dep_share_16 | st + bsec")),
             data = bs, cluster = ~bank, notes = FALSE)
  ib <- inference_battery(m, "shk_fc_dep_share_16", dk_lag = h + 1)
  res16[[length(res16) + 1]] <- data.table(h = h, b = ib$b, se_bank = ib$se_bank,
                                           p_bank = ib$p_bank, n = ib$n)
}
write_micro_csv(rbindlist(res16), "Micro_DesignB_Fixed2016_Exposure.csv")

# ---------------------------------------------------------------------------
# 3. Figures
# ---------------------------------------------------------------------------
d1 <- resB[spec == "single" & level == "bank_sector"]
d1[, Exposure := EXPOSURES[exposure]]
pB <- ggplot(d1, aes(h, b)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_ribbon(aes(ymin = b - 1.645 * se_bank, ymax = b + 1.645 * se_bank),
              fill = "darkorange3", alpha = 0.2) +
  geom_line(color = "darkorange4", linewidth = 0.9) +
  facet_wrap(~Exposure, scales = "free_y") +
  labs(title = "Design B: bank lending channel - shock x exposure (bank x sector level)",
       subtitle = "Sector-time FE absorb credit demand; effect per 1 SD exposure per 1 SD dollar shock; 90% bands",
       x = "Horizon (months)", y = "Coefficient (pp)") +
  theme_micro()
save_png(pB, "263_DesignB_Exposure_Interactions.png", w = 11, h = 7)

plot_dynamic_beta(d1[exposure == "fc_dep_share"],
  "Design B: FC deposit share x dollar shock (primary liability-dollarization test)",
  "Bank x sector panel; sector-time + bank-sector FE; 90% bands (cluster by bank)",
  "264_DesignB_FCDeposit_Dynamic.png")

write_rds_micro(list(bk = bk, bs = bs), "micro_designB_panels.rds")
cat(sprintf("=== 32 done in %.1f min ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
