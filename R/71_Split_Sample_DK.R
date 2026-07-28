# ============================================================================
# 71_Split_Sample_DK.R
# Driscoll-Kraay standard errors for the SPLIT-SAMPLE currency differential —
# the rows main-text Table 7 reports and Figure 3(b) displays.
#
# WHY THIS EXISTS. Both referee-side audits flagged an apparent clash between
# Figure 3's caption ("90% intervals under the declared primary convention,
# Driscoll-Kraay") and Table 7's note ("secondary rows for which the primary
# convention is not computed; brackets denote bank-clustered SEs"). Neither
# document is wrong about itself; the two exhibits display DIFFERENT OBJECTS:
#
#   Table 7 split rows  = SEPARATE-SAMPLE estimates. 33_Bank_Micro_DesignC_
#     Hedging.R:116-129 runs g_h ~ shk_usd + CTRLS | bst + bsc separately on the
#     hedged and unhedged subsamples and stores ONLY se_bank in
#     Micro_DesignC_Split_Sample.csv. Values quoted in Section 5.3: +2.06
#     (hedged, h=12, n=5,298) against +0.43 (unhedged, n=12,318). For these
#     rows the primary convention genuinely was NOT computed, so Table 7's note
#     is accurate.
#   Figure 3(b)          = the INTERACTED "sharp" variant from
#     Micro_DesignC_Hedging_Results.csv (b = +2.072, se_dk = 2.073), which does
#     carry Driscoll-Kraay SEs (61_Submission_Figure9_DesignAligned.R:27 uses
#     se_dk), so the caption is accurate too.
#
# The real defect is therefore neither the caption nor the note but the fact
# that the figure and the table show different estimates of the same idea while
# reading as though they show one. The fix adopted is to compute the primary
# convention for the object the TABLE reports, so that Table 7's split rows and
# Figure 3(b) can be placed on one convention and one object.
#
# Specification is copied verbatim from 33:116-129; the only addition is
# inference_battery(), which supplies two-way bank x month and
# Driscoll-Kraay(h+1) alongside the bank-clustered SE already published.
# The bank-clustered column is reproduced as a gate: it must match
# Micro_DesignC_Split_Sample.csv exactly.
#
# ADDITIVE: writes only Micro_DesignC_Split_DK.csv.
# ============================================================================

t0 <- Sys.time()
source("micro_helpers.R")
library(fixest)
set.seed(20260710)
cat("=== 71: Driscoll-Kraay SEs for the split-sample differential ===\n")

dA <- read_rds_micro("micro_designA_panel.rds")

CTRLS        <- c("glag_usd", "tier1_usd", "size_usd", "fcdep_usd")
FE_A         <- "bst + bsc"
SECT_TINY_FX <- c("CONSUMO", "VIVIENDA")        # 33:44

dA[, hedge_class := fifelse(sector %in% SECT_HEDGED,   "hedged",
                    fifelse(sector %in% SECT_UNHEDGED, "unhedged",
                    fifelse(sector %in% SECT_AMBIGUOUS, "ambiguous", "excluded")))]

res <- list()
for (grp in c("hedged", "unhedged")) {
  dg <- dA[hedge_class == grp & (!sector %in% SECT_TINY_FX | total > 5000)]
  for (h in HORIZONS) {
    ds <- dg[get(paste0("ok", h)) == TRUE]
    m <- tryCatch(
      feols(as.formula(paste0("g", h, " ~ shk_usd + ", paste(CTRLS, collapse = "+"),
                              " | ", FE_A)), data = ds, cluster = ~bank,
            notes = FALSE),
      error = function(e) NULL)
    if (is.null(m)) next
    ib <- inference_battery(m, "shk_usd", dk_lag = h + 1)
    res[[length(res) + 1]] <- data.table(
      group = grp, h = h, b = ib$b,
      se_bank = ib$se_bank, p_bank = ib$p_bank,
      se_twoway = ib$se_twoway, p_twoway = ib$p_twoway,
      se_dk = ib$se_dk, p_dk = ib$p_dk, n = ib$n)
  }
}
out <- rbindlist(res)
write_micro_csv(out, "Micro_DesignC_Split_DK.csv")

cat("\n--- Headline horizons (h = 6, 12, 18) ---\n")
print(out[h %in% c(6, 12, 18),
          .(group, h, b = round(b, 3),
            se_bank = round(se_bank, 3), p_bank = round(p_bank, 4),
            se_dk = round(se_dk, 3), p_dk = round(p_dk, 4), n)])

# ---- Gate: bank-clustered column must reproduce the published split file ----
pub <- fread(file.path(micro_paths$out_csv, "Micro_DesignC_Split_Sample.csv"))
chk <- merge(out[, .(group, h, b_new = b, se_new = se_bank)],
             pub[, .(group, h, b_pub = b, se_pub = se_bank)], by = c("group", "h"))
chk[, `:=`(db = abs(b_new - b_pub), dse = abs(se_new - se_pub))]
cat(sprintf("\nGATE vs Micro_DesignC_Split_Sample.csv: max |db| = %.3e, max |dse| = %.3e over %d rows\n",
            max(chk$db), max(chk$dse), nrow(chk)))
if (max(chk$db) < 1e-10 && max(chk$dse) < 1e-10) {
  cat("GATE PASSED: specification reproduced exactly; the DK columns are additions only.\n")
} else {
  cat("GATE FAILED: specification does not reproduce - do NOT quote the DK columns.\n")
}

cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
