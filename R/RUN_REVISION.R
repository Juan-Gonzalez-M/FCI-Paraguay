# ============================================================================
# RUN_REVISION.R — Master script for the revision-extras phase (scripts 37-44,
# coworker-advised exercises). Runs each script in an isolated R session.
# ADDITIVE: touches nothing outside output/revision/.
#
# Usage:  cd R/ && Rscript RUN_REVISION.R
# Env:    MICRO_WCB_B=99 Rscript RUN_REVISION.R   # fewer bootstrap reps (37)
# ============================================================================

scripts <- c(
  "37_Micro_Freeze_Checks.R",
  "38_CHR_Plausibly_Exogenous.R",
  "39_PostIT_ExpandingFCI_IV.R",
  "40_LagAugmented_LP.R",
  "41_Falsification_Placebos.R",
  "42_FCI_Component_Exclusion.R",
  "43_External_Data_Fetch.R",
  "44_Enhanced_Instrument_Robustness.R",
  "45_COVID_Reclassification_Check.R",
  "46_USD_ShiftShare_Decomposition.R"
)

t_all <- Sys.time()
log <- data.frame(script = scripts, status = NA_character_, minutes = NA_real_)

for (i in seq_along(scripts)) {
  s <- scripts[i]
  cat("\n############################################################\n")
  cat("##  Running", s, "\n")
  cat("############################################################\n")
  t1 <- Sys.time()
  rc <- system2("Rscript", shQuote(s))
  log$minutes[i] <- round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 2)
  log$status[i]  <- ifelse(rc == 0, "OK", paste0("FAILED (exit ", rc, ")"))
  if (rc != 0) cat("\n*** ", s, "failed — continuing with remaining scripts. ***\n")
}

cat("\n============================================================\n")
print(log, row.names = FALSE)
cat(sprintf("Total: %.1f min\n", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
write.csv(log, "../output/revision/csv/Rev_Run_Log.csv", row.names = FALSE)
