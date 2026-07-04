# ============================================================================
# RUN_MICRO.R — Master script for the bank-level micro phase (plan:
# Bank_Level_Micro_Plan_IREF_Revision.md). Runs scripts 30-35 in order, each
# in an isolated R session. NEW PHASE: does not touch the aggregate pipeline
# (RUN_ALL.R, scripts 01-23) or any existing output.
#
# Usage:  cd R/ && Rscript RUN_MICRO.R
# Env:    MICRO_WCB_B=199 Rscript RUN_MICRO.R   # quick pass (fewer bootstrap reps)
# Output: ../output/micro/{csv,png,rds}
# ============================================================================

scripts <- c(
  "30_Bank_Micro_Data_Construction.R",
  "31_Bank_Micro_DesignA_Currency.R",
  "32_Bank_Micro_DesignB_Exposure.R",
  "33_Bank_Micro_DesignC_Hedging.R",
  "34_Bank_Micro_DesignD_Mechanisms.R",
  "35_Bank_Micro_Robustness.R"
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
  if (rc != 0) {
    cat("\n*** ", s, "failed - stopping the pipeline. ***\n")
    break
  }
}

cat("\n============================================================\n")
print(log, row.names = FALSE)
cat(sprintf("Total: %.1f min\n", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
write.csv(log, "../output/micro/csv/Micro_Run_Log.csv", row.names = FALSE)
