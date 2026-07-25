# Replication README — Global Dollar Strength and Domestic Credit Conditions in a Partially Dollarized Economy (IREF-D-26-04765)

This file is the operational replication manifest for the manuscript and its online appendix. The online appendix documents definitions, specifications, and results; everything operational (scripts, output files, seeds, commands) lives here.

## 1. Software

- **R** (>= 4.2.0)
- Packages: `readxl`, `dplyr`, `tidyr`, `zoo`, `FactoMineR`, `vars`, `MARSS`, `ggplot2`, `gridExtra`, `lmtest`, `sandwich`, `tseries`, `urca`, `quantreg`, `ivreg`; Phases 2–3 additionally `data.table`, `fixest`, `modelsummary`, `ARDL`, `car`, `strucchange`.
- LaTeX (pdflatex) for the manuscript; pandoc + xelatex for the appendix PDF.

## 2. Running the pipeline

```bash
cd R/
Rscript RUN_ALL.R                             # Phase 1: aggregate pipeline (~8-10 min)
MICRO_WCB_B=9999 Rscript RUN_MICRO.R          # Phase 2: bank-level analysis (canonical bootstrap; ~4-5 h)
Rscript 36_Aggregate_FXAdjusted_Credit_LP.R   # constant-exchange-rate LPs (<1 min)
Rscript RUN_REVISION.R                        # Phase 3: identification/robustness extensions (~15 min)
```

Phase 1 writes to `output/png/` and `output/csv/`; Phase 2 to `output/micro/`; Phase 3 to `output/revision/`. Phases 2–3 are additive and require Phase 1 outputs. `MICRO_WCB_B` controls wild-cluster-bootstrap replications; all published bank-level bootstrap p-values use the canonical `MICRO_WCB_B=9999`. A quick pass (`MICRO_WCB_B=199`) reproduces coefficients exactly and bootstrap p-values to Monte Carlo error.

## 3. Seeds

All bootstrap and permutation results are seeded for exact reproducibility:
`set.seed(20260703)` (script 31, Design A battery — the canonical Table O.1 source), `set.seed(20260712)` (scripts 33/35/48), `set.seed(20260713)` (script 50), `set.seed(20260716)` (scripts 56/58/59). Timing permutations use 500 draws; the null-imposed joint moving-block reduced-form bootstrap uses 999 draws (block length 18 months).

## 4. Exhibit-to-script manifest (main text)

| Exhibit | Content | Script(s) | Key outputs |
|---|---|---|---|
| Table 1 | Institutional features | (hand-compiled from cited sources) | — |
| Table 2 | FCI variables and signs | `01_FCI_Complete.R` | `FCI_Complete_Results.csv` |
| Table 3 | FCI variant taxonomy (incl. one-sided and FCI_OS rows) | `01`, `52`, `56` | `FCI_Robustness_Versions.csv`, `Rev_Aligned_Construction.csv` |
| Table 4 | Identification strategy summary | — (summary of below) | — |
| Table 5 | Credit LP (one-sided primary; panels B–C) | `52`, `54`, `15`, `05` | `Rev_OneSided_Matched.csv`, `Rev_Canonical_Pipeline.csv`, `PostIT_LP_Credit.csv` |
| Table 6 | Candidate-instrument screening (composite) | `18_FCI_Improved_IV_LP.R` | `IV_First_Stage_Battery.csv` |
| Table 7 | Aligned reduced form, AR sets, bootstrap (FCI_OS) | `56_Aligned_Identification_Chain.R`, `59_Aligned_HAC_AR.R` | `Rev_Aligned_ReducedForm.csv`, `Rev_Aligned_HAC_AR.csv`, `Rev_Aligned_RF_Bootstrap.csv` |
| Table 8 | Bank-level within-cell estimates | `31`–`33`, `48` | `Micro_DesignA_Main.csv`, `Micro_Split_WCB.csv`, `Micro_Gradient_*.csv` |
| Table 9 | FCI x ToT interaction | `21_FCI_Commodity_Puzzle.R`, `58` (predetermined state) | `Commodity_Credit_Interaction_LP.csv`, `Rev_S3_ToT_Predetermined.csv` |
| Table 10 | Sectoral output responses + group tests | `27_FCI_Published_Quarterly_LP.R`, `57_Sectoral_Group_Tests.R` | `Output_Puzzle_Quarterly_Transmission.csv`, `Rev_Sector_GroupTests.csv`, `Rev_Sector_Family_FDR.csv` |
| Table 11 | Robustness summary | `42`, `58` (sign reversal, weights) | `Rev_FCI_LOO_Battery.csv`, `Rev_S5_*.csv` |
| Figure 1 | FCI history | `01` | — |
| Figure 2 | Credit IRFs (primary one-sided + composite) | `60_Submission_Figure3_TwoPanel.R` | `Figure_3.png` |
| Figure 3 | Post-IT IV-LP (expanding index) | `39_PostIT_ExpandingFCI_IV.R` | `Rev_PostIT_IV_ExpandingFCI.csv`, `Figure_8.png` |
| Figure 4 | Conley bounds (aligned) | `56` (breakdown), `49` Part 5 (joint bootstrap) | `Rev_Aligned_CHR_Breakdown.csv`, `Rev_CHR_Calibration_JointBootstrap.csv` |
| Figure 5 | Bank-level margins (design-aligned bands) | `61_Submission_Figure9_DesignAligned.R` | `Figure_9.png` |
| Figures 6–7 | ToT marginal effects; component-exclusion ladder | `51`, `42` | `Figure_5.png`, `Figure_10.png` |

## 5. Exhibit-to-script manifest (online appendix, selected)

- **A.5 timing conventions; one-sided co-baseline**: `52_OneSided_FCI_Cobaseline.R` → `Rev_OneSided_*.csv`
- **A.6 effective weights and sign sensitivity**: `58_Small_Robustness_Items.R` → `Rev_S5_Effective_Weights.csv`, `Rev_S5_SignReversal_LP.csv`
- **A.7 component-exclusion ladder**: `42_FCI_Component_Exclusion.R` → `Rev_FCI_LOO_Battery.csv`
- **E.0a cumulative log-level LP**: `58` → `Rev_S2_CumLogLevel_LP.csv`
- **E.0b lag-augmented LP**: `40_LagAugmented_LP.R` → `Rev_LagAugmented_LP*.csv`
- **E.3b currency-migration test**: `53_Currency_Migration_Test.R` → `Rev_Currency_Migration_*.csv`
- **F.2/F.6 falsification battery**: `41_Falsification_Placebos.R` → `Rev_Falsification_*.csv`
- **F.8 group-level sectoral tests**: `57_Sectoral_Group_Tests.R` → `Rev_Sector_*.csv`
- **K.1b–K.1c IV/AR scaling audit**: `49_IV_AR_Audit.R` → `Rev_IV_AR_Audit*.csv`
- **K.1d persistence battery (incl. Gregory-Hansen, lag grid, aligned rows, RF bootstrap)**: `50_Persistence_Valid_Aggregate.R`, `56`, `59` → `Rev_PV_*.csv`, `Rev_Aligned_*.csv`
- **K.6–K.7 instrument battery, CHR calibration detail**: `44_Enhanced_Instrument_Robustness.R`, `38`, `49` → `Rev_Instrument_Comparison.csv`, `Rev_CHR_*.csv`
- **N.4 matched ex-TCN exclusion**: `58` → `Rev_S4_ExTCN_Matched_*.csv`
- **N.6 rolling/expanding table; predetermined ToT**: `39`, `58`
- **O.1–O.7 bank-level construction, batteries, shift-share**: `30`–`36`, `45`–`48` → `Micro_*.csv`, `Rev_ShiftShare_*.csv`, `Rev_S6_Exposure_Timing.csv`
- **B.7 transmission schematic**: `docs/transmission_diagram.svg` → `output/png/transmission_diagram.png` (rsvg-convert)

## 6. Data

See Online Appendix Q (Data Sources and Access). `output/revision/external/` caches FRED downloads (fetched by script 43).
