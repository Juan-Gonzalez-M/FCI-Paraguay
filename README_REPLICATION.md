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

## 2b. Building the submission documents

**Manuscript** (run twice so cross-references resolve):

```bash
cd output/submission/
pdflatex -interaction=nonstopmode FCI_Paraguay_IREF_Submission.tex
pdflatex -interaction=nonstopmode FCI_Paraguay_IREF_Submission.tex
pdflatex -interaction=nonstopmode Cover_Letter_IREF.tex
```

**Online appendix and response letter** — both use the same pandoc invocation:

```bash
cd output/reports/
pandoc FCI_Paraguay_Online_Appendix.md \
  -o ../pdf/FCI_Paraguay_Online_Appendix.pdf \
  --pdf-engine=xelatex -V geometry:margin=1in -V fontsize=11pt

cd ../../private/
pandoc Response_to_Reviewers_DRAFT.md \
  -o Response_to_Reviewers.pdf \
  --pdf-engine=xelatex -V geometry:margin=1in -V fontsize=11pt
```

Four details are load-bearing and easy to get wrong:

- **Run pandoc from `output/reports/`, never from the repository root.** The appendix embeds 18 figures via `../../output/png/...` and `../submission/...` relative paths. From the root none of them resolve, and pandoc *succeeds silently*: the result is a 66-page, 289 KB PDF with every figure missing instead of the correct 79-page, 4.86 MB file. Check the page count **and** the file size after every appendix rebuild.
- **Do not pass `--toc`.** The source carries its own hand-written Contents table; `--toc` adds a second, generated one and inflates the document by ~3 pages.
- **`-V fontsize=11pt` is required.** Omitting it silently falls back to 10pt (leading 12.0 pt instead of 13.6 pt) and changes the pagination.
- **`-V geometry:margin=1in` is required.** Without it the `article` defaults give a much narrower text block.

The YAML `header-includes` block at the top of the `.md` is picked up automatically and must be preserved. It now carries `\usepackage{float}` and `\floatplacement{figure}{H}` as well as `\usepackage[labelformat=empty]{caption}`: pinning figures in place prevents a float from drifting across a section boundary (Figure M.1 previously migrated into Appendix N in the typeset output, splitting Table N.3). Note that `placeins` is *not* available in the TinyTeX distribution used here; `float` is.

Verification: the appendix should build to **letter size, ~92 pages and >4.5 MB** (a page count in the 60s with <300 KB means the figures did not resolve --- see the working-directory item above) and the response letter to **~9 pages**, with `pdfinfo` reporting `Creator: LaTeX via pandoc` / `Producer: xdvipdfmx`. Rebuilding the frozen V2 sources (`output/archive/IREF_V2_2026-07-17/`) with these commands reproduces both shipped V2 PDFs with byte-identical text layers (66 and 16 pages) — use that as a regression check if the settings are ever in doubt.

The stale `private/Response_to_Reviewers_DRAFT.tex` and `.log` are leftovers from an earlier build path and are not used; the `.md` is the single source.

Interval cells in appendix tables are written in math mode (`$[-4.15,\,-0.19]$`) because LaTeX otherwise line-breaks at the hyphen of a negative bound; a non-breaking space does not prevent this.

Note on the PDF format: pandoc/xelatex writes PDF 1.7 with a cross-reference *stream* rather than a classic `xref` table. Naive byte-level checks for `/Catalog`, `trailer`, or `%%EOF` will therefore report the file as malformed; it is not. Validate with `pdfinfo` and `pdftotext`, both of which parse it correctly.

## 3. Seeds

All bootstrap and permutation results are seeded for exact reproducibility:
`set.seed(20260703)` (script 31, Design A battery — the canonical Table O.1 source), `set.seed(20260712)` (scripts 33/35/48), `set.seed(20260713)` (script 50), `set.seed(20260716)` (scripts 56/58/59). Timing permutations use 500 draws; the null-imposed joint moving-block reduced-form bootstrap uses 999 draws (block length 18 months).

## 4. Exhibit-to-script manifest (main text)

Verified against the compiled manuscript on 29 July 2026: **8 tables, 4 figures**. (Earlier versions of this manifest listed 10 tables and 6 figures, and included a Conley--Hansen--Rossi bounds figure that has since moved to the online appendix; corrected here.) Note that the submission PNG filenames are not in exhibit order.

| Exhibit | Content | Script(s) | Key outputs |
|---|---|---|---|
| Table 1 | Institutional features | (hand-compiled from cited sources) | — |
| Table 2 | FCI variables and sign conventions | `01_FCI_Complete.R` | `FCI_Complete_Results.csv` |
| Table 3 | Identification strategy summary (three principal designs) | — (summary of below) | — |
| Table 4 | Credit LP (one-sided primary; Panels B–C) | `52`, `54`, `15`, `05` | `Rev_OneSided_Matched.csv`, `Rev_Canonical_Pipeline.csv`, `PostIT_LP_Credit.csv` |
| Table 5 | Candidate-instrument screening | `18_FCI_Improved_IV_LP.R` | `IV_First_Stage_Battery.csv` |
| Table 6 | Aligned reduced form, AR sets, bootstrap (FCI_IV) | `56_Aligned_Identification_Chain.R`, `59_Aligned_HAC_AR.R` | `Rev_Aligned_ReducedForm.csv`, `Rev_Aligned_HAC_AR.csv`, `Rev_Aligned_RF_Bootstrap.csv` |
| Table 7 | Bank-level within-cell estimates, exposure gradient, matched VIX placebo and difference test (**Driscoll–Kraay declared primary**) | `31`–`33`, `48`, `64` | `Micro_DesignA_Main.csv`, `Micro_Split_WCB.csv`, `Micro_Gradient_*.csv`, `Micro_VIXGradient_*.csv`, `Micro_Gradient_DXY_vs_VIX_Diff.csv` |
| Table 8 | Group-level sectoral tests (finance-dependent vs insulated) | `27_FCI_Published_Quarterly_LP.R`, `63_Sectoral_Group_Wald.R`, **`67_Sector_Group_Wald_OneSided.R`** | `Output_Puzzle_Quarterly_Transmission.csv`, `Rev_Sector_Group_Wald.csv`, `Rev_Sector_Group_Wald_OneSided.csv` |
| Figure 1 | Composite FCI with stress episodes | `62_Submission_Figure1.R` | `Figure_1.png` |
| Figure 2 | Credit IRFs (primary one-sided + composite) | `60_Submission_Figure3_TwoPanel.R` | `Figure_3.png` |
| Figure 3 | Bank-level currency-composition margin (Driscoll–Kraay primary bands) | `61_Submission_Figure9_DesignAligned.R` | `Figure_9.png` |
| Figure 4 | Component-exclusion ladder | `42_FCI_Component_Exclusion.R` | `Figure_10.png` |

**Relocated to the online appendix in the round-17 revision** (the manuscript no longer carries these as main-text exhibits):

| Appendix exhibit | Content | Script | Key outputs |
|---|---|---|---|
| Table N.6a | FCI x ToT interaction — predetermined state (was main-text Table 9) | `66_Submission_Figure5_Predetermined.R` | `Rev_ToT_Interaction_BothStates.csv` |
| Figure N.2 | ToT marginal effects — predetermined state (was main-text Figure 5) | `66_Submission_Figure5_Predetermined.R` | `Rev_Marginal_Effects_Predetermined.csv`, `Figure_5.png` |
| Table E.0d | Pointwise **and simultaneous (sup-$t$)** bands, h = 1--18, both samples | `70_OneSided_Band_SDs_CurrencyTest.R` | `Rev_OneSided_Profile_Band.csv` |
| Table E.3c | Full-sample stacked MN-vs-USD currency difference | `70_OneSided_Band_SDs_CurrencyTest.R` | `Rev_Currency_Equality_Test.csv` |
| Tables K.5b(i)--(ii) | Persistence-matched candidate screening (3 transformations x 2 lag conventions x 7 candidates) | `69_Screening_Transformation_Matched.R` | `Rev_Screening_RankSummary.csv`, `Rev_Screening_TransformationMatched.csv` |
| Table O.4a, DK column | Driscoll--Kraay p-values for the split rows (primary convention) | `71_Split_Sample_DK.R` | `Micro_DesignC_Split_DK.csv` |
| Table F.8e | Group-equality test under the one-sided timing standard | **`67_Sector_Group_Wald_OneSided.R`** | `Rev_Sector_Group_Wald_OneSided.csv` |
| Appendix A | FCI descriptive moments (were Table 2 columns) | `01_FCI_Complete.R` | `FCI_Complete_Results.csv` |
| Appendix F | Descriptive individual-sector rows (were Table 10 Panel A) | `27_FCI_Published_Quarterly_LP.R` | `Output_Puzzle_Quarterly_Transmission.csv` |

The FCI variant taxonomy (previously main-text Table 3, never cross-referenced) is now carried solely by the estimand dictionary at the front of the online appendix.

### 4b. Scripts added in the round-15 and round-17 revisions

| Script | Purpose | Outputs |
|---|---|---|
| `64_Gradient_VIX_Placebo.R` | Matched VIX x USD x 2016-FX-share exposure gradient (full inference battery); joint DXY-vs-VIX difference test; export of the 2016 baseline sector FX shares | `Micro_VIXGradient_TimeRobust.csv`, `_SectorRobust.csv`, `_FewCluster.csv`, `_LOSO.csv`, `Micro_Gradient_DXY_vs_VIX_Diff.csv`, `Micro_Sector_FX_Shares_2016.csv` |
| `65_Sector_Group_Stability.R` | Group-equality test under leave-one-component-out, two alternative classifications, and the post-IT subsample | `Rev_Sector_Group_Stability.csv` |
| `66_Submission_Figure5_Predetermined.R` | Predetermined-ToT interaction LP with the HAC covariance block; delta-method marginal effects; appendix Table N.6a / Figure N.2 and the contemporaneous counterpart | `Rev_ToT_Interaction_BothStates.csv`, `Rev_Marginal_Effects_Predetermined.csv`, `Figure_5.png`, `Figure_A_ToT_Contemporaneous.png` |
| **`67_Sector_Group_Wald_OneSided.R`** | Sectoral group-equality Wald test re-estimated on timing-clean one-sided indices (`FCI_COMP_ZSCORE`, `FCI_exCredit_ZSCORE`) alongside the published composite. Gated: the composite block must reproduce `Rev_Sector_Group_Wald.csv` to machine precision or the script stops. Reports per-SD as well as raw differences, since index scales differ by a factor of two | `Rev_Sector_Group_Wald_OneSided.csv` |
| **`count_words.sh`** | The manuscript word-count rule quoted in the cover and response letters: `pdftotext -layout`, every word before the `References` heading. Resolves paths from its own location, so it runs from `R/`, from the repository root, or by absolute path; prints words, total pages and main-text pages for the archived V1, the previous revision and the current submission | (stdout) |

Seeds: 64 reuses script 48's seeds (20260710 / 20260712) so the dollar and VIX gradients are compared under identical bootstrap draws. Canonical bootstrap size is `MICRO_WCB_B=9999`; the moving-block bootstrap is fixed at B = 499 and, being drawn after the wild bootstraps, is sensitive to `MICRO_WCB_B` through the RNG stream — reproduce it at the canonical setting.

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

**One series cannot be redistributed.** The DXY (ICE US Dollar Index) comes from Bloomberg under a terminal licence, monthly-averaged from daily observations, and lives in the `Global_Financial_Conditions` sheet of `data/FCI_data_1.xlsx`. Users without Bloomberg access can reproduce every dollar exercise on the public Federal Reserve broad dollar index (`DTWEXBGS` spliced with `DTWEXB`), which script 43 downloads automatically; Appendix K.6 reports the full instrument battery under that measure.
