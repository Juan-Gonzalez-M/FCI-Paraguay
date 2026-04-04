# Financial Conditions and Credit Transmission in a Dual-Engine Economy: Evidence from Paraguay

**Author:** Juan Manuel Gonzalez Masulli
**Affiliation:** Banco Central del Paraguay; Universidad Catolica Nuestra Senora de la Asuncion
**Submitted to:** Emerging Markets Review

## Abstract

In dollarized commodity-exporting economies, the dollar index (DXY) --- not the VIX --- drives domestic financial conditions, operating through balance sheet channels that risk-appetite frameworks overlook. We construct and validate the first Financial Conditions Index (FCI) for Paraguay (January 1995--August 2025), a partially dollarized commodity exporter. Three findings emerge. First, the DXY explains 23.6% of FCI variation versus near-zero for the VIX; identification-robust Anderson-Rubin inference confirms a negative credit channel peaking at 6--12 months. Second, transmission is state-dependent on the commodity cycle: pro-cyclical leverage during booms amplifies the marginal FCI effect on credit growth from -3.7 pp to -10.8 pp. Third, finance-dependent sectors contract sharply while agriculture remains orthogonal, a dual-engine structure reconciling strong credit effects with weak aggregate output responses.

## Data

All input data is in `data/`:

| File | Description |
|------|-------------|
| `FCI_data_1.xlsx` | Main database: 12 financial variables, monthly, Jan 1995--Aug 2025 (368 obs) |
| `Variables_exo.xlsx` | External variables: Global Financial Conditions + Main Commodity Prices |
| `Output_puzzle.xlsx` | Sectoral output data for dual-engine analysis |
| `legacy/Data_FCI_Paraguay.xlsx` | Archived Nov 2020 version (not used in current analysis) |

## Replication

### Requirements

- **R version:** 4.5.1 (tested; earlier 4.x versions likely compatible)
- **Required packages:** readxl, dplyr, tidyr, zoo, FactoMineR, vars, MARSS, ggplot2, gridExtra, lmtest, sandwich, tseries, urca, quantreg, ivreg

Missing packages are installed automatically by the pipeline.

### Running the Full Pipeline

```bash
cd R/
Rscript RUN_ALL.R
```

Runtime: approximately 8--10 minutes. Produces ~171 PNG charts and ~127 CSV files in `output/png/` and `output/csv/`.

### Running Individual Scripts

```bash
cd R/
Rscript 01_FCI_Complete.R     # FCI construction only (~1 min)
Rscript 05_FCI_Local_Projections.R  # Local projections only
```

Scripts must be run from the `R/` directory. Script 01 must run first; subsequent scripts depend on its outputs.

## Project Structure

```
FCI-Paraguay/
├── README.md               # This file (replication guide)
├── R/                      # Analysis scripts
│   ├── CLAUDE.md           # Development reference (not needed for replication)
│   ├── RUN_ALL.R           # Master pipeline script
│   ├── 01_FCI_Complete.R   # FCI construction (must run first)
│   ├── 02-23_*.R           # Analysis scripts (21 total)
│   └── ...
├── data/                   # Input data
│   ├── FCI_data_1.xlsx
│   ├── Variables_exo.xlsx
│   └── Output_puzzle.xlsx
├── output/                 # Generated outputs
│   ├── png/                # Charts (~171 files)
│   ├── csv/                # Data tables (~127 files)
│   ├── pdf/                # Compiled paper
│   └── reports/            # Paper source (LaTeX + supporting reports)
└── docs/                   # Template files and diagrams
```

## R Scripts

| Script | Purpose |
|--------|---------|
| `01_FCI_Complete.R` | FCI construction at all levels (comprehensive, endo/exo, channels, purified) |
| `02_FCI_Effects_Analysis.R` | Granger causality, IRF, out-of-sample forecasting |
| `03_FCI_Stationarity_Tests.R` | ADF, PP, KPSS unit root tests |
| `05_FCI_Local_Projections.R` | Local Projections for credit channel and macro effects |
| `06_FCI_Rolling_Stability.R` | Rolling/expanding window PCA loading stability |
| `07_FCI_Robustness_Endogeneity.R` | Endogeneity tests (FCI_exNPL, FCI_exCredit) |
| `08_FCI_Regime_TVP_Analysis.R` | Regime-switching (IT/COVID), time-varying parameters |
| `09_FCI_Monetary_Policy_Interaction.R` | FCI x TPM interaction, systematic forecasting |
| `10_FCI_Growth_at_Risk.R` | Growth-at-Risk quantile regressions |
| `11_FCI_Block_SVAR.R` | Block-exogenous SVAR for external spillovers |
| `12_FCI_Output_Puzzle_Investigation.R` | Output puzzle: IMAEP_SANB LP, VECM, transmission channels |
| `14_FCI_Output_Puzzle_Sectoral.R` | Sectoral credit/output decomposition |
| `15_FCI_PostIT_Subsample_Analysis.R` | Post-IT subsample re-estimation |
| `16_FCI_TCN_Reclassification.R` | TCN reclassification sensitivity |
| `17_FCI_New_External_Data.R` | GFC index, regional proxy, disaggregated commodities |
| `18_FCI_Improved_IV_LP.R` | Expanded IV-LP with Anderson-Rubin inference |
| `19_FCI_Proxy_SVAR.R` | Proxy-SVAR (Mertens & Ravn 2013) |
| `20_FCI_Enriched_Block_SVAR.R` | 7-variable external block SVAR, historical decomposition |
| `21_FCI_Commodity_Puzzle.R` | Disaggregated commodity effects, ToT interaction |
| `22_FCI_Regional_Spillovers.R` | Brazil Selic spillover analysis |
| `23_FCI_Identification_Triangulation.R` | DXY IV-LP, GFC-PCA IV, globally-purged FCI |

## Paper-to-Output Mapping

Key tables and figures in the paper and their generating scripts/output files:

| Paper Element | Output File(s) | Script |
|---|---|---|
| **Table 1** (Variables) | `csv/FCI_Summary_Statistics.csv` | `01` |
| **Table 2** Panel A (Post-IT LP) | `csv/PostIT_LP_Credit.csv` | `15` |
| **Table 2** Panel B (Specifications) | `csv/PostIT_LP_Credit.csv`, `csv/LP_Credit_Standard.csv`, `csv/Triangulation_Globally_Purged_LP.csv` | `05`, `15`, `23` |
| **Table 2** Panel C (Currency) | `csv/LP_Credit_Standard.csv`, `csv/PostIT_LP_Credit.csv` | `05`, `15` |
| **Table 3** (First-stage) | `csv/IV_First_Stage_Battery.csv` | `18` |
| **Table 4** (Anderson-Rubin) | `csv/Triangulation_DXY_IV_LP.csv` | `23` |
| **Table 5** (FCI x ToT) | `csv/Commodity_Credit_Interaction_LP.csv` | `21` |
| **Table 6** (Sectoral output) | `csv/LP_Sectoral_Output_Full.csv` | `14` |
| **Table 7** (Robustness) | Multiple CSVs (see notes below) | Multiple |
| **Figure 2** (FCI time series) | `png/01_FCI_Methods_Comparison.png` | `01` |
| **Figure 3** (Triangulation) | `png/251_Identification_Triangulation.png` | `23` |
| **Figure 4** (Marginal effects) | `png/234_Marginal_Effects_FCI_Credit.png` | `21` |
| **Figure 5** (Sectoral output) | `png/143_LP_Sectoral_Output_Full.png` | `14` |

**Table 7 sources:** `csv/LP_Credit_Standard.csv` (baseline), `csv/PostIT_LP_Credit.csv` (post-IT), `csv/LP_Method_by_Method.csv` (rates-only, ENDO variants), `csv/TCN_Reclass_LP_Comparison.csv` (TCN reclassified), `csv/Triangulation_Globally_Purged_LP.csv` (purged), `csv/LP_Placebo_Test.csv` (placebo), `csv/VAR_Ordering_Robustness.csv` (VAR ordering).

## Notes

- The paper source is `output/reports/FCI_Paraguay_EMR_Submission.tex`. Compile with `pdflatex` + `bibtex` using the Elsevier CAS single-column template (included in `docs/emr-template/`).
- `R/CLAUDE.md` is a development reference for AI-assisted coding and is not needed for replication.
- Positive FCI values indicate tighter financial conditions throughout.
- The FCI_CORE index (12 variables, full sample) is the primary index; FCI_FULL (13 variables, includes EMBI from 2013+) is used for robustness.
