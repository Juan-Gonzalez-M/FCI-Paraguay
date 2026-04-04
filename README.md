# Financial Conditions and Credit Transmission in a Dual-Engine Economy: Evidence from Paraguay

**Author:** Juan Manuel Gonzalez Masulli  
**Affiliation:** Banco Central del Paraguay; Universidad Catolica Nuestra Senora de la Asuncion  
**Contact:** jgonzalezm@bcp.gov.py

## Overview

This repository contains the data, code, and generated outputs for the construction and analysis of the first Financial Conditions Index (FCI) for Paraguay. The FCI combines 12 financial variables using four methodologies (Z-Score, PCA, VAR, Dynamic Factor Model) to produce a composite monthly index spanning January 1996 to December 2025 (360 observations).

The analysis provides evidence that:

1. **Dollar Channel Dominance:** The US Dollar Index (DXY), not US interest rates or global risk proxies (VIX), is the primary external driver of Paraguay's domestic financial conditions, operating through the balance sheet channel in a partially dollarized economy (~40% of bank credit in USD).
2. **Credit Transmission:** A one standard deviation FCI tightening is associated with a 5.6 pp decline in real credit growth at 12 months (post-IT sample), confirmed by convergent evidence from four identification strategies (OLS-LP, IV-LP with Anderson-Rubin inference, Proxy-SVAR, Block-SVAR).
3. **Commodity Amplification:** Transmission is state-dependent on the commodity cycle: pro-cyclical leverage during commodity booms amplifies the credit effect of financial tightening.
4. **Dual-Engine Structure:** Finance-dependent sectors (investment, construction, services) contract sharply under tightening, while agricultural output shows no persistent response --- helping reconcile strong credit effects with weak aggregate output responses commonly observed in commodity exporters.

## Data

All input data are contained in a single Excel file:

| File | Content | Observations |
|------|---------|--------------|
| `data/FCI_data_1.xlsx` | 6 sheets: FCI input variables, macro controls, global financial conditions, commodity prices, quarterly national accounts, monthly sectoral activity | 360 monthly (Jan 1996--Dec 2025); 127 quarterly (1994Q1--2025Q3) |

**Sheets:** Main_variables (12 FCI inputs + credit/deposit stocks), Datos_macro (IMAEP, CPI, credit), Global_Financial_Conditions (US 10Y, S&P 500, DXY, Selic, ToT), Main_Commodities_Prices (9 commodities), Quarterly_SA (GDP and expenditure components), Monthly_SA (sectoral activity indices).

## Replication

### Requirements

- **R version:** 4.3+ (tested on 4.5.1)
- **Required packages:** readxl, dplyr, tidyr, zoo, FactoMineR, vars, MARSS, ggplot2, gridExtra, lmtest, sandwich, tseries, urca, quantreg, ivreg

Missing packages are installed automatically by the pipeline.

### Running the Full Pipeline

```bash
cd R/
Rscript RUN_ALL.R
```

Runtime: approximately 8--10 minutes. Produces ~144 PNG charts and ~127 CSV data files in `output/png/` and `output/csv/`.

### Running Individual Scripts

```bash
cd R/
Rscript 01_FCI_Complete.R          # FCI construction only (~1 min)
Rscript 05_FCI_Local_Projections.R # Local projections only
```

Script `01_FCI_Complete.R` must run first; all subsequent scripts depend on its outputs. Scripts must be executed from the `R/` directory.

## Project Structure

```
FCI-Paraguay/
├── README.md
├── R/                          # Analysis scripts (28 files)
│   ├── RUN_ALL.R               # Master pipeline script
│   ├── 01_FCI_Complete.R       # FCI construction (must run first)
│   ├── 02--23_*.R              # Core analysis scripts (see table below)
│   ├── 24--28_*.R              # Verification and revision computations
│   └── CLAUDE.md               # Development reference (not needed for replication)
├── data/
│   └── FCI_data_1.xlsx         # Single input file (6 sheets)
├── output/                     # Generated outputs
│   ├── csv/                    # Data tables (~127 files)
│   └── png/                    # Charts (~144 files)
└── docs/
    └── emr-template/           # Elsevier CAS LaTeX template files
```

> **Note:** The manuscript source (LaTeX), compiled PDFs, and submission materials are maintained separately and are not included in this repository. The compiled paper and online appendix are available upon request from the corresponding author.

## Analysis Scripts

| Script | Purpose |
|--------|---------|
| `01_FCI_Complete.R` | FCI construction at all levels (comprehensive, endo/exo, channels, purified) |
| `02_FCI_Effects_Analysis.R` | Granger causality, impulse responses, out-of-sample forecasting |
| `03_FCI_Stationarity_Tests.R` | ADF, PP, KPSS unit root tests |
| `05_FCI_Local_Projections.R` | Local Projections for credit channel and macro effects |
| `06_FCI_Rolling_Stability.R` | Rolling/expanding window PCA loading stability |
| `07_FCI_Robustness_Endogeneity.R` | Endogeneity robustness (FCI_exNPL, FCI_exCredit) |
| `08_FCI_Regime_TVP_Analysis.R` | Regime-switching (IT/COVID), time-varying parameters |
| `09_FCI_Monetary_Policy_Interaction.R` | FCI x TPM interaction, systematic forecasting |
| `10_FCI_Growth_at_Risk.R` | Growth-at-Risk quantile regressions (Adrian et al. 2019) |
| `11_FCI_Block_SVAR.R` | Block-exogenous SVAR for external spillovers |
| `12_FCI_Output_Puzzle_Investigation.R` | Output puzzle: IMAEP_SANB LP, VECM, transmission channels |
| `14_FCI_Output_Puzzle_Sectoral.R` | Sectoral credit/output decomposition |
| `15_FCI_PostIT_Subsample_Analysis.R` | Post-IT subsample re-estimation (May 2011+) |
| `16_FCI_TCN_Reclassification.R` | Exchange rate reclassification sensitivity |
| `17_FCI_New_External_Data.R` | GFC index construction, regional proxy, disaggregated commodities |
| `18_FCI_Improved_IV_LP.R` | Expanded IV-LP with Anderson-Rubin weak-IV-robust inference |
| `19_FCI_Proxy_SVAR.R` | Proxy-SVAR identification (Mertens & Ravn 2013) |
| `20_FCI_Enriched_Block_SVAR.R` | 7-variable external block SVAR, historical decomposition |
| `21_FCI_Commodity_Puzzle.R` | Disaggregated commodity effects, ToT interaction LP |
| `22_FCI_Regional_Spillovers.R` | Brazil Selic spillover analysis |
| `23_FCI_Identification_Triangulation.R` | DXY-only IV-LP, GFC-PCA IV, globally-purged FCI, cross-method comparison |
| `27_FCI_Published_Quarterly_LP.R` | Quarterly sectoral LP using published national accounts |
| `28_FCI_Published_Monthly_LP.R` | Monthly sectoral LP using published activity indices |

## FCI Construction

The FCI combines 12 sign-adjusted variables (positive = tighter conditions) using four methods, then averages the normalized outputs:

| Variable | Channel | Sign |
|----------|---------|------|
| TPM (policy rate) | Rates | + |
| Lending-deposit spread | Rates | + |
| Market-policy spread | Rates | + |
| Credit growth (YoY) | Banking | - |
| Credit-to-deposit ratio | Banking | - |
| NPL ratio | Banking | + |
| ROE | Banking | - |
| Liquidity ratio | Banking | - |
| Nominal exchange rate (PYG/USD) | External/Domestic | + |
| Commodities (soybean + beef avg) | External | - |
| US Fed Funds rate | External | + |
| VIX | External | + |

All variables are standardized using a rolling 60-month window.

## Key Dates

| Event | Date |
|-------|------|
| Inflation Targeting adoption | May 2011 |
| COVID-19 period | March 2020 -- December 2021 |

## Citation

If you use this code or data, please cite:

> Gonzalez Masulli, J.M. (2026). Financial Conditions and Credit Transmission in a Dual-Engine Economy: Evidence from Paraguay. Working Paper.

## License

The views expressed are those of the author and do not represent the official position of the Banco Central del Paraguay.
