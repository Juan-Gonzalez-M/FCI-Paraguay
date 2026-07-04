# Dollar Shocks and Domestic Credit Conditions in a Partially Dollarized Economy: Evidence from Paraguay

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
| `data/Boletin_Bancos_May2026.xlsx` | Bank-level data from the BCP's published *Boletín de Bancos* (public source): credit by bank x sector x currency, balance-sheet items, ratios | 125 monthly (Jan 2016--May 2026), 20 banks |

**FCI_data_1 sheets:** Main_variables (12 FCI inputs + credit/deposit stocks), Datos_macro (IMAEP, CPI, credit), Global_Financial_Conditions (US 10Y, S&P 500, DXY, Selic, ToT), Main_Commodities_Prices (9 commodities), Quarterly_SA (GDP and expenditure components), Monthly_SA (sectoral activity indices).

External series used by the revision-extras scripts (Fed broad dollar index, BRL, EUR, GBP, TED spread) are downloaded from FRED automatically and cached in `output/revision/external/`.

## Replication

### Requirements

- **R version:** 4.3+ (tested on 4.5.1)
- **Required packages:** readxl, dplyr, tidyr, zoo, FactoMineR, vars, MARSS, ggplot2, gridExtra, lmtest, sandwich, tseries, urca, quantreg, ivreg

Missing packages are installed automatically by the pipeline.

### Running the Full Pipeline

The analysis is organized in three sequential phases, each with its own master script:

```bash
cd R/
Rscript RUN_ALL.R          # Phase 1: aggregate FCI pipeline (~8-10 min)
Rscript RUN_MICRO.R        # Phase 2: bank-level micro analysis (~12 min)
Rscript 36_Aggregate_FXAdjusted_Credit_LP.R   # FX-adjusted aggregate LPs (<1 min)
Rscript RUN_REVISION.R     # Phase 3: identification & robustness extensions (~10 min)
```

Phase 1 produces ~144 PNG charts and ~127 CSV files in `output/png|csv/`; Phase 2 writes to `output/micro/`; Phase 3 writes to `output/revision/` (and requires Phases 1-2 to have run once). Additional packages for Phases 2-3: `data.table`, `fixest`, `modelsummary` (installed automatically if missing).

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
├── R/                          # Analysis scripts
│   ├── RUN_ALL.R               # Phase 1 master: aggregate FCI pipeline
│   ├── 01--28_*.R              # Phase 1: FCI construction and aggregate analysis
│   ├── RUN_MICRO.R             # Phase 2 master: bank-level analysis
│   ├── micro_helpers.R         # Phase 2 shared utilities
│   ├── 30--36_*.R              # Phase 2: micro panels, designs A-D, FX-adjusted LPs
│   ├── README_MICRO.md         # Phase 2 documentation
│   ├── RUN_REVISION.R          # Phase 3 master: identification & robustness extensions
│   ├── revision_helpers.R      # Phase 3 shared utilities (LP/IV engines, FRED fetch)
│   ├── 37--47_*.R              # Phase 3: CHR bounds, post-IT IV, batteries, audits
│   ├── README_REVISION.md      # Phase 3 documentation
│   └── CLAUDE.md               # Development reference (not needed for replication)
├── data/
│   ├── FCI_data_1.xlsx         # Aggregate input file (6 sheets)
│   └── Boletin_Bancos_May2026.xlsx  # Bank-level input (public BCP bulletin)
├── output/                     # Generated outputs
│   ├── csv/  png/              # Phase 1 tables and charts
│   ├── micro/                  # Phase 2 outputs (csv, png, rds caches)
│   └── revision/               # Phase 3 outputs (csv, png, external data cache)
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

### Phase 2: Bank-Level Analysis (`RUN_MICRO.R`; details in `R/README_MICRO.md`)

| Script | Purpose |
|--------|---------|
| `30_Bank_Micro_Data_Construction.R` | Panels from the Boletín de Bancos; valuation adjustment; shocks; merger protocol; reconciliation |
| `31_Bank_Micro_DesignA_Currency.R` | Within-bank x sector currency-denomination LP; VIX and permutation placebos |
| `32_Bank_Micro_DesignB_Exposure.R` | Bank lending channel via exposure heterogeneity (no instrument needed) |
| `33_Bank_Micro_DesignC_Hedging.R` | Hedging heterogeneity: triple interaction, FX-share gradient, split samples |
| `34_Bank_Micro_DesignD_Mechanisms.R` | Deposits and NPLs by currency |
| `35_Bank_Micro_Robustness.R` | COVID, alt shocks, mergers, thresholds, valuation variants, inference battery |
| `36_Aggregate_FXAdjusted_Credit_LP.R` | Aggregate credit LPs with constant-exchange-rate outcomes (validated replication of script 05) |

### Phase 3: Identification & Robustness Extensions (`RUN_REVISION.R`; details in `R/README_REVISION.md`)

| Script | Purpose |
|--------|---------|
| `37_Micro_Freeze_Checks.R` | Gradient robustness (ex-COVID, tightening-only); TCN-timing audit; reconciliation decomposition |
| `38_CHR_Plausibly_Exogenous.R` | Conley-Hansen-Rossi plausibly-exogenous bounds with trade-calibrated direct effect |
| `39_PostIT_ExpandingFCI_IV.R` | Expanding-window FCI; post-IT IV-LP with Anderson-Rubin sets and effective F |
| `40_LagAugmented_LP.R` | Lag-augmented LPs (Montiel Olea & Plagborg-Moller 2021) |
| `41_Falsification_Placebos.R` | Structurally insulated placebo outcomes (quarterly) |
| `42_FCI_Component_Exclusion.R` | Leave-one-out FCI battery; quantities-at-t-1 variant; rates-only co-baseline |
| `43_External_Data_Fetch.R` | FRED downloads (broad dollar, BRL, EUR, GBP, TED), cached |
| `44_Enhanced_Instrument_Robustness.R` | Broad-dollar IV; EUR/GBP dollar-free placebo; enhanced purge; funding-stress horse race |
| `45_COVID_Reclassification_Check.R` | COVID loan-reclassification sensitivity for the micro panel |
| `46_USD_ShiftShare_Decomposition.R` | Within/between decomposition of the aggregate USD-PYG credit differential |
| `47_FXAdj_PostIT_Asymmetric_Check.R` | FX-adjustment audit: post-IT subsample and asymmetric LPs |

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
