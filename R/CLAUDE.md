# CLAUDE.md

Instructions for Claude Code when working with this repository.

## Project Summary

R-based Financial Conditions Index (FCI) for Paraguay. Calculates composite index using 4 methods (Z-Score, PCA, VAR, DFM) at multiple aggregation levels. Developed by the Banco Central del Paraguay for financial stability monitoring and monetary policy analysis.

**Single FCI**: 12 variables, full sample 1996-2025

**Sample**: January 1996 - December 2025 (360 monthly observations)

**Pipeline execution**: ~8-10 minutes for complete analysis (~298 output files: ~171 PNG, ~127 CSV)

## How to Run

```bash
cd R/
Rscript RUN_ALL.R             # Complete pipeline (all scripts)
# OR
Rscript 01_FCI_Complete.R     # Main FCI calculation only
```

## Project Structure

```
FCI-Paraguay/
├── R/                  # Analysis scripts + CLAUDE.md + RUN_ALL.R
├── data/               # Input data
│   └── FCI_data_1.xlsx         # Main database (all sheets: Main_variables, Datos_macro, Global_Financial_Conditions, Main_Commodities_Prices, Quarterly_SA, Monthly_SA)
├── output/             # Generated outputs (~171 PNG + ~127 CSV)
│   ├── png/            # PNG charts
│   ├── csv/            # CSV data files
│   ├── pdf/            # Compiled paper PDF
│   └── reports/        # Paper source (LaTeX) + supporting reports
└── docs/               # EMR template, SVG diagrams, verification log
```

## Scripts

| File | Purpose |
|------|---------|
| `01_FCI_Complete.R` | **Main script** - FCI at all levels (comprehensive, endo/exo, channels, purified) |
| `02_FCI_Effects_Analysis.R` | Predictive power, Granger causality, IRF, out-of-sample forecasting, output puzzle |
| `03_FCI_Stationarity_Tests.R` | Unit root tests (ADF, PP, KPSS) for all variables and indices |
| `05_FCI_Local_Projections.R` | Local Projections (Jordà 2005) for macro effects and credit channel |
| `06_FCI_Rolling_Stability.R` | Rolling/expanding window stability analysis for PCA loadings |
| `07_FCI_Robustness_Endogeneity.R` | Endogeneity tests with FCI_exNPL and FCI_exCredit |
| `08_FCI_Regime_TVP_Analysis.R` | Regime-switching (IT/COVID), TVP analysis, sign reversal investigation |
| `09_FCI_Monetary_Policy_Interaction.R` | FCI × TPM interaction analysis and systematic out-of-sample forecasting |
| `10_FCI_Growth_at_Risk.R` | Growth-at-Risk (Adrian et al. 2019) quantile regression analysis |
| `11_FCI_Block_SVAR.R` | Block-exogenous SVAR for external spillover analysis |
| `12_FCI_Output_Puzzle_Investigation.R` | Output puzzle deep dive: IMAEP_SANB LP, VECM cointegration, transmission channels |
| `14_FCI_Output_Puzzle_Sectoral.R` | Output puzzle sectoral extension: sectoral credit/output analysis |
| `15_FCI_PostIT_Subsample_Analysis.R` | Post-IT subsample (May 2011+): FCI re-estimation, credit/output LP, Granger, asymmetric effects, full-vs-post-IT comparison |
| `16_FCI_TCN_Reclassification.R` | TCN reclassification sensitivity: reclassifies exchange rate from domestic to external, variance decomposition, LP re-estimation, diagnostic tests |
| `17_FCI_New_External_Data.R` | New external data construction: GFC index (simple + PCA), regional FC proxy (Brazil Selic), disaggregated commodity YoY changes, correlation analysis |
| `18_FCI_Improved_IV_LP.R` | Improved IV-LP: expanded instruments (DXY, US10Y, SP500, Selic, ToT), first-stage diagnostics, Anderson-Rubin weak-IV robust inference |
| `19_FCI_Proxy_SVAR.R` | Proxy-SVAR (Mertens & Ravn 2013): structural FCI shock identification using external proxy innovations, wild bootstrap CIs, multi-proxy comparison |
| `20_FCI_Enriched_Block_SVAR.R` | Enriched Block-SVAR: 7-var external block (adds DXY, US10Y, SP500, Selic, ToT), historical decomposition of FCI by shock |
| `21_FCI_Commodity_Puzzle.R` | Commodity puzzle investigation: disaggregated commodity regressions, ToT vs aggregate, offsetting effects (agricultural vs oil), nonlinear effects, commodity-credit interaction LP |
| `22_FCI_Regional_Spillovers.R` | Regional spillover analysis: Brazil Selic → Paraguay FCI contagion vs common shocks, asymmetric spillovers, spillover LP, bivariate vs controlled VAR |
| `23_FCI_Identification_Triangulation.R` | Identification triangulation: DXY-only IV-LP (single strong instrument), GFC-PCA composite IV, globally-purged FCI, 4-panel cross-method comparison figure |
| `RUN_ALL.R` | Master script - runs complete pipeline |

## Data

- **Input**: `../data/FCI_data_1.xlsx` (single file, 6 sheets)
- **Output**: `../output/` (PNG charts and CSV files)

### Data Sheets in FCI_data_1.xlsx

| Sheet | Content | Rows | Period |
|-------|---------|------|--------|
| `Main_variables` | FCI input variables (12 vars + credit/deposit stocks) | 360 | Jan 1996 – Dec 2025 |
| `Datos_macro` | Macro controls (IMAEP, IMAEP_SANB, IPC, Creditos, Creditos_deflactados) | 360 | Jan 1996 – Dec 2025 |
| `Global_Financial_Conditions` | US_10Y, S&P 500, DXY, Selic_rate, IPE, Paraguay_ToT | 336 | Jan 1998 – Dec 2025 |
| `Main_Commodities_Prices` | Soybean, Soybean_oil, Soybean_flour, Cotton, Corn, Wheat, Meat, Sugar, Oil_Brent | 360 | Jan 1996 – Dec 2025 |
| `Quarterly_SA` | National accounts: PIB, Consumo Privado, FBCF, sectoral GDP | 127 | 1994Q1 – 2025Q3 |
| `Monthly_SA` | Sectoral activity: IMAEP, IMAEP sin agri, sectoral indices | 145 | Jan 2014 – Jan 2026 |

**Note**: FBKf (investment) and Consumo (consumption) are available at quarterly frequency in Quarterly_SA. Scripts 12 and 14 use quarterly data for transmission channel analysis.

## Credit Variable

The database contains **stock variables** for credit. The script automatically calculates **YoY growth**:

```r
Crecimiento_creditos = (Creditos_Sector_privado_totales / lag(Creditos_Sector_privado_totales, 12) - 1) * 100
```

## Sign Conventions

**Positive FCI = tighter conditions**

| Variable | Sign | Rationale |
|----------|------|-----------|
| TPM (policy rate) | +1 | Higher rates = tighter |
| Spread_activas_pasivas | +1 | Higher spread = tighter |
| Spread_mercado_TPM | +1 | Higher market spread = tighter |
| Crecimiento_creditos | -1 | Higher growth = looser |
| Ratio_Cred_Depo | -1 | Higher intermediation = looser |
| Morosidad (NPL) | +1 | Higher NPL = tighter |
| Rentabilidad (ROE) | -1 | Higher profitability = looser |
| Liquidez | -1 | Higher liquidity = looser |
| TCN (exchange rate) | +1 | Depreciation = tighter (net importer) |
| Commodities (soybean + beef avg) | -1 | Higher prices = looser (commodity exporter) |
| FFER (US Fed rate) | +1 | Higher US rates = tighter |
| VIX | +1 | Higher volatility = tighter |

## Output Structure

### Level 1: Comprehensive FCI
- `FCI_COMP_AVG` - 12 variables, full sample 1996-2025

### Level 2: Endogenous vs Exogenous
- `FCI_ENDO_AVG` - Domestic conditions (9 variables, policy-influenced)
- `FCI_EXO_AVG` - External conditions (3 variables: FFER, VIX, Commodities)
- `FCI_ENDO_ORTHO` - Domestic net of external spillovers

### Level 3: By Economic Channel
- `FCI_RATES_AVG` - Interest rate channel (3 variables)
- `FCI_BANKING_AVG` - Credit/banking channel (5 variables)
- `FCI_EXTERNAL_AVG` - External/global channel (4 variables)

### Robustness Indices (from 01_FCI_Complete.R, also computed in 07_FCI_Robustness_Endogeneity.R)
- `FCI_exNPL_AVG` - Excludes Morosidad (11 vars, for NPL local projections)
- `FCI_exCredit_AVG` - Excludes Crecimiento_creditos (11 vars, for credit-LHS regressions)
- `FCI_ENDO_exCredit_AVG` - Excludes Crecimiento_creditos from ENDO (8 vars, for credit-LHS decomposition)

### Level 4: Purified FCI (from 01_FCI_Complete.R)
- `FCI_exTPM_AVG` - Excludes TPM (simple exclusion)
- `FCI_PURIFIED_AVG` - Variables orthogonalized to TPM + lags (residual-based, ρ=0.359 with baseline)

### Regime Indicators (from 08_FCI_Regime_TVP_Analysis.R)
- `regime_IT` - 0 before May 2011, 1 from May 2011+ (Inflation Targeting)
- `regime_COVID` - 1 during Mar 2020 - Dec 2021, 0 otherwise
- `regime_IT_exCOVID` - IT period excluding COVID months

## Variable Groups

```r
# Level 3: By Economic Channel
rates         = c("TPM", "Spread_activas_pasivas", "Spread_mercado_TPM")           # 3 vars
banking       = c("Crecimiento_creditos", "Ratio_Cred_Depo", "Morosidad",
                  "Rentabilidad", "Liquidez")                                       # 5 vars
external = c("TCN", "Commodities", "FFER", "VIX")                                  # 4 vars

# Level 2: Endogenous vs Exogenous
endogenous = rates + banking + TCN                                                 # 9 vars
exogenous  = c("FFER", "VIX", "Commodities")                                       # 3 vars

# Comprehensive
all_vars = rates + banking + external                                              # 12 vars
```

## Methodologies

Each FCI is calculated using four methods, then averaged:

1. **Z-Score**: Equal-weighted average of standardized variables
2. **PCA**: First principal component (variance-weighted)
3. **VAR**: Average of VAR fitted values across equations
4. **DFM**: Dynamic factor model using MARSS state-space estimation

Rolling 60-month window used for standardization (with partial window for initial observations).

## Required Packages

```r
# Core packages
readxl, dplyr, tidyr, zoo, FactoMineR, vars, MARSS, ggplot2, gridExtra, lmtest, sandwich

# Stationarity tests (03_FCI_Stationarity_Tests.R)
tseries, urca

# Growth-at-Risk (10_FCI_Growth_at_Risk.R) and IV-LP (18_FCI_Improved_IV_LP.R)
quantreg, ivreg
```

## Output Files

Output files are organized into subfolders:
- `output/png/` — PNG charts (~171 files)
- `output/csv/` — CSV data files (~127 files)
- `output/pdf/` — Compiled paper PDF (1 file)
- `output/reports/` — Paper source (LaTeX + bib + CAS template) + supporting reports (15 files)

`RUN_ALL.R` automatically moves new PNG/CSV/PDF files into the appropriate subfolders after each pipeline run.

### CSV Files (output/csv/)
| File | Description |
|------|-------------|
| `FCI_Complete_Results.csv` | All FCI indices |
| `FCI_Summary_Statistics.csv` | Summary statistics for AVG indices |
| `FCI_Correlation_Matrix.csv` | Correlation matrix of AVG indices |
| `FCI_Internal_Granger_Causality.csv` | FCI_EXO → FCI_ENDO causality tests |
| `FCI_Extended_Granger_All.csv` | All FCI × macro variable Granger tests |
| `LP_Macro_Standard.csv` | Local projection results for macro variables |
| `LP_Credit_Standard.csv` | Local projection results by credit type |
| `LP_NPL_Standard.csv` | Local projection results for NPL |
| `LP_NPL_exNPL.csv` | NPL LP results using FCI_exNPL |
| `LP_IV_Results.csv` | IV-LP (2SLS) and OLS coefficients for credit regressions |
| `LP_IV_Diagnostics.csv` | IV diagnostics: first-stage F, Sargan, Wu-Hausman by horizon |
| `LP_Placebo_Test.csv` | Placebo (backward LP) and forward LP coefficients |
| `LP_Purged_FCI.csv` | LP results using macro-orthogonalized FCI (conservative lower bound) |
| `VAR_Ordering_Robustness.csv` | Credit IRF under 4 alternative Cholesky orderings |
| `LP_Credit_Endogeneity_Test.csv` | Comparison: Full FCI vs FCI_exCredit effects |
| `LP_Macro_Asymmetric.csv` | Asymmetric LP results (tightening vs easing) |
| `LP_Credit_Asymmetric.csv` | Asymmetric credit LP results |
| `LP_NPL_Asymmetric.csv` | Asymmetric NPL LP results |
| `FCI_Versions_Correlation.csv` | Correlation between FCI versions |
| `FCI_Credit_Currency_Correlations.csv` | FCI × credit by currency correlations |
| `FCI_Weight_Stability_Summary.csv` | PCA loading stability metrics |
| `FCI_Rolling_Loadings.csv` | Time-varying PCA loadings (rolling window) |
| `FCI_Expanding_Loadings.csv` | Time-varying PCA loadings (expanding window) |
| `FCI_Rolling_Correlations.csv` | Rolling correlations between FCI versions |
| `FCI_Robustness_Versions.csv` | FCI_exNPL and FCI_exCredit series |
| `FCI_Credit_Weight_Decomposition.csv` | Credit growth weight in FCI by methodology (Z-Score, PCA, VAR, DFM) |
| `FCI_Regime_Indicators.csv` | Regime dummy variables (IT, COVID) |
| `FCI_Regime_PCA_Comparison.csv` | Pre-IT vs Post-IT loading comparison |
| `LP_Regime_IT_Results.csv` | LP results with IT regime interaction |
| `LP_Regime_COVID_Results.csv` | LP results with COVID regime interaction |
| `FCI_TVP_Loadings.csv` | Time-varying PCA loadings (exponential weighting) |
| `FCI_TVP_vs_Rolling_Comparison.csv` | Comparison of TVP vs rolling loadings |
| `FCI_TVP_Sign_Consistency.csv` | Sign consistency analysis over time |
| `FCI_Threshold_Indicator.csv` | Crisis threshold indicator (FCI > 1 SD) |
| `LP_Threshold_Effects.csv` | Threshold LP results (crisis vs normal) |
| `FCI_TPM_Interaction_Amplification.csv` | dTPM × FCI interaction results (amplification) |
| `FCI_TPM_Interaction_Stance.csv` | FCI × TPM stance interaction results |
| `Forecast_RMSE_Comparison.csv` | RMSE by model/horizon/window/target |
| `Forecast_DM_Tests.csv` | Diebold-Mariano test results |
| `Forecast_Model_Rankings.csv` | Model rankings by target/horizon |
| `FCI_Sample_Verification.csv` | Execution date, sample range, n_observations |
| `FCI_Stationarity_Tests.csv` | ADF, PP, KPSS test results |
| `FCI_Sign_Flip_Diagnostic.csv` | Sign reversal analysis across regimes |
| `FCI_Inconsistent_Loadings.csv` | Variables with sign inconsistencies |
| `FCI_Endogeneity_Summary.csv` | Summary of endogeneity robustness |
| `FCI_Variance_Decomposition.csv` | Domestic vs external variance decomposition |
| `LP_Threshold_Sensitivity.csv` | Threshold effects at 1.0, 1.5, 2.0 SD |
| `FCI_Window_Sensitivity.csv` | Rolling window size comparison |
| `Growth_at_Risk_Results.csv` | Quantile regression coefficients |
| `GaR_Asymmetry_Tests.csv` | Tests for tail asymmetry |
| `Block_SVAR_IRF.csv` | Structural IRF results |
| `Block_SVAR_Spillover_Summary.csv` | External spillover analysis |
| `Block_SVAR_FEVD.csv` | Block-SVAR variance decomposition |
| `Output_Puzzle_Diagnostic.csv` | Output puzzle investigation |
| `Output_Puzzle_LP_Results.csv` | LP results for all output measures (IMAEP, IMAEP_SANB, FBKf, Consumo) |
| `Output_Puzzle_Summary.csv` | Summary comparison at key horizons |
| `Output_Puzzle_VECM_Summary.csv` | Johansen cointegration test results |
| `Credit_Output_LP_Results.csv` | LP: Credit growth → Output growth at all horizons |
| `Credit_Output_Granger.csv` | Granger causality tests: Credit ↔ Output |
| `Credit_Output_Mediation.csv` | Mediation analysis: Does credit mediate FCI → Output? |
| `LP_vs_SVAR_Comparison.csv` | Comparison of LP and SVAR impulse responses |
| `SVAR_IRF_Results.csv` | Structural VAR impulse response functions |
| `SVAR_FEVD_Results.csv` | SVAR forecast error variance decomposition |
| `PostIT_FCI_Diagnostics.csv` | Post-IT PCA loadings, variance explained, cross-method corr, var decomp |
| `PostIT_Granger.csv` | Post-IT Granger causality results |
| `PostIT_LP_Credit.csv` | Post-IT credit LP (all horizons, FCI types, credit types) |
| `PostIT_LP_Output.csv` | Post-IT output LP (all horizons, all variables) |
| `PostIT_LP_Asymmetric.csv` | Post-IT asymmetric LP with Wald tests |
| `PostIT_Comparison_Summary.csv` | Full-sample vs post-IT comparison with z-tests |
| `TCN_Reclass_FCI_Correlations.csv` | Correlations between baseline and alternative (TCN reclassified) FCI variants |
| `TCN_Reclass_Variance_Decomposition.csv` | Variance decomposition under baseline vs alternative TCN classification (4 rows) |
| `TCN_Reclass_LP_Results.csv` | Full LP results for all TCN reclassification specifications (all horizons) |
| `TCN_Reclass_LP_Comparison.csv` | LP comparison at h=6,12,18 with z-tests for equality |
| `TCN_Reclass_Diagnostics.csv` | TCN driver analysis (partial R-squared, Granger, PCA loadings) |
| `New_External_Variables.csv` | All constructed indices + commodity YoY changes (shared by scripts 18-22) |
| `GFC_Index_Diagnostics.csv` | GFC PCA loadings, variance explained, method correlations |
| `New_Variables_Correlation_Matrix.csv` | Full correlation matrix: new variables × FCI variants |
| `Commodity_Disaggregated_Correlations.csv` | Individual commodity correlations with FCI and credit |
| `New_External_Sample_Info.csv` | Sample ranges, missing data, overlap summary |
| `IV_First_Stage_Battery.csv` | All instrument diagnostics (t-stats, partial R², F-stats) |
| `IV_LP_Improved_Results.csv` | IV-LP coefficients at all horizons for all credit types |
| `IV_LP_Improved_Diagnostics.csv` | First-stage F, Hansen J, Wu-Hausman by horizon |
| `IV_Weak_Robust_Inference.csv` | Anderson-Rubin test results, AR confidence intervals |
| `ProxySVAR_IRF_Results.csv` | IRFs for all proxy×response×horizon combinations |
| `ProxySVAR_Diagnostics.csv` | Relevance F, exogeneity tests per proxy |
| `ProxySVAR_Credit_Comparison.csv` | Credit IRF at key horizons across all proxies + LP baseline |
| `ProxySVAR_Bootstrap_CI.csv` | Bootstrap confidence intervals for Proxy-SVAR |
| `Enriched_BSVAR_IRF.csv` | All IRFs (shock×response×horizon) for 7-var external block |
| `Enriched_BSVAR_FEVD.csv` | FEVD at all horizons by individual shock |
| `Enriched_BSVAR_Historical_Decomposition.csv` | Historical shock contributions to FCI by date |
| `Enriched_BSVAR_Comparison.csv` | 3-var vs 7-var FEVD comparison table |
| `Commodity_Disaggregated_Regressions.csv` | All individual + joint commodity regression results |
| `Commodity_ToT_Model_Comparison.csv` | R², AIC, coefficients for aggregate vs ToT vs disaggregated |
| `Commodity_Credit_Interaction_LP.csv` | Interaction LP coefficients (FCI × ToT) at all horizons |
| `Marginal_Effects_FCI_ToT.csv` | Delta-method marginal effects of FCI on credit by ToT regime (h=6,12,18) |
| `Commodity_Granger_Results.csv` | All commodity Granger causality results |
| `Commodity_Puzzle_Summary.csv` | Summary table of key commodity puzzle findings |
| `Regional_Spillover_Regressions.csv` | All spillover regression results (Selic + global controls) |
| `Regional_Granger_Results.csv` | Granger causality tests (bivariate + conditional) |
| `Regional_Spillover_LP.csv` | LP coefficients: ΔSelic → FCI and → Credit at all horizons |
| `Regional_Spillover_VAR_IRF.csv` | VAR IRFs: bivariate and controlled |
| `Triangulation_DXY_IV_LP.csv` | DXY-only 2SLS and OLS coefficients at all horizons with AR CIs |
| `Triangulation_GFC_PCA_IV.csv` | GFC-PCA 2SLS coefficients, AR CIs, first-stage F by horizon |
| `Triangulation_Globally_Purged_LP.csv` | LP results for baseline, macro-purged, and globally-purged FCI |
| `Triangulation_Summary_Comparison.csv` | Cross-method comparison at h=6,12,18,24 with all diagnostics |

### PNG Visualizations (output/png/)
| Series | Content |
|--------|---------|
| `01-11_*.png` | Core FCI visualizations (methods, levels, contributions) |
| `23-29_*.png` | Effects analysis (Granger, IRF, forecasting, FEVD, VAR ordering robustness) |
| `30_*.png` | Stationarity test summary |
| `50-59c_*.png` | Local projections (macro, credit, NPL, asymmetric, IV, placebo, purged FCI) |
| `60-64_*.png` | Rolling/expanding stability analysis, window comparison |
| `70-72_*.png` | Robustness/endogeneity comparisons |
| `80-89_*.png` | Regime-switching, TVP analysis, sign consistency |
| `90-99_*.png` | Monetary policy interaction and systematic forecasting |
| `100-102_*.png` | Growth-at-Risk fan charts and coefficient comparison |
| `110-112_*.png` | Block-SVAR spillovers and FEVD |
| `120-126_*.png` | Output puzzle investigation (LP vs SVAR, credit channel, mediation) |
| `140-165_*.png` | Sectoral credit/output LP analysis (long/short horizons, heatmaps, asymmetric, endo/exo) |
| `170-178_*.png` | Post-IT subsample analysis (FCI construction, credit/output LP, asymmetric, full-vs-post-IT comparison, PCA loadings, dashboard) |
| `179-185_*.png` | TCN reclassification sensitivity (FCI comparison, variance decomposition, LP overlays ENDO/EXO, diagnostics, sensitivity table, dashboard) |
| `186-191_*.png` | New external data (GFC index, regional FC proxy, correlation heatmap, disaggregated commodities, dashboard) |
| `196-200_*.png` | Improved IV-LP (first-stage diagnostics, credit IRFs, OLS vs 2SLS, Anderson-Rubin test, dashboard) |
| `206-212_*.png` | Proxy-SVAR (individual proxy IRFs, credit comparison across proxies, diagnostics, dashboard) |
| `216-222_*.png` | Enriched Block-SVAR (7-var external IRFs to FCI/Credit, FEVD, FEVD comparison, historical decomposition, episodes, dashboard) |
| `226-234_*.png` | Commodity puzzle (disaggregated coefficients, ToT comparison, offsetting effects, lagged/nonlinear, interaction LP, Granger heatmap, dashboard, marginal effects) |
| `241-248_*.png` | Regional spillovers (direct spillover, contagion vs common shocks, Granger, asymmetric, LP FCI/Credit, VAR IRF, dashboard) |
| `251-255_*.png` | Identification triangulation (4-panel cross-method figure, DXY-only IV-LP, GFC-PCA diagnostics, globally-purged LP, dashboard) |

### Reports (output/reports/)
Paper source and supporting documentation (15 files: LaTeX source, bib, CAS template files, 8 supporting markdown reports).

## Key Analytical Features

1. **Orthogonalization**: Regresses FCI_ENDO on FCI_EXO to extract domestic conditions net of external spillovers
2. **Asymmetric effects**: Local projections allow different responses to tightening vs easing
3. **Credit channel decomposition**: Tests effects on Total, MN (local currency), and USD (dollarized) credit
4. **Endogeneity robustness**: All credit-LHS regressions use FCI_exCredit (excludes credit growth); NPL-LHS uses FCI_exNPL. Both computed in script 01 and available to all downstream scripts. Weight decomposition quantifies credit's contribution by methodology.
5. **Parameter stability**: Rolling window analysis of PCA loadings over time
6. **Regime-dependent analysis**: Tests if FCI effects differ pre/post Inflation Targeting (May 2011) and during COVID
7. **Time-varying parameters**: Exponentially-weighted PCA for smooth time-varying loadings (λ=0.97)
8. **Threshold effects**: Tests if FCI effects are stronger in crisis periods (multiple thresholds: 1.0, 1.5, 2.0 SD)
9. **Monetary policy interaction**: Tests if FCI amplifies/dampens TPM transmission and if TPM stance affects FCI transmission
10. **Systematic forecasting**: Out-of-sample RMSE comparison across 7 models with Diebold-Mariano tests
11. **Stationarity verification**: ADF, PP, and KPSS tests on all variables and indices
12. **Sign reversal investigation**: Analyzes PCA loading sign changes across regimes with economic interpretation
13. **Purified FCI**: FCI orthogonalized to monetary policy (FCI_exTPM and FCI_PURIFIED)
14. **Growth-at-Risk**: Quantile regressions testing if FCI affects left tail more than median (Adrian et al. 2019)
15. **Block-exogenous SVAR**: Structural analysis of external spillovers with small-open-economy identification
16. **Output puzzle analysis**: Investigation of why FCI predicts credit but not output
17. **Agricultural contamination test**: LP comparison of IMAEP (aggregate) vs IMAEP_SANB (non-agricultural) to test if agriculture masks credit channel
18. **VECM cointegration**: Johansen test for long-run equilibrium between FCI, credit, and output
19. **Transmission channel decomposition**: Separate LP for investment (FBKf) and consumption to identify FCI transmission mechanisms
20. **Credit channel mediation**: Tests if credit mediates FCI → Output relationship (Baron & Kenny approach)
21. **LP vs SVAR robustness**: Compares Local Projections with Structural VAR (Cholesky) for credit-output transmission
22. **Sequential timing analysis**: Identifies optimal lags in FCI → Credit → Output chain
23. **Real vs nominal credit comparison**: Replicates LP, Granger, VAR, SVAR, and mediation using deflated credit growth
24. **Controlled LP**: LP with IMAEP and IPC controls to isolate direct financial channel (vs unconditional LP)
25. **Full-system VAR comparison**: Tests whether swapping nominal for real credit changes IRFs to output and prices
26. **IV-LP identification**: Uses external variables (VIX, FFER, Commodities) as instruments for domestic FCI under small-open-economy exclusion restriction; 2SLS with Sargan overidentification test
27. **Placebo test**: Backward LP (regressing past credit on current FCI) as falsification — confirms forward-looking content
28. **Purged FCI**: Macro-orthogonalized FCI (residuals after removing IMAEP, IPC + 2 lags) provides conservative lower-bound credit effect
29. **VAR ordering robustness**: Tests credit IRF stability across 4 Cholesky orderings (FCI first/last, output first, credit first)
30. **Post-IT subsample analysis**: Re-estimates FCI and core LP specifications on post-IT data (May 2011+, ~175 obs) to confirm results hold under current regime with theory-consistent PCA loadings
31. **TCN reclassification sensitivity**: Reclassifies exchange rate from domestic (9-var ENDO) to external (4-var EXO with TCN), re-estimates variance decomposition and credit LP on both full and post-IT samples, with z-tests for equality and diagnostic tests (TCN drivers, Granger, loading stability)
32. **Global Financial Conditions index**: Constructs GFC index from DXY, US10Y, SP500, VIX, FFER using simple average and PCA methods; regional FC proxy from Brazil Selic
33. **Improved IV-LP**: Expanded instrument set (Selic, DXY, SP500, US10Y, ToT) with first-stage diagnostics, Anderson-Rubin weak-IV robust test, and comparison to original IV-LP
34. **Proxy-SVAR**: Mertens & Ravn (2013) external instruments SVAR identifying structural FCI shock using proxy innovations (Selic, DXY, VIX, GFC); wild bootstrap CIs; multi-proxy robustness
35. **Enriched Block-SVAR**: 7-variable external block (VIX, SP500, DXY, US10Y, FFER, Selic, ToT) with historical decomposition of FCI into shock contributions across crisis episodes
36. **Commodity puzzle diagnosis**: Disaggregated commodity price effects, Paraguay-specific Terms of Trade, agricultural export vs oil import offsetting effects, nonlinear/threshold tests, commodity-credit interaction LP
37. **Regional spillover analysis**: Brazil Selic → Paraguay FCI contagion tests controlling for global shocks, asymmetric spillovers, spillover LP (direct + indirect credit channel), bivariate vs controlled VAR
38. **Identification triangulation**: DXY-only IV-LP (single strong instrument, F≈33), GFC-PCA composite IV (Barnichon & Mesters 2025), globally-purged FCI (removes all external variation), 4-panel cross-method comparison figure

## Key Dates for Paraguay

| Event | Date | Notes |
|-------|------|-------|
| Inflation Targeting | May 2011 | Monetary policy regime change |
| COVID-19 Start | March 2020 | Pandemic onset |
| COVID-19 End | December 2021 | End of acute pandemic period |

## Key Empirical Findings

Results from the comprehensive analysis (March 2026, 360-observation sample Jan 1996–Dec 2025):

### Credit Transmission
- **1 SD FCI tightening → -7.66 pp real credit growth** at 12 months (p=0.006), full-sample FCI
- Full-sample domestic (FCI_ENDO, h=12): -7.95 pp, p=0.005
- Post-IT domestic (FCI_ENDO_exCredit, h=12): -5.60 pp, p=0.004
- Post-IT USD credit: -12.84 pp (p=0.0001); MN credit: -1.69 pp (p=0.285)
- Credit weight in FCI: Z-Score 8.3%, PCA 14.6%, VAR ~8.3%, DFM 13.3% (avg ~11.1%)
- FCI_exCredit correlation with full FCI: 0.975

### Variance Decomposition
- **94.8% domestic factors**, 5.2% external factors
- Confirms FCI primarily reflects domestic financial conditions

### External Spillovers (Block-SVAR)
| Shock | → FCI_ENDO | Significant? |
|-------|------------|--------------|
| VIX | +0.079 (h=10) | **Yes** |
| FFER | -0.038 (h=24) | No |
| Commodities | +0.005 (h=1) | No |

**Key insight**: Global risk aversion (VIX) significantly tightens domestic conditions; US monetary policy does not.

### Growth-at-Risk
- **No significant asymmetric effects** detected
- FCI affects growth distribution symmetrically
- Linear models appropriate for Paraguay

### Structural Breaks
- **3 sign flips** between pre-IT and post-IT periods (Spread, Commodities, TCN)
- Explained by inflation targeting regime change (May 2011)
- Post-IT loadings more consistent with economic theory

### Output Puzzle
- FCI strongly predicts credit but **not output (IMAEP)**
- Explanation: ~45% of GDP structurally insensitive (crops ~7%, manufacturing ~20% — partly due to agro-processing component, binational electricity ~9%, taxes ~8%; all at constant 2014 prices); informal credit channels
- Contribution: "Credit channel without output effects" in commodity-exporting SOEs

### Stationarity
- 2/12 variables stationary in levels (TPM, VIX)
- FCI indices show mixed stationarity (appropriate for financial conditions)

## Report Files

| File | Description |
|------|-------------|
| `output/reports/FCI_Paraguay_EMR_Submission.tex` | **PRIMARY: journal submission LaTeX source** |
| `output/pdf/FCI_Paraguay_EMR_Submission.pdf` | Compiled submission PDF |
| `output/reports/FCI_Paraguay_EMR_Submission.md` | Markdown mirror of final paper |
| `output/reports/FCI_Paraguay_Online_Appendix.md` | Online appendix (referenced by paper) |
| `output/reports/FCI_Paraguay_Technical_Report.md` | Technical robustness documentation |
| `output/reports/fci-paraguay.bib` | Bibliography database |
| `docs/Paper_Corrections_Verification.md` | Verification log |

## References

- Brave & Butters (2011) - Chicago Fed NFCI
- Hatzius et al. (2010) - Goldman Sachs FCI
- Arrigoni, Bobasu & Venditti (2022) - FCIs for Latin America
- Jordà (2005) - Local Projections methodology
- Diebold & Mariano (1995) - Comparing Predictive Accuracy
- Adrian, Boyarchenko & Giannone (2019) - Vulnerable Growth / Growth-at-Risk
- Cushman & Zha (1997) - Block-exogenous SVAR
- Dickey & Fuller (1979) - ADF test
- Phillips & Perron (1988) - PP test
- Kwiatkowski et al. (1992) - KPSS test
- Johansen (1988, 1991) - Cointegration Tests and VECM
- Lütkepohl (2005) - New Introduction to Multiple Time Series Analysis

## Target Journals

With current implementation, suitable for:
- **Field journals**: Journal of International Money and Finance, Journal of Banking & Finance, Economic Modelling
- **Policy outlets**: BIS Working Papers, IMF Working Papers, CEMLA publications

## Last Updated

March 2026 - Journal submission cleanup: removed orphaned script 13 outputs, outdated reports/PDFs/docs, LaTeX build artifacts; created README.md for referee replicability; verified all paper numbers against CSV outputs.

March 2026 - Output subfolder split: `output/png/`, `output/csv/`, `output/pdf/` for type-based organization; `RUN_ALL.R` auto-sorts new files.

March 2026 - Data consolidation: EMBI removed, commodity variable changed to soybean+beef average levels, sample updated to Jan 1996–Dec 2025, all external data sheets consolidated into FCI_data_1.xlsx. Scripts 12/14 adapted for quarterly transmission analysis.

March 2026 - Project folder reorganization: `docs/` for working documents, `output/reports/` for markdown reports; script 14 added to pipeline.

January 2026 - Complete implementation of reviewer feedback including Growth-at-Risk, Block-SVAR, endogeneity robustness, and comprehensive empirical paper.
