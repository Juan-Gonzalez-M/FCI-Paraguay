# Bank-Level Micro Phase

Bank-level extension of the FCI-Paraguay analysis: within-bank identification
of the dollar-shock credit channel using publicly available data from the
Central Bank of Paraguay's published Boletín de Bancos. **Fully additive**:
nothing in the aggregate pipeline (scripts 01–23, `RUN_ALL.R`,
`output/csv|png|pdf`) is touched. All outputs go to `output/micro/{csv,png,rds}`.

## How to run

```bash
cd R/
Rscript RUN_MICRO.R                    # full pipeline (999 wild-bootstrap reps)
MICRO_WCB_B=99 Rscript RUN_MICRO.R     # quick pass (fewer bootstrap reps)
```

Individual scripts can be run standalone in order (30 must run first).

## Scripts

| File | Plan phase | Content |
|---|---|---|
| `micro_helpers.R` | — | Paths, constants, valuation/growth utilities, native restricted wild-cluster-bootstrap-t (Rademacher, cluster = bank), inference battery (cluster / two-way / Driscoll–Kraay), plot helpers |
| `30_Bank_Micro_Data_Construction.R` | Phase 0 (§1) | Panels P1–P4 + deposits; shock series (purged DXY innovation = primary, raw DXY, VIX, FCI shock); valuation adjustment (§1.2); entry/exit & merger protocol (§1.3); aggregate reconciliation; exposure variables (§1.5) |
| `31_Bank_Micro_DesignA_Currency.R` | Phase 1 (§2) | Within-bank×sector×month currency-denomination LP (h=1–18), bank×sector×time + bank×sector×currency FE; VIX placebo; 500-draw random-timing permutation placebo |
| `32_Bank_Micro_DesignB_Exposure.R` | Phase 2 (§3) | Bank lending channel: shock × exposure (FC deposit share, external funding, Tier 1, size, liquidity), bank level and bank×sector level with sector×time FE (Degryse et al. 2019 demand absorption); 2016-fixed exposure robustness |
| `33_Bank_Micro_DesignC_Hedging.R` | Phase 3 (§4) | Triple interaction Shock × USD × Unhedged; sharp test, mid-dollarization primary contrast, continuous FX-share gradient, month-of-year×sector×currency FE seasonality check; split-sample hedged vs unhedged |
| `34_Bank_Micro_DesignD_Mechanisms.R` | §3.3 | Deposit flows by currency (bank×time FE); NPL response in pp by currency and hedging status (answers R1 #5) |
| `35_Bank_Micro_Robustness.R` | Phase 4 (§5) | COVID; alt shocks; merger treatments (no-drop / balanced-13 / pro-forma); cell thresholds 1/5/10 bn; deflation variants incl. deliberately-unadjusted valuation-bias demo; asymmetry; cell trends; winsorization; extensive margin LPM; forest plot |
| `RUN_MICRO.R` | — | Master runner (isolated sessions, run log) |

## Key data-construction decisions (documented deviations / notes)

0. **Input workbook:** `data/Boletin_Bancos_May2026.xlsx` (BCP Boletín de
   Bancos, public; monthly, Jan 2016 – May 2026, 20 entities, 9 sheets).
1. **TCN & CPI end Dec 2025** while the bank panel runs to May 2026, and the
   workbook's `TC` sheet is credit-card data, not exchange rates. FC-book
   valuation adjustment therefore stops at Dec 2025. There is a marked
   `TCN_2026` / `IPC_2026` slot at the top of script 30 — paste the five 2026
   monthly values there to unlock the extra outcome months.
2. **Size** = ln(risk-weighted assets) from the Ratios sheet ("ACTIVOS Y
   CONTINGENTES PONDERADOS"); total accounting assets are not a single line in
   the bulletin.
3. **Liquidity** = "Disponible + Inversiones Temporales/Depósitos" (Ratios)
   instead of hand-built EEFF aggregate.
4. **External funding share** = Externo / (Externo + total deposits) (total
   liabilities not a single line).
5. **US-CPI deflation variant** (plan §5.5) not available in the project
   database; the FC book is nominal USD (baseline) with a PYG-nominal variant.
   As the plan notes, this does not affect within-month currency-contrast
   identification.
6. **Design C continuous treatment** uses the 2016 sector FX share as the
   dollar-earning-intensity gradient (national-accounts export intensity is
   not available at the bulletin's 13-sector breakdown).
7. **Fed broad dollar index / FOMC-window instrument** not in the database;
   alternative shocks implemented: raw DXY innovation, purged DXY, VIX, FCI
   shock.
8. **`fwildclusterboot` is archived on CRAN**; the wild cluster bootstrap
   (restricted, Rademacher, bootstrap-t) is implemented natively in
   `micro_helpers.R` and applied with cluster = bank throughout.

## Verified facts (script 30 output, matching plan §0)

- 20 entities, 13 continuous over 125 months; exits 2022m6 (1007→1004),
  2023m6 (1028→1008), 2024m5 (1039→1046) classified as mergers, 2025m7 (1040)
  wind-down; entries 2022m11/2023m12/2024m3.
- 23,798 bank×sector×month cells > 1 bn PYG; 84.9% dual-currency.
- System FX credit share 51.6% (2016m1) → 43.7% (2026m5).
- Sector FX shares (May 2026): Agricultura 0.89, Comercio Mayorista 0.53,
  Industria 0.51, Servicios 0.47, Ganadería 0.45, Consumo/Vivienda 0.08.
- Aggregate reconciliation: micro/aggregate level ratio ≈ 0.95 (financieras
  wedge), YoY growth correlation 0.926.

## Output files (`output/micro/csv/`)

| File | Content |
|---|---|
| `Micro_Shock_Series.csv` | All shock series (standardized, 2016–25 window) |
| `Micro_Aggregate_Reconciliation.csv` | Micro vs aggregate credit reconciliation |
| `Micro_Exit_Events.csv`, `Micro_Entry_Events.csv` | Merger/entry protocol |
| `Micro_Sector_FX_Shares.csv`, `Micro_Sample_Info.csv` | Descriptives |
| `Micro_DesignA_Main.csv` | Headline β_h + full inference battery (incl. VIX placebo) |
| `Micro_DesignA_Placebo_Permutation.csv` | 500-draw permutation distribution |
| `Micro_DesignB_Exposure_Results.csv` | All exposure interactions, both levels, single + joint |
| `Micro_DesignB_Fixed2016_Exposure.csv` | 2016-fixed exposure robustness |
| `Micro_DesignC_Hedging_Results.csv` | Triple-interaction variants |
| `Micro_DesignC_Split_Sample.csv` | Hedged vs unhedged split β_h |
| `Micro_DesignD_Mechanisms.csv` | Deposits and NPL responses |
| `Micro_Robustness_Battery.csv` | All §5 variants at h = 6/12/18 |
| `Micro_Extensive_Margin_LPM.csv` | FX-lending activity LPM |
| `Micro_Run_Log.csv` | Pipeline run log |

PNGs `260–269` in `output/micro/png/` (dynamic β_h, placebos, exposure
facets, hedged-vs-unhedged, deposits/NPL, robustness forest plot).
