# Revision-Extras Phase (IREF Revision, Coworker-Advised Exercises)

Implements the empirical/econometric exercises recommended in the coworker
advice on the IREF-D-26-04765 decision (see `COWORKER ADVISES ON JOURNAL
DECISION.docx` and `COWORKER ASSESMENT MICRO_EXERCISES.docx`), as filtered
through the expert assessment. **Fully additive**: nothing in the aggregate
pipeline (01–23), the micro phase (30–36), the data, or the paper is touched.
All outputs go to `output/revision/{csv,png,external}`.

Scope decisions (user-approved): external series fetched automatically from
FRED (cached in `output/revision/external/`, `data/` untouched); the
FOMC-window instrument and the cross-country dollarization gradient are
**deliberately skipped** (first-to-cut items; name in limitations).

## How to run

```bash
cd R/
Rscript RUN_REVISION.R                  # full phase (~10 min)
MICRO_WCB_B=99 Rscript RUN_REVISION.R   # quick pass (fewer bootstrap reps in 37)
```

Requires the micro phase (`RUN_MICRO.R`) and script 36 to have run once
(reads their rds/csv caches). Script 43 must precede 44 (handled by the
runner).

## Scripts and provenance

| Script | Exercise | Advice source |
|---|---|---|
| `revision_helpers.R` | Shared LP/IV engines (script-05/23 conventions), AR sets, MOP effective F, FRED fetch, rolling/expanding standardization | — |
| `37_Micro_Freeze_Checks.R` | Design C gradient & split ex-COVID / tightening-only; h=1–3 TCN-timing audit; reconciliation wedge decomposition | Coworker 1 (micro assessment) |
| `38_CHR_Plausibly_Exogenous.R` | Conley–Hansen–Rossi (2012) plausibly-exogenous bounds with trade-calibrated γ | Coworker 1 §1 / WP1-3; Coworker 3 C2-4 |
| `39_PostIT_ExpandingFCI_IV.R` | Expanding-window FCI reconstruction (validated vs appendix K.5/N.6 pattern); rolling-vs-expanding table; **new** post-IT IV-LP with AR sets + effective F | Coworker 1 §2 / WP1-5; Coworkers 2, 3 |
| `40_LagAugmented_LP.R` | Lag-augmented LP with EHW SEs (MOP-PM 2021) vs NW baseline | Coworker 1 WP6 |
| `41_Falsification_Placebos.R` | Insulated-placebo battery (binational electricity, public consumption, agriculture) vs responding outcomes | Coworker 3 C1-2 |
| `42_FCI_Component_Exclusion.R` | Leave-one-out FCI battery + quantities-at-t−1 FCI; rates-only co-baseline row | Coworker 3 C5-2; Coworker 1 WP2 |
| `43_External_Data_Fetch.R` | FRED fetch: broad dollar (DTWEXB/DTWEXBGS splice), BRL, EUR, GBP, TED; Wu–Xia best-effort | Coworker 1 WP1-1/2, WP6 |
| `44_Enhanced_Instrument_Robustness.R` | Enhanced purge (+Selic, BRL); broad-dollar IV; **EUR/GBP dollar-free placebo**; TED funding-stress horse race | Coworker 1 WP1-1/2, WP6 |
| `45_COVID_Reclassification_Check.R` | COVID-reclassified balances (Medida Excepcional/Medidas transitorias): bias signing, bank×currency add-back LP, apportioned gradient sensitivity | Colleague 1 (feedback on exercises) — "the one remaining mandatory exercise" |
| `46_USD_ShiftShare_Decomposition.R` | Within/between shift-share decomposition reconciling the aggregate USD contraction with the positive within-cell differential | Colleague 1 (feedback on changes & response) — "mandatory" |
| `RUN_REVISION.R` | Master runner + `Rev_Run_Log.csv` | — |

## Documented deviations and judgment calls

1. **EUR/USD is not used as the placebo** (Coworker 1's Step-2 suggestion):
   EUR/USD is itself a dollar price and shows the mirror-image first stage
   (t ≈ −2.9), which is *confirmatory*, not falsifying. The clean dollar-free
   placebo is the **EUR/GBP cross** (first stage t = 0.2, F = 0.1).
2. **Levels vs innovations**: the credit-relevant dollar variation is the
   persistent (levels) component — monthly innovation instruments (raw or
   purged) have strong first stages but near-zero IV coefficients at h ≤ 12.
   Exogeneity of the levels instrument is defended by the CHR bounds (38),
   the ER controls, and the placebo architecture, not by innovation-purging,
   which strips the persistent identifying variation together with the
   confounds. This nuance belongs in the response letter.
3. **Expanding-window validation is qualitative**: the appendix table's
   producing script is not in the repo; coefficient signs/magnitudes and the
   rolling-collapse / expanding-strength pattern reproduce under both
   classical and HAC conventions (gates printed in 39).
4. **LOO variants are z-score-method** approximations (anchored: full-12
   corr 0.96 with `FCI_COMP_AVG`; published 4-method `FCI_exCredit_AVG` and
   `FCI_RATES_AVG` rows included in the same table for direct comparison).
5. **Wu–Xia shadow rate** download failed at run time (Atlanta Fed URL);
   purge v2 runs without it (slot in 43 retries on next run). TED ends
   Jan 2022 (discontinued).
6. **Lee et al. (2022) tF** deliberately not implemented (applies narrowly to
   just-identified 2SLS; AR sets + effective F already cover weak-IV
   inference — avoid stacking procedures). "Monotonicity" framing (Coworker 3)
   also deliberately dropped.

## Key outputs (`output/revision/csv/`)

`Rev_DesignC_Gradient_Robustness.csv`, `Rev_DesignC_Split_Robustness.csv`,
`Rev_TCN_Timing_Audit.csv`, `Rev_Reconciliation_*.csv`, `Rev_CHR_Bounds*.csv`,
`Rev_ExpandingFCI_*.csv`, `Rev_PostIT_IV_ExpandingFCI.csv`,
`Rev_LagAugmented_LP*.csv`, `Rev_Falsification_*.csv`, `Rev_FCI_LOO_*.csv`,
`Rev_External_Series.csv`, `Rev_Instrument_Comparison.csv`,
`Rev_BroadDollar_IV_LP.csv`, `Rev_FundingStress_HorseRace.csv`,
`Rev_Run_Log.csv`. Figures `280–288` in `output/revision/png/`.

Results narrative: `output/revision/Revision_Exercises_Report.md`.
