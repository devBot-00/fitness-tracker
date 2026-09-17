# Kalman-filter maintenance estimator for HDH Log

Files
- `maintenance-kf.js` — the estimator. No dependencies. Loads as a `<script>` (global `MaintKF`) or via `require`.
- `ols_ar1.js` — a faithful port of the estimator currently in `index.html` (reproduces 2,698 ± 810 on Sep 5–16 exactly). Used only for benchmarking.
- `bench.js` — Monte Carlo benchmark; `node bench.js 400` reruns it. `bench_results.txt` is the last run.
- `real_data.js` — the Sep 5–16 rows used in the checks below.

## What it is

A 3-state Kalman filter run one step per calendar day:

| state | meaning |
|---|---|
| `W` | tissue weight, lb (glycogen/water excluded) |
| `M` | maintenance, kcal/day |
| `U` | AR(1) water/gut noise, lb (same carryover idea as the AR(1) band in the app) |

Dynamics: `W[t+1] = W[t] + (I[t] − M[t])/ρ`, `M[t+1] = M[t] + β·ΔW + drift`, `U[t+1] = φ·U[t] + noise`.
Observation: `scale[t] = W + U + G + v`, where `G` is a glycogen/water offset computed from logged carbs.

Same energy-balance idea as the OLS, but recursive and Bayesian, which is what buys the four things below.

## What changed vs. the current OLS, and what it's worth

Benchmark: 400 simulated cuts per scenario. The simulator is richer than the filter's model on purpose — adaptive thermogenesis, AR(1) water plus carb glycogen with different constants than the filter assumes, 4% outlier weigh-ins, 8% missing weigh-ins, weekends eaten bigger and unlogged 30% of the time, ±75 kcal daily logging error, true energy density 3,100–3,350. Error is estimate − true *current* maintenance.

Steady cut (glycogen already settled):

| day | OLS (app) RMSE / band | KF + formula prior RMSE / band |
|---|---|---|
| 12 | 316 / ±1,431 | 194 / ±461 |
| 21 | 182 / ±564 | 121 / ±268 |
| 28 | 123 / ±353 | 86 / ±193 |
| 42 | 89 / ±163 (89% coverage) | 63 / ±125 (96%) |
| 60 | 89 / ±162 (89% coverage) | 50 / ±96 (94%) |

Cut starting day 0 (2–3 lb water leaves in week one): OLS reports nothing until day 20 (10-day skip + 10-day span); KF reports from day 1 (prior-dominated), and by day 28 is at RMSE 147 vs 227, day 42 at 83 vs 113.

Where the gain comes from (ablations in `bench_results.txt`):
1. **Energy density 3,200 in a cut instead of 3,500.** Most of the late-window gain and the reason OLS coverage drops to 89% after 6 weeks (its band is right-sized but the estimate is biased +35–50). Pure fat is 3,500 kcal/lb; a high-protein cut loses 10–15% lean, ~3,200.
2. **A prior on M** (Katch–McArdle, or the previous phase's estimate). Day-12 RMSE 239 → 194, band ±600 → ±460, and the phase-start case goes from "building" to a usable number. Cost: a prior that's 300 kcal wrong biases day 12 by ~150 and is washed out by day 28–42.
3. **Robust gate + carb glycogen term.** Day 21–42 RMSE ~15% lower and bias near zero. A weigh-in more than 2.5σ from expectation is downweighted, not discarded (9/15's 149.8 is the example). Carb swings are subtracted before they can pose as tissue.
4. **Honest band.** Coverage is 94–97% at every horizon; the OLS band is 3× too wide early (98–99% coverage at ±1,400) and too narrow late.

Also: unlogged days widen the band instead of silently averaging; there's a daily time series of `M ± se` (`series`) for a chart; and `daysTo100` projects how many more daily weigh-ins until the 95% band is ±100.

What it does **not** fix: under-logging. If unlogged days really run +250 (weekends), every estimator reads ~35 kcal low. Set `unloggedBump` if that's the pattern.

## Sep 5–16, your data

| | maintenance |
|---|---|
| OLS (app) | 2,698 ± 810 |
| KF, no prior | 2,603 ± 543 |
| KF, Katch–McArdle prior (2,488 ± 300) | **2,538 ± 398** |
| KF, prior, carb term off | 2,471 ± 537 |

9/15 gated (innovation −2.85 lb). Projected 33 more daily weigh-ins to reach ±100.

## Integration (index.html)

Paste the body of `maintenance-kf.js` into a `<script>` before the app script (it defines `MaintKF`). Run it **once over the whole log** — it carries information across phase boundaries itself, so the carried-forward / settle-skip / window-chip plumbing in `estimateMaint()` isn't needed for it:

```js
function estimateMaintKF(){
  const data=allData(), rows=MaintKF.adapterFromLog(data, dayTotals);          // {k, weight, cals, carbs} per day
  const first=rows.find(r=>r.weight); if(!first) return { status:'building' };
  const bf0=(data[first.k].body||{}).bf || 18;
  return MaintKF.estimateMaintKF(rows, {
    end: todayKey(),
    prior: MaintKF.katchPrior(first.weight, bf0, 1.55),                          // anchors only the start of the log; data takes over
    phaseStarts: phaseList().map(p=>({ k:p.start, phase: p.type==='maint'?'maintain':p.type })),
    phase: 'cut',                                                                // phase before the first logged phase
    // unloggedBump: 250,                                                        // if unlogged days are usually bigger
  });
}
// K.tdee, K.ci, K.trendW, K.lbWk, K.pctWk, K.balance, K.avgIn, K.daysTo100,
// K.series[{d, tdee, se, weight, gated}]  → chart of maintenance ± band over time
// K.warnings[]                            → same slot as the OLS warnings
```

Per phase it switches energy density (`kcalPerLbByPhase`: cut 3,200 / maintain 3,500 / bulk 3,000), widens the maintenance band by `phaseInflateM` (75 kcal — activity usually changes with the phase), and for `settleDays` (10) after each phase start lets weigh-ins anchor weight but teach nothing about maintenance — the 10-day skip without discarding the data or the prior. The filter's memory is set by `sigmaM` (6 kcal/day per day of maintenance drift), so older data fades on its own; the 14/28/42 window chips can go, or stay as a display filter for the chart.

To show both estimators side by side for a while (recommended: they should agree within the bands once you're 4+ weeks into a steady phase), keep `estimateMaint()` as is and add the KF card next to it.

Parameters worth exposing in the UI: activity multiplier for the prior (1.4–1.7), `unloggedBump`, and a `kcalPerLb` override. Everything else has a defensible default in `DEFAULTS`.
