# Clinical Trial Statistics Demo (R)

This repository is a reproducible, simulated clinical-trial statistics example in R.
Designed to show implementation of typical clinical trial analyses.

## What this demo includes

Typical analyses in clinical trial statistics include:

- Baseline characteristics tables
- Continuous endpoint analysis with ANCOVA
- Binary endpoint analysis with logistic regression
- Time-to-event analysis with Kaplan-Meier curves and Cox regression
- Repeated-measures analysis for longitudinal outcomes
- Subgroup analyses with a simple forest plot
- Adverse event summaries

The data in this repo are simulated, so it is safe to run and share.

## Files

- `R/analysis.R` — generates the dummy data, runs the analyses, and saves tables/figures
- `results/` — CSV and text outputs created by the script
- `figures/` — PNG figures created by the script

## How to run

From the repository root:

```r
Rscript R/analysis.R
```

After running, you should see outputs such as:

- `results/baseline_summary.csv`

| Variable           | Placebo     | Active      |
| ------------------ | ----------- | ----------- |
| Age, mean (SD)     | 58.4 (11.3) | 59.2 (11.4) |
| Female, n (%)      | 50 (41.7%)  | 59 (49.2%)  |
| High strata, n (%) | 51 (42.5%)  | 52 (43.3%)  |
  
- `results/continuous_endpoint.csv`

| Endpoint                     | Estimator                  | Estimate | CI_low | CI_high | P_value  |
| ---------------------------- | -------------------------- | -------- | ------ | ------- | -------- |
| Week 12 change from baseline | ANCOVA (Active vs Placebo) | -8.16    | -10.5  | -5.81   | 7.63e-11 |

- `results/binary_endpoint.csv`

| Endpoint             | Estimator                                  | OR   | CI_low | CI_high | P_value  |
| -------------------- | ------------------------------------------ | ---- | ------ | ------- | -------- |
| Responder at week 12 | Logistic regression OR (Active vs Placebo) | 6.52 | 3.59   | 11.83   | 6.87e-10 |

- `results/time_to_event.csv`

| Endpoint            | Estimator                     | HR   | CI_low | CI_high | P_value |
| ------------------- | ----------------------------- | ---- | ------ | ------- | ------- |
| Time to first event | Cox PH HR (Active vs Placebo) | 0.59 | 0.42   | 0.81    | 0.00112 |

- `results/mmrm_summary.txt`

- `results/subgroup_forest.csv`

| Subgroup | Level  | Estimate          | Low               | High              |
| -------- | ------ | ----------------- | ----------------- | ----------------- |
| strata   | High   | -9.27891183076301 | -12.7920968602426 | -5.76572680128337 |
| strata   | Low    | -7.06562206506464 | -10.3324723227326 | -3.7987718073967  |
| sex      | Female | -7.79284266278517 | -11.4304741680476 | -4.15521115752276 |
| sex      | Male   | -8.50758852792628 | -11.67145662151   | -5.34372043434259 |

- `results/adverse_events_summary.csv`


| trt     | subject_count | subjects_with_any_ae | events_total | subjects_with_serious_ae | percent_any_ae | percent_serious_ae |
| ------- | ------------- | -------------------- | ------------ | ------------------------ | -------------- | ------------------ |
| Placebo | 120           | 65                   | 100          | 6                        | 54.2           | 5                  |
| Active  | 120           | 84                   | 146          | 15                       | 70             | 12.5               |
  

### Figures

- **Figure 1.** Kaplan-Meier curve for time to first event by treatment arm.
![png](figures/km_curve.png)
- **Figure 2.** Mean biomarker trajectory over time with SD bars.
![png](figures/longitudinal_means.png)
- **Figure 3.** Forest plot of subgroup treatment effects for the continuous endpoint.
![png](figures/forest_plot.png)
## Notes on the analyses

The code uses a fairly standard trial-statistics workflow:

- `lm()` for the primary continuous endpoint
- `glm(..., family = binomial())` for the responder analysis
- `survival::survfit()` and `survival::coxph()` for time-to-event analyses
- `nlme::lme()` for a repeated-measures model
