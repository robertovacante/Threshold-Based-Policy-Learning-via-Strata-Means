# Threshold Policy Learning via Strata Means (Stata)

This repository provides a Stata implementation to learn simple threshold policies from a stratified randomized experiment, evaluate their policy value with inverse-probability weighting (Horvitz–Thompson), construct stratified bootstrap bands along a percentile path, and visualize group composition and mismatch diagnostics. If you include a short methods note (`methods.pdf`) in the repo, this README will refer to it for background.

## What you need

Use a single `.dta` file that contains identifiers and strata, the experimental assignment, and midline outcomes. Strata are defined as `female × location`. Treatment is the experimental assignment. The code expects: `id`, `female`, `location`, `vti_treated`, and the outcomes `m_has_job`, `m_jobactive`, `m_start_business`, `m_has_savings`, `m_risk_general`, `m_risk_job`, `m_tot_savings`, `m_monthly_income`, `m_life_satis_now`. You can rename or replace these at the top of the do-file.

## How to run

Open Stata in the repo folder and edit the first lines of `tpl_strata_means.do` to point to your data and output folders:

```stata
local DATADIR  "path/to/your/data"
local OUTDIR   "path/to/output"
local DATAFILE "base_midline.dta"
```

Then run:

```stata
do tpl_strata_means.do
```

The script reads `${DATADIR}/${DATAFILE}`, builds a simple strata-mean score, traces a percentile path (default `p = 5,10,…,95`), computes the IPW value and diagnostics at each cutoff, bootstraps stratified percentile bands (default `reps = 500`), and writes plots and logs to `${OUTDIR}/OPL_results/`.

## What the script does

For each outcome, the code computes a score at the stratum level and applies threshold rules given by the empirical percentiles of that score. At each percentile it estimates the policy value via Horvitz–Thompson using the design propensities within strata, forms the difference between the policy value and the RCT at the median, and records the four joint cells between the learned rule and the experimental assignment. From these cells it reports the mismatch rate (reallocations relative to the RCT mix) and the treated share under the rule (coverage). A stratified bootstrap resamples within `female × location` to produce pointwise percentile confidence bands for the policy value and the RCT.

## Outputs

For each outcome, the script writes: (i)	combined_<outcome>.png and .pdf: a two-panel figure with (top) composition bars and mismatch curve and (bottom) value-vs-RCT with a 95% stratified bootstrap band; (ii)	combined_<outcome>.gph: the Stata graph file.

An example of the output is included in this repository. See Figure 1 (figure1.pdf) and Figure 2 (figure2.pdf) for representative panels produced by the script; they illustrate the composition and mismatch diagnostics and the value path with its bootstrap band.
## Customization

Edit the treatment and outcome variable lists near the top of the do-file:

```stata
local treatvar vti_treated
local bin_outcomes  m_has_job m_jobactive m_start_business m_has_savings
local cont_outcomes m_risk_general m_risk_job m_tot_savings m_monthly_income m_life_satis_now
```

Adjust the percentile grid by changing `forvalues p = 5(5)95`. The number of bootstrap replications is controlled by `local reps = 500`. Randomness is set by `set seed 13579`. Aesthetics (scheme and fonts) are optional; remove those lines if the fonts are not installed.

## Reproducibility and interpretation

Bootstrap results are reproducible given the fixed seed. The IPW step uses stratum propensities estimated within each bootstrap draw, consistent with design-based evaluation. With discrete, stratum-level scores, the value curve changes in steps because a small threshold shift can toggle an entire stratum; mismatch tends to grow as the learned rule departs from the experimental allocation, which increases variance under IPW.

## Troubleshooting

If bands look very wide or flat, check for small or near-deterministic strata; these inflate variance. If Stata reports missing variables, align your dataset’s names with those referenced above. For methodological details and notation, see `methods.pdf`.
