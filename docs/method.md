# Independent specification and interpretation

Source: Harp et al. (2026), GeoHealth 10(2), e2025GH001657,
https://doi.org/10.1029/2025GH001657, Sections 2.1–2.6, Tables 1–2, and published
Supporting Information S1, Texts S1–S3 and figure captions. The supplement was
obtained from Europe PMC's `supplementaryFiles` endpoint; only its PDF was
extracted. No deposited author code, data, forecasts, or figure assets were used.

## Observations and spatial support

Use county annual WNND counts for 2005–2022. Assign each county to a NEON domain
using its Census 2020 population center. Add CDC's separately published 2022
Connecticut file; discard its planning-region placeholder rows. Map old Shannon
County 46113 to Oglala Lakota 46102, and former Bedford city 51515 into Bedford
County 51019 to match 2020 boundaries. Counts absent from CDC's activity table
are treated as zero reported cases, not as absence of infection. Explicit missing
values are rejected, rather than zero-filled. Exclude AK, HI and territories.

The paper does not state the vintage of population centers or boundaries.
Use Census 2020 centers, 2020 1:500,000 cartographic county polygons, and October
2024 NEON polygons. Population centers outside a domain polygon (coastal islands
and boundary differences) are assigned to the nearest polygon in EPSG:5070;
distances are saved. These assignments can change regional cases and weights.
The 11 named regions are retained even though the current data give Mid Atlantic
987 cases rather than the stated minimum 1,000. The caption's 50 mean annual
cases and the text's 1,000 total cases are also different thresholds. The 11 regions
cover 92.7% of assigned continental county cases, compared with the reported >93%.

Compute original-grid county means using grid cells whose centers lie in a
county, weighted by cosine latitude. Use the nearest valid cell at the population
center for counties with no valid cell centers, recording counts per grid.
This is a numerical approximation to fractional-cell polygon integration; it
does not claim exact reproduction of the unspecified original zonal statistic.
Only display geometries are simplified, by 0.002 degrees. Actual means use the
unsimplified downloaded county polygons. Displayed domain polygons are clipped
to 125–66 degrees west and 24–50 degrees north to omit noncontinental fragments.

For county c, let H_c be total reported WNND cases in 2005–2022. Regional climate
is sum(log(H_c) * x_c) / sum(log(H_c)), excluding H_c=0. H_c=1 contributes zero
weight. Weights use all 18 years, following the written historical weighting.
This uses information from the held-out year and is not a prospective procedure.

## Climate inputs and selection

Use PRISM monthly tmean, tdmean, ppt, vpdmax at 4 km, downloaded from the NC State
Climate Office public THREDDS service. These grids identify Oregon State PRISM
as their source and are general public climate inputs unrelated to the target
paper. Temperature and dew point are degrees Celsius, precipitation is mm/month,
and maximum VPD is hPa, subject to each file's metadata. These are contemporary
reprocessed grids, not the paper's frozen 2023 input vintage.

Use gridMET's public pentad PDSI, averaging available pentad endpoint values
within each calendar month. The paper does not specify its monthly reduction.
The accessible gridMET metadata identify a 1980–2019 reference period. No
published correlation values are used to construct predictor data.

NLDAS is specified by the paper, but its land-surface model and product version
are not. A future authenticated run can use NLDAS_NOAH0125_M.2.0 SoilM_0_10cm
(kg/m^2), interpreting 50 mm as the midpoint of the top 0–100 mm layer. This is
an explicit interpretation, not an established fact about the paper's input.
NASA HTTPS returned 401 and its old GrADS endpoint returned 410 during this run.
No alternate soil product is substituted. Files from that chosen NLDAS product
can be placed under data/raw/nldas/ and read by the processing code.

For each region, compute Pearson correlations between ln(annual cases) and each
climate variable in previous October–December and concurrent January–July.
Select the largest absolute correlation, then the largest from a different group:
(tmean, tdmean, vpdmax), (PDSI, soil moisture), or (precipitation). Ties use lexical
predictor order. The primary analysis performs this screening once over all years,
as described in the paper. It is not nested cross-validation. With missing soil
moisture, this produces a clearly labeled restricted-candidate analysis.
In parallel, fit the exact published Table 1 predictor sets wherever available.
Table 1 is a permitted written model specification, not derived input data.

## Likelihood, priors, and computation

For region r and year t, Y_rt ~ NB2(mu_rt, phi_r), with
log(mu_rt) = alpha_r + sum_j beta_rj z_rtj and
Var(Y_rt) = mu_rt + mu_rt^2 / phi_r. Each model is fit separately by region.
Predictors z are centered and scaled using only the training fold's mean and
sample standard deviation. Coefficients on this scale have independent N(0,2.5^2)
priors, including the centered intercept. phi has an Exponential(1) prior.
These implement the documented rstanarm default scales for a non-Gaussian GLM.
The supplement does not print prior numerical values. Its linear mean expression
omits the default log link and would permit negative means; we use the documented
stan_glm.nb log link. Raw integer cases, not logged counts, enter the likelihood.
The county extension's prose also conflicts with its NB equation about logging.

Inference uses self-normalized importance sampling with a multivariate Student-t
proposal (7 degrees of freedom; scale 1.2 times inverse MAP Hessian) and exact
posterior/proposal weights. An analytic gradient locates the posterior mode, and
a finite-difference Hessian estimates proposal curvature;
it does not replace the posterior with a Gaussian. Four independent batches start
with 8,192 proposals each; sizes double up to four attempts when effective sample
size is below 4,000 or Pareto k is at least 0.7. Save raw-weight ESS, Pareto k,
optimizer status, and between-batch Monte Carlo error. Numerical pass also requires
maximum mean MCSE / posterior SD below 0.05. Importance sampling replaces the
paper's 4 Stan chains of 500 warmup + 500 retained draws and is checked against
four independent random-walk Metropolis chains in a known-parameter experiment.
The same comparison is made for the observed fold with largest Pareto k and for
the Desert Southwest 2021 published-predictor holdout. These chains use 4,000
burn-in iterations and 10,000 retained iterations each, with proposal scale
2.38/sqrt(d) times a Cholesky factor of the MAP covariance. Pass requires R-hat
at most 1.01, bulk ESS at least 400, and posterior mean differences below 0.15 SD.

Use 30,000 posterior predictive integer draws per held-out case to estimate the
23 specified quantiles (1%, 2.5%, 5%, 10%, 15%, ..., 95%, 97.5%, 99%), with an
inverted empirical CDF. Fixed random seeds are specified in src/settings.py.
The benchmark is the identical Bayesian model with only an intercept and phi.
All 18 leave-one-year-out fits train on the remaining 17 years, including future
years relative to the target year. This estimates retrospective association-based
forecast performance and does not establish prospective forecast skill.

## Scoring and uncertainty

Transform both observed counts and forecast quantiles by ln(1+x). For the 11
central prediction intervals, alpha_k=2q_lower,k, w_k=alpha_k/2, w_0=1/2, and
WIS = [w_0 |y-m| + sum_k w_k IS_alpha,k] / 11.5. IS is interval width plus
2/alpha times the amount outside the interval. The median term is assigned to
over- or underprediction according to sign. Components sum exactly to WIS.

Resample 18 paired years with replacement 10,000 times for regional scores.
For national uncertainty, also resample 11 regions, retaining historical case
weights and normalizing after resampling. Year draws are paired between competing
models but independent between regions. This does not preserve synchronous
regional shocks and can affect national interval widths. The paper's wording does not uniquely
specify whether regions are resampled or only years; this choice is documented.
National score is a weighted average of regional WIS, not the WIS of summed
national count forecasts. Percent improvement is 100*(WIS_base-WIS_model)/WIS_base.
Report 2.5%, 50%, 97.5% bootstrap percentiles and the fraction of nonpositive
improvements, not a posterior probability of a causal effect.

PIT uses the midpoint of the predictive CDF jump at the observed integer count.
The supplement does not specify treatment of discrete ties. Interval coverage
uses inclusive endpoints. Save largest absolute and signed forecast errors.

## Scope limits

The county mixed model and competition comparisons are not reproduced in this
run. The public competition outputs are excluded by the task's independence
constraint. No numerical competition forecast or published figure is substituted.
The full regional code supports all six covariates, but a complete all-input
empirical reproduction requires the missing NLDAS input. Restricted-candidate
results must never be reported as an exact recreation of the paper's full model.
