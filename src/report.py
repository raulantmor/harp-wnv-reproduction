"""Build the panel registry, numerical comparison tables, and LaTeX report."""
import json
from pathlib import Path
import pandas as pd
import numpy as np
from settings import *

def tex(s):
    s=str(s)
    return ''.join({'&':r'\&','%':r'\%','_':r'\_','#':r'\#','$':r'\$','{':r'\{','}':r'\}'}.get(c,c) for c in s)

def registry():
    rows=[]
    def add(panel,quantity,documentation,data,parameters,command,plot,status,reason):
        rows.append(dict(panel=panel,quantity=quantity,documentation=documentation,data=data,parameters=parameters,
                         data_command=command,plot_command='Not run' if status=='not reproduced' else '.venv/bin/python src/plot.py',plot=plot,status=status,reason=reason))
    for letter,level,stat in [('a','regional','mean'),('b','regional','median'),('c','county','mean'),('d','county','median')]:
        add('1'+letter,f'{level} annual WNND {stat}','Section 2.2; Figure 1 caption','CDC; Census 2020; NEON Oct 2024','2005-2022; zeros included',
            '.venv/bin/python src/prepare.py geography','figure1.pdf','independently reproduced','Descriptive observed-data reconstruction, not a model result; revised inputs and geography differ')
    for n,var in [(2,'PDSI'),(3,'temperature')]:
        for letter,month in zip('abcd',['previous October','January','April','July']):
            add(str(n)+letter,f'Pearson r: ln(WNND) vs {month} {var}',f'Section 2.3; Figure {n} caption',f'CDC + '+('gridMET' if n==2 else 'PRISM'),
                '18 years; log-total-case county weights','.venv/bin/python src/prepare.py climate; .venv/bin/python src/forecast.py',f'figure{n}.pdf','independently reproduced','Independent correlations; numerical differences expected from input vintage and spatial processing')
    add('4','Prairie Peninsula annual probabilistic forecasts','Sections 2.4-2.5; Table 1; Text S2','CDC + available PRISM/gridMET','17 training years; 23 quantiles; NB2',
        '.venv/bin/python src/forecast.py','figure4.pdf','partially reproduced','Published primary soil moisture unavailable; generated restricted-candidate forecasts and benchmark')
    add('5','Regional and national percent WIS improvement','Section 2.5; Text S1; Figure 5','Saved independent forecasts','Log1p WIS; 10000 paired bootstrap samples',
        '.venv/bin/python src/forecast.py','figure5.pdf','partially reproduced','National and regional restricted-candidate analysis omits soil moisture; published-pair scores saved where possible')
    for n,r in enumerate(REGIONS,1):
        add(f'S{n}',f'{REGIONS[r]} monthly climate correlations',f'Section 2.3; Figure S{n} caption','CDC; PRISM; gridMET; soil row missing','Oct(-1)-Jul; mask abs(r)<0.32',
            '.venv/bin/python src/forecast.py',f'figureS{n}.pdf','partially reproduced','Five of six climate rows independently computed; soil moisture row explicitly missing')
    ts=pd.read_csv(RESULTS/'time_series_panels.csv')
    for _,z in ts[ts.figure!='figure4'].iterrows():
        n=z.figure.replace('figure','')
        add(n,f'{REGIONS[z.region]} annual probabilistic forecasts',f'Table 1; Text S2; Figure {n} caption','CDC; PRISM/gridMET','18 LOYO forecasts; NB2; 23 quantiles',
            '.venv/bin/python src/forecast.py',z.figure+'.pdf','independently reproduced' if z.published_predictors else 'partially reproduced',
            'Both published predictor sets independently fitted; input/prior/numerical assumptions documented' if z.published_predictors else 'Soil moisture required by published model unavailable; restricted-candidate curves and benchmark shown')
    add('S22','Absolute WIS improvement','Section 2.5; Figure S22','Saved independent scores','10000 paired bootstraps',
        '.venv/bin/python src/forecast.py','figureS22.pdf','partially reproduced','Restricted-candidate model comparison')
    for n,m in [(23,'univariable'),(24,'bivariable')]:
        add(f'S{n}',f'{m} PIT calibration histogram',f'Figure S{n}; Text S2','Saved posterior predictive draws','198 held-out predictions; 10 equal bins; midpoint PIT',
            '.venv/bin/python src/forecast.py',f'figureS{n}.pdf','partially reproduced','Restricted-candidate models; discrete PIT convention assumed')
    add('Table S3 / Section 3.2 competition comparisons','County mixed-model and regional ensemble forecast skill','Section 2.6; Text S3; Table S3',
        'Competition forecast submissions excluded','Not evaluated','Not run: prohibited comparison inputs','none','not reproduced','County extension not run; competition outputs excluded by independence constraint')
    pd.DataFrame(rows).to_csv(ROOT/'figure_registry.csv',index=False)
    return rows

def figure(name,caption,width=r'\textwidth'):
    return r'\begin{figure}[H]\centering\includegraphics[width='+width+']{../figures/'+name+'.pdf}\caption{'+caption+r'}\end{figure}'+'\n'

def main():
    rows=registry();f=pd.read_csv(RESULTS/'forecasts.csv');skills=pd.read_csv(RESULTS/'skill.csv')
    published=pd.read_csv(ROOT/'docs/published_skill.csv').set_index('region')
    selected=pd.read_csv(RESULTS/'selected_predictors.csv').set_index('region')
    v=json.loads((RESULTS/'validation.json').read_text());d=json.loads((RESULTS/'observed_validation.json').read_text())
    comp=[]
    for r in ORDER+['National']:
        for m in ['uni','bi']:
            s=skills[(skills.region==r)&(skills.model=='available_'+m)].iloc[0]
            comp.append(dict(region=r,model=m,published=published.loc[r,m],published_lower=published.loc[r,m+'_lower'],published_upper=published.loc[r,m+'_upper'],independent=s['median'],difference_pp=s['median']-published.loc[r,m],lower=s.lower,upper=s.upper))
    comp=pd.DataFrame(comp);comp.to_csv(RESULTS/'published_comparison.csv',index=False)
    # Compare published-predictor fits separately from candidate-screening fits.
    exact=[];rng=np.random.default_rng(SEED);ix=rng.integers(18,size=(10000,18))
    for r in ORDER:
        for m in ['uni','bi']:
            a=f[(f.region==r)&(f.model=='published_'+m)].sort_values('year')
            if len(a)!=18:continue
            b=f[(f.region==r)&(f.model=='benchmark')].sort_values('year')
            bb=b.wis.to_numpy()[ix].mean(axis=1);aa=a.wis.to_numpy()[ix].mean(axis=1)
            q=np.quantile(100*(bb-aa)/bb,[.025,.5,.975]);exact.append(dict(region=r,model=m,lower=q[0],median=q[1],upper=q[2],published=published.loc[r,m]))
    exact=pd.DataFrame(exact);exact.to_csv(RESULTS/'published_predictor_skill.csv',index=False)
    national=comp[comp.region=='National'].set_index('model')
    bm=pd.read_csv(ROOT/'docs/published_benchmark.csv').set_index('region').rename(columns={'wis':'published'})
    bm['independent']=f[f.model=='benchmark'].groupby('region').wis.mean()
    bm['difference']=bm.independent-bm.published
    bm.to_csv(RESULTS/'benchmark_comparison.csv')
    pc=pd.read_csv(ROOT/'docs/published_coefficients.csv').set_index('region')
    params=pd.read_csv(RESULTS/'parameters.csv')
    ours=params[(params.model=='published_uni')&(params.parameter=='beta1')].groupby('region')[['mean','lower','upper']].mean().add_prefix('independent_')
    pc=pc.join(ours)
    pc['same_sign']=np.sign(pc['mean'])==np.sign(pc.independent_mean)
    pc.to_csv(RESULTS/'coefficient_comparison.csv')
    npanels=sum(z['status']=='independently reproduced' for z in rows)
    partial=', '.join(z['panel'] for z in rows if z['status']=='partially reproduced')
    reproduced=', '.join(z['panel'] for z in rows if z['status']=='independently reproduced')
    text=r'''\documentclass[11pt]{article}
\usepackage[margin=0.85in]{geometry}
\usepackage{amsmath,graphicx,booktabs,longtable,float,xcolor,caption,subcaption}
\usepackage[authoryear,round]{natbib}
\usepackage[colorlinks=true,linkcolor=blue,urlcolor=blue,citecolor=blue]{hyperref}
\setlength{\parindent}{0pt}\setlength{\parskip}{6pt}
\renewcommand{\thefigure}{R\arabic{figure}}
\title{Independent reconstruction of regional climate-informed West Nile neuroinvasive disease forecasts}
\author{Reproducible research implementation for Raul Moreno}
\date{8 October 2026}
\begin{document}\maketitle
\begin{abstract}
We implement the regional negative-binomial Bayesian regression and retrospective
evaluation described by Harp et al. (2026), using independently acquired public
county disease counts and gridded climate inputs. We reconstruct county and regional
case summaries, climate correlations, probabilistic forecasts, and score uncertainty.
The executed analysis is incomplete: NASA NLDAS soil moisture required by several
published predictor sets returned HTTP 401. We therefore distinguish exact published
predictor fits where inputs are available from an explicitly restricted-candidate
analysis using PRISM and gridMET. Synthetic recovery, independent MCMC, predictive
checks, and scoring identities test the implementation. No original author code,
deposited data, submitted competition forecast, or digitized plot enters the analysis.
\end{abstract}
\section{Scope and evidence}
The target is \citet{harp}, GeoHealth 10(2), e2025GH001657. Written Sections 2.1--2.6,
Tables 1--2, and Supporting Information S1 provide the specification. All plotted
numbers in this report were calculated by this implementation from public inputs.
No published image appears in the report. Reproduced means that we generated the
specified quantity independently; it does not imply numerical identity with the
publication. Input vintages and explicit interpretation choices remain material.
Descriptive case maps and correlation maps are identified separately from model
forecasts. The study-region locator is a schematic and is not a reproduced model result.
Report figures carry an R prefix; source-paper figure identifiers are stated in
their captions and in the panel registry.
'''
    text+=f'Independently reproduced panels: {tex(reproduced)}. Partially reproduced panels: {tex(partial)}. The competition comparisons and county mixed-model extension are not reproduced.\n'
    text+=r'''\section{Inputs and spatial processing}
CDC's public ArboNET dashboard supplies annual county WNND counts. We use its
historical release finalized 16 September 2026 and subset 2005--2022. A separate
CDC file supplies 2022 Connecticut counties, whose rows in the main file are
unavailable planning-region placeholders. Missing county-year activity rows are
interpreted as zero reported cases; explicit missing case values are rejected.
Legacy Shannon/Oglala Lakota and Bedford identifiers are harmonized to Census 2020
boundaries. Alaska, Hawaii, and territories are excluded. County totals match
independently downloaded contiguous-state totals in 17 of 18 years. Two Iowa cases
in 2010 use unknown county code 19999; we retain this documented discrepancy
rather than assign cases to a region without evidence.

Counties are assigned to NEON October 2024 polygons using Census 2020 population
centers. Coastal centers outside a polygon use the nearest polygon in EPSG:5070;
the exception log gives all distances. We use 2020 cartographic county polygons at
1:500,000. Each county climate mean weights grid cells by cosine latitude, including
cells whose centers lie inside the original polygon. Counties without valid cell
centers use the nearest valid cell at their population center. This approximates
fractional-cell integration; display geometries alone are simplified. Domain maps
are clipped to the continental viewing window to omit noncontinental fragments.

The paper's text selects regions with at least 1,000 historical cases, whereas
Figure 1's caption describes a 50-case annual-mean threshold. Current inputs give
Mid Atlantic 987 cases; we retain the 11 explicitly named regions to preserve the
comparison. These 11 regions contain 92.7 percent of the independently assigned
continental county cases, compared with the paper's statement of more than
93 percent. Regional climate weights are $w_c=\log H_c$, where $H_c$ is each county's
total 2005--2022 cases. Counties with $H_c=0$ are excluded and $H_c=1$ has zero weight.
The regional covariate is $x_{rt}=\sum_{c\in r}w_c x_{ct}/\sum_{c\in r}w_c$.

We download 720 PRISM monthly grids from the NC State Climate Office's public
THREDDS service, with Oregon State PRISM identified in their metadata. Variables
are mean temperature, mean dew point (${}^\circ\mathrm{C}$), precipitation
($\mathrm{mm}$ per month), and maximum vapor pressure deficit ($\mathrm{hPa}$).
The seasonal window is previous October--December and concurrent January--July.
These contemporary reprocessed inputs need not equal the paper's 2023 input vintage.
gridMET pentad PDSI is averaged by calendar month before county reduction; the
accessible metadata describe a 1980--2019 reference period. The paper does not
specify how pentads were reduced to months.

NLDAS soil moisture could not be downloaded anonymously (HTTP 401); its legacy
GrADS endpoint returned HTTP 410. The paper does not identify the land-surface
model or version. A proposed, unexecuted interpretation is Noah monthly v2.0,
0--100 mm layer moisture, whose midpoint is 50 mm. No alternate soil product is
substituted here. Our explicit missing soil row limits the six-variable screening
and any model requiring soil moisture. Data versions, URL hashes, access times,
transformations, and failed attempts are recorded in the repository's provenance
files. Large raw downloads are excluded from Git and downloadable by script.
'''
    text+=figure('region_schematic','Study-region locator, redrawn from official NEON polygons. Labels identify the 11 modeled domains; unlabelled domains are excluded. This is a schematic, not a numerical reproduction.')
    text+=r'''\section{Regression and retrospective evaluation}
For each region and year, we model the raw count $Y_t$ by
\begin{equation}
Y_t\mid\mu_t,\phi\sim\mathrm{NB2}(\mu_t,\phi),\qquad
\log\mu_t=\alpha+\sum_j\beta_j z_{tj},\qquad
\mathrm{Var}(Y_t)=\mu_t+\mu_t^2/\phi.
\end{equation}
Here $z$ is the climate predictor centered and scaled by the training fold's mean
and sample standard deviation. The univariable and bivariable models have one and
two climate predictors. The historical benchmark has an intercept alone, with
the same dispersion and prior treatment. All models are independent across regions.
The supplement writes a linear expression for the mean but names
\texttt{stan\_glm.nb}, whose default link is logarithmic. We interpret the displayed
linear expression as the linear predictor and use the log link. Counts are not
logged before fitting the count likelihood.

We infer numerical prior settings from the official rstanarm defaults
\citep{rstanarm}: independent $\alpha,\beta_j\sim N(0,2.5^2)$ on the centered,
standardized scale and $\phi\sim\mathrm{Exponential}(1)$. These are weakly informative,
not flat priors. The written supplement does not enumerate numerical prior scales.
Our sampler includes the Jacobian for $\log\phi$. It uses a 7-degree-of-freedom
multivariate Student proposal centered at the posterior mode, found with an
analytic gradient, with scale 1.2 times the inverse finite-difference Hessian.
Exact posterior-to-proposal weights correct
the proposal. Thus a Gaussian curvature calculation does not replace the posterior.
Four independent batches begin with 8,192 proposals each; sample sizes double if
importance effective sample size is below 4,000 or Pareto $k\geq0.7$. Numerical pass
also requires mean Monte Carlo error below 0.05 posterior standard deviations.
This replaces the paper's four Stan chains of 500 warmup and 500 retained draws.
Independent Metropolis chains test the numerical substitution.

Predictor selection ranks absolute Pearson correlations between $\log Y_t$ and
the 60 potential month-variable combinations. Secondary variables must come from
a different group: air conditions (temperature, dew point, VPD), antecedent soil
conditions (PDSI, soil moisture), or precipitation. The executed restricted analysis
screens 50 available candidates, then uses the same two-group constraint. We also
fit the exact Table 1 predictor sets when their inputs are present. Screening and
county climate weights use the full 18-year record, as described in the paper.
They are not recomputed inside each validation fold; this can favor retrospective
skill and precludes interpreting these scores as a prospective validation.

Each held-out year is predicted from a fit to the other 17 years. Thus early years
are predicted using later observations. We generate 30,000 posterior predictive
integer draws and save 23 empirical quantiles at 1, 2.5, 5, 10, 15,\ldots,95,97.5,
and 99 percent. Seeds, parameter summaries, fold predictions, and diagnostics are
saved. The full-model input pathway supports soil moisture when independent raw
files become available, but those results are absent from this run.

\subsection{Skill and uncertainty}
We apply $\log(1+x)$ to observed counts and forecast quantiles. For $K=11$ central
intervals with lower and upper bounds $l_k,u_k$ and miscoverage $a_k$, the score is
\begin{align}
\mathrm{IS}_{a_k}&=(u_k-l_k)+\frac{2}{a_k}(l_k-y)\mathbf{1}_{y<l_k}
 +\frac{2}{a_k}(y-u_k)\mathbf{1}_{y>u_k},\\
\mathrm{WIS}&=\frac{\frac12|y-m|+\sum_{k=1}^{11}\frac{a_k}{2}\mathrm{IS}_{a_k}}{11.5}.
\end{align}
The median is $m$; the three components are width, overprediction, and
underprediction \citep{bracher}. We allocate the median term to the appropriate
signed component, so the components sum exactly. Lower WIS is better.
Percent improvement is $100(\mathrm{WIS}_{B}-\mathrm{WIS}_{M})/\mathrm{WIS}_{B}$.

Regional uncertainty uses 10,000 paired resamples of 18 years. National scores
weight regional scores by historical cases. National uncertainty resamples both
regions and years and renormalizes those weights; this is an explicit interpretation
of incompletely specified bootstrap details. Year draws are paired between models
but independent between regions, so synchronous regional shocks are not preserved;
this choice can affect national interval widths. The national score is not the score
of a national total-count forecast. Intervals are bootstrap 2.5th and 97.5th
percentiles. PIT uses the midpoint of the predictive CDF jump at an integer
observation; interval coverage includes endpoints.
\section{Numerical results and comparison}
'''
    for m in ['uni','bi']:
        z=national.loc[m]
        text+=f'The restricted-candidate {"univariable" if m=="uni" else "bivariable"} model gives a national median improvement of {z.independent:.1f}\\% ({z.lower:.1f} to {z.upper:.1f}\\%), compared with the published {z.published:.1f}\\% ({z.published_lower:.1f} to {z.published_upper:.1f}\\%). The difference is {z.difference_pp:+.1f} percentage points. These estimates compare different predictor-selection spaces because soil moisture is missing.\n\n'
    text+=f'As a benchmark sanity check, the independent regional mean historical-model WIS differs from the rounded published values in Supplementary Table S5 by a mean absolute {bm.difference.abs().mean():.3f} and a maximum absolute {bm.difference.abs().max():.3f}. The full regional comparison is saved in benchmark\\_comparison.csv. This check isolates the disease data and count-distribution model from the missing climate input.\n'
    text+=r'\begin{table}[H]\centering\small\begin{tabular}{lrrrr}\toprule Region & Paper uni (\%) & Ours uni (\%) & Paper bi (\%) & Ours bi (\%)\\\midrule'+'\n'
    for r in ORDER+['National']:
        z=comp[comp.region==r].set_index('model');text+=f'{tex(REGIONS.get(r,r))} & {published.loc[r,"uni"]:.1f} & {z.loc["uni","independent"]:.1f} & {published.loc[r,"bi"]:.1f} & {z.loc["bi","independent"]:.1f} \\\\\n'
    text+=r'\bottomrule\end{tabular}\caption{Published Table 2 medians versus independent restricted-candidate screening. Published numbers are transcribed text references and do not enter fitting or generated plots. Full intervals and differences are saved in results/published\_comparison.csv.}\end{table}'+'\n'
    text+=r'\begin{table}[H]\centering\small\begin{tabular}{llrr}\toprule Region & Published predictor model & Paper (\%) & Independent median [95\% interval]\\\midrule'+'\n'
    for _,z in exact.iterrows():text+=f'{tex(REGIONS[z.region])} & {z.model} & {z.published:.1f} & {z["median"]:.1f} [{z.lower:.1f}, {z.upper:.1f}] \\\\\n'
    text+=r'\bottomrule\end{tabular}\caption{Separate fits using the exact written Table 1 predictor sets wherever available. Differences can still reflect data vintage, geometry, monthly aggregation, prior interpretation, and numerical integration.}\end{table}'+'\n'
    text+=r'\begin{table}[H]\centering\scriptsize\begin{tabular}{lllll}\toprule Region & Paper primary & Selected primary & Paper secondary & Selected secondary\\\midrule'+'\n'
    for r in ORDER:
        z=selected.loc[r];text+=' & '.join(map(tex,[REGIONS[r],z.paper_primary,z.primary,z.paper_secondary,z.secondary]))+r' \\'+'\n'
    text+=r'\bottomrule\end{tabular}\caption{Predictor names encode variable and calendar month; months 10--12 are previous-year values. Selection uses all years and five available variables. Soil moisture candidates are absent, so equality with published selection is not expected.}\end{table}'+'\n'
    both=selected.primary_match & selected.secondary_match
    text+=f'Independent screening agrees with {int(selected.primary_match.sum())} of 11 published primary predictors and {int(selected.secondary_match.sum())} of 11 secondary predictors. Both match in {int(both.sum())} regions: {tex(", ".join(REGIONS[r] for r in selected.index[both]))}. In several soil-dependent regions, the available-input screen promotes the published secondary variable; this is a change in model specification rather than evidence that the original soil effect is absent.\n'
    text+=r'\begin{table}[H]\centering\small\begin{tabular}{lrrr}\toprule Region & Table S1 slope & Independent slope & Independent bounds\\\midrule'+'\n'
    for r,z in pc.iterrows():text+=f'{tex(REGIONS[r])} & {z["mean"]:+.3f} & {z.independent_mean:+.3f} & [{z.independent_lower:+.3f}, {z.independent_upper:+.3f}] \\\\\n'
    text+=r'\bottomrule\end{tabular}\caption{Univariable coefficients using the exact written predictors, in their raw input units. Each entry averages the posterior mean or 95 percent bound over the 18 held-out fits, matching the stated Table S1 summary. These averaged bounds are not a pooled posterior interval.}\end{table}'+'\n'
    text+=r'''The Northern Plains temperature and Central Plains dew-point slopes have
opposite signs to the entries printed in Supplementary Table S1. Other obtainable
slopes agree closely in sign and magnitude. Our signs follow the independently
calculated input associations and fitted likelihood. The written descriptions do
not resolve whether the disagreement reflects a table error, an undocumented
transformation, or different inputs. We retain the discrepancy; no sign was changed
to match the reference. All fold-specific coefficients and dispersion estimates
are saved in results/parameters.csv.
'''
    text+=figure('figure1','Independent Figure 1 reconstruction: (a,b) regional mean and median annual CDC WNND cases; (c,d) corresponding county quantities, 2005--2022. Grey county polygons have zero cases; hatching indicates excluded domains. Regional colors are linear and county colors logarithmic. These are descriptive observational summaries, not fitted model results.')
    text+=figure('figure2','Independent Figure 2 reconstruction: correlation between log annual WNND and gridMET PDSI in (a) previous October, (b) January, (c) April, and (d) July. Positive PDSI denotes wetter conditions. Red means positive correlation and blue negative; grey domains are excluded. Each correlation uses 18 years and log-case county weights.')
    text+=figure('figure3','Independent Figure 3 reconstruction: the same calculation as Figure 2 for PRISM monthly mean temperature. The data are independently processed public grids, with no values extracted from published maps.')
    text+=figure('figure4','Partial Figure 4 reconstruction for Prairie Peninsula. Black points are CDC observations. Boxes show 25th, 50th and 75th predictive percentiles; whiskers span 2.5--97.5 percent. Grey is the historical benchmark. Green curves come from our restricted-candidate univariable and bivariable models because the published primary soil predictor is unavailable. They are independent calculations but do not complete reproduction of the published model.')
    text+=figure('figure5','Partial Figure 5 reconstruction: relative WIS improvement for the independently selected available-input models. Boxes show bootstrap quartiles and whiskers span 95 percent intervals. National scores weight regional scores by case totals. The absent soil input limits comparability to the published full model.')
    text+=r'\subsection{Timing and magnitude of forecast errors}'+'\n'
    errors=pd.read_csv(RESULTS/'largest_absolute_errors.csv');errors=errors[errors.model=='available_bi'].set_index('region')
    text+=r'\begin{table}[H]\centering\small\begin{tabular}{lrrrr}\toprule Region & Year & Observed & Forecast median & Signed error\\\midrule'+'\n'
    for r in ORDER:
        z=errors.loc[r];text+=f'{tex(REGIONS[r])} & {int(z.year)} & {int(z.observed)} & {z.q11:.0f} & {z.median_error:+.0f} \\\\\n'
    text+=r'\bottomrule\end{tabular}\caption{Largest absolute median forecast error in each region for the available-input bivariable model. Positive values mean overprediction. Units are WNND cases. These errors are calculated from saved held-out forecasts.}\end{table}'+'\n'
    az=f[(f.region=='D14')&(f.year==2021)&(f.model=='available_bi')].iloc[0]
    text+=f'The publication identifies the Desert Southwest 2021 outbreak as an influential failure case. Our 2021 observation is {int(az.observed)} cases; the available-input bivariable median is {az.q11:.0f}, with a 95\\% predictive interval [{az.q01:.0f}, {az.q21:.0f}]. Exact published year-by-year predictions are not reported numerically in the text and were not digitized, so a numerical forecast-to-forecast error comparison is unavailable.\n'
    text+=r'''\section{Validation and limitations}
The likelihood agrees with an independent SciPy negative-binomial parameterization.
The WIS calculation agrees with its 23-quantile pinball-loss identity; a point
forecast equal to the observation scores exactly zero. Identical models remain
identical under the paired bootstrap. Training-fold standardization is isolated
from held-out covariates, and secondary predictor selection respects grouping.
'''
    text+=f'Across {d["fits"]} leave-one-year-out fits, the smallest importance effective sample size is {d["min_importance_ess"]:.0f}, largest Pareto $k$ is {d["max_pareto_k"]:.3f}, and maximum mean Monte Carlo error is {d["max_mcse_over_sd"]:.3f} posterior standard deviations. {d["failed_importance_checks"]} fits fail the specified importance-sampling checks. The optimizer reports precision/status warnings for {d["map_optimizer_warnings"]} proposals; these are retained in diagnostics. Since importance weights correct the proposal, these warnings alone do not establish a biased posterior, but poor weight diagnostics would fail the run.\n'
    text+=f'The known-parameter experiment uses 240 synthetic counts with intercept 3, slopes 0.55 and -0.4, and NB shape 4. All four true values fall inside their 95\\% posterior intervals. Over 24 repeated simulations, coverage is {tex(", ".join(f"{100*z:.1f}%" for z in v["recovery_24_simulations"]))} for intercept, slope 1, slope 2, and shape. Four independent Metropolis chains give maximum R-hat {v["mcmc"]["rhat_max"]:.2f}, minimum bulk ESS {v["mcmc"]["ess_bulk_min"]:.0f}, and posterior mean differences below {v["mcmc"]["max_mean_difference_in_posterior_sd"]:.3f} posterior SD relative to importance sampling.\n'
    text+=figure('validation','Implementation validation using synthetic data, separate from manuscript reproduction. Left: four independent Metropolis traces for a standardized slope. Right: observed synthetic counts and posterior predictive distribution. Known parameters and intervals are saved in synthetic\_recovery.csv.')
    mc=pd.read_csv(RESULTS/'observed_mcmc_checks.csv')
    text+=f'We also compare four Metropolis chains with importance sampling for the observed fold with largest Pareto $k$ and for the Desert Southwest 2021 holdout using the published predictor pair. Across these two checks, maximum R-hat is {mc.rhat_max.max():.3f}, minimum bulk ESS is {mc.ess_bulk_min.min():.0f}, and the largest posterior mean difference is {mc.max_mean_difference_in_posterior_sd.max():.3f} posterior SD. All chain draws and check summaries are saved.\n'
    ppc=pd.read_csv(RESULTS/'posterior_predictive_checks.csv'); nout=int(((ppc.observed<ppc.lower)|(ppc.observed>ppc.upper)).sum())
    text+=f'Full-data posterior predictive checks compare mean, variance, and maximum counts for each of three available-input models in all 11 regions. {nout} of {len(ppc)} observed summaries lie outside their 95\\% posterior predictive intervals. These checks diagnose gross distributional mismatch; they are not independent tests of predictive skill. Held-out coverage, PIT values, and signed errors are saved separately.\n'
    text+=r'''The principal limitation is missing NLDAS soil moisture. Additional discrepancies
may arise from CDC revisions, spatial vintages, coastal assignments, grid-cell
approximation, PDSI monthly averaging, and inferred defaults. Full-record predictor
selection and climate weights allow information from held-out outcomes to enter
the workflow, and are therefore a limitation of the retrospective design itself.
No national result from the restricted-candidate run establishes or refutes the
paper's full-model improvement. The climate associations are predictive regressions
and do not isolate causal biological mechanisms.

The regional and county-level comparisons against the 2022 Forecasting Challenge
are not reproduced because external forecast outputs are excluded by the project
instructions. The county random-intercept extension is not implemented in this
regional pipeline. The main text's reference to logged county counts conflicts
with its supplement's count likelihood and also needs an explicit resolution
before that extension can be claimed reproduced.

\section{Build procedure and outcome}
The environment uses project-local Python 3.11.9 created by
\texttt{/usr/bin/python3.11 -m venv .venv}. All Python commands use the explicit
interpreter path. Exact installed dependencies are pinned in requirements.txt.
The clean procedure is:
\begin{verbatim}
bash scripts/setup.sh
bash scripts/rebuild.sh
\end{verbatim}
The first command recreates only .venv. The second downloads raw public inputs,
reconstructs processed data, fits models, validates numerical methods, creates all
plots, and compiles this report. For this delivered run, input downloads, processing,
fitting, validation, plotting, and PDF compilation were executed; build logs and
the final build outcome are stored in results. Cached raw files make repeat builds
faster. Fold caches are checked against a hash of inputs, predictor specifications,
random seeds, and modeling source before reuse. They are not independent input data.
The repository includes saved
processed inputs and computed outputs, so figures can also be regenerated offline.

TeX Live 2018 packages were checked with kpsewhich before writing the report:
article, geometry, amsmath, graphicx, booktabs, longtable, float, xcolor, caption,
subcaption, natbib, and hyperref. All were present. No system package was installed.
The unavailable siunitx package is not used; units are typeset with standard
mathematical roman text. The requested minimal-research and academic writing
guidance is linked and summarized in docs/guidance.md; no project-specific
vdv-damage content was used.

\section{Figure-by-figure record}
The machine-readable figure\_registry.csv contains one row per panel, source
documentation, input data, parameters, data command, plot command, and status.
For compactness, command keys below mean:
\begin{verbatim}
G: .venv/bin/python src/prepare.py geography
C: .venv/bin/python src/prepare.py climate
F: .venv/bin/python src/forecast.py
P: .venv/bin/python src/plot.py
\end{verbatim}
All numerical plots use saved independent output. R means independently reproduced,
P means partially reproduced, and N means not reproduced. R does not assert exact
agreement. Full source-specific assumptions are given above and in docs/method.md.
\small
\begin{longtable}{p{.09\textwidth}p{.24\textwidth}p{.32\textwidth}p{.21\textwidth}}
\toprule Panel & Quantity & Documentation / inputs / parameters & Status; data / plot command\\\midrule\endhead
'''
    for z in rows:
        command='G / P' if z['panel'].startswith('1') and len(z['panel'])==2 else ('C, F / P' if z['panel'][0] in '23' else 'F / P')
        if z['status']=='not reproduced':command='not run'
        state={'independently reproduced':'R','partially reproduced':'P','not reproduced':'N'}[z['status']]
        text+=tex(z['panel'])+' & '+tex(z['quantity'])+' & '+tex(z['documentation']+'; '+z['data']+'; '+z['parameters'])+' & '+tex(state+'; '+command+'; '+z['reason'])+r' \\'+'\n'
    text+=r'\bottomrule\end{longtable}\normalsize\clearpage\appendix\section{Supplementary figure attempts}'+'\n'
    for n in range(1,25):
        caption=next(z for z in rows if z['panel']==f'S{n}')
        if n<=11: detail='Correlations use 18 years; white cells mask absolute r below 0.32, NA denotes unavailable soil moisture, and black/grey borders mark selected primary/secondary available predictors.'
        elif n<=21:detail='Box/whisker definitions match Figure 4. Titles distinguish published predictor sets from restricted-candidate models. Black points are CDC observations; all forecasts are independently generated.'
        elif n==22:detail='Absolute log1p WIS improvement from 10,000 paired bootstrap samples; positive values favor climate information. Green shades distinguish one and two predictors.'
        else:detail='Proportion of 198 held-out observations in each of ten midpoint-PIT bins. The dashed line is 0.10; nonuniformity indicates calibration departures under the stated discrete convention.'
        text+=figure(f'figureS{n}',tex(f'Figure S{n}: '+caption['quantity']+'. '+caption['status']+'. '+detail))
        if n>=12 and n%2==1:text+=r'\clearpage'+'\n'
    text+=r'''\begin{thebibliography}{9}
\bibitem[Harp et al.(2026)]{harp} Harp, R. D., et al. (2026).
A Regionally Determined Climate-Informed West Nile Virus Forecast Technique.
\emph{GeoHealth}, 10(2), e2025GH001657.
\url{https://doi.org/10.1029/2025GH001657}. Supporting Information S1, Texts S1--S3,
Figures S1--S24, Tables S1--S7.
\bibitem[Bracher et al.(2021)]{bracher} Bracher, J., Ray, E. L., Gneiting, T., and Reich, N. G. (2021).
Evaluating epidemic forecasts in an interval format. \emph{PLoS Computational Biology},
17, e1008618. \url{https://doi.org/10.1371/journal.pcbi.1008618}.
\bibitem[Stan Development Team(2026)]{rstanarm} Stan Development Team.
Prior distributions for rstanarm models. Accessed 8 October 2026.
\url{https://mc-stan.org/rstanarm/articles/priors.html}.
\bibitem[CDC(2026)]{cdc} Centers for Disease Control and Prevention.
Historic West Nile virus data. Finalized 16 September 2026; accessed 8 October 2026.
\url{https://www.cdc.gov/west-nile-virus/data-maps/historic-data.html}.
\bibitem[NEON(2024)]{neon} National Ecological Observatory Network.
Domain polygons, October 2024. \url{https://www.neonscience.org/field-sites/spatial-data-maps}.
\bibitem[Abatzoglou(2013)]{gridmet} Abatzoglou, J. T. (2013). Development of gridded
surface meteorological data for ecological applications and modelling.
\emph{International Journal of Climatology}, 33, 121--131.
\url{https://doi.org/10.1002/joc.3413}.
\end{thebibliography}
\end{document}
'''
    (ROOT/'report/report.tex').write_text(text)

if __name__=='__main__':main()
