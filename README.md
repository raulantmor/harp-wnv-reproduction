# Independent Harp et al. (2026) reproduction

This repository reconstructs the regional Bayesian WNND forecasting method in
Harp et al., *A Regionally Determined Climate-Informed West Nile Virus Forecast
Technique*, GeoHealth 10(2), e2025GH001657,
[doi:10.1029/2025GH001657](https://doi.org/10.1029/2025GH001657).

**This is a partial empirical reproduction.** The Python regional pipeline is
implemented and exercised with independent CDC, Census, NEON, PRISM, and gridMET
inputs. NLDAS soil moisture requires authentication and was not obtained. The user
approved documenting dependent panels as partial. No substitute soil product or
published figure is used. The county mixed-model extension and competition
comparisons are not reproduced. This is retrospective association-based modeling,
not an operational prospective forecast.

- [Compiled PDF report](report/report.pdf)
- [LaTeX source](report/report.tex)
- [Panel-by-panel registry](figure_registry.csv)
- [Method, equations, priors, assumptions](docs/method.md)
- [Data provenance and failed attempts](docs/provenance.md)
- [Validation results](results/validation.json)
- [Observed-fit diagnostics](results/observed_validation.json)
- [Published vs independently calculated skill](results/published_comparison.csv)

The documented rebuild completed on 2026-10-08 with exit code 0 and produced a
30-page report. It used Python 3.11.9; all 792 forecast fits passed importance
sampling checks with zero optimizer warnings, six tests passed, and synthetic
and observed-data MCMC comparisons passed. The final restricted-input national
skill improvements were 15.4% (95% interval 4.6–25.5%) and 18.7% (5.6–30.3%).
These are partial comparisons with the published 18.8% and 21.8% improvements.
See [build status](results/build_status.json), [full rebuild log](results/rebuild_run.txt),
and [final figure/report rebuild log](results/presentation_rebuild.txt).
Files labeled `initial` preserve earlier setup and diagnostic runs.

## Rebuild

Linux with `/usr/bin/python3.11`, network access, and existing TeX Live 2018
(`pdflatex`, `latexmk`, `kpsewhich`) is required. No system package installation is
performed. The system Python 3.6.8 is never used. All dependencies are pinned.
Run the following commands from the repository root.

```bash
bash scripts/setup.sh
bash scripts/rebuild.sh
```

`setup.sh` recreates the project's `.venv` and reports its Python version.
Alternatively `bash scripts/rebuild.sh --fresh-env` performs both steps. Raw
downloads occupy roughly a few GB and are cached under `data/raw/`, excluded
from Git. The download scripts record URL, access time, response metadata, size,
and SHA-256 hashes. The PDSI service can be slow; its downloader uses three
independent processes, each with small validated requests. Re-run after an
interrupted download; incomplete `.part` files are not processed.

Individual steps:

```bash
.venv/bin/python src/download.py base
.venv/bin/python src/download.py prism
.venv/bin/python scripts/pdsi_parallel.py
.venv/bin/python src/prepare.py geography
.venv/bin/python src/prepare.py climate
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/harp-matplotlib .venv/bin/python src/forecast.py
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/harp-matplotlib .venv/bin/python src/postcheck.py
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/harp-matplotlib .venv/bin/python src/validate.py
.venv/bin/python -m pytest -q tests
.venv/bin/python src/plot.py
make report
```

`make all` runs the same sequence. `PYTHON` defaults to `.venv/bin/python`.
The report uses installed standard LaTeX packages; it does **not** use `siunitx`.
All units use ordinary math typography. Dependencies and interpreter paths are
contained entirely within the project.

To regenerate plots and the report offline from the saved independent results:

```bash
.venv/bin/python src/plot.py
make report
```

Fold caches validate a hash of covariates, observations, predictor names, seeds,
and modeling source before reuse. The posterior predictive draws in
`results/predictive_draws.npz` and quantiles in `results/forecasts.csv` are independent
outputs, not copied publication data. Full raw climate grids are never committed.
Upstream live datasets can change; recorded hashes and the committed processed
inputs identify the executed version.

## Scientific boundaries and outputs

Written paper and supplement descriptions are specifications. Links and an account
of the two applied guidance documents are in `docs/guidance.md`. No original author repository,
institutional deposit, author-derived dataset, forecast submission, plotting
script, or digitized plot was consulted. Published Table 2 values are transcribed
in `docs/published_skill.csv` and benchmark scores from Supplementary Table S5 in
`docs/published_benchmark.csv` solely for explicitly labeled numerical comparisons.
`docs/published_coefficients.csv` similarly records six slopes from Table S1;
two sign disagreements are discussed in the report.

The independent observed-data summaries and correlation maps reconstruct Figures
1–3. NLDAS-dependent forecasts and all six-variable supplementary heatmaps remain
partial. Where inputs permit, the exact Table 1 predictor pairs are fitted in
parallel with independent screening; the registry distinguishes those results.
Figure 4, Figure 5, and the national restricted-candidate skill estimates must not
be interpreted as reproductions of the complete six-variable model.

`available_uni` and `available_bi` screen the obtainable candidates.
`published_uni` and `published_bi` use the corresponding Table 1 specification,
only when every required predictor exists. `benchmark` is intercept-only NB2.
Every region/year/model has quantiles, signed error, WIS components, interval
coverage, and a midpoint-PIT value. Bootstrap draws, parameter summaries, selected
predictors, posterior predictive checks, and convergence diagnostics are saved.

An authenticated, independent NLDAS_NOAH0125_M.2.0 file collection may be placed in
`data/raw/nldas/` before processing. The reader expects `SoilM_0_10cm`, `time`,
`lat`, and `lon`. This chosen product interprets the unspecified 50-mm depth as
the 0–100-mm layer midpoint and requires a documented sensitivity analysis; it
does not resolve the paper's missing product detail automatically. This pathway
has not been validated against real NLDAS files in the delivered run.

## Upload destination

The only authorized remote is
`git@github.com:raulantmor/harp-wnv-reproduction.git`.
GitHub SSH commands use:

```bash
GIT_SSH_COMMAND='ssh -F /dev/null -i ~/.ssh/id_ed25519 -o IdentitiesOnly=yes' git ls-remote git@github.com:raulantmor/harp-wnv-reproduction.git
```

Credentials, raw downloads, virtual environments, and temporary build files are
excluded from the repository.
