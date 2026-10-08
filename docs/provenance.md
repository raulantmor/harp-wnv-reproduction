# Data provenance and permitted-source boundary

Access date for this run: **2026-10-08 UTC**. Exact per-download timestamps,
response metadata and SHA-256 hashes are stored in `data/download_manifest.jsonl`.
Files downloaded during initial source discovery are indexed in
`data/initial_download_manifest.json`. Large raw inputs are excluded from Git.
`docs/published_skill.csv` is a hand transcription of **Table 2**,
`docs/published_benchmark.csv` contains benchmark scores from **Table S5**, and
`docs/published_coefficients.csv` contains the six obtainable univariable slopes
from **Table S1**. These are comparison references only; they are not inputs to
any fit or generated manuscript panel.

| Input | Version / scope | Independent source | Processing |
|---|---|---|---|
| WNND county counts | CDC live historical release finalized 2026-09-16; subset 2005–2022 | CDC ArboNET historical dashboard's downloadable CSV | FIPS normalization, CT repair, zero reporting assumptions, county-to-NEON sums |
| CT 2022 correction | CDC separate legacy-county CSV | CDC dashboard `ct-2022.json` | Replace unavailable planning-region placeholders; preserve actual 090xx county IDs |
| State WNND totals | CDC historical dashboard release | CDC `wnv_hum_historic - Historic Sta.csv` | Independent national conservation check; exclude AK, HI, PR |
| County population centers | 2020 decennial Census | Census `CenPop2020_Mean_CO.txt` | Point-in-NEON-polygon assignment; nearest-polygon exception log |
| County polygons | 2020 cartographic boundaries, 1:500,000 | Census `cb_2020_us_county_500k.zip` | EPSG:4326 grid cell-center rasterization; area weights; separate simplified display geometries |
| NEON regions | October 2024 | NEON `NEONDomains_2024.zip` | Continental 17 domains for maps, explicit 11 for regression |
| tmean, tdmean, ppt, vpdmax | PRISM 4 km monthly; NC State files created 2025-09-12 where recorded; 2004-10 to 2022-07 seasonal subset | NC State Climate Office THREDDS; source metadata identify Oregon State PRISM | Decode fill values; use calendar date in filename because `months since` metadata are not CF-decodable with the stated calendar; county means, log-case weighted regional means |
| PDSI | gridMET pentad, endpoint metadata dated 2026-10-05, baseline 1980–2019 | University of Idaho / Northwest Knowledge THREDDS | Decode packed unsigned data and scale/offset; monthly arithmetic mean of pentad endpoints; county and regional weighting |
| Soil moisture | Intended NLDAS Noah monthly v2.0; not obtained | NASA GES DISC | None in executed run; no surrogate used |

The climate window contains 180 month-year combinations per variable. PRISM
contributes 720 downloaded grids. Source hashes identify mutable upstream
versions; a future clean download may differ if a provider revises its files.
Saved processed inputs and numerical outputs preserve the executed analysis.

Some NC State mirror fields have blank units and `UNKNOWN` long names (tdmean and
vpdmax). Units are assigned from the PRISM product definitions: degrees Celsius
and hPa respectively, rather than inferred from the mirror's uninformative fields.
The mirror's `max` attributes are not used to validate or rescale grids; actual
decoded grid values are used. Sample global and variable attributes are saved in
`data/climate_metadata.json`. Standardizing covariates makes a constant unit-scale
choice irrelevant to forecasts, but it would affect coefficients in raw units.

## Failed access and quality-control attempts

- Initial sandbox DNS access was blocked; authorized downloads were repeated with
  network escalation.
- PMC's direct supplement URL returned browser-check HTML, not a PDF. The publisher
  supplement endpoint returned HTTP 403. Europe PMC's supplement package succeeded;
  only the published S1 PDF was extracted. Packaged manuscript images were not used.
- NASA's NLDAS HTTPS file returned HTTP 401 Unauthorized. The legacy GrADS metadata
  endpoint returned HTTP 410 Gone. Authentication was not bypassed.
- The gridMET NCSS year request returned HTTP 500. Its OPeNDAP service succeeded,
  but one large response contained spurious constant -15 fields (packed zeroes).
  That cache is renamed `.rejected.nc` and excluded from processing. Smaller
  six-pentad requests, per-grid variance checks, and single-grid retries are used.
- The first PRISM xarray inspection failed on nonstandard `months since` time
  encoding. `decode_times=False` and explicit filename dates resolve the issue.
- The region-threshold check found Mid Atlantic at 987 cases. This is reported and
  the explicitly listed 11 domains are retained for comparability.
- County totals agree with the independent contiguous-state totals in 17 of 18
  years. Two 2010 Iowa cases have CDC county code 19999 (unknown county), which
  cannot be assigned to a NEON region; they are explicitly excluded rather than
  allocated arbitrarily. The exclusion and state-conservation tables are saved.

No institutional entry for the target paper, original author code, archived
author inputs, forecast submissions, plotting script, or digitized graph was
used. Written Table 1 predictor names are scientific specifications; written
Tables 2, S1 and S5 provide numerical comparison references only.
