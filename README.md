# CMIP6 Equatorial-Wave Diagnostics in Python

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22884295.svg)](https://doi.org/10.5281/zenodo.22884295)

A Python workflow for Wheeler–Kiladis precipitation spectra, Kelvin and equatorial-Rossby filtered precipitation, skill-score NetCDF files, individual figures, and multi-model comparison panels. Precipitation-based equatorial-wave skill metrics are converted from Qiang's NCL scripts.

Data are not distributed with this repository.

## Example output

![Five-model equatorial-wave comparison](examples/Panel_EQ-Wave_OBS-ITM_ACCESS-CM2_multi_python.png)

The upper panels compare the observed and ACCESS-CM2 symmetric wavenumber–frequency spectra. The lower panel compares ten precipitation-based equatorial-wave skill metrics for ACCESS-CM2, ACCESS-ESM1-5, HadGEM3-GC31-LL, KACE-1-0-G, and UKESM1-0-LL.


### Bottom-panel metrics

| Number | Label | Diagnostic |
|---:|---|---|
| 01 | Eastward | Agreement with the observed eastward-propagating symmetric precipitation spectrum (positive zonal wavenumbers). |
| 02 | Westward | Agreement with the observed westward-propagating symmetric precipitation spectrum (negative zonal wavenumbers). |
| 03 | EK [MAM] | Equatorial Kelvin-wave precipitation-variance skill for March–May. |
| 04 | EK [JJA] | Equatorial Kelvin-wave precipitation-variance skill for June–August. |
| 05 | EK [SON] | Equatorial Kelvin-wave precipitation-variance skill for September–November. |
| 06 | EK [DJF] | Equatorial Kelvin-wave precipitation-variance skill for December–February. |
| 07 | ER [MAM] | Equatorial Rossby-wave precipitation-variance skill for March–May. |
| 08 | ER [JJA] | Equatorial Rossby-wave precipitation-variance skill for June–August. |
| 09 | ER [SON] | Equatorial Rossby-wave precipitation-variance skill for September–November. |
| 10 | ER [DJF] | Equatorial Rossby-wave precipitation-variance skill for December–February. |

The Eastward and Westward metrics compare selected regions of the normalized wavenumber–frequency spectrum. The EK and ER metrics compare the spatial patterns and amplitudes of seasonally grouped, wave-filtered precipitation variance between each model and the observational reference.

Each metric uses the skill score:

```text
score = (1 + pattern_correlation)^4 / [4 (SDR + 1/SDR)^2]
```

Here, `SDR` is the model spatial standard deviation divided by the observed spatial standard deviation. A score approaches 1 when both spatial pattern and amplitude agree with observations; values approaching 0 indicate poor agreement. The horizontal 0.5 line is a visual reference, not a universal pass/fail threshold.

## Required data

Each model needs global daily precipitation with:

- variable `pr` (also recognizes `PREC`, `precip`, and `PRECT`);
- dimensions/coordinates `time`, `lat`, and `lon`;
- coverage of 1980–2014;
- units `kg m-2 s-1` or `mm/day`;
- one or more consecutive files from the same ensemble and grid.

Two observational reference files are required:

```text
DATA_ROOT/obs/SpaceTime.obs.PRECT.nc
DATA_ROOT/obs/kf_filter_GPCP.nc
```

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On NCI Gadi, `conda/analysis3` provides the dependencies:

```bash
module use /g/data/xp65/public/modules
module load conda/analysis3
export PYTHONNOUSERSITE=1
```

## Run one model from raw precipitation

```bash
python eqwave_pipeline.py all \
  --model ACCESS-CM2 \
  --data-root /path/to/work-data \
  --raw-files /path/to/pr_1980_1999.nc /path/to/pr_2000_2014.nc \
  --output Panel_EQ-Wave_OBS-ACCESS-CM2.png
```

The `data-root` must contain the `obs` directory shown above. Model results are written beneath that root.

Available stages:

```bash
python eqwave_pipeline.py spectrum ...
python eqwave_pipeline.py filter ...
python eqwave_pipeline.py scores ...
python eqwave_pipeline.py plot ...
```

To reuse existing expensive spectrum/filter intermediates:

```bash
python eqwave_pipeline.py all ... --reuse-existing
```

## Discover CMIP6 model files

First make a lightweight filename index:

```bash
find /path/to/CMIP6/CMIP -type f -name 'pr_day_*_historical_*.nc' \
  > cmip6_pr_files.txt
```

Then build the catalogue:

```bash
python build_cmip6_catalog.py \
  --file-list cmip6_pr_files.txt \
  --output cmip6_daily_pr_1980_2014.csv
```

The selector prefers `r1i1p1f1`, the native `gn` grid, the latest version, and complete 1980–2014 coverage.

## Make a multi-model panel

The top panels contain OBS and one selected model. Repeated `--result` arguments add models to the bottom skill-score panel:

```bash
python plot_multi_model.py \
  --obs-spectrum /path/to/obs/SpaceTime.obs.PRECT.nc \
  --top-model ACCESS-CM2 \
  --result ACCESS-CM2 /results/ACCESS-CM2 \
  --result ACCESS-ESM1-5 /results/ACCESS-ESM1-5 \
  --result HadGEM3-GC31-LL /results/HadGEM3-GC31-LL \
  --result KACE-1-0-G /results/KACE-1-0-G \
  --result UKESM1-0-LL /results/UKESM1-0-LL \
  --output multi_model_panel.png
```

Each result directory must contain the `CMIP6/SpaceTime`, `CMIP6/STScore`, and `CMIP6/EQ-WaveScore` subdirectories produced by `eqwave_pipeline.py`.

## PBS on Gadi

Copy and edit the included template:

```bash
cp submit_model.pbs.example submit_model.pbs
```

Replace `PROJECT`, prepare a text file containing one raw precipitation path per line, then submit:

```bash
qsub -v MODEL=ACCESS-CM2,RAW_FILES_FILE=/path/raw_files.txt,OBS_ROOT=/path/data,OUTPUT_ROOT=/scratch/PROJECT/user/ACCESS-CM2 submit_model.pbs
```

## Outputs

```text
CMIP6/SpaceTime/SpaceTime.MODEL.PRECT.nc
CMIP6/Kr-filter/kf_filter_Pre_MODEL.nc
CMIP6/STScore/Skill-Score_MODEL_Pre_ST_1980-2014.nc
CMIP6/EQ-WaveScore/Skill-Score_EQ-Wave_Pre_MODEL.nc
Panel_EQ-Wave_OBS-MODEL.png
```

## Scientific scope

These diagnostics measure tropical space–time spectral power and Kelvin/ER wave skill. They are relevant to MJO evaluation but do not alone measure every aspect of MJO simulation. See `NOTICE.md` for provenance and validation notes.
