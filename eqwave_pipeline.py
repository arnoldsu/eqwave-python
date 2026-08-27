#!/usr/bin/env python3
"""ACCESS-CM2 equatorial-wave diagnostics: preprocess, score, and plot."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date
from scipy.stats import pearsonr

from kf_filter import kf_filter
import wavenumber_frequency_functions as wf


SEASONS = ("MAM", "JJA", "SON", "DJF")
LEVELS = np.array([.6, .7, .8, .9, 1., 1.1, 1.15, 1.2, 1.25,
                   1.3, 1.35, 1.4, 1.45, 1.5, 1.6])


def _lat_slice(coord: xr.DataArray, south: float, north: float) -> slice:
    return slice(south, north) if coord[0] < coord[-1] else slice(north, south)


def _precip(paths: list[Path], south: float, north: float) -> xr.DataArray:
    ds = (xr.open_dataset(paths[0]) if len(paths) == 1 else
          xr.open_mfdataset(paths, combine="by_coords"))
    candidates = [v for v in ("pr", "PREC", "precip", "PRECT") if v in ds]
    if not candidates:
        raise ValueError(f"No precipitation variable found in {paths}")
    da = ds[candidates[0]].sel(lat=_lat_slice(ds.lat, south, north))
    da = da.where((da.time.dt.year >= 1980) & (da.time.dt.year <= 2014), drop=True)
    units = da.attrs.get("units", "")
    if units in {"kg m-2 s-1", "kg/m2/s", "kg m**-2 s**-1"}:
        da = da * 86400.0
        da.attrs["units"] = "mm/day"
    return da


def make_spectrum(raw: list[Path], output: Path) -> None:
    """Reproduce NCL wkSpaceTime settings (96-day windows, 30-day stride)."""
    pr = _precip(raw, -15, 15).load()
    power = wf.spacetime_power(
        pr, segsize=96, noverlap=66, spd=1, latitude_bounds=(-15, 15),
        dosymmetries=True, rmvLowFrq=True,
    )
    power.loc[{"frequency": 0}] = np.nan
    mean_power = power.mean("component")
    kernel = np.array([[0., 1., 0.], [1., 4., 1.], [0., 1., 0.]]) / 8.
    background = wf.smooth_wavefreq(
        mean_power, kern=kernel, nsmooth=50, freq_name="frequency"
    )

    def region(da: xr.DataArray) -> xr.DataArray:
        da = da.sel(frequency=slice(0, .5), wavenumber=slice(-15, 15))
        return da.transpose("frequency", "wavenumber").rename(
            frequency="freq", wavenumber="wave"
        )

    sym = region(power.sel(component="symmetric", drop=True))
    asym = region(power.sel(component="antisymmetric", drop=True))
    back = region(background)
    out = xr.Dataset({
        "FIG_1_SYM": np.log10(sym),
        "FIG_1_ASYM": np.log10(asym),
        "FIG_3_BACK": np.log10(back),
        "FIG_3_SYM": sym / back,
        "FIG_3_ASYM": asym / back,
    }).astype("float32")
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(output)
    print(f"Wrote {output}")


def make_filtered(raw: list[Path], output: Path) -> None:
    pr = _precip(raw, -20, 20).load()
    dims = ("time", "lat", "lon")
    pr = pr.transpose(*dims)
    result = {}
    settings = {
        "Kelvin": (2.5, 30., 1, 14),
        "ER": (6., 90., -8, -1),
        "MRG": (3., 10., -10, -1),
    }
    for wave, (tmin, tmax, kmin, kmax) in settings.items():
        vals = np.empty(pr.shape, dtype="float32")
        for j in range(pr.sizes["lat"]):
            print(f"{wave}: latitude {j + 1}/{pr.sizes['lat']}", flush=True)
            vals[:, j, :] = kf_filter(
                pr.isel(lat=j), 1, tmin, tmax, kmin, kmax, 8., 90., wave
            ).values
        result[wave] = (dims, vals)
    output.parent.mkdir(parents=True, exist_ok=True)
    xr.Dataset(result, coords={d: pr[d] for d in dims}).to_netcdf(output)
    print(f"Wrote {output}")


def _skill(obs: np.ndarray, model: np.ndarray, weights=None) -> float:
    valid = np.isfinite(obs) & np.isfinite(model)
    a, b = obs[valid], model[valid]
    if weights is None:
        pc = pearsonr(a, b).statistic
    else:
        w = weights[valid]
        am, bm = np.average(a, weights=w), np.average(b, weights=w)
        pc = np.average((a-am)*(b-bm), weights=w) / np.sqrt(
            np.average((a-am)**2, weights=w) * np.average((b-bm)**2, weights=w)
        )
    sdr = np.std(b, ddof=1) / np.std(a, ddof=1)
    return float((1 + pc)**4 / (4 * (sdr + 1/sdr)**2))


def score_spectrum(obs_file: Path, model_file: Path, output: Path, model: str) -> None:
    obs, mod = xr.open_dataset(obs_file).FIG_3_SYM, xr.open_dataset(model_file).FIG_3_SYM
    east_o = obs.sel(freq=slice(0, .4375), wave=slice(0, 15))
    west_o = obs.sel(freq=slice(0, .1875), wave=slice(-15, 0))
    east_m = mod.sel(freq=slice(0, .4375), wave=slice(0, 15))
    west_m = mod.sel(freq=slice(0, .1875), wave=slice(-15, 0))
    scores = np.array([_skill(east_o.values, east_m.values),
                       _skill(west_o.values, west_m.values)])[:, None]
    ds = xr.Dataset({"SE": (("SpaceTime", "Model_Name"), scores)},
                    coords={"SpaceTime": ["East", "West"], "Model_Name": [model]})
    output.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(output)
    print(f"Wrote {output}")


def _season_mask(time: xr.DataArray, season: str) -> xr.DataArray:
    months = {"MAM": (3, 4, 5), "JJA": (6, 7, 8),
              "SON": (9, 10, 11), "DJF": (12, 1, 2)}[season]
    return time.dt.month.isin(months)


def _seasonal_variances(path: Path, variable: str, chunk: int = 128):
    """Compute four temporal variance maps without loading a large cube."""
    with Dataset(path) as nc:
        var, time = nc.variables[variable], nc.variables["time"]
        lat = np.asarray(nc.variables["lat"][:])
        lon = np.asarray(nc.variables["lon"][:])
        count = np.zeros(4, dtype=np.int64)
        total = np.zeros((4, lat.size, lon.size), dtype=np.float64)
        total2 = np.zeros_like(total)
        season_index = {3: 0, 4: 0, 5: 0, 6: 1, 7: 1, 8: 1,
                        9: 2, 10: 2, 11: 2, 12: 3, 1: 3, 2: 3}
        for start in range(0, time.size, chunk):
            stop = min(start + chunk, time.size)
            dates = num2date(time[start:stop], time.units,
                             getattr(time, "calendar", "standard"))
            data = np.ma.filled(var[start:stop], np.nan).astype(np.float64)
            for season in range(4):
                idx = np.array([season_index[d.month] == season for d in dates])
                if idx.any():
                    block = data[idx]
                    count[season] += block.shape[0]
                    total[season] += np.nansum(block, axis=0)
                    total2[season] += np.nansum(block * block, axis=0)
        variance = total2 / count[:, None, None] - (total / count[:, None, None])**2
    return variance, lat, lon


def score_waves(obs_file: Path, model_file: Path, output: Path, model: str) -> None:
    scores = {"Kelvin": [], "ER": []}
    for wave in scores:
        ovar, olat, olon = _seasonal_variances(obs_file, wave)
        mvar, mlat, mlon = _seasonal_variances(model_file, wave)
        for season in range(4):
            ov = xr.DataArray(ovar[season], dims=("lat", "lon"), coords={"lat": olat, "lon": olon})
            mv = xr.DataArray(mvar[season], dims=("lat", "lon"), coords={"lat": mlat, "lon": mlon})
            mv = mv.interp(lat=olat, lon=olon)
            w = np.broadcast_to(np.cos(np.deg2rad(olat))[:, None], ov.shape)
            scores[wave].append(_skill(ov.values, mv.values, w))
    ds = xr.Dataset({
        "Score_Kelvin": (("season", "model"), np.array(scores["Kelvin"])[:, None]),
        "Score_ER": (("season", "model"), np.array(scores["ER"])[:, None]),
    }, coords={"season": list(SEASONS), "model": [model]})
    output.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(output)
    print(f"Wrote {output}")


def plot_panel(obs_spectrum: Path, model_spectrum: Path, st_score: Path,
               wave_score: Path, output: Path, model: str) -> None:
    spectra = [xr.open_dataset(obs_spectrum).FIG_3_SYM,
               xr.open_dataset(model_spectrum).FIG_3_SYM]
    fig = plt.figure(figsize=(14, 15), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=(1, .07, .72))
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    cmap = plt.get_cmap("RdYlBu_r", len(LEVELS) + 1)
    norm = BoundaryNorm(LEVELS, cmap.N)
    contour = None
    for i, (ax, da, title) in enumerate(zip(axes, spectra, ("OBS", model))):
        contour = ax.contourf(da.wave, da.freq, da, levels=LEVELS,
                              cmap=cmap, norm=norm, extend="both")
        ax.set(xlim=(-15, 15), ylim=(0, .8), xlabel="Zonal Wave Number",
               ylabel="Frequency (cpd)", title=f"({chr(97+i)}) West     {title} [Symmetric]     East")
        for period in (3, 6, 30):
            ax.axhline(1/period, color="k", ls="--", lw=.9)
            ax.text(-14.7, 1/period+.008, f"{period} days", fontsize=9)
        ax.axvline(0, color="k", ls="--", lw=.9)
        freq, wave = wf.genDispersionCurves(Ahe=[50, 25, 12])
        for wi in (3, 4, 5):
            for hi in range(3):
                ax.plot(wave[wi, hi], freq[wi, hi], "k", lw=1)
        ax.text(11, .40, "Kelvin"); ax.text(-10.7, .07, "n=1 ER")
        ax.text(-3, .45, "n=1 IG")
    fig.colorbar(contour, cax=fig.add_subplot(gs[1, :]), orientation="horizontal")

    ax = fig.add_subplot(gs[2, :])
    st, eq = xr.open_dataset(st_score), xr.open_dataset(wave_score)
    vals = np.r_[st.SE.values[:, 0], eq.Score_Kelvin.values[:, 0], eq.Score_ER.values[:, 0]]
    x = np.arange(1, 11)
    ax.scatter(x, vals, s=100, marker="P", color="purple", label=model)
    ax.set(xlim=(.4, 10.6), ylim=(0, 1), ylabel="Worse <== Skill Score ==> Better",
           xticks=x, xticklabels=[f"({i:02d})" for i in x], title="(c)")
    ax.axhline(.5, color="k"); ax.axvline(2.5, color="k"); ax.axvline(6.5, color="k")
    ax.grid(axis="y", ls=":"); ax.legend(loc="lower center")
    labels = ["Eastward", "Westward", "EK[MAM]", "EK[JJA]", "EK[SON]",
              "EK[DJF]", "ER[MAM]", "ER[JJA]", "ER[SON]", "ER[DJF]"]
    for xi, label in zip(x, labels):
        ax.text(xi, -.09, label, rotation=55, ha="right", va="top",
                transform=ax.get_xaxis_transform())
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output}")


def paths(root: Path, model: str, raw_files: list[Path] | None = None) -> dict[str, object]:
    raw = raw_files or [root / "CMIP6/year_daily" /
                        f"pr_day_{model}_historical_r1i1p1f1_gn_19800101-20141231.nc"]
    return {
        "raw": [path.resolve() for path in raw],
        "obs_spec": root / "obs/SpaceTime.obs.PRECT.nc",
        "model_spec": root / "CMIP6/SpaceTime" / f"SpaceTime.{model}.PRECT.nc",
        "obs_kf": root / "obs/kf_filter_GPCP.nc",
        "model_kf": root / "CMIP6/Kr-filter" / f"kf_filter_Pre_{model}.nc",
        "st": root / "CMIP6/STScore" / f"Skill-Score_{model}_Pre_ST_1980-2014.nc",
        "eq": root / "CMIP6/EQ-WaveScore" / f"Skill-Score_EQ-Wave_Pre_{model}.nc",
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=("all", "spectrum", "filter", "scores", "plot"))
    p.add_argument("--model", default="ACCESS-CM2")
    p.add_argument("--data-root", type=Path, default=Path("../data"))
    p.add_argument("--raw-files", type=Path, nargs="+",
                   help="one or more consecutive historical daily-pr files")
    p.add_argument("--output", type=Path)
    p.add_argument("--reuse-existing", action="store_true",
                   help="for 'all', retain existing expensive intermediate files")
    a = p.parse_args(); q = paths(a.data_root.resolve(), a.model, a.raw_files)
    if a.stage in ("all", "spectrum") and not (a.reuse_existing and q["model_spec"].exists()):
        make_spectrum(q["raw"], q["model_spec"])
    if a.stage in ("all", "filter") and not (a.reuse_existing and q["model_kf"].exists()):
        make_filtered(q["raw"], q["model_kf"])
    if a.stage in ("all", "scores"):
        score_spectrum(q["obs_spec"], q["model_spec"], q["st"], a.model)
        score_waves(q["obs_kf"], q["model_kf"], q["eq"], a.model)
    if a.stage in ("all", "plot"):
        out = a.output or Path(f"Panel_EQ-Wave_OBS-ITM_{a.model}_python.png")
        plot_panel(q["obs_spec"], q["model_spec"], q["st"], q["eq"], out, a.model)


if __name__ == "__main__":
    main()
