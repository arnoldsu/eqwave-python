#!/usr/bin/env python3
"""Plot OBS/one-model spectra above a multi-model skill comparison."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import xarray as xr

from eqwave_pipeline import LEVELS
import wavenumber_frequency_functions as wf


def files(directory: Path, model: str) -> dict[str, Path]:
    return {
        "spectrum": directory / "CMIP6/SpaceTime" / f"SpaceTime.{model}.PRECT.nc",
        "st": directory / "CMIP6/STScore" / f"Skill-Score_{model}_Pre_ST_1980-2014.nc",
        "eq": directory / "CMIP6/EQ-WaveScore" / f"Skill-Score_EQ-Wave_Pre_{model}.nc",
    }


def scores(paths: dict[str, Path]) -> np.ndarray:
    st, eq = xr.open_dataset(paths["st"]), xr.open_dataset(paths["eq"])
    return np.r_[st.SE.values[:, 0], eq.Score_Kelvin.values[:, 0],
                 eq.Score_ER.values[:, 0]]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs-spectrum", type=Path, required=True)
    parser.add_argument("--top-model", default="ACCESS-CM2")
    parser.add_argument("--result", nargs=2, action="append", metavar=("MODEL", "DIR"),
                        required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    results = {model: files(Path(directory), model)
               for model, directory in args.result}
    missing = [str(path) for item in results.values() for path in item.values()
               if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required outputs:\n" + "\n".join(missing))
    if args.top_model not in results:
        raise ValueError("--top-model must also be supplied with --result")

    obs = xr.open_dataset(args.obs_spectrum).FIG_3_SYM
    model = xr.open_dataset(results[args.top_model]["spectrum"]).FIG_3_SYM
    fig = plt.figure(figsize=(14, 15), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=(1, .07, .72))
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    cmap = plt.get_cmap("RdYlBu_r", len(LEVELS) + 1)
    norm = BoundaryNorm(LEVELS, cmap.N)
    contour = None
    freq_curve, wave_curve = wf.genDispersionCurves(Ahe=[50, 25, 12])
    for i, (axis, data, title) in enumerate(zip(axes, (obs, model),
                                                 ("OBS", args.top_model))):
        contour = axis.contourf(data.wave, data.freq, data, levels=LEVELS,
                                cmap=cmap, norm=norm, extend="both")
        axis.set(xlim=(-15, 15), ylim=(0, .8), xlabel="Zonal Wave Number",
                 ylabel="Frequency (cpd)",
                 title=f"({chr(97+i)}) West     {title} [Symmetric]     East")
        for period in (3, 6, 30):
            axis.axhline(1 / period, color="k", ls="--", lw=.9)
            axis.text(-14.7, 1 / period + .008, f"{period} days", fontsize=9)
        axis.axvline(0, color="k", ls="--", lw=.9)
        for wave_type in (3, 4, 5):
            for depth in range(3):
                axis.plot(wave_curve[wave_type, depth],
                          freq_curve[wave_type, depth], "k", lw=1)
        axis.text(11, .40, "Kelvin"); axis.text(-10.7, .07, "n=1 ER")
        axis.text(-3, .45, "n=1 IG")
    fig.colorbar(contour, cax=fig.add_subplot(grid[1, :]), orientation="horizontal")

    axis = fig.add_subplot(grid[2, :])
    markers = ("P", "o", "s", "D", "^", "X", "v", "*")
    colors = ("purple", "blue", "red", "cyan", "seagreen", "brown", "orange", "black")
    x = np.arange(1, 11)
    for index, (name, paths) in enumerate(results.items()):
        axis.scatter(x, scores(paths), s=85, marker=markers[index % len(markers)],
                     color=colors[index % len(colors)], label=name)
    axis.set(xlim=(.4, 10.6), ylim=(0, 1),
             ylabel="Worse <== Skill Score ==> Better", xticks=x,
             xticklabels=[f"({i:02d})" for i in x], title="(c)")
    axis.axhline(.5, color="k"); axis.axvline(2.5, color="k"); axis.axvline(6.5, color="k")
    axis.grid(axis="y", ls=":"); axis.legend(loc="upper center", ncol=3)
    labels = ("Eastward", "Westward", "EK[MAM]", "EK[JJA]", "EK[SON]",
              "EK[DJF]", "ER[MAM]", "ER[JJA]", "ER[SON]", "ER[DJF]")
    for position, label in zip(x, labels):
        axis.text(position, -.09, label, rotation=55, ha="right", va="top",
                  transform=axis.get_xaxis_transform())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
