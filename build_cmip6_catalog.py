#!/usr/bin/env python3
"""Select one historical daily-pr realization/grid/version per CMIP6 model."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re


PATTERN = re.compile(
    r"pr_day_(?P<model>.+?)_historical_(?P<member>r\d+i\d+p\d+f\d+)_"
    r"(?P<grid>[^_]+)_(?P<start>\d{8})-(?P<end>\d{8})\.nc$"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path,
                        default=Path("/g/data/oi10/replicas/CMIP6/CMIP"))
    parser.add_argument("--output", type=Path, default=Path("cmip6_daily_pr_1980_2014.csv"))
    parser.add_argument("--file-list", type=Path,
                        help="newline-delimited paths from find(1)")
    args = parser.parse_args()

    records = []
    paths = ((Path(line) for line in args.file_list.read_text().splitlines())
             if args.file_list else args.root.rglob("pr_day_*_historical_*.nc"))
    for path in paths:
        match = PATTERN.match(path.name)
        if not match:
            continue
        item = match.groupdict()
        if item["end"] < "19800101" or item["start"] > "20141231":
            continue
        item["path"] = str(path)
        item["version"] = path.parent.name
        records.append(item)

    selected = []
    for model in sorted({row["model"] for row in records}):
        rows = [row for row in records if row["model"] == model]
        members = sorted({row["member"] for row in rows},
                         key=lambda x: (x != "r1i1p1f1", int(re.match(r"r(\d+)", x).group(1)), x))
        rows = [row for row in rows if row["member"] == members[0]]
        grids = sorted({row["grid"] for row in rows}, key=lambda x: (x != "gn", x))
        rows = [row for row in rows if row["grid"] == grids[0]]
        versions = sorted({row["version"] for row in rows}, reverse=True)
        rows = [row for row in rows if row["version"] == versions[0]]
        rows.sort(key=lambda row: row["start"])
        if rows[0]["start"] <= "19800101" and rows[-1]["end"] >= "20141230":
            selected.extend(rows)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=(
            "model", "member", "grid", "version", "start", "end", "path"
        ))
        writer.writeheader(); writer.writerows(selected)
    print(f"Wrote {args.output}: {len(set(r['model'] for r in selected))} models, "
          f"{len(selected)} files")


if __name__ == "__main__":
    main()
