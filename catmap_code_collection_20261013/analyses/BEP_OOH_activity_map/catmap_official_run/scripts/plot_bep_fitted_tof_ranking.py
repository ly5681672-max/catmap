#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read archived CatMAP BEP-fitted TOF data and redraw Figure 04.
Plotting/export only: no DFT, BEP refit, CatMAP or input modifications.
"""
from __future__ import annotations

import csv
import html
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "13site_BEP_fitted_TOF.csv"
FIGURE = ROOT / "figures" / "04_13site_BEP_fitted_TOF_ranking.svg"
RANKING = ROOT / "data" / "13site_BEP_fitted_TOF_ranking.csv"
COLUMNS = (
    "rank", "surface_name", "catalyst", "TOF_C6H12O",
    "log10_TOF_C6H12O", "independent_TS_above_both_endpoints",
    "catmap_endpoint_clamp_applied",
)


def scientific(value: str) -> str:
    mantissa, exponent = f"{float(value):.3e}".split("e")
    return f"{mantissa}e{int(exponent):+03d}"


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 13:
        raise ValueError(f"Expected 13 catalysts, got {len(rows)}")
    for row in rows:
        tof = float(row["TOF_C6H12O"])
        lg = float(row["log10_TOF_C6H12O"])
        if row["status"] != "CONVERGED" or not (tof > 0) or abs(math.log10(tof) - lg) > 1e-8:
            raise ValueError(f"Invalid or inconsistent source row: {row['surface_name']}")
    rows.sort(key=lambda row: float(row["TOF_C6H12O"]), reverse=True)
    with RANKING.open("w", newline="", encoding="utf-8-sig") as out:
        writer = csv.writer(out)
        writer.writerow(COLUMNS)
        for rank, row in enumerate(rows, 1):
            writer.writerow([rank] + [row[name] for name in COLUMNS[1:]])

    width, height, left, right = 1440, 820, 320, 1158
    x_min, x_max, y0, gap, bar_h = -28, 1, 160, 40, 21

    def xpos(v: float) -> float:
        return left + (v - x_min) / (x_max - x_min) * (right - left)

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" '
        'aria-label="CatMAP 13-site BEP-fitted turnover frequency ranking">',
        f'<rect width="{width}" height="{height}" fill="#fff"/>',
        '<text x="72" y="50" font-family="Arial, Helvetica, sans-serif" '
        'font-weight="700" font-size="29" fill="#22252B">CatMAP BEP-fitted TOF ranking</text>',
        '<text x="73" y="82" font-family="Arial, Helvetica, sans-serif" font-size="17" fill="#54606F">'
        '13 sites  ·  Four-parameter OOH-extended BEP (in-sample)  ·  CatMAP 0.3.1  ·  333.15 K</text>',
        '<text x="73" y="117" font-family="Arial, Helvetica, sans-serif" font-size="13" fill="#626D7C">'
        'Ranked by CatMAP-calculated epoxide TOF; only TS1 is replaced by the BEP-fitted value for each site.</text>',
        '<text x="73" y="145" font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#707988">RANK</text>',
        '<text x="137" y="145" font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#707988">CATALYST</text>',
        f'<text x="{right+30}" y="145" font-family="Arial, Helvetica, sans-serif" '
        'font-size="12" fill="#707988">TOF (s&#8315;&#185;)</text>',
    ]
    for tick in (-28, -25, -20, -15, -10, -5, 0, 1):
        px = xpos(tick)
        stroke = "#B3BBC5" if tick == 0 else "#E6E9ED"
        sw = 1.6 if tick == 0 else 1
        svg += [
            f'<line x1="{px:.2f}" y1="153" x2="{px:.2f}" y2="690" stroke="{stroke}" stroke-width="{sw}"/>',
            f'<text x="{px:.2f}" y="714" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" '
            f'font-size="13" fill="#5C6674">{tick}</text>',
        ]
    for idx, row in enumerate(rows):
        y = y0 + idx * gap
        co = row["surface_name"] == "tico"
        limited = row["catmap_endpoint_clamp_applied"].lower() == "true"
        px = xpos(float(row["log10_TOF_C6H12O"]))
        name = html.escape(row["catalyst"]) + ("†" if limited else "")
        color = "#9E4266" if co else "#728BA7"
        svg += [
            f'<rect x="60" y="{y-14}" width="1290" height="35" fill="{"#FFFFFF" if idx % 2 else "#F8F9FB"}"/>',
            f'<text x="80" y="{y+8}" text-anchor="end" font-family="Arial, Helvetica, sans-serif" '
            f'font-size="17" fill="#7D8795">{idx+1:02d}</text>',
            f'<text x="137" y="{y+8}" font-family="Arial, Helvetica, sans-serif" '
            f'font-weight="{700 if co else 500}" font-size="18" fill="#29313C">{name}</text>',
            f'<rect x="{left}" y="{y-9}" width="{max(0,px-left):.2f}" height="{bar_h}" '
            f'rx="3" fill="{color}"/>',
            f'<text x="{right+30}" y="{y+8}" font-family="Consolas, Menlo, monospace" '
            f'font-size="16" fill="#26303C">{scientific(row["TOF_C6H12O"])}</text>',
        ]
    svg += [
        f'<line x1="{left}" y1="698" x2="{right}" y2="698" stroke="#9EA7B4" stroke-width="1.2"/>',
        f'<text x="{(left+right)/2}" y="745" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" '
        'font-size="16" fill="#27313D">log&#8321;&#8320;[TOF / s&#8315;&#185;]    '
        '(bar starts at display limit: &#8722;28)</text>',
        '<text x="73" y="777" font-family="Arial, Helvetica, sans-serif" font-size="13" fill="#5C6674">'
        '† CatMAP endpoint-limited TS1; no independent barrier above the initial state.</text>',
        '<text x="73" y="801" font-family="Arial, Helvetica, sans-serif" font-size="13" fill="#5C6674">'
        'In-sample fitting; NOT leave-one-out prediction. Source: data/13site_BEP_fitted_TOF.csv.</text>',
        '</svg>',
    ]
    FIGURE.write_text("\n".join(svg) + "\n", encoding="utf-8")
    print(f"Saved: {FIGURE}")
    print(f"Saved: {RANKING}")
    for idx, row in enumerate(rows, 1):
        print(f"{idx:2d}  {row['catalyst']:<8s}  {float(row['TOF_C6H12O']):.6e}")


if __name__ == "__main__":
    main()
