"""Read-only, stdlib-only audit of CatMAP BEP OOH six-step net-rate consistency.
Requires the existing grid_61x61.csv and grid_61x61_net_rates.csv.
No CatMAP solver, DFT, file writes, or reruns.
"""
from pathlib import Path
import csv
import json

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data"
N = 61

def load(name):
    with (DATA / name).open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))

def main():
    grid = load("grid_61x61.csv")
    rates = load("grid_61x61_net_rates.csv")
    assert len(grid) == N*N and len(rates) == N*N, "Expected 3721 grid and net-rate rows"
    lookup = {(int(r["grid_i"]), int(r["grid_j"])): r for r in grid}
    assert len(lookup) == N*N, "Duplicate or missing grid coordinate"
    checked = []
    for j in range(N):
        for k in range(N):
            i = k if j % 2 == 0 else N-1-k  # Native run_grid serpentine order
            g = lookup[(i, j)]
            record = rates[j*N+k]
            values = [float(record[f"step_{n}_net_rate"]) for n in range(1, 7)]
            spread = (max(values)-min(values)) / max(max(map(abs,values)), 1e-300)
            reported = float(g["flux_spread_relative"])
            assert abs(spread-reported) < 1e-12, (
                f"Grid/rate row order mismatch at {(i, j)}: {spread} != {reported}")
            checked.append({"i": i, "j": j, "logtof": float(g["log10_TOF_C6H12O"]),
                            "spread": spread, "residual": float(g["residual"]),
                            "status": g["status"]})
    high = [r for r in checked if r["spread"] > 0.1]
    active = [r for r in checked if r["logtof"] >= -30]
    summary = {
        "grid_count": len(checked),
        "converged_count": sum(r["status"] == "CONVERGED" for r in checked),
        "spread_gt_0.01": sum(r["spread"] > 0.01 for r in checked),
        "spread_gt_0.1": len(high),
        "spread_gt_0.5": sum(r["spread"] > 0.5 for r in checked),
        "all_high_spread_below_logTOF_minus30": all(r["logtof"] < -30 for r in high),
        "logTOF_ge_minus30_count": len(active),
        "logTOF_ge_minus30_max_spread": max(r["spread"] for r in active),
        "max_solver_residual": max(r["residual"] for r in checked),
        "most_inconsistent": sorted(checked, key=lambda r: -r["spread"])[:8],
    }
    print(json.dumps(summary, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
