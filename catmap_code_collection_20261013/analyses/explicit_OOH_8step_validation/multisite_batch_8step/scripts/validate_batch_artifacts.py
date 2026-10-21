"""Deterministic integrity checks and SHA256 manifest for the provisional batch."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

BATCH = Path(__file__).resolve().parents[1]
MANIFEST = BATCH / "SHA256SUMS.csv"


def read_rows(name: str) -> list[dict[str, str]]:
    with (BATCH / name).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def validate() -> None:
    status = read_rows("MULTISITE_CALCULATION_STATUS.csv")
    status_by_site = {row["surface_name"]: row for row in status}
    outputs = read_rows("MULTISITE_TOF_COVERAGE.csv")
    barriers = read_rows("MULTISITE_REACTION_BARRIERS.csv")
    gates = read_rows("NUMERICAL_QUALITY_GATE.csv")
    drc = read_rows("MULTISITE_DRC_PROVISIONAL.csv")
    assert len(status) == 14 and len({r["surface_name"] for r in status}) == 14
    out_names = {r["surface_name"] for r in outputs}
    assert len(outputs) == 10
    assert out_names == {r["surface_name"] for r in status if r["calculation_status"] in {"COMPUTED_PROVISIONAL", "ALREADY_COMPUTED_PROVISIONAL"}}
    assert len(barriers) == 28 and len(gates) == 14 and len(drc) == 20
    for row in outputs:
        assert row["calculation_status"] in {"COMPUTED_PROVISIONAL", "ALREADY_COMPUTED_PROVISIONAL"}
        assert status_by_site[row["surface_name"]]["direction_semantics"] == "REVERSIBLE_NATIVE_CATMAP"
        assert row["physical_gate"].startswith("PROVISIONAL_")
        assert row["TOF_C6H12O_s-1"] and float(row["TOF_C6H12O_s-1"]) > 0
        for step in range(1, 10):
            assert row[f"forward_flux_step_{step}"]
            assert row[f"reverse_flux_step_{step}"]
            assert row[f"net_flux_step_{step}"]
        assert float(row["reverse_flux_step_2"]) > 0
    signed = {(r["surface_name"], r["path"]): float(r["signed_TS_minus_IS_eV"]) for r in barriers}
    for key in (("tiw", "NEB1"), ("tita", "NEB1"), ("timo", "NEB2"), ("tiv", "NEB2"), ("ticr", "NEB2")):
        assert signed[key] < 0, (key, signed[key])
    assert all(r["kind"].startswith(("CONDITIONAL_", "HISTORICAL_CONDITIONAL_")) for r in drc)
    assert all(r["physical_gate"].startswith("PROVISIONAL_") for r in gates if r["surface_name"] not in {"tico", "tife"})


def write_manifest() -> int:
    files = sorted(p for p in BATCH.rglob("*") if p.is_file() and p != MANIFEST and "__pycache__" not in p.parts and p.suffix != ".pyc")
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["relative_path", "sha256", "size_bytes"])
        for path in files:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            writer.writerow([path.relative_to(BATCH).as_posix(), digest, path.stat().st_size])
    return len(files)


if __name__ == "__main__":
    validate()
    print(f"VALIDATION_PASS; files_hashed={write_manifest()}; manifest={MANIFEST}")
