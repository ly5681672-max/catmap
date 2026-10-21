"""Write reproducible input-provenance and output SHA256 manifests."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[2]
SOURCE_MANIFEST = PROJECT / "analyses" / "NEB_coadsorption_TiCo_reference_audit" / "data" / "tico_source_manifest.json"

INPUTS = [
    ROOT / "candidate_tico.mkm",
    ROOT / "candidate_energies.txt",
    ROOT / "candidate_reaction_energy_barriers.csv",
    ROOT / "candidate_steady_state_coverages.csv",
    ROOT / "candidate_tof_baseline.csv",
    ROOT / "candidate_solver_attempts_baseline.json",
    ROOT / "run_candidate_zero_inventory_strict.py",
    ROOT / "run_conditional_sensitivity.py",
    ROOT / "run_rate_audit.py",
    SOURCE_MANIFEST,
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


missing = [str(p) for p in INPUTS if not p.is_file()]
if missing:
    raise SystemExit("missing input(s):\n" + "\n".join(missing))

input_lines = ["# Input SHA256 (analysis copies and audit provenance)"]
for p in INPUTS:
    input_lines.append(f"{sha256(p)}  {p.relative_to(ROOT) if p.is_relative_to(ROOT) else p.relative_to(PROJECT)}")
source_data = json.loads(SOURCE_MANIFEST.read_text(encoding="utf-8"))
input_lines.append("\n# Original source SHA256 values recorded by the read-only source audit")
input_lines.extend(f"{digest}  {path}  [source audit manifest]" for path, digest in sorted(source_data.get("sources_sha256", {}).items()))
(ROOT / "INPUT_SHA256.txt").write_text("\n".join(input_lines) + "\n", encoding="utf-8")

output_lines = ["# SHA256 for files in this analysis package; manifests are excluded"]
for p in sorted(ROOT.rglob("*")):
    if not p.is_file() or p.name in {"OUTPUT_SHA256.txt", "INPUT_SHA256.txt"}:
        continue
    output_lines.append(f"{sha256(p)}  {p.relative_to(ROOT).as_posix()}")
(ROOT / "OUTPUT_SHA256.txt").write_text("\n".join(output_lines) + "\n", encoding="utf-8")
print(f"input entries={len(INPUTS) + len(source_data.get('sources_sha256', {}))}; output files={len(output_lines) - 1}")
