"""Idempotently split or verify conditional DRC and case-level tables."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
source = ROOT / "TICO_PROVISIONAL_DRC.csv"
with source.open(encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))
summary = [r for r in rows if r.get("drc_central")]
cases = [r for r in rows if not r.get("drc_central")]
cases_path = ROOT / "EXPLICIT_TS_PERTURBATION_CASES.csv"

# The generator may already have emitted the two files separately. In that
# case, verify the existing case table instead of trying to split the summary
# file a second time.
if len(summary) == 4 and len(cases) == 0 and cases_path.exists():
    with cases_path.open(encoding="utf-8-sig", newline="") as f:
        cases = list(csv.DictReader(f))
elif len(summary) == 4 and len(cases) == 8:
    fields = list(dict.fromkeys(k for row in summary for k in row))
    with source.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(summary)
    fields = list(dict.fromkeys(k for row in cases for k in row))
    with cases_path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(cases)
else:
    raise SystemExit(f"unexpected split counts: DRC={len(summary)}, cases={len(cases)}")

if len(summary) != 4 or len(cases) != 8:
    raise SystemExit(f"invalid finalized tables: DRC={len(summary)}, cases={len(cases)}")
print(f"verified DRC={len(summary)}, explicit TS cases={len(cases)}")
