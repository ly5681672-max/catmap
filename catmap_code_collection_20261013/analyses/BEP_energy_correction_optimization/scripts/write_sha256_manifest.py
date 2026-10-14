from hashlib import sha256
from pathlib import Path
import csv

package = Path(__file__).resolve().parents[1]
target = package / "reports" / "SHA256SUMS.csv"
rows = []
for path in sorted(package.rglob("*")):
    if path.is_file() and path != target:
        rows.append({
            "relative_path": path.relative_to(package).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path.read_bytes()).hexdigest(),
        })
with target.open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=["relative_path", "size_bytes", "sha256"])
    writer.writeheader()
    writer.writerows(rows)
print(f"hashed {len(rows)} package files; manifest={target}")
