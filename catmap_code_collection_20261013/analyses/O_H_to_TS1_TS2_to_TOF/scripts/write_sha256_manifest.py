from __future__ import annotations

import csv
import hashlib
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
MANIFEST = PACKAGE / "SHA256SUMS.csv"
SOURCE_FILES = {
    r"H:\catmap\catmap\catmap_code_collection_20261013\catmap_official_run\energies.txt": "strict model energy input",
    r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\empirical_O_H_descriptor_analysis\data\matched_data.csv": "historical O/Ha/Hb + strict baseline TOF source",
    r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\empirical_O_H_descriptor_analysis\data\reference_slab_comparison.csv": "historical vs strict slab-reference comparison",
    r"H:\catmap\catmap\catmap_code_collection_20261013\catmap_official_run\TS1_TS2_descriptor_grid\run_official_grid.py": "prior official CatMAP grid runner",
    r"H:\catmap\catmap\catmap_code_collection_20261013\catmap_official_run\TS1_TS2_descriptor_grid\runs\grid_50x50\catmap_ts1_ts2_grid_50x50.log": "prior 50x50 solver log",
    r"H:\catmap\catmap\catmap_code_collection_20261013\catmap_official_run\TS1_TS2_descriptor_grid\runs\grid_50x50\grid_summary.json": "prior 50x50 grid summary",
    r"H:\catmap\catmap\catmap_code_collection_20261013\catmap_official_run\TS1_TS2_descriptor_grid\runs\grid_50x50\catmap_VectorMap_interpolated_400x400_xyz.csv": "prior official VectorMap interpolation source",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    records = []
    for path in sorted(p for p in PACKAGE.rglob("*") if p.is_file() and p != MANIFEST):
        records.append({"scope": "package", "path": path.relative_to(PACKAGE).as_posix(),
                        "sha256": sha256(path), "bytes": path.stat().st_size,
                        "note": "generated/copy snapshot in this package"})
    for raw_path, note in SOURCE_FILES.items():
        path = Path(raw_path)
        if not path.is_file():
            raise FileNotFoundError(f"Missing source for hash verification: {path}")
        records.append({"scope": "original_source_read_only", "path": str(path),
                        "sha256": sha256(path), "bytes": path.stat().st_size,
                        "note": note})
    with MANIFEST.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["scope", "path", "sha256", "bytes", "note"])
        writer.writeheader()
        writer.writerows(records)
    print(f"manifest={MANIFEST}")
    print(f"package_files={sum(r['scope'] == 'package' for r in records)}")
    print(f"source_files={sum(r['scope'] == 'original_source_read_only' for r in records)}")


if __name__ == "__main__":
    main()
