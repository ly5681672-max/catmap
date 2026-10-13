"""Run V2 mapping tests and the frozen V1 regression suite; write logs only in V2."""

from __future__ import annotations

import contextlib
import io
import sys
import unittest
from pathlib import Path


V2_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = V2_ROOT.parent
LOG_PATH = V2_ROOT / "reports" / "V2_MAPPING_TEST_LOG.txt"


def run_suite(start_dir: Path, pattern: str, top_level_dir: Path | None = None) -> tuple[bool, int, str]:
    suite = unittest.TestLoader().discover(
        str(start_dir), pattern=pattern,
        top_level_dir=str(top_level_dir) if top_level_dir else None,
    )
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        result = unittest.TextTestRunner(stream=buffer, verbosity=2).run(suite)
    return result.wasSuccessful(), result.testsRun, buffer.getvalue()


def main() -> int:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    v2_ok, v2_count, v2_log = run_suite(V2_ROOT / "tests", "test_v2_source_row_mapping.py")
    v1_ok, v1_count, v1_log = run_suite(PROJECT_ROOT / "tests", "test_*.py")
    output = (
        "V2 source-row mapping tests\n===========================\n"
        f"tests_run={v2_count}; passed={v2_ok}\n{v2_log}\n"
        "V1 frozen regression tests\n===========================\n"
        f"tests_run={v1_count}; passed={v1_ok}\n{v1_log}"
    )
    LOG_PATH.write_text(output, encoding="utf-8")
    print(f"v2_tests={v2_count}; v2_passed={v2_ok}; v1_tests={v1_count}; v1_passed={v1_ok}; log={LOG_PATH}")
    print(output)
    return 0 if v2_ok and v1_ok and v2_count == 8 and v1_count == 35 else 1


if __name__ == "__main__":
    sys.exit(main())
