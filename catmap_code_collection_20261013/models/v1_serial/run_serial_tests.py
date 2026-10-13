"""Run the complete project unittest suite and save a reviewable text log."""

from __future__ import annotations

import contextlib
import io
import sys
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LOG_FILE = ROOT / "reports" / "SERIAL_V1_TEST_LOG.txt"


def run() -> tuple[bool, int, str]:
    buffer = io.StringIO()
    loader = unittest.TestLoader()
    suite = loader.discover(str(ROOT / "tests"), pattern="test_*.py")
    runner = unittest.TextTestRunner(stream=buffer, verbosity=2)
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        result = runner.run(suite)
    log = buffer.getvalue()
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    LOG_FILE.write_text(log, encoding="utf-8")
    return result.wasSuccessful(), result.testsRun, log


if __name__ == "__main__":
    passed, count, log = run()
    print(f"tests_run={count}; passed={passed}; log={LOG_FILE}")
    print(log)
    sys.exit(0 if passed else 1)
