#!/usr/bin/env python3
"""Standalone audit: check the bundled data without the original model."""
from run_model import validate

if __name__ == "__main__":
    validate()
    print("Source workbook and all calculated energies belong to this directory.")
    print("Full NEB image-by-image physical validation remains a separate task.")
