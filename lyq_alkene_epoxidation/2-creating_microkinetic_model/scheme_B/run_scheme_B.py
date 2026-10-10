"""Run Scheme B then audit grid and solve exact catalyst descriptor points."""
import os
import subprocess
import sys
from pathlib import Path
d=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(d/"build_scheme_B.py")],cwd=d,check=True)
try:
    from catmap import ReactionModel
except ImportError:
    raise SystemExit("Activate the Python environment with CatMAP installed.")
os.chdir(d)
model=ReactionModel(setup_file="alkene_epoxidation.mkm")
model.output_variables=["coverage","production_rate","rate"]
model.run()
print("CatMAP finished. Now auditing 400 grid points and solving 14 exact catalysts.")
from analyze_scheme_B import analyze
analyze(model=model,tolerance=.01,min_signal=0)
print("Scheme B results saved in diagnostics/ (UTF-8 CSV).")
