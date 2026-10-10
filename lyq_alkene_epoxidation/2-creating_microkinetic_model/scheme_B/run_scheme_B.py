"""Run Scheme B in its own directory; do not touch the baseline A."""
import os, subprocess, sys
from pathlib import Path
d=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(d/'build_scheme_B.py')],cwd=d,check=True)
try:
    from catmap import ReactionModel
except ImportError:
    raise SystemExit('Cannot import CatMAP. Activate the Python environment where CatMAP is installed.')
os.chdir(d)
model=ReactionModel(setup_file='alkene_epoxidation.mkm')
model.output_variables=['coverage','production_rate','rate']
model.run()
print('Scheme B finished; inspect scheme_B.log / scheme_B.pkl')
