from pathlib import Path
import hashlib, json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
D,O,F,R=ROOT/"data",ROOT/"origin_exports",ROOT/"figures",ROOT/"reports"
HROOT=Path(r"H:\catmap\catmap\catmap_code_collection_20261013")
HTARGET=HROOT/"analyses"/"BEP_scaling_and_TS_reconstruction"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda:f.read(1024*1024),b""): h.update(c)
    return h.hexdigest()

sources={
"source_strict_energies.txt":HROOT/"catmap_official_run/energies.txt",
"source_matched_O_H_TOF.csv":HROOT/"analyses/TS1_TS2_descriptor_volcano/data/matched_data.csv",
"source_prior_BEP_points.csv":HROOT/"analyses/O_H_to_TS1_TS2_to_TOF/data/BEP_catalyst_values_and_predictions.csv",
"source_prior_BEP_summary.csv":HROOT/"analyses/O_H_to_TS1_TS2_to_TOF/data/BEP_model_summary.csv",
"source_prior_chain_metrics.csv":HROOT/"analyses/O_H_to_TS1_TS2_to_TOF/data/O_H_to_BEP_chain_model_metrics.csv",
"source_prior_chain_predictions.csv":HROOT/"analyses/O_H_to_TS1_TS2_to_TOF/data/O_H_to_BEP_chain_TS_predictions.csv",
"source_carbon_O_H_TS1.csv":HROOT/"analyses/Carbon_containing_descriptor_to_TS1/data/matched_carbon_O_H_TS1_data.csv",
"source_descriptor_inventory.csv":HROOT/"analyses/Carbon_containing_descriptor_to_TS1/data/descriptor_source_inventory.csv",
"source_catmap_grid_50x50.csv":HROOT/"analyses/TS1_TS2_descriptor_volcano/origin_exports/catmap/catmap_TS1_TS2_50x50_xyz.csv",
"source_catmap_VectorMap_400x400.csv":HROOT/"analyses/TS1_TS2_descriptor_volcano/origin_exports/catmap/catmap_VectorMap_interpolated_400x400_xyz.csv",
"source_validated_catalyst_points.csv":HROOT/"analyses/TS1_TS2_descriptor_volcano/origin_exports/catmap/validated_catalyst_points.csv",
"source_grid_summary.json":HROOT/"analyses/TS1_TS2_descriptor_volcano/origin_exports/catmap/grid_summary.json",
"source_grid_setup.mkm":HROOT/"catmap_official_run/TS1_TS2_descriptor_grid/ts1_ts2_descriptor_grid.mkm",
"source_strict_config_tire.mkm":HROOT/"catmap_official_run/results/official_strict_tolerance/tire/alkene_epoxidation.mkm",
"source_convergence_13sites.csv":HROOT/"analyses/microkinetic_mechanism_13sites/data/convergence_13sites.csv",
"source_mechanism_report.md":HROOT/"analyses/microkinetic_mechanism_13sites/reports/MECHANISM_ANALYSIS_REPORT.md",
"source_prior_BEP_report.md":HROOT/"analyses/O_H_to_TS1_TS2_to_TOF/reports/O_H_to_TS1_TS2_TOF_REPORT.md"}
hash_rows=[]; failures=[]
for name,source in sources.items():
    staged=D/name
    if not staged.exists() or not source.exists():
        failures.append(f"Missing staged/source file: {name} / {source}"); continue
    a,b=sha(staged),sha(source)
    hash_rows.append({"staged_name":name,"source_path":str(source),"staged_sha256":a,"current_source_sha256":b,"match":a==b})
    if a!=b: failures.append(f"Source snapshot mismatch: {name}")
pd.DataFrame(hash_rows).to_csv(R/"source_input_sha256.csv",index=False)

energy=pd.read_csv(D/"source_strict_energies.txt",sep="\t")
surfaces=sorted(set(energy.surface_name.dropna()))
tof=pd.read_csv(D/"source_matched_O_H_TOF.csv")
conv=tof.loc[tof.status.eq("CONVERGED")&~tof.surface_name.eq("tizr")]
if len(conv)!=13: failures.append(f"Expected 13 valid catalysts, found {len(conv)}")
if len(surfaces)!=14: failures.append(f"Expected 14 surfaces in energy table incl. Ti-Zr, found {len(surfaces)}")

for p in (O/"BEP1_13sites.csv",O/"BEP2_13sites.csv"):
    if not p.exists(): failures.append(f"Missing {p.name}"); continue
    x=pd.read_csv(p)
    if len(x)!=13: failures.append(f"{p.name} rows={len(x)} not 13")
    if x[["reaction_energy_eV","signed_TS_minus_IS_eV"]].isna().any().any(): failures.append(f"NaN in {p.name}")
    if set(x.surface_name)!=set(conv.surface_name): failures.append(f"Catalyst key mismatch in {p.name}")

co=pd.read_csv(O/"official_BEP_comparison.csv")
if len(co)!=2 or co[["delta_alpha","delta_beta"]].abs().to_numpy().max()>1e-10:
    failures.append("Official CatMAP vs OLS coefficient mismatch")
meta=json.loads((D/"official_BEP_scaler_metadata.json").read_text(encoding="utf-8"))
if meta.get("reaction_model_run_called") is not False: failures.append("Unexpected CatMAP run flag")
mapping=pd.read_csv(O/"BEP_TOF_mapping_summary.csv")
if len(mapping)!=13: failures.append(f"TOF mapping rows={len(mapping)} not 13")
if mapping.nearest_saved_interpolated_log10_TOF.notna().sum()!=mapping.existing_grid_bounds_supported.sum(): failures.append("Grid mapping missingness does not match bounds flag")
for stem in ["NEB1_path_BEP","NEB2_path_BEP","BEP_residuals_by_catalyst","BEP_extended_model_comparison","BEP_TS_reconstruction_parity","BEP_chain_coordinates_on_existing_CatMAP_surface"]:
    for ext in ("png","pdf"):
        if not (F/f"{stem}.{ext}").exists(): failures.append(f"Missing figure {stem}.{ext}")
archive_match = None
if HTARGET.exists():
    marker=HTARGET/"README.md"
    manifest=HTARGET/"reports"/"SHA256SUMS.csv"
    archive_match=(marker.exists() and manifest.exists() and "BEP Scaling and TS Reconstruction" in marker.read_text(encoding="utf-8"))
    if not archive_match: failures.append("Existing H target does not carry this package identity; refuse overwrite")

log=["BEP package validation",f"valid catalyst rows: {len(conv)}",f"energy surfaces incl. unresolved: {len(surfaces)}",
     f"staged input snapshots checked: {len(hash_rows)}",f"snapshot SHA matches: {sum(r['match'] for r in hash_rows)}",
     f"official CatMAP coefficient rows: {len(co)}",f"saved VectorMap coordinate mappings: {int(mapping.nearest_saved_interpolated_log10_TOF.notna().sum())}/13",
     "CatMAP ReactionModel.run(): NOT CALLED","new steady-state/TOF calculation: NOT RUN","raw inputs modified: NO",
     f"H archive target absent or expected package identity confirmed: {not HTARGET.exists() or archive_match}","RESULT: "+("PASS" if not failures else "FAIL")]
if failures: log+= ["FAILURES:"]+["- "+e for e in failures]
(R/"VALIDATION_LOG.txt").write_text("\n".join(log)+"\n",encoding="utf-8")
print("\n".join(log))
if failures: raise SystemExit(1)
