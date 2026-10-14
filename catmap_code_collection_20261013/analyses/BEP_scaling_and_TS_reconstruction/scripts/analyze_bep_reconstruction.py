from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

ROOT = Path(__file__).resolve().parents[1]
D, O, R = ROOT / "data", ROOT / "origin_exports", ROOT / "reports"
OUT = ROOT / "figures"
O.mkdir(exist_ok=True, parents=True)
R.mkdir(exist_ok=True, parents=True)
ENERGY = pd.read_csv(D / "source_strict_energies.txt", sep="\t")
GAS = ENERGY.loc[ENERGY.surface_name.isna()].set_index("species_name").formation_energy
TOF = pd.read_csv(D / "source_matched_O_H_TOF.csv")
CARBON = pd.read_csv(D / "source_carbon_O_H_TS1.csv")
PRIOR1 = pd.read_csv(D / "source_prior_BEP_points.csv")
CHAIN = pd.read_csv(D / "source_prior_chain_predictions.csv")
GRID = pd.read_csv(D / "source_catmap_VectorMap_400x400.csv")

VALID = TOF.loc[TOF.status.eq("CONVERGED"), ["surface_name", "catalyst", "TOF", "log10_TOF"]].copy()
VALID = VALID.loc[~VALID.surface_name.eq("tizr")].reset_index(drop=True)
ENERGY = ENERGY.loc[ENERGY.surface_name.isin(VALID.surface_name)].copy()
PIV = ENERGY.pivot(index="surface_name", columns="species_name", values="formation_energy")
BASE = VALID.merge(PIV.reset_index(), on="surface_name", validate="one_to_one")
BASE["NEB1_dG_rxn_eV"] = BASE.EpoxideHHO - BASE.Rpair
BASE["NEB1_signed_TS_minus_IS_eV"] = BASE.TS1 - BASE.Rpair
BASE["NEB2_dG_rxn_eV"] = BASE.OHHa - BASE.HHO
BASE["NEB2_signed_TS_minus_IS_eV"] = BASE.TS2 - BASE.HHO
BASE["overall_coadsorption_dG1_eV"] = BASE.Rpair - GAS["H2O2"] - GAS["C6H12"]
BASE.to_csv(O / "strict_energy_bep_13sites.csv", index=False)

def fit_line(x, y):
    x = np.asarray(x, float).reshape(-1, 1)
    y = np.asarray(y, float)
    model = LinearRegression().fit(x, y)
    pred = model.predict(x)
    cv = cross_val_predict(LinearRegression(), x, y, cv=LeaveOneOut())
    return {"alpha": float(model.coef_[0]), "beta": float(model.intercept_),
            "R2": float(r2_score(y, pred)), "RMSE_eV": float(mean_squared_error(y, pred)**0.5),
            "MAE_eV": float(mean_absolute_error(y, pred)),
            "LOOCV_RMSE_eV": float(mean_squared_error(y, cv)**0.5),
            "LOOCV_R2": float(r2_score(y, cv)), "fitted": pred, "loocv": cv}

def pair_table(path, initial, ts, final, num):
    x = BASE[["surface_name", "catalyst", initial, ts, final, "TOF", "log10_TOF"]].copy()
    x["path"] = path
    x["reaction_energy_eV"] = x[final] - x[initial]
    x["signed_TS_minus_IS_eV"] = x[ts] - x[initial]
    fit = fit_line(x.reaction_energy_eV, x.signed_TS_minus_IS_eV)
    x["BEP_fitted_in_sample_eV"] = fit["fitted"]
    x["BEP_residual_eV"] = x.signed_TS_minus_IS_eV - fit["fitted"]
    x["BEP_predicted_LOOCV_eV"] = fit["loocv"]
    x["LOOCV_residual_eV"] = x.signed_TS_minus_IS_eV - fit["loocv"]
    x.to_csv(O / f"BEP{num}_13sites.csv", index=False)
    return x, fit

B1, F1 = pair_table("NEB1_path2", "Rpair", "TS1", "EpoxideHHO", 1)
B2, F2 = pair_table("NEB2_path4", "HHO", "TS2", "OHHa", 2)
official = pd.read_csv(D / "official_BEP_scaler_coefficients.csv")
comp = []
for name, fit in [("NEB1_path2", F1), ("NEB2_path4", F2)]:
    row = official.loc[official.path.eq(name)]
    if row.empty:
        raise ValueError(f"No official CatMAP BEP coefficient row for {name}")
    cat_alpha, cat_beta = float(row.iloc[0].official_alpha), float(row.iloc[0].official_beta_eV)
    comp.append({"path": name, "numpy_alpha": fit["alpha"], "numpy_beta": fit["beta"],
                 "CatMAP_alpha": cat_alpha, "CatMAP_beta": cat_beta,
                 "delta_alpha": fit["alpha"]-cat_alpha,
                 "delta_beta": fit["beta"]-cat_beta, "R2": fit["R2"],
                 "RMSE_eV": fit["RMSE_eV"], "LOOCV_RMSE_eV": fit["LOOCV_RMSE_eV"],
                 "LOOCV_R2": fit["LOOCV_R2"]})
pd.DataFrame(comp).to_csv(O / "official_BEP_comparison.csv", index=False)

# Retain negative signed barriers; extend the conventional BEP by one candidate term at a time.
CARBON = CARBON.merge(VALID[["surface_name"]], on="surface_name", validate="one_to_one")
features = {
    "overall_coadsorption_dG1_strict": BASE.set_index("surface_name").overall_coadsorption_dG1_eV,
    "HHO_strict_initial_state_overlap": BASE.set_index("surface_name").HHO,
    "OHHa_strict_final_state_overlap": BASE.set_index("surface_name").OHHa,
    "C6H12_star_legacy": CARBON.set_index("surface_name").C6H12_star_legacy,
    "C6H12_adsorption_legacy": CARBON.set_index("surface_name").C6H12_adsorption_energy_legacy,
    "C6H12O_star_legacy_postreaction": CARBON.set_index("surface_name").C6H12O_star_legacy,
    "O_star_legacy": CARBON.set_index("surface_name").O_star_legacy,
    "Ha_star_legacy": CARBON.set_index("surface_name").Ha_star_legacy,
    "Hb_star_legacy": CARBON.set_index("surface_name").Hb_star_legacy,
    "OOH_star_legacy": CARBON.set_index("surface_name").OOH_star_legacy,
    "OH_star_legacy": CARBON.set_index("surface_name").OH_star_legacy,
}
refs = {"overall_coadsorption_dG1_strict":"strict six-step gas and Rpair energy; not pure interaction energy",
        "HHO_strict_initial_state_overlap":"strict BEP2 initial-state energy; overlaps the reaction-energy coordinate",
        "OHHa_strict_final_state_overlap":"strict BEP2 final-state energy; overlaps the reaction-energy coordinate",
        "C6H12_star_legacy":"historical single-adsorbate reference; not closed against strict slab",
        "C6H12_adsorption_legacy":"historical gas-to-adsorbate difference; reference conditional",
        "C6H12O_star_legacy_postreaction":"historical product-side adsorbate; post-reaction diagnostic",
        "O_star_legacy":"historical single-adsorbate reference; slab closure not validated",
        "Ha_star_legacy":"historical single-adsorbate reference; slab closure not validated",
        "Hb_star_legacy":"historical single-adsorbate reference; slab closure not validated",
        "OOH_star_legacy":"historical descriptor; not an independent strict-network OOH state",
        "OH_star_legacy":"historical single-adsorbate reference; slab closure not validated"}

def extended(path, data, fit, ycol):
    rows=[]
    for name, vals in features.items():
        m=data.set_index("surface_name")[["reaction_energy_eV", ycol]].join(vals.rename("X"), how="inner").dropna()
        if len(m)<5: continue
        X=np.column_stack([m.reaction_energy_eV, m.X])
        mod=LinearRegression().fit(X,m[ycol]); pred=mod.predict(X)
        cv=cross_val_predict(LinearRegression(),X,m[ycol],cv=LeaveOneOut())
        rows.append({"path":path,"added_descriptor":name,"N":len(m),"alpha_dG":mod.coef_[0],"gamma_descriptor":mod.coef_[1],
                     "beta":mod.intercept_,"R2":r2_score(m[ycol],pred),"RMSE_eV":mean_squared_error(m[ycol],pred)**.5,
                     "MAE_eV":mean_absolute_error(m[ycol],pred),"LOOCV_RMSE_eV":mean_squared_error(m[ycol],cv)**.5,
                     "LOOCV_R2":r2_score(m[ycol],cv),"descriptor_reference_note":refs[name]})
    return rows

ext = extended("NEB1_path2",B1,F1,"signed_TS_minus_IS_eV") + extended("NEB2_path4",B2,F2,"signed_TS_minus_IS_eV")
pd.DataFrame(ext).sort_values(["path","RMSE_eV"]).to_csv(O / "BEP_extended_model_ranking.csv",index=False)

correlations=[]
for path, data, ycol in [("NEB1_path2",B1,"BEP_residual_eV"),("NEB2_path4",B2,"BEP_residual_eV")]:
    data_idx=data.set_index("surface_name")
    for fname, values in features.items():
        joined=data_idx[ycol].to_frame().join(values.rename("X"),how="inner").dropna()
        correlations.append({"path":path,"candidate_feature":fname,"N":len(joined),
                             "pearson_r":joined[ycol].corr(joined.X,method="pearson"),
                             "spearman_r":joined[ycol].corr(joined.X,method="spearman"),
                             "reference_note":refs[fname]})
pd.DataFrame(correlations).to_csv(O/"BEP_residual_feature_correlations.csv",index=False)

res = pd.concat([B1.assign(path="NEB1_path2"),B2.assign(path="NEB2_path4")],ignore_index=True)
res["abs_residual_eV"] = res.BEP_residual_eV.abs()
res.to_csv(O / "BEP_residual_analysis.csv",index=False)

# Separate true-state BEP reconstruction from descriptor-chain reconstruction.
recon=[]
for _, row in B1.iterrows():
    recon.append({"surface_name":row.surface_name,"catalyst":row.catalyst,"path":"NEB1_path2",
                  "actual_TS_eV":row.TS1,"actual_IS_eV":row.Rpair,"actual_dG_rxn_eV":row.reaction_energy_eV,
                  "BEP_in_sample_TS_reconstruction_eV":row.Rpair+row.BEP_fitted_in_sample_eV,
                  "BEP_LOOCV_TS_reconstruction_eV":row.Rpair+row.BEP_predicted_LOOCV_eV,
                  "actual_signed_TS_minus_IS_eV":row.signed_TS_minus_IS_eV})
for _, row in B2.iterrows():
    recon.append({"surface_name":row.surface_name,"catalyst":row.catalyst,"path":"NEB2_path4",
                  "actual_TS_eV":row.TS2,"actual_IS_eV":row.HHO,"actual_dG_rxn_eV":row.reaction_energy_eV,
                  "BEP_in_sample_TS_reconstruction_eV":row.HHO+row.BEP_fitted_in_sample_eV,
                  "BEP_LOOCV_TS_reconstruction_eV":row.HHO+row.BEP_predicted_LOOCV_eV,
                  "actual_signed_TS_minus_IS_eV":row.signed_TS_minus_IS_eV})
pd.DataFrame(recon).to_csv(O / "BEP_TS1_TS2_reconstruction.csv",index=False)

# Reuse the pre-existing O/H chain predictions (NEB1: O; NEB2: O/H entries) and map to saved VectorMap only.
# Existing model metrics identify Ha for NEB1 and Hb for NEB2 as best LOOCV chain predictors.
best = CHAIN.loc[((CHAIN.path.eq("NEB1_path2")) & CHAIN.descriptor_set.eq("Ha")) |
                 ((CHAIN.path.eq("NEB2_path4")) & CHAIN.descriptor_set.eq("Hb"))].copy()
grid_cols = list(GRID.columns)
xcol = next(c for c in grid_cols if "TS1" in c and "eV" in c)
ycol = next(c for c in grid_cols if "TS2" in c and "eV" in c)
zcol = next(c for c in grid_cols if "log10_TOF" in c)
gxy=GRID[[xcol,ycol,zcol]].dropna().to_numpy(float)
maps=[]
for surf in VALID.surface_name:
    q=best.loc[best.surface_name.eq(surf)]
    if len(q)!=2: continue
    r1=q.loc[q.path.eq("NEB1_path2")].iloc[0]; r2=q.loc[q.path.eq("NEB2_path4")].iloc[0]
    x=float(r1.predicted_TS_LOOCV_eV); y=float(r2.predicted_TS_LOOCV_eV)
    inrange=(gxy[:,0].min()<=x<=gxy[:,0].max() and gxy[:,1].min()<=y<=gxy[:,1].max())
    idx=int(np.argmin((gxy[:,0]-x)**2+(gxy[:,1]-y)**2))
    maps.append({"surface_name":surf,"catalyst":r1.catalyst,"NEB1_chain_predictor":"Ha","NEB2_chain_predictor":"Hb",
                 "TS1_LOOCV_prediction_eV":x,"TS2_LOOCV_prediction_eV":y,"existing_grid_bounds_supported":inrange,
                 "nearest_saved_VectorMap_TS1_eV":gxy[idx,0] if inrange else np.nan,
                 "nearest_saved_VectorMap_TS2_eV":gxy[idx,1] if inrange else np.nan,
                 "nearest_saved_interpolated_log10_TOF":gxy[idx,2] if inrange else np.nan,
                 "mapping_method":"nearest node in previously saved 400x400 CatMAP VectorMap interpolation; not a new solve"})
pd.DataFrame(maps).to_csv(O / "BEP_TOF_mapping_summary.csv",index=False)

# Detailed, reproducible audit of signed values and exceptional negative barriers.
neg1=B1.loc[B1.signed_TS_minus_IS_eV<0,"catalyst"].tolist()
neg2=B2.loc[B2.signed_TS_minus_IS_eV<0,"catalyst"].tolist()
summary={"n_catalysts":len(VALID),"surface_names":VALID.surface_name.tolist(),
         "NEB1":{"fit":{k:v for k,v in F1.items() if k not in ("fitted","loocv")},"negative_signed_TS_minus_IS":neg1},
         "NEB2":{"fit":{k:v for k,v in F2.items() if k not in ("fitted","loocv")},"negative_signed_TS_minus_IS":neg2},
         "CatMAP_run_called":False,"steady_state_rerun":False,"grid_reused":True}
for path, data in [("NEB1_path2",B1),("NEB2_path4",B2)]:
    d=data.loc[~data.catalyst.eq("Ti-Fe")]
    f=fit_line(d.reaction_energy_eV,d.signed_TS_minus_IS_eV)
    summary[path+"_TiFe_excluded_sensitivity"]={k:v for k,v in f.items() if k not in ("fitted","loocv")}
    best_ext=sorted([r for r in ext if r["path"]==path],key=lambda r:r["RMSE_eV"])[0]
    summary[path+"_best_one_extra_descriptor_in_sample"]={k:best_ext[k] for k in ("added_descriptor","R2","RMSE_eV","LOOCV_RMSE_eV","LOOCV_R2","descriptor_reference_note")}
(D/"bep_analysis_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
