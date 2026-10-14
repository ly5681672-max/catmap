"""Explore O/H historical formation energies -> TS energies -> saved CatMAP TOF.

This script never calls ReactionModel.run() or changes any CatMAP input.
All joins are keyed by surface_name. Empirical and saved-grid interpolation
products are exported separately, and points outside the saved TS grid remain
missing rather than being clamped to its boundary.
"""
from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
OUT = ROOT / "origin_exports"
for p in (FIG / "scaling", FIG / "bep", FIG / "tof_maps", OUT / "maps"):
    p.mkdir(parents=True, exist_ok=True)

STATE_NAMES = ["Rpair", "TS1", "EpoxideHHO", "HHO", "TS2", "OHHa", "H2O"]
STATE_COL = {s: f"{s}_s" for s in STATE_NAMES}
FEATURES = {"O": "O", "Ha": "Ha", "Hb": "Hb"}
FEATURE_SETS = [(x,) for x in FEATURES] + list(itertools.combinations(FEATURES, 2))
FEATURE_SETS_WITH_TRIPLE = FEATURE_SETS + [("O", "Ha", "Hb")]
SURFACE_ALIAS = {
    "tife": "Ti-Fe", "timn": "Ti-Mn", "tihf": "Ti-Hf", "tire": "Ti-Re",
    "timo": "Ti-Mo", "tiv": "Ti-V", "tico": "Ti-Co", "titi": "Ti-Ti",
    "tiw": "Ti-W", "tita": "Ti-Ta", "ticr": "Ti-Cr", "ti": "Ti", "tini": "Ti-Ni",
}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.titlesize": 12, "figure.dpi": 130,
                     "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42})


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            clean = {k: ("" if isinstance(v, float) and not math.isfinite(v) else v)
                     for k, v in row.items()}
            w.writerow(clean)


def read_energy_table(path: Path) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            surface = (row.get("surface_name") or "").strip()
            species = (row.get("species_name") or "").strip()
            if surface not in SURFACE_ALIAS or species not in STATE_NAMES:
                continue
            if species in result.setdefault(surface, {}):
                raise ValueError(f"Duplicate energy row: {surface}/{species}")
            result[surface][species] = float(row["formation_energy"])
    return result


def build_master() -> list[dict]:
    matched_rows = read_csv(DATA / "matched_data.csv")
    slab_rows = read_csv(DATA / "reference_slab_comparison.csv")
    energies = read_energy_table(DATA / "energies.txt")
    matched = {}
    for r in matched_rows:
        key = r["surface_name"].strip().lower()
        if key in matched:
            raise ValueError(f"Duplicate matched surface_name {key}")
        matched[key] = r
    slabs = {r["surface_name"].strip().lower(): r for r in slab_rows}
    expected = set(SURFACE_ALIAS)
    if set(matched) != expected:
        raise ValueError(f"TOF/OH surface keys differ: missing={expected-set(matched)}, extra={set(matched)-expected}")
    if set(energies) != expected:
        raise ValueError(f"Network energy surface keys differ: missing={expected-set(energies)}, extra={set(energies)-expected}")
    if set(slabs) != expected:
        raise ValueError(f"Slab comparison keys differ: missing={expected-set(slabs)}, extra={set(slabs)-expected}")

    rows = []
    for surface in SURFACE_ALIAS:
        m, e, slab = matched[surface], energies[surface], slabs[surface]
        if m["catalyst"] != SURFACE_ALIAS[surface]:
            raise ValueError(f"Catalyst alias mismatch for {surface}: {m['catalyst']}")
        if m["status"] != "CONVERGED" or abs(float(m["temperature_K"]) - 333.15) > 1e-8:
            raise ValueError(f"Unexpected baseline status/temperature for {surface}")
        missing_states = set(STATE_NAMES) - set(e)
        if missing_states:
            raise ValueError(f"Missing state energies for {surface}: {missing_states}")
        row = {"surface_name": surface, "catalyst": SURFACE_ALIAS[surface],
               "temperature_K": 333.15, "status": m["status"],
               "TOF_C6H12O": float(m["TOF"]), "log10_TOF_C6H12O": float(m["log10_TOF"]),
               "O_historical_eV": float(m["O"]), "Ha_historical_eV": float(m["Ha"]),
               "Hb_historical_eV": float(m["Hb"]),
               "legacy_slab_cell": slab["legacy_clean_slab_cell"],
               "legacy_clean_slab_G_eV": float(slab["legacy_clean_slab_G_eV"]),
               "strict_slab_cell": slab["strict_clean_slab_cell"],
               "strict_clean_slab_G_eV": float(slab["strict_clean_slab_G_eV"]),
               "strict_minus_legacy_slab_eV": float(slab["strict_minus_legacy_slab_eV"]),
               "legacy_descriptor_reference_status": "historical-reference exploratory; no slab correction applied"}
        for state in STATE_NAMES:
            row[STATE_COL[state]] = e[state]
        row["dG2_EpoxideHHO_minus_Rpair_eV"] = e["EpoxideHHO"] - e["Rpair"]
        row["dG2_act_TS1_minus_Rpair_eV"] = e["TS1"] - e["Rpair"]
        row["dG4_OHHa_minus_HHO_eV"] = e["OHHa"] - e["HHO"]
        row["dG4_act_TS2_minus_HHO_eV"] = e["TS2"] - e["HHO"]
        rows.append(row)
    return rows


def design(x: np.ndarray, means: np.ndarray, scales: np.ndarray) -> np.ndarray:
    return (x - means) / scales


def fit_linear(X: np.ndarray, y: np.ndarray) -> dict:
    means = X.mean(axis=0)
    scales = X.std(axis=0)
    scales[scales == 0] = 1.0
    Z = design(X, means, scales)
    A = np.column_stack([np.ones(len(y)), Z])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return {"means": means, "scales": scales, "beta": beta}


def predict_linear(model: dict, X: np.ndarray) -> np.ndarray:
    return np.column_stack([np.ones(len(X)), design(X, model["means"], model["scales"])]) @ model["beta"]


def metric_row(y: np.ndarray, pred: np.ndarray) -> dict:
    residual = pred - y
    denom = float(np.sum((y - y.mean()) ** 2))
    return {"R2_in_sample": 1.0 - float(np.sum(residual**2)) / denom if denom else float("nan"),
            "RMSE_in_sample": float(np.sqrt(np.mean(residual**2))),
            "MAE_in_sample": float(np.mean(np.abs(residual))),
            "max_abs_error_in_sample": float(np.max(np.abs(residual)))}


def loo_predict(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.empty(len(y), dtype=float)
    for i in range(len(y)):
        mask = np.arange(len(y)) != i
        out[i] = predict_linear(fit_linear(X[mask], y[mask]), X[i:i+1])[0]
    return out


def extended_metrics(y: np.ndarray, pred: np.ndarray, loo: np.ndarray) -> dict:
    denom = float(np.sum((y - y.mean())**2))
    return {**metric_row(y, pred),
            "LOOCV_RMSE": float(np.sqrt(np.mean((loo-y)**2))),
            "LOOCV_MAE": float(np.mean(np.abs(loo-y))),
            "LOOCV_R2": 1.0-float(np.sum((loo-y)**2))/denom if denom else float("nan")}


def feature_matrix(rows: list[dict], features: tuple[str, ...]) -> np.ndarray:
    return np.array([[float(r[f"{f}_historical_eV"]) for f in features] for r in rows])


def analyze_models(rows: list[dict]) -> tuple[list[dict], list[dict], dict]:
    metrics, predictions, selected = [], [], {}
    y_targets = {STATE_COL[s]: np.array([r[STATE_COL[s]] for r in rows], float) for s in STATE_NAMES}
    y_targets.update({
        "dG2_eV": np.array([r["dG2_EpoxideHHO_minus_Rpair_eV"] for r in rows], float),
        "dG4_eV": np.array([r["dG4_OHHa_minus_HHO_eV"] for r in rows], float),
    })
    for target, y in y_targets.items():
        target_models = []
        for fs in FEATURE_SETS_WITH_TRIPLE:
            X = feature_matrix(rows, fs)
            model = fit_linear(X, y)
            pred = predict_linear(model, X)
            loo = loo_predict(X, y)
            names = [FEATURES[f] for f in fs]
            row = {"target": target, "descriptor_set": "+".join(fs), "n_predictors": len(fs),
                   **extended_metrics(y, pred, loo), "intercept_centered": float(model["beta"][0])}
            for idx, f in enumerate(fs):
                row[f"mean_{f}"] = float(model["means"][idx])
                row[f"scale_{f}"] = float(model["scales"][idx])
                row[f"coef_{f}_standardized"] = float(model["beta"][idx+1])
            metrics.append(row)
            candidate = {"features": fs, "model": model, "metrics": row}
            target_models.append(candidate)
            for i, catalyst in enumerate(rows):
                predictions.append({"surface_name": catalyst["surface_name"], "catalyst": catalyst["catalyst"],
                                    "target": target, "descriptor_set": "+".join(fs),
                                    "actual_eV": float(y[i]), "predicted_in_sample_eV": float(pred[i]),
                                    "residual_pred_minus_actual_eV": float(pred[i]-y[i]),
                                    "predicted_LOOCV_eV": float(loo[i]),
                                    "LOOCV_residual_eV": float(loo[i]-y[i])})
        # Choose among 1- and 2-descriptor models only; prefer fewer terms only on exact ties.
        simple = [x for x in target_models if len(x["features"]) <= 2]
        selected[target] = min(simple, key=lambda x: (x["metrics"]["RMSE_in_sample"],
                                                      len(x["features"]), x["metrics"]["max_abs_error_in_sample"]))
    metrics.sort(key=lambda r: (r["target"], r["n_predictors"], r["RMSE_in_sample"]))
    return metrics, predictions, selected


def bep_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, np.ndarray]:
    # Return physical alpha/beta in y = alpha*x + beta, plus fitted y.
    model = fit_linear(x[:, None], y)
    alpha = float(model["beta"][1] / model["scales"][0])
    beta = float(model["beta"][0] - alpha * model["means"][0])
    return alpha, beta, alpha*x+beta


def analyze_bep(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    definitions = [
        ("NEB1_path2", "dG2_EpoxideHHO_minus_Rpair_eV", "dG2_act_TS1_minus_Rpair_eV", "Rpair", "TS1"),
        ("NEB2_path4", "dG4_OHHa_minus_HHO_eV", "dG4_act_TS2_minus_HHO_eV", "HHO", "TS2"),
    ]
    summaries, records = [], []
    for name, rxn_col, act_col, initial, ts in definitions:
        x = np.array([r[rxn_col] for r in rows], float)
        y = np.array([r[act_col] for r in rows], float)
        alpha, beta, pred = bep_fit(x, y)
        loo = np.empty(len(y))
        for i in range(len(y)):
            mask = np.arange(len(y)) != i
            a_i, b_i, _ = bep_fit(x[mask], y[mask])
            loo[i] = a_i*x[i]+b_i
        summary = {"path": name, "equation": f"{act_col} = alpha*{rxn_col} + beta",
                   "alpha": alpha, "beta_eV": beta, **extended_metrics(y, pred, loo),
                   "negative_TS_minus_IS_count": int(np.sum(y < 0)),
                   "negative_TS_minus_IS_catalysts": ";".join(r["catalyst"] for i,r in enumerate(rows) if y[i] < 0)}
        summaries.append(summary)
        for i, r in enumerate(rows):
            records.append({"surface_name": r["surface_name"], "catalyst": r["catalyst"], "path": name,
                            "initial_state": initial, "transition_state": ts,
                            "reaction_energy_eV": float(x[i]), "raw_TS_minus_IS_eV": float(y[i]),
                            "BEP_fitted_in_sample_eV": float(pred[i]),
                            "BEP_residual_eV": float(pred[i]-y[i]),
                            "BEP_predicted_LOOCV_eV": float(loo[i]),
                            "LOOCV_residual_eV": float(loo[i]-y[i])})
    return summaries, records


def chain_route(rows: list[dict], features: tuple[str, ...], which: str) -> tuple[dict, list[dict]]:
    if which == "NEB1_path2":
        init_col, dg_col, barrier_col, ts_name = "Rpair_s", "dG2_EpoxideHHO_minus_Rpair_eV", "dG2_act_TS1_minus_Rpair_eV", "TS1_s"
    else:
        init_col, dg_col, barrier_col, ts_name = "HHO_s", "dG4_OHHa_minus_HHO_eV", "dG4_act_TS2_minus_HHO_eV", "TS2_s"
    X = feature_matrix(rows, features)
    y_init = np.array([r[init_col] for r in rows], float)
    y_dg = np.array([r[dg_col] for r in rows], float)
    y_bar = np.array([r[barrier_col] for r in rows], float)
    y_ts = np.array([r[ts_name] for r in rows], float)
    alpha, beta, _ = bep_fit(y_dg, y_bar)
    init_model, dg_model = fit_linear(X, y_init), fit_linear(X, y_dg)
    init_pred, dg_pred = predict_linear(init_model, X), predict_linear(dg_model, X)
    pred_ts = init_pred + alpha*dg_pred + beta
    loo = np.empty(len(rows))
    for i in range(len(rows)):
        mask = np.arange(len(rows)) != i
        a_i, b_i, _ = bep_fit(y_dg[mask], y_bar[mask])
        mi, md = fit_linear(X[mask], y_init[mask]), fit_linear(X[mask], y_dg[mask])
        i_pred = predict_linear(mi, X[i:i+1])[0]
        d_pred = predict_linear(md, X[i:i+1])[0]
        loo[i] = i_pred + a_i*d_pred + b_i
    summary = {"path": which, "descriptor_set": "+".join(features), "n_predictors_each": len(features),
               "BEP_alpha": alpha, "BEP_beta_eV": beta, **extended_metrics(y_ts, pred_ts, loo),
               "initial_state_fit_RMSE_eV": metric_row(y_init, init_pred)["RMSE_in_sample"],
               "reaction_energy_fit_RMSE_eV": metric_row(y_dg, dg_pred)["RMSE_in_sample"],
               "chain_definition": f"predicted {init_col} + alpha*predicted {dg_col} + beta"}
    details = [{"surface_name": r["surface_name"], "catalyst": r["catalyst"], "path": which,
                "descriptor_set": "+".join(features), "actual_TS_eV": float(y_ts[i]),
                "predicted_initial_state_eV": float(init_pred[i]), "predicted_reaction_energy_eV": float(dg_pred[i]),
                "BEP_predicted_barrier_eV": float(alpha*dg_pred[i]+beta),
                "predicted_TS_chain_eV": float(pred_ts[i]), "chain_residual_eV": float(pred_ts[i]-y_ts[i]),
                "predicted_TS_LOOCV_eV": float(loo[i]), "LOOCV_residual_eV": float(loo[i]-y_ts[i])}
               for i,r in enumerate(rows)]
    return summary, details


def direct_best_models(rows: list[dict], selected: dict) -> list[dict]:
    out = []
    for target in ("TS1_s", "TS2_s"):
        cand = selected[target]
        X = feature_matrix(rows, cand["features"])
        y = np.array([r[target] for r in rows], float)
        pred = predict_linear(cand["model"], X)
        loo = loo_predict(X, y)
        for i, r in enumerate(rows):
            out.append({"surface_name": r["surface_name"], "catalyst": r["catalyst"], "target": target,
                        "descriptor_set": "+".join(cand["features"]), "actual_TS_eV": float(y[i]),
                        "predicted_TS_in_sample_eV": float(pred[i]), "residual_eV": float(pred[i]-y[i]),
                        "predicted_TS_LOOCV_eV": float(loo[i]), "LOOCV_residual_eV": float(loo[i]-y[i])})
    return out


def quadratic_tof(ts1: np.ndarray, ts2: np.ndarray, coef_data: dict) -> np.ndarray:
    c = np.asarray(coef_data["coefficients_standardized"], float)
    sc = coef_data["scales"]
    u = (ts1-float(sc["TS1_mean"])) / float(sc["TS1_std"])
    v = (ts2-float(sc["TS2_mean"])) / float(sc["TS2_std"])
    return c[0]+c[1]*u+c[2]*v+c[3]*u*u+c[4]*u*v+c[5]*v*v


def write_matrix(path: Path, xs: np.ndarray, ys: np.ndarray, z: np.ndarray) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x\\y", *[f"{v:.12g}" for v in ys]])
        for i, x in enumerate(xs):
            w.writerow([f"{x:.12g}", *["" if not np.isfinite(v) else f"{v:.15g}" for v in z[i,:]]])


def fit_oh_pair(rows: list[dict], pair: tuple[str, str], target: str) -> tuple[dict, np.ndarray, np.ndarray]:
    X = feature_matrix(rows, pair)
    y = np.array([r[target] for r in rows], float)
    m = fit_linear(X, y)
    return m, predict_linear(m, X), loo_predict(X, y)


def model_value(model: dict, pair: tuple[str, str], X: np.ndarray) -> np.ndarray:
    return predict_linear(model, X)


def closest_labels(ax, xs, ys, labels, fontsize=7.5, offset_by_label=None):
    # Deterministic, tiny stagger to keep point markers and names visible.
    offsets = [(5,6),(6,-10),(-34,8),(7,9),(-32,-10),(6,7),(-30,6),(7,-11),
               (6,8),(-32,9),(6,-10),(-20,9),(7,7)]
    offset_by_label = offset_by_label or {}
    for i, (x,y,label) in enumerate(zip(xs,ys,labels)):
        ax.annotate(label, (x,y), xytext=offset_by_label.get(label, offsets[i % len(offsets)]), textcoords="offset points",
                    fontsize=fontsize, color="#172b4d", annotation_clip=False)


def savefig(fig, stem: Path):
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)


def plot_model_parity(rows: list[dict], selected: dict):
    fig, axes = plt.subplots(1,2,figsize=(11,5),layout="constrained")
    for ax,target,title in zip(axes,["TS1_s","TS2_s"],["TS1* direct O/H mapping","TS2* direct O/H mapping"]):
        c=selected[target]; X=feature_matrix(rows,c["features"]); y=np.array([r[target] for r in rows],float); pred=predict_linear(c["model"],X)
        ax.scatter(y,pred,s=48,facecolors="white",edgecolors="#1769aa",linewidth=1.2,zorder=3)
        closest_labels(ax,y,pred,[r["catalyst"] for r in rows],7)
        lo=min(y.min(),pred.min()); hi=max(y.max(),pred.max()); ax.plot([lo,hi],[lo,hi],"--",color="0.4",lw=1)
        ax.set(xlabel="Actual strict-model formation energy (eV)",ylabel="Fitted formation energy (eV)",title=f"{title}\nfeatures={'+'.join(c['features'])}; RMSE={c['metrics']['RMSE_in_sample']:.3f} eV")
        ax.grid(alpha=.18)
    savefig(fig,FIG/"scaling"/"best_direct_TS_mapping_parity")


def plot_bep(bep_summary: list[dict], bep_rows: list[dict]):
    fig,axes=plt.subplots(1,2,figsize=(11,5),layout="constrained")
    for ax,path,title in zip(axes,["NEB1_path2","NEB2_path4"],["Path 2: Rpair* → TS1* → EpoxideHHO*","Path 4: HHO* → TS2* → OHHa*"]):
        rr=[r for r in bep_rows if r["path"]==path]; s=next(r for r in bep_summary if r["path"]==path)
        x=np.array([r["reaction_energy_eV"] for r in rr]); y=np.array([r["raw_TS_minus_IS_eV"] for r in rr]); yp=s["alpha"]*x+s["beta_eV"]
        ax.scatter(x,y,s=45,facecolors="white",edgecolors="#bd3b30",linewidth=1.2,zorder=3)
        closest_labels(ax,x,y,[r["catalyst"] for r in rr],7)
        order=np.argsort(x); ax.plot(x[order],yp[order],color="#255a91",lw=2,label="OLS BEP fit")
        ax.axhline(0,color="0.45",lw=.8,ls=":")
        ax.set(xlabel="Reaction energy ΔG (eV)",ylabel="Raw NEB TS−IS difference ΔG‡ (eV)",title=f"{title}\nα={s['alpha']:.3f}, β={s['beta_eV']:.3f} eV, R²={s['R2_in_sample']:.3f}, LOOCV RMSE={s['LOOCV_RMSE']:.3f} eV")
        ax.grid(alpha=.18); ax.legend(frameon=False)
    savefig(fig,FIG/"bep"/"NEB2_NEb4_BEP_relationships")


def plot_direct_chain(rows: list[dict], direct_rows: list[dict], chain_rows: list[dict]):
    fig,axes=plt.subplots(2,2,figsize=(11,9),layout="constrained")
    for rowi,target in enumerate(["TS1_s","TS2_s"]):
        for coli,route in enumerate(["direct","chain"]):
            ax=axes[rowi,coli]
            if route=="direct":
                data=[r for r in direct_rows if r["target"]==target]
                actual=np.array([r["actual_TS_eV"] for r in data]); pred=np.array([r["predicted_TS_in_sample_eV"] for r in data]); fs=data[0]["descriptor_set"]
            else:
                p="NEB1_path2" if target=="TS1_s" else "NEB2_path4"
                data=[r for r in chain_rows if r["path"]==p and r["descriptor_set"] in ("O+Ha","O+Hb")]
                fs="O+Ha and O+Hb"
                actual=np.array([r["actual_TS_eV"] for r in data]); pred=np.array([r["predicted_TS_chain_eV"] for r in data])
            ax.scatter(actual,pred,s=42,facecolors="white",edgecolors="#26734d" if route=="direct" else "#bf6b21",zorder=3)
            if route=="direct": labels=[r["catalyst"] for r in data]
            else: labels=[r["catalyst"]+"/"+r["descriptor_set"] for r in data]
            closest_labels(ax,actual,pred,labels,6.3)
            lo=min(actual.min(),pred.min()); hi=max(actual.max(),pred.max()); ax.plot([lo,hi],[lo,hi],"--",color="0.4",lw=1)
            err=float(np.sqrt(np.mean((pred-actual)**2)))
            ax.set(xlabel="Actual strict-model state energy (eV)",ylabel="Reconstructed state energy (eV)",title=f"{target}: {route} route; {fs}\nin-sample RMSE={err:.3f} eV")
            ax.grid(alpha=.18)
    savefig(fig,FIG/"scaling"/"direct_vs_BEP_chain_TS_reconstruction")


def plot_ts_coordinates(rows: list[dict], chain_rows: list[dict]):
    fig,ax=plt.subplots(figsize=(8.2,7),layout="constrained")
    colors={"O+Ha":"#2b6cb0","O+Hb":"#d97706"}; markers={"O+Ha":"o","O+Hb":"s"}
    actual_x=np.array([r["TS1_s"] for r in rows]); actual_y=np.array([r["TS2_s"] for r in rows])
    ax.scatter(actual_x,actual_y,marker="D",s=70,facecolor="white",edgecolor="black",label="Actual 13-site TS coordinates",zorder=5)
    closest_labels(ax,actual_x,actual_y,[r["catalyst"] for r in rows],7.2)
    for pair in ("O+Ha","O+Hb"):
        p1="NEB1_path2"; p2="NEB2_path4"
        a={r["surface_name"]:r for r in chain_rows if r["path"]==p1 and r["descriptor_set"]==pair}
        b={r["surface_name"]:r for r in chain_rows if r["path"]==p2 and r["descriptor_set"]==pair}
        xx=np.array([a[r["surface_name"]]["predicted_TS_chain_eV"] for r in rows]); yy=np.array([b[r["surface_name"]]["predicted_TS_chain_eV"] for r in rows])
        ax.scatter(xx,yy,marker=markers[pair],s=42,color=colors[pair],alpha=.85,label=f"O/H→reaction ΔG→BEP→TS ({pair})",zorder=4)
        for i,r in enumerate(rows):
            ax.plot([actual_x[i],xx[i]],[actual_y[i],yy[i]],color=colors[pair],lw=.55,alpha=.28,zorder=1)
    ax.set(xlabel="TS1_s formation energy (eV)",ylabel="TS2_s formation energy (eV)",title="Actual vs O/H→BEP reconstructed TS1–TS2 coordinates")
    ax.legend(frameon=True,fontsize=8); ax.grid(alpha=.18)
    savefig(fig,FIG/"scaling"/"actual_vs_OH_chain_TS_coordinates")


def load_saved_catmap_interpolant() -> tuple[RegularGridInterpolator, tuple[float,float,float,float]]:
    path=DATA/"catmap_VectorMap_interpolated_400x400_xyz.csv"
    arr=np.genfromtxt(path,delimiter=",",names=True,encoding="utf-8-sig")
    names=arr.dtype.names
    xname,yname=names[0],names[1]
    zname=next(n for n in names if "log10_TOF" in n)
    x=np.unique(arr[xname]); y=np.unique(arr[yname])
    if len(arr)!=len(x)*len(y): raise ValueError("Saved CatMAP VectorMap file is not a complete regular grid")
    z=arr[zname].reshape((len(x),len(y)))
    interp=RegularGridInterpolator((x,y),z,method="linear",bounds_error=False,fill_value=np.nan)
    return interp,(float(x.min()),float(x.max()),float(y.min()),float(y.max()))


def map_peak(z: np.ndarray, xs: np.ndarray, ys: np.ndarray) -> dict:
    valid=np.isfinite(z)
    if not valid.any(): return {"available":False}
    masked=np.where(valid,z,-np.inf); ij=np.unravel_index(np.argmax(masked),masked.shape); i,j=ij
    boundary=(i==0 or j==0 or i==len(xs)-1 or j==len(ys)-1)
    neighborhood=[]
    for di in (-1,0,1):
        for dj in (-1,0,1):
            ii,jj=i+di,j+dj
            if 0<=ii<len(xs) and 0<=jj<len(ys) and (di or dj): neighborhood.append(z[ii,jj])
    local=(not boundary and all(np.isfinite(v) and z[i,j]>v for v in neighborhood))
    return {"available":True,"x_eV":float(xs[i]),"y_eV":float(ys[j]),"log10_TOF":float(z[i,j]),
            "TOF":float(10**z[i,j]),"valid_grid_points":int(valid.sum()),
            "boundary_of_rectangular_grid":bool(boundary),"eight_neighbor_local_maximum":bool(local)}


def plot_map(xs: np.ndarray, ys: np.ndarray, z: np.ndarray, rows: list[dict], pair: tuple[str,str],
             activity_type: str, peak: dict):
    X,Y=np.meshgrid(xs,ys,indexing="ij")
    fig,ax=plt.subplots(figsize=(10,8),layout="constrained")
    cmap=plt.get_cmap("viridis").copy(); cmap.set_bad("#d9d9d9")
    cf=ax.contourf(X,Y,np.ma.masked_invalid(z),levels=28,cmap=cmap)
    cb=fig.colorbar(cf,ax=ax); cb.set_label(f"{activity_type}: predicted log10(TOF of C6H12O)")
    hx=f"{pair[0]}_historical_eV"; hy=f"{pair[1]}_historical_eV"
    px=np.array([r[hx] for r in rows]); py=np.array([r[hy] for r in rows])
    ax.scatter(px,py,marker="D",s=48,facecolor="white",edgecolor="black",lw=.8,label="13 strict-baseline catalyst points",zorder=6)
    label_offsets = {
        "Ti-Fe": (-42, 12), "Ti-Ti": (9, -17), "Ti-Mn": (8, -14),
        "Ti-Hf": (-34, 10), "Ti-Mo": (-40, -12), "Ti-Ta": (-38, 12),
        "Ti-Re": (8, 11), "Ti-Ni": (9, 11), "Ti-Co": (-40, 10),
        "Ti-Cr": (8, -14), "Ti-V": (8, 10), "Ti-W": (8, 10), "Ti": (8, 11),
    }
    closest_labels(ax,px,py,[r["catalyst"] for r in rows],7.3,label_offsets)
    if peak.get("available"):
        peak_label = ("Boundary maximum in mapped range" if peak.get("boundary_of_rectangular_grid")
                      else "Mapped-grid local maximum")
        ax.scatter([peak["x_eV"]],[peak["y_eV"]],marker="*",s=250,color="#ff3b30",edgecolor="black",lw=.8,label=peak_label,zorder=9,clip_on=False)
        # Place edge maxima inward so the star and its values remain visible.
        xmid=(xs[0]+xs[-1])/2; ymid=(ys[0]+ys[-1])/2
        xoff=-118 if peak["x_eV"] >= xmid else 10
        yoff=12 if peak["y_eV"] <= ymid else -66
        ax.annotate(f"{peak_label}\n({peak['x_eV']:.3f}, {peak['y_eV']:.3f}) eV\nlog10TOF={peak['log10_TOF']:.3f}\n{peak['valid_grid_points']}/160000 valid",
                    (peak["x_eV"],peak["y_eV"]),xytext=(xoff,yoff),textcoords="offset points",fontsize=8,
                    bbox={"boxstyle":"round,pad=.3","fc":"white","alpha":.92},arrowprops={"arrowstyle":"-","color":"#333","lw":.8},zorder=10,annotation_clip=False)
    ax.set(xlabel=f"{pair[0]} historical formation energy (eV)",ylabel=f"{pair[1]} historical formation energy (eV)",
           title=f"{pair[0]} + {pair[1]} → TS1/TS2 → {activity_type}\nO/H reference historical; exploratory; not new CatMAP solves")
    # A small display margin preserves the full rectangular grid while preventing
    # boundary stars, markers, or labels from being clipped by the axes frame.
    xpad=max(float(xs[-1]-xs[0])*.025,1e-6); ypad=max(float(ys[-1]-ys[0])*.025,1e-6)
    ax.set_xlim(xs[0]-xpad,xs[-1]+xpad); ax.set_ylim(ys[0]-ypad,ys[-1]+ypad)
    ax.legend(loc="best",fontsize=8,framealpha=.9); ax.grid(alpha=.12)
    savefig(fig,FIG/"tof_maps"/f"{pair[0]}_{pair[1]}_{activity_type}_mapped_activity")


def make_maps(rows: list[dict], selected: dict) -> tuple[list[dict], list[dict]]:
    coef_data=json.loads((DATA/"empirical_quadratic_coefficients.json").read_text(encoding="utf-8"))
    interp, ts_bounds=load_saved_catmap_interpolant()
    peak_records=[]; point_records=[]
    for pair in [("O","Ha"),("O","Hb")]:
        m1,pred1,loo1=fit_oh_pair(rows,pair,"TS1_s")
        m2,pred2,loo2=fit_oh_pair(rows,pair,"TS2_s")
        pair_name="+".join(pair)
        for i,r in enumerate(rows):
            point_records.append({"surface_name":r["surface_name"],"catalyst":r["catalyst"],"descriptor_set":pair_name,
                                  "O_historical_eV":r["O_historical_eV"],f"{pair[1]}_historical_eV":r[f"{pair[1]}_historical_eV"],
                                  "actual_TS1_eV":r["TS1_s"],"predicted_TS1_eV":float(pred1[i]),"TS1_LOOCV_predicted_eV":float(loo1[i]),
                                  "actual_TS2_eV":r["TS2_s"],"predicted_TS2_eV":float(pred2[i]),"TS2_LOOCV_predicted_eV":float(loo2[i]),
                                  "actual_TOF_C6H12O":r["TOF_C6H12O"],"actual_log10_TOF_C6H12O":r["log10_TOF_C6H12O"]})
        xvals=np.array([r[f"{pair[0]}_historical_eV"] for r in rows]); yvals=np.array([r[f"{pair[1]}_historical_eV"] for r in rows])
        dx=max(float(np.ptp(xvals))*.05,1e-6); dy=max(float(np.ptp(yvals))*.05,1e-6)
        xs=np.linspace(float(xvals.min()-dx),float(xvals.max()+dx),400); ys=np.linspace(float(yvals.min()-dy),float(yvals.max()+dy),400)
        X,Y=np.meshgrid(xs,ys,indexing="ij"); OH=np.column_stack([X.ravel(),Y.ravel()])
        ts1=model_value(m1,pair,OH).reshape(X.shape); ts2=model_value(m2,pair,OH).reshape(X.shape)
        empirical=quadratic_tof(ts1,ts2,coef_data)
        points=np.column_stack([ts1.ravel(),ts2.ravel()])
        inside=((points[:,0]>=ts_bounds[0])&(points[:,0]<=ts_bounds[1])&(points[:,1]>=ts_bounds[2])&(points[:,1]<=ts_bounds[3]))
        catmap=np.full(len(points),np.nan,dtype=float)
        catmap[inside]=interp(points[inside])
        catmap=catmap.reshape(X.shape)
        empirical_peak=map_peak(empirical,xs,ys); catmap_peak=map_peak(catmap,xs,ys)
        for activity,z,peak in [("empirical_TS_quadratic",empirical,empirical_peak),("saved_CatMAP_VectorMap_interpolation",catmap,catmap_peak)]:
            full_xyz=[]
            for i in range(400):
                for j in range(400):
                    full_xyz.append({f"{pair[0]}_historical_eV":float(xs[i]),f"{pair[1]}_historical_eV":float(ys[j]),
                                     "predicted_TS1_eV":float(ts1[i,j]),"predicted_TS2_eV":float(ts2[i,j]),
                                     "activity_type":activity,"predicted_log10_TOF_C6H12O":float(z[i,j]) if np.isfinite(z[i,j]) else float("nan"),
                                     "inside_saved_TS_grid":bool(inside[i*400+j]),
                                     "inside_observed_TS_training_rectangle":bool(ts1[i,j]>=-34.566395 and ts1[i,j]<=-32.451521 and ts2[i,j]>=-7.921799 and ts2[i,j]<=-4.880121)})
            slug=pair[0]+"_"+pair[1]+"_"+activity
            write_csv(OUT/"maps"/f"{slug}_400x400_xyz.csv",full_xyz)
            write_matrix(OUT/"maps"/f"{slug}_400x400_matrix.csv",xs,ys,z)
            write_csv(OUT/"maps"/f"{slug}_peak.csv",[{"descriptor_1":peak.get("x_eV"),"descriptor_2":peak.get("y_eV"),
                       "predicted_log10_TOF":peak.get("log10_TOF"),"predicted_TOF":peak.get("TOF"),
                       "valid_grid_points":peak.get("valid_grid_points"),"boundary_of_rectangular_grid":peak.get("boundary_of_rectangular_grid"),
                       "eight_neighbor_local_maximum":peak.get("eight_neighbor_local_maximum"),"activity_source":activity}])
            plot_map(xs,ys,z,rows,pair,activity,peak)
            peak_records.append({"descriptor_set":pair[0]+"+"+pair[1],"activity_source":activity,
                                 "x_range_eV":[float(xs.min()),float(xs.max())],"y_range_eV":[float(ys.min()),float(ys.max())],
                                 "TS_grid_bounds_eV":{"TS1":[ts_bounds[0],ts_bounds[1]],"TS2":[ts_bounds[2],ts_bounds[3]]},**peak})
    write_csv(OUT/"maps"/"O_H_mapped_TS_and_actual_catalyst_points.csv",point_records)
    return peak_records,point_records


def main():
    rows=build_master()
    write_csv(DATA/"matched_analysis_table.csv",rows)
    long=[]
    for r in rows:
        for state in STATE_NAMES:
            long.append({"surface_name":r["surface_name"],"catalyst":r["catalyst"],"state":state,
                         "formation_energy_eV":r[STATE_COL[state]],"source_file":"catmap_official_run/energies.txt",
                         "reference_note":"conditional strict-model input; underlying global-reference evidence remains per source audit"})
    write_csv(DATA/"network_state_energy_long.csv",long)
    metrics,predictions,selected=analyze_models(rows)
    write_csv(DATA/"linear_descriptor_model_metrics.csv",metrics)
    write_csv(DATA/"linear_descriptor_model_predictions.csv",predictions)
    bep_summary,bep_rows=analyze_bep(rows)
    write_csv(DATA/"BEP_model_summary.csv",bep_summary)
    write_csv(DATA/"BEP_catalyst_values_and_predictions.csv",bep_rows)
    chain_metrics=[]; chain_rows=[]
    for path in ["NEB1_path2","NEB2_path4"]:
        for fs in FEATURE_SETS:
            s,d=chain_route(rows,fs,path); chain_metrics.append(s); chain_rows.extend(d)
    write_csv(DATA/"O_H_to_BEP_chain_model_metrics.csv",chain_metrics)
    write_csv(DATA/"O_H_to_BEP_chain_TS_predictions.csv",chain_rows)
    direct_rows=direct_best_models(rows,selected)
    write_csv(DATA/"best_direct_TS_predictions.csv",direct_rows)
    plot_model_parity(rows,selected); plot_bep(bep_summary,bep_rows); plot_direct_chain(rows,direct_rows,chain_rows); plot_ts_coordinates(rows,chain_rows)
    peak_records,point_records=make_maps(rows,selected)
    selected_json={}
    for target,c in selected.items():
        selected_json[target]={"features":list(c["features"]),"model":{
            "means":c["model"]["means"].tolist(),"scales":c["model"]["scales"].tolist(),"beta_standardized":c["model"]["beta"].tolist()},
            "metrics":c["metrics"]}
    (DATA/"selected_linear_models.json").write_text(json.dumps(selected_json,indent=2,ensure_ascii=False),encoding="utf-8")
    (DATA/"mapped_activity_peaks.json").write_text(json.dumps(peak_records,indent=2,ensure_ascii=False),encoding="utf-8")
    summary={"n_catalysts":len(rows),"catalysts":[r["catalyst"] for r in rows],"temperature_K":333.15,
             "join_key":"surface_name","oh_reference":"historical; no slab correction applied",
             "excluded":"Ti-Zr (unresolved)","linear_models":metrics,"BEP":bep_summary,
             "chain_models":chain_metrics,"mapped_peaks":peak_records,
             "saved_CatMAP_grid":"CatMAP 0.3.1 VectorMap 400x400 interpolation of prior 50x50 converged solver grid; no new CatMAP solve in this script"}
    (DATA/"analysis_summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"n":len(rows),"best_TS1":selected["TS1_s"]["metrics"],"best_TS1_features":selected["TS1_s"]["features"],
                      "best_TS2":selected["TS2_s"]["metrics"],"best_TS2_features":selected["TS2_s"]["features"],
                      "BEP":bep_summary,"peaks":peak_records},indent=2,ensure_ascii=False))


if __name__ == "__main__":
    main()
