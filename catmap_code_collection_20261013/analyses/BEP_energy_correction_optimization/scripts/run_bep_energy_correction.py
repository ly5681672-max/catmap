from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PKG = Path(__file__).resolve().parents[1]
DATA = PKG / "data"
FIG = PKG / "figures"
EXP = PKG / "origin_exports"
REPORTS = PKG / "reports"
INPUT = DATA / "source_carbon_O_H_TS1.csv"
ENERGIES = DATA / "source_strict_energies.txt"
CATALYST_ORDER = ["tife", "timn", "tihf", "tire", "timo", "tiv", "tico", "titi", "tiw", "tita", "ticr", "ti", "tini"]
SPECIAL = {"tiw", "tita"}
FOCUS = {"tiw", "tita", "tico", "tife"}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})
COLORS = {"traditional_BEP": "#4C78A8", "model_A_OH_Ha_OOH": "#F58518", "model_B_OH_Ha_OHHa": "#54A24B", "variable_slope_OH": "#B279A2"}
MARKERS = {"Ti-W": "X", "Ti-Ta": "P"}
LABEL_OFFSETS = {
    "Ti-Fe": (-27, 8), "Ti-Mn": (-23, 10), "Ti-Hf": (6, -12), "Ti-Re": (8, -13),
    "Ti-Mo": (7, 10), "Ti-V": (-23, -13), "Ti-Co": (7, 9), "Ti-Ti": (-24, 10),
    "Ti-W": (7, 9), "Ti-Ta": (7, 8), "Ti-Cr": (6, 10), "Ti": (-14, -17), "Ti-Ni": (7, -13),
}
RIGHT_PANEL_OFFSETS = {
    "Ti-Fe": (-36, -17), "Ti-Mn": (-30, 14), "Ti-Hf": (-34, -12), "Ti-Re": (12, -25),
    "Ti-Mo": (12, -17), "Ti-V": (9, -19), "Ti-Co": (8, 8), "Ti-Ti": (-38, 15),
    "Ti-W": (8, 9), "Ti-Ta": (8, 8), "Ti-Cr": (-7, 18), "Ti": (-15, -20), "Ti-Ni": (9, 15),
}


def metric(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    residual = np.asarray(y) - np.asarray(pred)
    ss_res = float(np.dot(residual, residual))
    centered = np.asarray(y) - np.mean(y)
    ss_tot = float(np.dot(centered, centered))
    return {
        "R2": 1.0 - ss_res / ss_tot if ss_tot else float("nan"),
        "RMSE_eV": float(np.sqrt(np.mean(residual**2))),
        "MAE_eV": float(np.mean(np.abs(residual))),
        "MaxAbsError_eV": float(np.max(np.abs(residual))),
    }


def design(df: pd.DataFrame, model: str) -> tuple[np.ndarray, list[str]]:
    x = df["deltaG2_eV"].to_numpy(float)
    if model == "traditional_BEP":
        cols = [x]
        names = ["deltaG2_eV"]
    elif model == "model_A_OH_Ha_OOH":
        names = ["deltaG2_eV", "OH_star_legacy_eV", "Ha_star_legacy_eV", "OOH_star_legacy_eV"]
        cols = [df[n].to_numpy(float) for n in names]
    elif model == "model_B_OH_Ha_OHHa":
        names = ["deltaG2_eV", "OH_star_legacy_eV", "Ha_star_legacy_eV", "OHHa_strict_eV"]
        cols = [df[n].to_numpy(float) for n in names]
    elif model == "variable_slope_OH":
        z = df["OH_centered_eV"].to_numpy(float)
        names = ["deltaG2_eV", "deltaG2_x_OH_centered", "OH_centered_eV"]
        cols = [x, x * z, z]
    else:
        raise ValueError(model)
    return np.column_stack([*cols, np.ones(len(df))]), [*names, "intercept_eV"]


def fit(df: pd.DataFrame, model: str) -> tuple[np.ndarray, list[str], np.ndarray]:
    matrix, names = design(df, model)
    coefs, _, rank, _ = np.linalg.lstsq(matrix, df["Y_TS1_minus_Rpair_eV"].to_numpy(float), rcond=None)
    if rank < matrix.shape[1]:
        raise ValueError(f"rank deficient: {model}, rank={rank}, columns={matrix.shape[1]}")
    return coefs, names, matrix @ coefs


def loocv(df: pd.DataFrame, model: str) -> tuple[np.ndarray, float, float]:
    observed = df["Y_TS1_minus_Rpair_eV"].to_numpy(float)
    pred = np.empty_like(observed)
    for i in range(len(df)):
        train = df.drop(df.index[i])
        coefs, _, _ = fit(train, model)
        matrix, _ = design(df.iloc[[i]], model)
        pred[i] = float((matrix @ coefs).item())
    m = metric(observed, pred)
    return pred, m["RMSE_eV"], m["R2"]


def get_data() -> pd.DataFrame:
    raw = pd.read_csv(INPUT)
    energies = pd.read_csv(ENERGIES, sep="\t")
    ohha = energies.loc[(energies["species_name"] == "OHHa") & energies["surface_name"].notna(), ["surface_name", "formation_energy"]].copy()
    ohha = ohha.rename(columns={"formation_energy": "OHHa_strict_eV"})
    if ohha["surface_name"].duplicated().any():
        raise ValueError("strict energies contain duplicate OHHa surface rows")
    df = raw.merge(ohha, on="surface_name", how="left", validate="one_to_one")
    df = df[df["status"].astype(str).str.upper().eq("CONVERGED")].copy()
    df = df[~df["surface_name"].eq("tizr")].copy()
    df["surface_name"] = df["surface_name"].astype(str).str.lower()
    df = df.set_index("surface_name", drop=False).reindex(CATALYST_ORDER).reset_index(drop=True)
    required = ["catalyst", "Rpair_strict", "TS1_strict", "EpoxideHHO_strict", "OH_star_legacy", "Ha_star_legacy", "OOH_star_legacy", "OHHa_strict_eV"]
    missing = [c for c in required if c not in df or df[c].isna().any()]
    if missing:
        raise ValueError(f"missing required columns/values: {missing}")
    if len(df) != 13 or df["surface_name"].duplicated().any():
        raise ValueError(f"expected 13 unique catalyst rows, got {len(df)}")
    df["Y_TS1_minus_Rpair_eV"] = df["TS1_strict"].astype(float) - df["Rpair_strict"].astype(float)
    df["deltaG2_eV"] = df["EpoxideHHO_strict"].astype(float) - df["Rpair_strict"].astype(float)
    if not np.allclose(df["Y_TS1_minus_Rpair_eV"], df["NEB1_forward_energy_difference_TS1_minus_Rpair"].astype(float), atol=1e-10):
        raise ValueError("Y target differs from stored NEB1 signed difference")
    if not np.allclose(df["deltaG2_eV"], df["NEB1_reaction_energy_EpoxideHHO_minus_Rpair"].astype(float), atol=1e-10):
        raise ValueError("deltaG2 differs from stored NEB1 reaction energy")
    df["OH_star_legacy_eV"] = df["OH_star_legacy"].astype(float)
    df["Ha_star_legacy_eV"] = df["Ha_star_legacy"].astype(float)
    df["OOH_star_legacy_eV"] = df["OOH_star_legacy"].astype(float)
    df["OH_centered_eV"] = df["OH_star_legacy_eV"] - df["OH_star_legacy_eV"].mean()
    return df


MODELS = ["traditional_BEP", "model_A_OH_Ha_OOH", "model_B_OH_Ha_OHHa", "variable_slope_OH"]
MODEL_LABELS = {
    "traditional_BEP": "Traditional BEP",
    "model_A_OH_Ha_OOH": "BEP + OH + Ha + OOH",
    "model_B_OH_Ha_OHHa": "BEP + OH + Ha + OHHa",
    "variable_slope_OH": "OH-modulated slope",
}


def records_and_diagnostics(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metrics_rows, coef_rows, residual_rows, cv_rows, influence_rows = [], [], [], [], []
    coef_lookup: dict[tuple[str, str], tuple[np.ndarray, list[str]]] = {}
    for subset_name, subset in [("13site_refit", df), ("11site_refit", df[~df.surface_name.isin(SPECIAL)].copy())]:
        for model in MODELS:
            coefs, names, pred = fit(subset, model)
            cv_pred, cv_rmse, cv_r2 = loocv(subset, model)
            met = metric(subset["Y_TS1_minus_Rpair_eV"].to_numpy(float), pred)
            metrics_rows.append({"subset": subset_name, "model": model, "n": len(subset), **met, "LOOCV_RMSE_eV": cv_rmse, "LOOCV_R2": cv_r2})
            coef_lookup[(subset_name, model)] = (coefs, names)
            for term, value in zip(names, coefs):
                coef_rows.append({"subset": subset_name, "model": model, "term": term, "coefficient": float(value), "units": "dimensionless" if term == "deltaG2_eV" else "eV/eV" if term not in {"intercept_eV"} else "eV"})
            for i, (_, row) in enumerate(subset.iterrows()):
                residual_rows.append({"subset": subset_name, "model": model, "surface_name": row.surface_name, "catalyst": row.catalyst,
                                      "actual_Y_eV": row.Y_TS1_minus_Rpair_eV, "predicted_Y_eV": pred[i], "residual_actual_minus_predicted_eV": row.Y_TS1_minus_Rpair_eV-pred[i],
                                      "LOOCV_predicted_Y_eV": cv_pred[i], "LOOCV_residual_eV": row.Y_TS1_minus_Rpair_eV-cv_pred[i]})
            cv_rows.append({"subset": subset_name, "model": model, "surface_name": "ALL", "LOOCV_RMSE_eV": cv_rmse, "LOOCV_R2": cv_r2})
            for held in subset.surface_name.tolist():
                train = subset[subset.surface_name.ne(held)]
                c_hold, n_hold, _ = fit(train, model)
                full, full_names = coef_lookup[(subset_name, model)]
                if n_hold != full_names:
                    raise AssertionError("coefficient order mismatch")
                for term, fullval, holdval in zip(full_names, full, c_hold):
                    influence_rows.append({"subset": subset_name, "model": model, "held_out_catalyst": held, "term": term,
                                           "full_fit_coefficient": float(fullval), "leave_one_out_coefficient": float(holdval),
                                           "coefficient_change": float(holdval-fullval)})
    # Apply the 13-site coefficients to the same 11-site subset; this is not an 11-site refit.
    sub11 = df[~df.surface_name.isin(SPECIAL)].copy()
    for model in MODELS:
        coef, names = coef_lookup[("13site_refit", model)]
        mat, names2 = design(sub11, model)
        if names != names2:
            raise AssertionError("coefficient order mismatch in transfer evaluation")
        met = metric(sub11["Y_TS1_minus_Rpair_eV"].to_numpy(float), mat @ coef)
        metrics_rows.append({"subset": "13site_coefficients_applied_to_11", "model": model, "n": len(sub11), **met, "LOOCV_RMSE_eV": np.nan, "LOOCV_R2": np.nan})
        for i, (_, row) in enumerate(sub11.iterrows()):
            residual_rows.append({"subset": "13site_coefficients_applied_to_11", "model": model, "surface_name": row.surface_name, "catalyst": row.catalyst,
                                  "actual_Y_eV": row.Y_TS1_minus_Rpair_eV, "predicted_Y_eV": (mat@coef)[i], "residual_actual_minus_predicted_eV": row.Y_TS1_minus_Rpair_eV-(mat@coef)[i],
                                  "LOOCV_predicted_Y_eV": np.nan, "LOOCV_residual_eV": np.nan})
    return tuple(pd.DataFrame(x) for x in [metrics_rows, coef_rows, residual_rows, cv_rows, influence_rows])


def ablations(df: pd.DataFrame) -> pd.DataFrame:
    specs = {
        "A_OH_Ha_OOH": ["deltaG2_eV", "OH_star_legacy_eV", "Ha_star_legacy_eV", "OOH_star_legacy_eV"],
        "B_OH_Ha_OHHa": ["deltaG2_eV", "OH_star_legacy_eV", "Ha_star_legacy_eV", "OHHa_strict_eV"],
    }
    rows = []
    for subset_name, subset in [("13site", df), ("11site", df[~df.surface_name.isin(SPECIAL)].copy())]:
        for family, terms in specs.items():
            for n_terms in range(1, len(terms)+1):
                included = terms[:n_terms]
                mat = np.column_stack([subset[t].to_numpy(float) for t in included] + [np.ones(len(subset))])
                coefs, _, rank, _ = np.linalg.lstsq(mat, subset.Y_TS1_minus_Rpair_eV.to_numpy(float), rcond=None)
                pred = mat @ coefs
                met = metric(subset.Y_TS1_minus_Rpair_eV.to_numpy(float), pred)
                rows.append({"subset": subset_name, "family": family, "added_term": included[-1], "terms": ";".join(included), "n_terms": len(included), "rank": rank, **met})
    return pd.DataFrame(rows)


def contribution_diagnostics(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    correlations = []
    names_map = {
        "model_A_OH_Ha_OOH": ["deltaG2_eV", "OH_star_legacy_eV", "Ha_star_legacy_eV", "OOH_star_legacy_eV"],
        "model_B_OH_Ha_OHHa": ["deltaG2_eV", "OH_star_legacy_eV", "Ha_star_legacy_eV", "OHHa_strict_eV"],
    }
    base = fit(df, "traditional_BEP")[0][0]
    for model, predictors in names_map.items():
        mat = np.column_stack([df[predictors].to_numpy(float), np.ones(len(df))])
        y = df.Y_TS1_minus_Rpair_eV.to_numpy(float)
        coef, _, rank, _ = np.linalg.lstsq(mat, y, rcond=None)
        full_r2 = metric(y, mat@coef)["R2"]
        for removed in [None, *predictors]:
            keep = [p for p in predictors if p != removed]
            m = np.column_stack([df[keep].to_numpy(float), np.ones(len(df))])
            c, _, _, _ = np.linalg.lstsq(m, y, rcond=None)
            r2 = metric(y, m@c)["R2"]
            rows.append({"model": model, "diagnostic": "full" if removed is None else f"drop_{removed}", "removed_term": removed or "none", "R2": r2, "delta_R2_vs_full": r2-full_r2, "full_model_R2": full_r2})
        # Hold traditional BEP slope fixed while refitting intercept and correction terms.
        correction = predictors[1:]
        fixed_target = y - base * df.deltaG2_eV.to_numpy(float)
        m = np.column_stack([df[correction].to_numpy(float), np.ones(len(df))])
        c, _, _, _ = np.linalg.lstsq(m, fixed_target, rcond=None)
        pred = base*df.deltaG2_eV.to_numpy(float) + m@c
        rows.append({"model": model, "diagnostic": "traditional_slope_fixed", "removed_term": "none", "R2": metric(y,pred)["R2"], "delta_R2_vs_full": metric(y,pred)["R2"]-full_r2, "full_model_R2": full_r2})
        # Partial explanatory increment of each additional term, conditional on all other listed terms.
        for term in predictors:
            reduced = [p for p in predictors if p != term]
            mr = np.column_stack([df[reduced].to_numpy(float), np.ones(len(df))])
            cr, _, _, _ = np.linalg.lstsq(mr, y, rcond=None)
            rows.append({"model": model, "diagnostic": f"partial_increment_{term}", "removed_term": term, "R2": full_r2-metric(y,mr@cr)["R2"], "delta_R2_vs_full": np.nan, "full_model_R2": full_r2})
        # VIF and correlations among physical input columns; unstandardized condition number is reported with caution.
        corr = df[predictors].corr()
        for i, a in enumerate(predictors):
            for b in predictors[i+1:]:
                correlations.append({"model": model, "predictor_1": a, "predictor_2": b, "Pearson_r": corr.loc[a,b]})
        xmat = df[predictors].to_numpy(float)
        standardized = (xmat-xmat.mean(axis=0))/np.where(xmat.std(axis=0)==0,1,xmat.std(axis=0))
        cond = float(np.linalg.cond(standardized))
        for i, term in enumerate(predictors):
            other = np.delete(standardized, i, axis=1)
            target = standardized[:,i]
            pred_i = np.column_stack([other,np.ones(len(other))]) @ np.linalg.lstsq(np.column_stack([other,np.ones(len(other))]), target, rcond=None)[0]
            r2_i = metric(target,pred_i)["R2"]
            correlations.append({"model": model, "predictor_1": term, "predictor_2": "VIF", "Pearson_r": 1/(1-r2_i) if r2_i < 1 else np.inf, "standardized_design_condition_number": cond})
    return pd.DataFrame(rows), pd.DataFrame(correlations)


def effective_slope(df: pd.DataFrame, coef_df: pd.DataFrame) -> pd.DataFrame:
    c = coef_df[(coef_df.subset=="13site_refit") & (coef_df.model=="variable_slope_OH")].set_index("term").coefficient
    alpha0 = float(c["deltaG2_eV"])
    alpha1 = float(c["deltaG2_x_OH_centered"])
    out = df[["surface_name","catalyst","OH_star_legacy_eV","OH_centered_eV","deltaG2_eV","Y_TS1_minus_Rpair_eV"]].copy()
    out["alpha0"] = alpha0
    out["alpha1"] = alpha1
    out["alpha_eff"] = alpha0 + alpha1*out.OH_centered_eV
    return out


def leave_one_prediction_sensitivity(df: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for model in MODELS:
        full_coef,_,_=fit(df,model)
        full_matrix,_=design(df,model)
        full_pred=full_matrix@full_coef
        for held in df.surface_name.tolist():
            train=df[df.surface_name.ne(held)]
            coef,_,_=fit(train,model)
            pred=design(df,model)[0]@coef
            shift=pred-full_pred
            rows.append({"model":model,"held_out_catalyst":held,"rms_prediction_shift_all_13_eV":float(np.sqrt(np.mean(shift**2))),"max_abs_prediction_shift_all_13_eV":float(np.max(np.abs(shift)))})
    return pd.DataFrame(rows)


def ts_reconstruction(df: pd.DataFrame, residual_df: pd.DataFrame) -> pd.DataFrame:
    r = residual_df[(residual_df.subset=="13site_refit")].copy()
    r["actual_TS1_formation_energy_eV"] = r.surface_name.map(df.set_index("surface_name").TS1_strict.astype(float))
    r["Rpair_strict_eV"] = r.surface_name.map(df.set_index("surface_name").Rpair_strict.astype(float))
    r["predicted_TS1_formation_energy_eV"] = r.Rpair_strict_eV + r.predicted_Y_eV
    r["TS1_reconstruction_residual_eV"] = r.actual_TS1_formation_energy_eV-r.predicted_TS1_formation_energy_eV
    r["LOOCV_predicted_TS1_formation_energy_eV"] = r.Rpair_strict_eV + r.LOOCV_predicted_Y_eV
    return r


def savefig(name: str, fig: plt.Figure) -> None:
    fig.tight_layout()
    fig.savefig(FIG/f"{name}.png", bbox_inches="tight")
    fig.savefig(FIG/f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def annotate_points(ax, df: pd.DataFrame, x: np.ndarray, y: np.ndarray, *, size=26, alpha=0.92) -> None:
    for i, row in df.reset_index(drop=True).iterrows():
        special = row.catalyst in SPECIAL or row.surface_name in SPECIAL
        color = "#D62728" if special else "#1B4965"
        marker = "X" if row.surface_name == "tiw" else "P" if row.surface_name == "tita" else "o"
        ax.scatter(x[i], y[i], s=size+9 if special else size, color=color, marker=marker, edgecolor="white", linewidth=0.6, zorder=4, alpha=alpha)
        ax.annotate(row.catalyst, (x[i],y[i]), xytext=LABEL_OFFSETS.get(row.catalyst,(4,4)), textcoords="offset points", fontsize=7, color="#222222", zorder=5)


def make_figures(df: pd.DataFrame, metrics_df: pd.DataFrame, coef_df: pd.DataFrame, residual_df: pd.DataFrame, ablation_df: pd.DataFrame, slope_df: pd.DataFrame) -> pd.DataFrame:
    curve_rows = []
    # 1. Traditional BEP.
    b, _, _ = fit(df,"traditional_BEP")
    xx = np.linspace(df.deltaG2_eV.min()-0.08,df.deltaG2_eV.max()+0.08,240)
    yy = b[0]*xx+b[1]
    fig,ax=plt.subplots(figsize=(7.3,5.2))
    for i,row in df.iterrows():
        marker = "X" if row.surface_name=="tiw" else "P" if row.surface_name=="tita" else "o"
        c="#D62728" if row.surface_name in SPECIAL else "#1B4965"
        ax.scatter(row.deltaG2_eV,row.Y_TS1_minus_Rpair_eV,marker=marker,color=c,s=58 if row.surface_name in SPECIAL else 40,edgecolor="white",linewidth=.6,zorder=3)
        ax.annotate(row.catalyst,(row.deltaG2_eV,row.Y_TS1_minus_Rpair_eV),xytext=LABEL_OFFSETS.get(row.catalyst,(4,4)),textcoords="offset points",fontsize=7)
    legend_handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor="#1B4965",markeredgecolor="white",label="Other catalysts",markersize=7),Line2D([0],[0],marker="X",color="none",markerfacecolor="#D62728",markeredgecolor="white",label="Ti-W",markersize=8),Line2D([0],[0],marker="P",color="none",markerfacecolor="#D62728",markeredgecolor="white",label="Ti-Ta",markersize=8),Line2D([0],[0],color=COLORS["traditional_BEP"],lw=2,label=f"OLS: Y={b[0]:.3f} ΔG₂ {b[1]:+.3f}; R²={metric(df.Y_TS1_minus_Rpair_eV.to_numpy(),fit(df,'traditional_BEP')[2])['R2']:.3f}")]
    ax.axhline(0,color="#777",lw=.7); ax.set(xlabel="NEB1 reaction energy, ΔG₂ = G(EpoxideHHO*) − G(Rpair*) (eV)",ylabel="Signed NEB1 difference, Y = G(TS1*) − G(Rpair*) (eV)",title="Traditional BEP baseline (13 catalysts)"); ax.legend(handles=legend_handles,frameon=False,ncol=2)
    savefig("01_traditional_BEP",fig)
    for xval,yval in zip(xx,yy): curve_rows.append({"plot":"traditional_BEP","x_deltaG2_eV":xval,"fitted_Y_eV":yval})
    # 2/3. Actual vs predicted A and B.
    for idx,model in enumerate(["model_A_OH_Ha_OOH","model_B_OH_Ha_OHHa"],start=2):
        part=residual_df[(residual_df.subset=="13site_refit")&(residual_df.model==model)]
        fig,ax=plt.subplots(figsize=(6.4,5.5))
        ax.scatter(part.actual_Y_eV,part.predicted_Y_eV,color=COLORS[model],s=42,edgecolor="white",linewidth=.6,zorder=3)
        for _,r in part.iterrows(): ax.annotate(r.catalyst,(r.actual_Y_eV,r.predicted_Y_eV),xytext=LABEL_OFFSETS.get(r.catalyst,(4,4)),textcoords="offset points",fontsize=7)
        lo=min(part.actual_Y_eV.min(),part.predicted_Y_eV.min())-.1; hi=max(part.actual_Y_eV.max(),part.predicted_Y_eV.max())+.1
        ax.plot([lo,hi],[lo,hi],"--",color="#555",lw=1,label="y = x")
        met=metrics_df[(metrics_df.subset=="13site_refit")&(metrics_df.model==model)].iloc[0]
        ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel="Observed signed TS1−Rpair (eV)",ylabel="Predicted signed TS1−Rpair (eV)",title=MODEL_LABELS[model]+f"\nR²={met.R2:.3f}; RMSE={met.RMSE_eV:.3f} eV; LOOCV RMSE={met.LOOCV_RMSE_eV:.3f} eV")
        ax.legend(frameon=False); savefig(f"0{idx}_{'model_A_parity' if idx==2 else 'model_B_parity'}",fig)
    # 4. Residual comparison.
    fig,ax=plt.subplots(figsize=(9.2,5.3))
    subset=residual_df[residual_df.subset=="13site_refit"]
    names=df.catalyst.tolist(); xpos=np.arange(len(names)); width=.19
    for j,model in enumerate(MODELS):
        r=subset[subset.model==model].set_index("catalyst").reindex(names)
        ax.bar(xpos+(j-1.5)*width,r.residual_actual_minus_predicted_eV,width,label=MODEL_LABELS[model],color=COLORS[model])
    ax.axhline(0,color="#333",lw=.8); ax.set_xticks(xpos,names,rotation=35,ha="right"); ax.set_ylabel("Residual: observed − fitted Y (eV)"); ax.set_title("In-sample residuals by catalyst (signed NEB1 difference)"); ax.legend(ncol=2,frameon=False)
    savefig("04_residual_comparison",fig)
    # 5. OH-dependent slope diagnostic, with actual and model Y-X relation.
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(11,4.8))
    for _,r in slope_df.iterrows():
        special=r.surface_name in SPECIAL
        marker="X" if r.surface_name=="tiw" else "P" if r.surface_name=="tita" else "o"
        ax1.scatter(r.OH_star_legacy_eV,r.alpha_eff,marker=marker,color="#D62728" if special else COLORS["variable_slope_OH"],s=37,zorder=3)
        ax1.annotate(r.catalyst,(r.OH_star_legacy_eV,r.alpha_eff),xytext=LABEL_OFFSETS.get(r.catalyst,(3,3)),textcoords="offset points",fontsize=6.5)
    ax1.plot(slope_df.sort_values("OH_star_legacy_eV").OH_star_legacy_eV,slope_df.sort_values("OH_star_legacy_eV").alpha_eff,color=COLORS["variable_slope_OH"],lw=1.2,zorder=2)
    ax1.set(xlabel="Historical OH* energy (eV)",ylabel="Effective slope α_eff",title="OH-dependent BEP slope")
    sc=None
    for _,r in df.iterrows():
        marker="X" if r.surface_name=="tiw" else "P" if r.surface_name=="tita" else "o"
        sc=ax2.scatter(r.deltaG2_eV,r.Y_TS1_minus_Rpair_eV,c=[r.OH_star_legacy_eV],cmap="viridis",vmin=df.OH_star_legacy_eV.min(),vmax=df.OH_star_legacy_eV.max(),marker=marker,s=49,edgecolor="white",linewidth=.5,zorder=3)
    # Fitted values across the measured OH-energy range for selected ΔG2 quantiles.
    co=coef_df[(coef_df.subset=="13site_refit")&(coef_df.model=="variable_slope_OH")].set_index("term").coefficient
    xgrid=np.linspace(df.deltaG2_eV.min()-.05,df.deltaG2_eV.max()+.05,150)
    for zval,ls in [(df.OH_centered_eV.quantile(.2),"--"),(0,"-"),(df.OH_centered_eV.quantile(.8),":")]:
        zgrid=np.full_like(xgrid,zval)
        ygrid=co["deltaG2_eV"]*xgrid+co["deltaG2_x_OH_centered"]*xgrid*zgrid+co["OH_centered_eV"]*zgrid+co["intercept_eV"]
        ax2.plot(xgrid,ygrid,ls=ls,lw=1.6,label=f"OH centered={zval:.2f} eV")
    for _,r in df.iterrows(): ax2.annotate(r.catalyst,(r.deltaG2_eV,r.Y_TS1_minus_Rpair_eV),xytext=RIGHT_PANEL_OFFSETS.get(r.catalyst,(3,3)),textcoords="offset points",fontsize=6.2)
    ax2.set(xlabel="ΔG₂ (eV)",ylabel="Signed TS1−Rpair (eV)",title="Observed points and variable-slope model"); ax2.legend(frameon=False,fontsize=7)
    fig.colorbar(sc,ax=ax2,label="Historical OH* energy (eV)"); savefig("05_variable_slope_mechanism",fig)
    for zval in np.linspace(df.OH_centered_eV.min(),df.OH_centered_eV.max(),100):
        curve_rows.append({"plot":"variable_slope_alpha_eff","OH_centered_eV":zval,"alpha_eff":co["deltaG2_eV"]+co["deltaG2_x_OH_centered"]*zval})
    # 6. Fit/LOOCV performance across all model families and both refit subsets.
    p=metrics_df[metrics_df.subset.isin(["13site_refit","11site_refit"])].copy()
    fig,axes=plt.subplots(1,3,figsize=(13,4.5))
    for ax,metric_name,title in zip(axes,["R2","RMSE_eV","LOOCV_RMSE_eV"],["In-sample R²","In-sample RMSE (eV)","LOOCV RMSE (eV)"]):
        x=np.arange(len(MODELS)); w=.36
        for j,(sub,color) in enumerate([("13site_refit","#4C78A8"),("11site_refit","#F58518")]):
            v=p[p.subset==sub].set_index("model").reindex(MODELS)[metric_name].to_numpy(float)
            ax.bar(x+(j-.5)*w,v,w,label=sub.replace("_refit"," refit"),color=color)
        ax.set_xticks(x,["BEP","A","B","Variable\nslope"]); ax.set_title(title); ax.grid(axis="y",alpha=.2)
    axes[0].legend(frameon=False); fig.suptitle("Separately refitted subsets; 11-site group excludes Ti-W and Ti-Ta",y=1.03); savefig("06_model_performance",fig)
    # 7. Ablation by nested additions.
    fig,ax=plt.subplots(figsize=(9,5.1))
    for subset_name,color,marker in [("13site","#4C78A8","o"),("11site","#F58518","s")]:
        for family,ls in [("A_OH_Ha_OOH","-"),("B_OH_Ha_OHHa","--")]:
            q=ablation_df[(ablation_df.subset==subset_name)&(ablation_df.family==family)]
            label=f"{family.split('_')[0]} / {subset_name}"
            ax.plot(np.arange(1,len(q)+1),q.R2,marker=marker,ls=ls,color=color,label=label)
            for _,r in q.iterrows(): curve_rows.append({"plot":"ablation","subset":subset_name,"family":family,"n_terms":r.n_terms,"R2":r.R2,"added_term":r.added_term})
    ax.set_xticks([1,2,3,4],["ΔG₂","+ OH","+ Ha","+ OOH or OHHa"]); ax.set_ylabel("In-sample R²"); ax.set_title("Nested variable ablation (each subset refitted)"); ax.grid(alpha=.2); ax.legend(ncol=2,frameon=False)
    savefig("07_variable_ablation",fig)
    return pd.DataFrame(curve_rows)


def main() -> None:
    for folder in [FIG,EXP,REPORTS]: folder.mkdir(parents=True,exist_ok=True)
    df=get_data()
    metrics_df,coef_df,residual_df,cv_df,influence_df=records_and_diagnostics(df)
    ablation_df=ablations(df)
    contribution_df,corr_df=contribution_diagnostics(df)
    slope_df=effective_slope(df,coef_df)
    prediction_sensitivity_df=leave_one_prediction_sensitivity(df)
    ts_df=ts_reconstruction(df,residual_df)
    # Core tables.
    metrics_df.to_csv(DATA/"BEP_MODEL_COMPARISON.csv",index=False,float_format="%.12g")
    coef_df.to_csv(DATA/"BEP_MODEL_COEFFICIENTS.csv",index=False,float_format="%.12g")
    metrics_df[metrics_df.subset.isin(["13site_refit","11site_refit","13site_coefficients_applied_to_11"])].to_csv(DATA/"BEP_13SITE_11SITE_COMPARISON.csv",index=False,float_format="%.12g")
    ablation_df.to_csv(DATA/"BEP_VARIABLE_ABLATION.csv",index=False,float_format="%.12g")
    residual_df.to_csv(DATA/"BEP_CATALYST_RESIDUALS.csv",index=False,float_format="%.12g")
    slope_df.to_csv(DATA/"BEP_EFFECTIVE_SLOPE.csv",index=False,float_format="%.12g")
    ts_df.to_csv(DATA/"BEP_TS1_RECONSTRUCTION.csv",index=False,float_format="%.12g")
    contribution_df.to_csv(DATA/"BEP_CONTRIBUTION_DIAGNOSTICS.csv",index=False,float_format="%.12g")
    corr_df.to_csv(DATA/"BEP_COLLINEARITY_DIAGNOSTICS.csv",index=False,float_format="%.12g")
    influence_df.to_csv(DATA/"BEP_LEAVE_ONE_COEFFICIENT_SENSITIVITY.csv",index=False,float_format="%.12g")
    prediction_sensitivity_df.to_csv(DATA/"BEP_LEAVE_ONE_PREDICTION_SENSITIVITY.csv",index=False,float_format="%.12g")
    cv_df.to_csv(DATA/"BEP_LOOCV_SUMMARY.csv",index=False,float_format="%.12g")
    # Origin exports: joined raw data, fitted results, response curves, model coefficients and diagnostics.
    df.to_csv(EXP/"BEP_source_points_13catalysts.csv",index=False,float_format="%.12g")
    residual_df.to_csv(EXP/"BEP_all_model_predictions_residuals.csv",index=False,float_format="%.12g")
    ts_df.to_csv(EXP/"BEP_TS1_reconstruction.csv",index=False,float_format="%.12g")
    coef_df.to_csv(EXP/"BEP_model_coefficients.csv",index=False,float_format="%.12g")
    ablation_df.to_csv(EXP/"BEP_variable_ablation.csv",index=False,float_format="%.12g")
    slope_df.to_csv(EXP/"BEP_effective_slope.csv",index=False,float_format="%.12g")
    prediction_sensitivity_df.to_csv(EXP/"BEP_leave_one_prediction_sensitivity.csv",index=False,float_format="%.12g")
    curve_df=make_figures(df,metrics_df,coef_df,residual_df,ablation_df,slope_df)
    curve_df.to_csv(EXP/"BEP_fitted_curves_and_ablation.csv",index=False,float_format="%.12g")
    # Report is deliberately generated from the saved estimates.
    m13=metrics_df[metrics_df.subset=="13site_refit"].set_index("model")
    m11=metrics_df[metrics_df.subset=="11site_refit"].set_index("model")
    mtrans=metrics_df[metrics_df.subset=="13site_coefficients_applied_to_11"].set_index("model")
    coeff_lookup={(r.subset,r.model,r.term):r.coefficient for r in coef_df.itertuples()}
    def formula(model, subset="13site_refit"):
        pairs=[(term,coeff_lookup[(subset,model,term)]) for term in coef_df[(coef_df.subset==subset)&(coef_df.model==model)].term]
        beta = next(v for t,v in pairs if t == "intercept_eV")
        terms = [f"({v:+.6g})·{t}" for t,v in pairs if t != "intercept_eV"]
        return " + ".join(terms) + f" + ({beta:+.6g} eV)"
    top_res= residual_df[(residual_df.subset=="13site_refit")&(residual_df.model.isin(MODELS))].copy()
    worst=top_res.loc[top_res.groupby("model").residual_actual_minus_predicted_eV.apply(lambda s:s.abs().idxmax())]
    alpha0=coeff_lookup[("13site_refit","variable_slope_OH","deltaG2_eV")]
    alpha1=coeff_lookup[("13site_refit","variable_slope_OH","deltaG2_x_OH_centered")]
    influence_summary=[]
    for model in MODELS:
        top=prediction_sensitivity_df[prediction_sensitivity_df.model==model].sort_values("rms_prediction_shift_all_13_eV",ascending=False).iloc[0]
        influence_summary.append(f"{MODEL_LABELS[model]}: {top.held_out_catalyst}（RMS预测变化 {top.rms_prediction_shift_all_13_eV:.3f} eV）")
    report=[
        "# Ti–M氧转移NEB1扩展BEP与能量修正分析",
        "",
        "## 研究边界与数据",
        "本分析仅使用已归档的13个位点严格CatMAP形成能、历史单吸附态OH*/Ha*/OOH*数值及严格模型OHHa态；未运行DFT、CatMAP稳态/TOF或条件扫描。正式目标为有符号 `Y = G(TS1*) − G(Rpair*)`，传统BEP自变量为 `ΔG₂ = G(EpoxideHHO*) − G(Rpair*)`。负的TS1−Rpair值保留。Ti–Zr不纳入。历史OH*/Ha*/OOH*参考的slab闭合尚未完成，相关扩展系数均属探索性统计关系，不解释为已验证因果效应。",
        "",
        "样本按surface_name一对一匹配；13个正式位点为："+", ".join(df.catalyst)+"。11位点对照仅排除Ti-W、Ti-Ta，未替代正式13位点结果。",
        "",
        "## 模型与复现指标",
        "OLS系数按eV输入与eV输出；LOOCV每次仅用训练折重新估计系数。11-site refit与13-site系数用于11点是两类不同评估。",
        "",
        "|模型|13-site R²|RMSE (eV)|MAE (eV)|LOOCV RMSE (eV)|11-site refit R²|11-site refit RMSE|13-site系数用于11点 R²|",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for model in MODELS:
        a=m13.loc[model]; b=m11.loc[model]; c=mtrans.loc[model]
        report.append(f"|{MODEL_LABELS[model]}|{a.R2:.4f}|{a.RMSE_eV:.4f}|{a.MAE_eV:.4f}|{a.LOOCV_RMSE_eV:.4f}|{b.R2:.4f}|{b.RMSE_eV:.4f}|{c.R2:.4f}|")
    report += ["", "### 13位点拟合公式", ""]
    for model in MODELS: report.append(f"- **{MODEL_LABELS[model]}**：`Y = {formula(model)}`")
    report += [
        "",
        "## 结果要点",
        f"13位点样本内最高R²为OH变斜率模型（{m13.loc['variable_slope_OH'].R2:.4f}），其LOOCV RMSE={m13.loc['variable_slope_OH'].LOOCV_RMSE_eV:.3f} eV，较传统BEP低；但11位点重拟合R²仅{m11.loc['variable_slope_OH'].R2:.3f}、LOOCV RMSE={m11.loc['variable_slope_OH'].LOOCV_RMSE_eV:.3f} eV，显示模型对样本构成敏感。方案B在11位点的LOOCV RMSE最低（{m11.loc['model_B_OH_Ha_OHHa'].LOOCV_RMSE_eV:.3f} eV），但将13点系数直接应用于11点时R²为{mtrans.loc['model_B_OH_Ha_OHHa'].R2:.3f}，说明其跨子集迁移也没有得到支持。",
        "嵌套消融显示OH与Ha加入后明显改善样本内拟合，A中的OOH追加贡献近零；B中OHHa有一定增量。完整模型中ΔG₂系数反号且删去ΔG₂后的R²损失较小，因此改善主要不是传统正斜率BEP项得到强化，而是能量修正项/交互结构对当前样本的解释增加。历史OH/Ha/OOH的参考不确定性、变量共线性和小样本仍限制物理解读。",
        "",
        "## 变斜率定义",
        f"OH中心值为 `{df.OH_star_legacy_eV.mean():.8g} eV`。模型为 `Y = α₀ΔG₂ + α₁ΔG₂·Z + γZ + β`，其中 `Z = E_OH − mean(E_OH)`；13位点 `α₀={alpha0:.6g}`，`α₁={alpha1:.6g}`，有效斜率 `α_eff=α₀+α₁Z` 的样本范围为 `{slope_df.alpha_eff.min():.6g}` 至 `{slope_df.alpha_eff.max():.6g}`。此为样本内交互模型，不等同已验证的机制斜率定律。",
        "",
        "## 贡献与稳健性诊断",
        "完整项与逐项删除、ΔG₂删除、传统BEP斜率固定后的结果见 `BEP_CONTRIBUTION_DIAGNOSTICS.csv`；变量相关、VIF与标准化设计矩阵条件数见 `BEP_COLLINEARITY_DIAGNOSTICS.csv`；逐一留出一个催化剂重拟合的系数变化见 `BEP_LEAVE_ONE_COEFFICIENT_SENSITIVITY.csv`。这些量用于诊断，不将多变量OLS系数当作独立因果贡献。",
        f"按删除项后的R²差衡量，方案A去掉ΔG₂仅下降{float(contribution_df[(contribution_df.model=='model_A_OH_Ha_OOH')&(contribution_df.diagnostic=='drop_deltaG2_eV')].delta_R2_vs_full.iloc[0]):.4f}；方案B去掉ΔG₂下降{float(contribution_df[(contribution_df.model=='model_B_OH_Ha_OHHa')&(contribution_df.diagnostic=='drop_deltaG2_eV')].delta_R2_vs_full.iloc[0]):.4f}。相反，把ΔG₂斜率固定为传统BEP的13点斜率后，A/B的R²分别为{float(contribution_df[(contribution_df.model=='model_A_OH_Ha_OOH')&(contribution_df.diagnostic=='traditional_slope_fixed')].R2.iloc[0]):.4f}/{float(contribution_df[(contribution_df.model=='model_B_OH_Ha_OHHa')&(contribution_df.diagnostic=='traditional_slope_fixed')].R2.iloc[0]):.4f}。A中OOH项的条件增量R²仅约{float(contribution_df[(contribution_df.model=='model_A_OH_Ha_OOH')&(contribution_df.diagnostic=='partial_increment_OOH_star_legacy_eV')].R2.iloc[0]):.5f}；因此不能把A的提升归因于OOH单项。",
        "",
        "最大绝对样本内残差（按模型）："+"；".join(f"{r.model}: {r.catalyst} ({r.residual_actual_minus_predicted_eV:+.4f} eV)" for r in worst.itertuples()),
        "",
        "Ti–W、Ti–Ta、Ti–Co、Ti–Fe的单点残差在 `BEP_CATALYST_RESIDUALS.csv` 保留。Ti–W/Ti–Ta的有符号NEB1差值为负值，不做截断。",
        "",
        "|催化剂|传统BEP残差|方案A残差|方案B残差|OH变斜率残差|",
        "|---|---:|---:|---:|---:|",
        *["|"+r.catalyst+"|"+"|".join(f"{float(residual_df[(residual_df.subset=='13site_refit')&(residual_df.model==model)&(residual_df.catalyst==r.catalyst)].residual_actual_minus_predicted_eV.iloc[0]):+.4f}" for model in MODELS)+"|" for r in df[df.surface_name.isin(FOCUS)].itertuples()],
        "逐一移除一个位点后重拟合，并比较该系数集对全部13位点的预测变化，影响最大的留出位点分别为："+"；".join(influence_summary)+"。这是小样本影响诊断，不等价于删除这些点。",
        "",
        "### 逐项消融（13-site）",
        "",
        "|模型序列|新增项|R²|RMSE (eV)|",
        "|---|---|---:|---:|",
    ]
    for r in ablation_df[ablation_df.subset=="13site"].itertuples(): report.append(f"|{r.family}|{r.added_term}|{r.R2:.4f}|{r.RMSE_eV:.4f}|")
    report += [
        "",
        "## TS1重建",
        "每种模型以 `G(TS1)_pred = G(Rpair)_actual + Y_pred` 重建。`BEP_TS1_RECONSTRUCTION.csv` 同时包含样本内与LOOCV预测及残差。这里仅诊断能量重建，不重算TOF。",
        "",
        "## 图件与可复现文件",
        "7类图均输出PNG/PDF：传统BEP、方案A预测一致性、方案B预测一致性、四模型残差、OH调制斜率机制图、13/11拟合性能对照、变量消融。图中保留13个位点名称，Ti–W/Ti–Ta使用不同标记。Origin所需原始点、拟合值、残差、系数和曲线数据位于 `origin_exports/`。",
        "",
        "复现：在分析包根目录运行 `python scripts/run_bep_energy_correction.py`。脚本只读取 `data/` 快照并写入本包，不调用CatMAP运行接口。",
        "",
        "## 解释限制",
        "训练R²提高不保证外推；13个样本而多变量项数较多，LOOCV及逐位点系数变化需一并考虑。历史单吸附态数据与严格共吸附模型能量的slab/热修正参考尚未闭合，因此方案A中的OH*/Ha*/OOH*及交叉比较具有参考不确定性。方案B使用严格能量表OHHa态，但仍是跨反应网络能量项的统计耦合，不能据此断言OHHa直接控制NEB1。",
        "",
        "## 参考与软件",
        "数值拟合使用本机NumPy OLS（`numpy.linalg.lstsq`）、逐点留一重拟合及Matplotlib；CatMAP官方结果仅作为只读数据源，本分析未调用 `ReactionModel.run()`。",
        "图件设计采用原始数据与派生分析分离、图例/符号双重编码、保留全部点与有符号值等原则；参考：Kassis T, Agarwal V, He Y, Patel D, Brueckner AM. Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065v2 (2026). https://arxiv.org/abs/2609.00065",
    ]
    (REPORTS/"BEP_ENERGY_CORRECTION_REPORT.md").write_text("\n".join(report)+"\n",encoding="utf-8")
    # Hash all package files except the manifest itself.
    manifest=[]
    for p in sorted(PKG.rglob("*")):
        if p.is_file() and p.name!="SHA256SUMS.csv":
            manifest.append({"relative_path":str(p.relative_to(PKG)).replace("\\","/"),"size_bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()})
    pd.DataFrame(manifest).to_csv(REPORTS/"SHA256SUMS.csv",index=False)
    print(json.dumps({"n":len(df),"models_13":metrics_df[metrics_df.subset=="13site_refit"].to_dict(orient="records"),"models_11":metrics_df[metrics_df.subset=="11site_refit"].to_dict(orient="records"),"out":str(PKG)},ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
