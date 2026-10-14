from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import ConnectionPatch


PROJECT = Path(r"H:\catmap\catmap\catmap_code_collection_20261013")
LEGACY = PROJECT / "analyses/C_O_H_to_TS1_scaling/data/source_DFT_formation_catmap.txt"
RAW_LEGACY = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\DFT_formation_catmap.txt")
CONVERTER = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\1-generating_input_file\convert_dft_to_catmap.py")
STRICT = PROJECT / "catmap_official_run/energies.txt"
MATCHED = PROJECT / "analyses/empirical_O_H_descriptor_analysis/data/matched_data.csv"
SLAB_AUDIT = PROJECT / "analyses/C_O_H_to_TS1_scaling/data/source_reference_slab_comparison.csv"
WORKBOOK = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\己烯环氧化.xlsx")
OUT = Path(__file__).resolve().parents[1]

SURFACES = ["tife", "timn", "tihf", "tire", "timo", "tiv", "tico",
            "titi", "tiw", "tita", "ticr", "ti", "tini"]
DISPLAY = {
    "tife": "Ti-Fe", "timn": "Ti-Mn", "tihf": "Ti-Hf", "tire": "Ti-Re",
    "timo": "Ti-Mo", "tiv": "Ti-V", "tico": "Ti-Co", "titi": "Ti-Ti",
    "tiw": "Ti-W", "tita": "Ti-Ta", "ticr": "Ti-Cr", "ti": "Ti",
    "tini": "Ti-Ni",
}
SCATTER_OFFSETS = {
    "Ti-Fe": (10, -15), "Ti-Mn": (10, -20), "Ti-Hf": (-24, 10), "Ti-Re": (10, 13),
    "Ti-Mo": (10, -14), "Ti-V": (10, -17), "Ti-Co": (-28, 12), "Ti-Ti": (-30, -16),
    "Ti-W": (10, 9), "Ti-Ta": (-29, 13), "Ti-Cr": (-28, 10), "Ti": (-16, -18), "Ti-Ni": (10, 8),
}
MAP_OFFSETS = {
    "tife": (10, 12), "timn": (10, -20), "tihf": (10, 9), "tire": (-36, 17),
    "timo": (11, -17), "tiv": (11, -23), "tico": (10, 10), "titi": (-36, -17),
    "tiw": (-34, 10), "tita": (10, 18), "ticr": (-36, 8), "ti": (10, -7), "tini": (10, 8),
}
WORKBOOK_LABELS = {"tife": "TiFe", "timn": "TiMn", "tihf": "TiHf", "tire": "TiRe", "timo": "TiMo",
                   "tiv": "TiV", "tico": "TiCo", "titi": "TiTi", "tiw": "TiW", "tita": "TiTa",
                   "ticr": "TiCr", "ti": "单核Ti", "tini": "TiNi", "tinb": "TiNb", "tizr": "TiZr"}
LEGACY_SPECIES = {
    "C6H12_star": "C6H12", "C6H12O_star": "C6H12O", "O_star": "O",
    "Ha_star": "Ha", "Hb_star": "Hb", "OOH_star": "OOH", "OH_star": "OH",
    "OOH_C6H12_candidate": "OOH-C6H12", "Hb_O_candidate": "Hb-O",
}
PRETTY = {
    "C6H12_star_legacy": "C6H12* (legacy formation energy)",
    "C6H12O_star_legacy": "C6H12O* (legacy formation energy)",
    "O_star_legacy": "O* (legacy formation energy)", "Ha_star_legacy": "Ha* (legacy formation energy)",
    "Hb_star_legacy": "Hb* (legacy formation energy)", "OOH_star_legacy": "OOH* (legacy formation energy)",
    "OH_star_legacy": "OH* (legacy formation energy)",
    "C6H12O_minus_C6H12_legacy": "C6H12O* - C6H12* (legacy state-energy difference)",
    "Rpair_strict": "Rpair* (current strict-model formation energy)",
}
CARBON_FEATURES = ["C6H12_star_legacy", "C6H12O_star_legacy", "C6H12O_minus_C6H12_legacy"]
O_H_FEATURES = ["O_star_legacy", "Ha_star_legacy", "Hb_star_legacy", "OOH_star_legacy", "OH_star_legacy"]


def mkdirs() -> None:
    for folder in ["data", "figures", "origin_exports", "reports"]:
        (OUT / folder).mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def unique_energy_map(df: pd.DataFrame, label: str) -> pd.DataFrame:
    if df.duplicated(["surface_name", "species_name"]).any():
        dup = df.loc[df.duplicated(["surface_name", "species_name"], keep=False), ["surface_name", "species_name"]]
        raise ValueError(f"{label}: duplicate surface/species rows: {dup.to_dict('records')}")
    return df


def fit_ols(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xm = x.mean(axis=0)
    ym = y.mean()
    slopes, *_ = np.linalg.lstsq(x - xm, y - ym, rcond=None)
    intercept = float(ym - xm @ slopes)
    pred = intercept + x @ slopes
    return np.r_[slopes, intercept], pred


def metrics(x: np.ndarray, y: np.ndarray) -> tuple[dict, np.ndarray, np.ndarray]:
    coef, pred = fit_ols(x, y)
    residual = pred - y
    rmse = float(np.sqrt(np.mean(residual**2)))
    mae = float(np.mean(np.abs(residual)))
    ss_tot = float(np.sum((y - y.mean())**2))
    r2 = float(1 - np.sum(residual**2) / ss_tot) if ss_tot else float("nan")
    loo = []
    for idx in range(len(y)):
        keep = np.arange(len(y)) != idx
        c, _ = fit_ols(x[keep], y[keep])
        loo.append(float(c[-1] + x[idx] @ c[:-1] - y[idx]))
    cv_rmse = float(np.sqrt(np.mean(np.square(loo))))
    return {"N": len(y), "R2_in_sample": r2, "RMSE_eV": rmse, "MAE_eV": mae,
            "Max_abs_error_eV": float(np.max(np.abs(residual))), "LOOCV_RMSE_eV": cv_rmse}, coef, pred


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    legacy_all = pd.read_csv(LEGACY, sep="\t", keep_default_na=False)
    strict_all = pd.read_csv(STRICT, sep="\t", keep_default_na=False)
    legacy = legacy_all.copy()
    strict = strict_all.copy()
    matched = pd.read_csv(MATCHED)
    slab = pd.read_csv(SLAB_AUDIT)
    for name, df in [("legacy", legacy), ("strict", strict)]:
        unique_energy_map(df, name)
    legacy = legacy[legacy.surface_name.isin(SURFACES)].copy()
    strict = strict[strict.surface_name.isin(SURFACES)].copy()
    if set(matched.surface_name) != set(SURFACES):
        raise ValueError("matched TOF surfaces do not equal the declared 13-site cohort")
    if len(matched) != 13 or matched.surface_name.duplicated().any():
        raise ValueError("TOF source is not 13 unique catalyst surfaces")

    rows = pd.DataFrame({"surface_name": SURFACES})
    for feature, species in LEGACY_SPECIES.items():
        data = legacy[(legacy.species_name == species) & legacy.surface_name.isin(SURFACES)]
        if len(data) != 13 or data.surface_name.duplicated().any():
            raise ValueError(f"legacy {species}: expected exactly 13 unique entries, found {len(data)}")
        rows = rows.merge(data[["surface_name", "formation_energy", "reference"]].rename(
            columns={"formation_energy": f"{feature}_legacy", "reference": f"{feature}_reference"}),
            on="surface_name", validate="one_to_one")

    gases = legacy_all[legacy_all.site_name == "gas"].set_index("species_name")
    old_gas = {sp: float(gases.loc[sp, "formation_energy"]) for sp in ["C6H12", "C6H12O"]}

    strict_pivot = strict.pivot(index="surface_name", columns="species_name", values="formation_energy")
    for state in ["Rpair", "TS1", "EpoxideHHO"]:
        if state not in strict_pivot or strict_pivot[state].isna().any():
            raise ValueError(f"current strict model missing {state}")
        rows[f"{state}_strict"] = rows.surface_name.map(strict_pivot[state])
    rows["NEB1_forward_energy_difference_TS1_minus_Rpair"] = rows.TS1_strict - rows.Rpair_strict
    rows["NEB1_reaction_energy_EpoxideHHO_minus_Rpair"] = rows.EpoxideHHO_strict - rows.Rpair_strict

    rows = rows.merge(slab[["surface_name", "strict_minus_legacy_slab_eV"]], on="surface_name", validate="one_to_one")
    for feature in LEGACY_SPECIES:
        col = f"{feature}_legacy"
        rows[f"{feature}_slabshift_conditional"] = rows[col] - rows.strict_minus_legacy_slab_eV
    rows["C6H12_adsorption_energy_legacy"] = rows.C6H12_star_legacy - old_gas["C6H12"]
    rows["C6H12O_adsorption_energy_legacy"] = rows.C6H12O_star_legacy - old_gas["C6H12O"]
    rows["C6H12O_minus_C6H12_legacy"] = rows.C6H12O_star_legacy - rows.C6H12_star_legacy
    rows["C6H12O_minus_C6H12_slabshift_conditional"] = rows.C6H12O_star_slabshift_conditional - rows.C6H12_star_slabshift_conditional
    rows = rows.merge(matched[["surface_name", "catalyst", "TOF", "log10_TOF", "status", "temperature_K"]], on="surface_name", validate="one_to_one")
    rows["catalyst"] = rows.surface_name.map(DISPLAY)
    if rows.isna().any().any():
        bad = rows.columns[rows.isna().any()].tolist()
        raise ValueError(f"missing values after keyed joins: {bad}")

    source_meta = {
        "old_gas_C6H12_formation_energy_eV": old_gas["C6H12"],
        "old_gas_C6H12O_formation_energy_eV": old_gas["C6H12O"],
        "current_model_gas_C6H12_formation_energy_eV": float(strict_all[(strict_all.site_name == "gas") & (strict_all.species_name == "C6H12")].formation_energy.iloc[0]),
        "current_model_gas_C6H12O_formation_energy_eV": float(strict_all[(strict_all.site_name == "gas") & (strict_all.species_name == "C6H12O")].formation_energy.iloc[0]),
        "source_hashes": {str(p): sha256(p) for p in [LEGACY, RAW_LEGACY, CONVERTER, STRICT, MATCHED, SLAB_AUDIT]},
        "raw_legacy_matches_archived_source_copy": sha256(RAW_LEGACY) == sha256(LEGACY),
    }
    return rows, strict, source_meta


def export_workbook_crosswalk(data: pd.DataFrame) -> dict:
    """Read-only source-cell provenance using the legacy converter's row definitions."""
    from openpyxl import load_workbook
    row_map = {"C6H12_star": ("catmap数据集", 20, "C6H12"),
               "C6H12O_star": ("catmap数据集", 24, "C6H12O"),
               "O_star": ("catmap数据集", 29, "O"), "Ha_star": ("catmap数据集", 8, "Ha"),
               "Hb_star": ("catmap数据集", 12, "Hb"), "OOH_star": ("catmap数据集", 16, "OOH"),
               "OH_star": ("catmap数据集", 34, "OH"),
               "OOH_C6H12_candidate": ("catmap数据集", 43, "OOH-C6H12"),
               "Hb_O_candidate": ("catmap数据集", 47, "Hb-O")}
    wb_f = load_workbook(WORKBOOK, read_only=True, data_only=False)
    wb_v = load_workbook(WORKBOOK, read_only=True, data_only=True)
    records = []
    try:
        ws_f, ws_v = wb_f["catmap数据集"], wb_v["catmap数据集"]
        starts = {}
        for col in range(1, ws_f.max_column + 1):
            val = ws_f.cell(1, col).value
            for surface, label in WORKBOOK_LABELS.items():
                if val == label:
                    starts[surface] = col
        if not set(SURFACES).issubset(set(starts)):
            raise ValueError(f"workbook catalyst headers mismatch: {sorted(starts)}")
        ordered = sorted((col, surf) for surf, col in starts.items())
        gcols = {}
        for i, (start, surf) in enumerate(ordered):
            stop = ordered[i+1][0] if i + 1 < len(ordered) else start + 9
            candidates = [col for col in range(start, stop) if ws_f.cell(3, col).value == "G"]
            if len(candidates) != 1:
                raise ValueError(f"{surf}: expected one G column in catalyst block, got {candidates}")
            gcols[surf] = candidates[0]
        from openpyxl.utils import get_column_letter
        formation = pd.read_csv(LEGACY, sep="\t", keep_default_na=False)
        lookup = {(r.surface_name, r.species_name): float(r.formation_energy) for r in formation.itertuples() if r.surface_name in SURFACES}
        for surface in SURFACES:
            col = gcols[surface]
            for feature, (sheet_name, row, species) in row_map.items():
                cell = f"{get_column_letter(col)}{row}"
                formula = ws_f[cell].value
                cached = ws_v[cell].value
                records.append({"surface_name": surface, "catalyst": DISPLAY[surface], "species_or_state": species,
                                "source_workbook": str(WORKBOOK), "source_workbook_sha256": sha256(WORKBOOK),
                                "source_sheet": sheet_name, "source_cell": cell, "source_formula": formula,
                                "cached_G_eV": cached, "legacy_formation_energy_eV": lookup.get((surface, species)),
                                "formation_energy_table_file": str(RAW_LEGACY), "formation_energy_table_sha256": sha256(RAW_LEGACY),
                                "converter_file": str(CONVERTER), "converter_sha256": sha256(CONVERTER),
                                "legacy_reference_tag": next((r.reference for r in formation.itertuples() if r.surface_name == surface and r.species_name == species), None)})
        # Summary sheet has a different column layout; independently reproduce the converter's per-sheet lookup.
        ws2f, ws2v = wb_f["吉布斯自由能汇总"], wb_v["吉布斯自由能汇总"]
        starts2 = {}
        for col in range(1, ws2f.max_column + 1):
            val = ws2f.cell(1, col).value
            for surface, label in WORKBOOK_LABELS.items():
                if val == label:
                    starts2[surface] = col
        if not set(WORKBOOK_LABELS).issubset(set(starts2)):
            raise ValueError("summary worksheet catalyst headers incomplete")
        ordered2 = sorted((col, surf) for surf, col in starts2.items())
        tife_col = None
        for i, (start, surf) in enumerate(ordered2):
            stop = ordered2[i+1][0] if i + 1 < len(ordered2) else start + 9
            cols = [col for col in range(start, stop) if ws2f.cell(3, col).value == "G"]
            if len(cols) != 1:
                raise ValueError(f"summary {surf}: expected one G column, got {cols}")
            if surf == "tife":
                tife_col = cols[0]
        if tife_col is None:
            raise ValueError("summary sheet TiFe Gibbs column was not found")
        for species, row in [("C6H12_g", 58), ("C6H12O_g", 60)]:
            cell = f"{get_column_letter(tife_col)}{row}"
            records.append({"surface_name": "gas_reference", "catalyst": "shared-gas", "species_or_state": species,
                            "source_workbook": str(WORKBOOK), "source_workbook_sha256": sha256(WORKBOOK),
                            "source_sheet": "吉布斯自由能汇总", "source_cell": cell, "source_formula": ws2f[cell].value,
                            "cached_G_eV": ws2v[cell].value, "legacy_formation_energy_eV": None,
                            "formation_energy_table_file": str(RAW_LEGACY), "formation_energy_table_sha256": sha256(RAW_LEGACY),
                            "converter_file": str(CONVERTER), "converter_sha256": sha256(CONVERTER),
                            "legacy_reference_tag": "DFT_Gibbs_gas; converter GAS_ROWS"})
    finally:
        wb_f.close(); wb_v.close()
    pd.DataFrame(records).to_csv(OUT / "data/source_workbook_crosswalk.csv", index=False, encoding="utf-8-sig")
    return {"workbook_path": str(WORKBOOK), "workbook_sha256": sha256(WORKBOOK), "crosswalk_rows": len(records),
            "gas_reference_source": "catmap数据集 and 吉布斯自由能汇总; cell rows follow converter_dft_to_catmap.py"}


def feature_values(df: pd.DataFrame, alignment: str) -> dict[str, np.ndarray]:
    suffix = "_legacy" if alignment == "legacy_reference" else "_slabshift_conditional"
    features = {
        "C6H12_star" + suffix: df.C6H12_star_legacy.to_numpy() if alignment == "legacy_reference" else df.C6H12_star_slabshift_conditional.to_numpy(),
        "C6H12O_star" + suffix: df.C6H12O_star_legacy.to_numpy() if alignment == "legacy_reference" else df.C6H12O_star_slabshift_conditional.to_numpy(),
        "O_star" + suffix: df.O_star_legacy.to_numpy() if alignment == "legacy_reference" else df.O_star_slabshift_conditional.to_numpy(),
        "Ha_star" + suffix: df.Ha_star_legacy.to_numpy() if alignment == "legacy_reference" else df.Ha_star_slabshift_conditional.to_numpy(),
        "Hb_star" + suffix: df.Hb_star_legacy.to_numpy() if alignment == "legacy_reference" else df.Hb_star_slabshift_conditional.to_numpy(),
        "OOH_star" + suffix: df.OOH_star_legacy.to_numpy() if alignment == "legacy_reference" else df.OOH_star_slabshift_conditional.to_numpy(),
        "OH_star" + suffix: df.OH_star_legacy.to_numpy() if alignment == "legacy_reference" else df.OH_star_slabshift_conditional.to_numpy(),
        "C6H12O_minus_C6H12" + suffix: df.C6H12O_minus_C6H12_legacy.to_numpy() if alignment == "legacy_reference" else df.C6H12O_minus_C6H12_slabshift_conditional.to_numpy(),
        # Network-derived input: physically meaningful initial state but directly enters ΔG‡ = TS1−Rpair.
        "Rpair_strict": df.Rpair_strict.to_numpy(),
    }
    return features


def build_candidate_fits(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ytargets = {
        "TS1_formation_energy": df.TS1_strict.to_numpy(),
        "NEB1_forward_energy_difference_TS1_minus_Rpair": df.NEB1_forward_energy_difference_TS1_minus_Rpair.to_numpy(),
        "NEB1_reaction_energy_EpoxideHHO_minus_Rpair": df.NEB1_reaction_energy_EpoxideHHO_minus_Rpair.to_numpy(),
    }
    metric_rows, coef_rows, pred_rows = [], [], []
    specs = []
    for alignment in ["legacy_reference", "slabshift_conditional"]:
        feat = feature_values(df, alignment)
        for key, values in feat.items():
            if key == "Rpair_strict":
                group, reference_note = "network_initial_state", "contains Rpair; interpret target coupling, especially for TS1−Rpair"
            elif key.startswith("C6H12O_minus_C6H12"):
                group, reference_note = "derived_carbon_state_difference", "slab and common atom-reference terms cancel algebraically; thermal/source consistency remains conditional"
            else:
                group, reference_note = ("carbon_adsorbate" if key.startswith(("C6H12_star", "C6H12O_star")) else "O_H_adsorbate"), "historical standalone-state formation energy; source reference not fully validated against current strict model"
            specs.append((key, values, alignment, group, reference_note))
        # Candidate TS-like states are preserved as leakage/overlap diagnostics only.
        diag = df.OOH_C6H12_candidate_legacy.to_numpy()
        specs.append(("OOH_C6H12_candidate_overlap_diagnostic", diag, alignment,
                      "transition_state_candidate_diagnostic", "excluded from primary ranking: source calls it a TS candidate; structural identity/independence from NEB1 TS1 is unverified"))

    # Primary descriptors must not duplicate by constant gas-phase shifts. The two reference tracks are reported separately.
    for target, y in ytargets.items():
        for alignment in ["legacy_reference", "slabshift_conditional"]:
            features = feature_values(df, alignment)
            eligible = [(k, v) for k, v in features.items()]
            # Single fits.
            for names in [(k,) for k, _ in eligible]:
                xs = np.column_stack([dict(eligible)[names[0]]])
                m, c, p = metrics(xs, y)
                info = next(s for s in specs if s[0] == names[0] and s[2] == alignment)
                status = "diagnostic_target_coupling" if "Rpair" in names[0] else "primary_candidate"
                metric_rows.append({"target": target, "alignment": alignment, "model": "single_linear", "features": " + ".join(names), **m,
                                    "candidate_group": info[3], "interpretation_flag": status, "reference_note": info[4]})
                coef_rows.append({"target": target, "alignment": alignment, "model": "single_linear", "features": " + ".join(names),
                                  "feature_coefficients": json.dumps({names[0]: float(c[0])}), "intercept": float(c[-1])})
                for i, (_, row) in enumerate(df.iterrows()):
                    pred_rows.append({"target": target, "alignment": alignment, "model": "single_linear", "features": " + ".join(names),
                                      "surface_name": row.surface_name, "catalyst": row.catalyst, "actual_eV": y[i], "predicted_eV": p[i], "residual_pred_minus_actual_eV": p[i]-y[i]})
            # All unique pairs in the nine-feature physical inventory. Aliases are not included twice.
            for (n1, x1), (n2, x2) in itertools.combinations(eligible, 2):
                names = (n1, n2)
                xs = np.column_stack([x1, x2])
                m, c, p = metrics(xs, y)
                coupled = any("Rpair" in name for name in names)
                groups = [next(s[3] for s in specs if s[0] == n and s[2] == alignment) for n in names]
                flag = "diagnostic_target_coupling" if coupled else "primary_candidate"
                note = "; ".join(next(s[4] for s in specs if s[0] == n and s[2] == alignment) for n in names)
                metric_rows.append({"target": target, "alignment": alignment, "model": "two_feature_linear", "features": " + ".join(names), **m,
                                    "candidate_group": " + ".join(groups), "interpretation_flag": flag, "reference_note": note})
                coef_rows.append({"target": target, "alignment": alignment, "model": "two_feature_linear", "features": " + ".join(names),
                                  "feature_coefficients": json.dumps({n1: float(c[0]), n2: float(c[1])}), "intercept": float(c[-1])})
                for i, (_, row) in enumerate(df.iterrows()):
                    pred_rows.append({"target": target, "alignment": alignment, "model": "two_feature_linear", "features": " + ".join(names),
                                      "surface_name": row.surface_name, "catalyst": row.catalyst, "actual_eV": y[i], "predicted_eV": p[i], "residual_pred_minus_actual_eV": p[i]-y[i]})
        # Leakage diagnostic separate from normal candidate ranking.
        x = df.OOH_C6H12_candidate_legacy.to_numpy().reshape(-1, 1)
        m, c, p = metrics(x, y)
        metric_rows.append({"target": target, "alignment": "legacy_reference", "model": "single_linear_diagnostic_only", "features": "OOH-C6H12* candidate", **m,
                            "candidate_group": "transition_state_candidate_diagnostic", "interpretation_flag": "EXCLUDED_FROM_PRIMARY_RANKING", "reference_note": "possible transition-state/target information overlap; structural relation to the actual NEB1 TS1 is unverified"})
    # Conventional one-variable BEP check, reported separately from adsorbate descriptor ranking.
    y = ytargets["NEB1_forward_energy_difference_TS1_minus_Rpair"]
    x = df.NEB1_reaction_energy_EpoxideHHO_minus_Rpair.to_numpy().reshape(-1, 1)
    m, c, p = metrics(x, y)
    metric_rows.append({"target": "NEB1_forward_energy_difference_TS1_minus_Rpair", "alignment": "current_strict_network",
                        "model": "BEP_linear", "features": "NEB1_reaction_energy_EpoxideHHO_minus_Rpair", **m,
                        "candidate_group": "reaction-energy_BEP", "interpretation_flag": "separate_physical_relation",
                        "reference_note": "BEP-style fit; both signed targets are derived from the same current network energies; not an independent descriptor validation"})
    coef_rows.append({"target": "NEB1_forward_energy_difference_TS1_minus_Rpair", "alignment": "current_strict_network",
                      "model": "BEP_linear", "features": "NEB1_reaction_energy_EpoxideHHO_minus_Rpair",
                      "feature_coefficients": json.dumps({"NEB1_reaction_energy_EpoxideHHO_minus_Rpair": float(c[0])}), "intercept": float(c[-1])})
    for i, (_, row) in enumerate(df.iterrows()):
        pred_rows.append({"target": "NEB1_forward_energy_difference_TS1_minus_Rpair", "alignment": "current_strict_network",
                          "model": "BEP_linear", "features": "NEB1_reaction_energy_EpoxideHHO_minus_Rpair",
                          "surface_name": row.surface_name, "catalyst": row.catalyst, "actual_eV": y[i], "predicted_eV": p[i],
                          "residual_pred_minus_actual_eV": p[i]-y[i]})
    metrics_df = pd.DataFrame(metric_rows)
    coeff_df = pd.DataFrame(coef_rows)
    pred_df = pd.DataFrame(pred_rows)
    metrics_df.to_csv(OUT / "data/fit_metrics_all_candidates.csv", index=False, encoding="utf-8-sig")
    coeff_df.to_csv(OUT / "data/fit_coefficients_all_candidates.csv", index=False, encoding="utf-8-sig")
    pred_df.to_csv(OUT / "data/predictions_residuals_all_candidates.csv", index=False, encoding="utf-8-sig")
    return metrics_df, coeff_df, pred_df


def plot_performance(m: pd.DataFrame) -> None:
    data = m[(m.alignment == "legacy_reference") & (m.interpretation_flag == "primary_candidate")]
    fig, axes = plt.subplots(2, 1, figsize=(13, 10), sharex=True)
    for ax, col, label in [(axes[0], "R2_in_sample", "In-sample R² (higher is better)"), (axes[1], "RMSE_eV", "In-sample RMSE (eV; lower is better)")]:
        top = data[data.target == "TS1_formation_energy"].sort_values(col, ascending=(col != "R2_in_sample")).head(12)
        names = top.features.str.replace("_legacy", "", regex=False).str.replace("_strict", "", regex=False)
        ax.barh(names.iloc[::-1], top[col].iloc[::-1], color="#4472C4")
        ax.set_xlabel(label)
        ax.grid(axis="x", alpha=.25)
    fig.suptitle("TS1 formation energy reconstruction: top descriptor models\n13 matched catalysts; ordinary linear fits; legacy descriptor reference")
    fig.tight_layout()
    fig.savefig(OUT / "figures/TS1_top_descriptor_performance.png", dpi=220)
    fig.savefig(OUT / "figures/TS1_top_descriptor_performance.pdf")
    plt.close(fig)


def annotate_points(ax, x, y, labels):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    labels = list(labels)
    ax.scatter(x, y, s=48, color="#1f77b4", edgecolors="white", linewidths=.6, zorder=3, label="13 validated catalyst data")
    order = np.argsort(np.asarray(y, dtype=float))
    stack = np.linspace(.05, .95, len(labels))
    for pos, idx in zip(stack, order):
        ax.annotate(labels[idx], (x[idx], y[idx]), xytext=(1.04, pos), textcoords="axes fraction", ha="left", va="center",
                    fontsize=8, zorder=4, annotation_clip=False,
                    arrowprops={"arrowstyle": "-", "color": "#777777", "lw": .55, "shrinkA": 1, "shrinkB": 3},
                    bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "#bbbbbb", "lw": .35, "alpha": .92})


def plot_actual_predicted(metrics_df: pd.DataFrame, pred: pd.DataFrame, target: str, output_stem: str, title: str, feature_filter=None):
    candidates = metrics_df[(metrics_df.target == target) & (metrics_df.alignment == "legacy_reference") &
                            (metrics_df.model == "two_feature_linear") & (metrics_df.interpretation_flag == "primary_candidate")]
    if feature_filter:
        candidates = candidates[candidates.features.map(feature_filter)]
    candidates = candidates.sort_values(["RMSE_eV", "LOOCV_RMSE_eV"]).head(2)
    if candidates.empty:
        return []
    paths = []
    for _, cand in candidates.iterrows():
        d = pred[(pred.target == target) & (pred.alignment == "legacy_reference") & (pred.model == cand.model) & (pred.features == cand.features)].copy()
        fig, ax = plt.subplots(figsize=(9.3, 6.2))
        fig.subplots_adjust(right=.71)
        annotate_points(ax, d.actual_eV, d.predicted_eV, d.catalyst)
        low = float(min(d.actual_eV.min(), d.predicted_eV.min()))
        high = float(max(d.actual_eV.max(), d.predicted_eV.max()))
        ax.plot([low, high], [low, high], "k--", lw=1.2, label="y = x")
        ax.set_xlabel(f"Actual {target.replace('_', ' ')} (eV)")
        ax.set_ylabel("Predicted (eV)")
        ax.set_title(f"{title}\n{cand.features}; R²={cand.R2_in_sample:.3f}, RMSE={cand.RMSE_eV:.3f} eV")
        ax.legend(frameon=False)
        ax.grid(alpha=.2)
        fig.tight_layout()
        stem = f"{output_stem}_{cand.features.replace(' + ', '__').replace('*', 'star')}"
        p1, p2 = OUT / f"figures/{stem}.png", OUT / f"figures/{stem}.pdf"
        fig.savefig(p1, dpi=220); fig.savefig(p2); plt.close(fig)
        paths.extend([str(p1), str(p2)])
    return paths


def export_response_map(df: pd.DataFrame, metrics_df: pd.DataFrame, coeff_df: pd.DataFrame, target: str, alignment: str, features: tuple[str, str], stem: str, nx=200, ny=200):
    m = metrics_df[(metrics_df.target == target) & (metrics_df.alignment == alignment) &
                   (metrics_df.model == "two_feature_linear") &
                   (metrics_df.features == " + ".join(features))]
    c = coeff_df[(coeff_df.target == target) & (coeff_df.alignment == alignment) &
                 (coeff_df.model == "two_feature_linear") & (coeff_df.features == " + ".join(features))]
    if m.empty or c.empty:
        return None
    coef_json = json.loads(c.iloc[0].feature_coefficients)
    intercept = float(c.iloc[0].intercept)
    b1, b2 = [coef_json[f] for f in features]
    col_suffix = "_legacy" if alignment == "legacy_reference" else "_slabshift_conditional"
    def get_real(feat):
        # Descriptor names in the feature inventory map back to frame column names.
        base = feat.removesuffix(col_suffix)
        if base == "C6H12O_minus_C6H12":
            col = base + col_suffix
        elif base == "Rpair_strict":
            col = base
        else:
            col = base + col_suffix
        return df[col].to_numpy(), col
    x_real, x_col = get_real(features[0]); y_real, y_col = get_real(features[1])
    def bounds(vals):
        lo, hi = float(np.min(vals)), float(np.max(vals))
        pad = .05 * (hi - lo) if hi > lo else .05
        return lo - pad, hi + pad
    xmin, xmax = bounds(x_real); ymin, ymax = bounds(y_real)
    gx, gy = np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny)
    xx, yy = np.meshgrid(gx, gy)
    zz = intercept + b1*xx + b2*yy
    xyz = pd.DataFrame({"x": xx.ravel(), "y": yy.ravel(), "predicted_target_eV": zz.ravel()})
    xyz.to_csv(OUT / f"origin_exports/{stem}_{nx}x{ny}_xyz.csv", index=False, encoding="utf-8-sig")
    matrix = pd.DataFrame(zz, index=gy, columns=gx)
    matrix.index.name = "y\\x"
    matrix.to_csv(OUT / f"origin_exports/{stem}_{nx}x{ny}_matrix.csv", encoding="utf-8-sig")
    points = pd.DataFrame({"surface_name": df.surface_name, "catalyst_name": df.catalyst, "x": x_real, "y": y_real,
                           "actual_target_eV": df.TS1_strict if target == "TS1_formation_energy" else df.NEB1_forward_energy_difference_TS1_minus_Rpair,
                           "TOF_C6H12O": df.TOF, "log10_TOF": df.log10_TOF, "point_type": "validated catalyst; empirical fit observation"})
    points.to_csv(OUT / f"origin_exports/{stem}_catalyst_points.csv", index=False, encoding="utf-8-sig")
    fig = plt.figure(figsize=(12.8, 8.0))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.0, .22], height_ratios=[1.0, .065], wspace=.08, hspace=.18)
    ax = fig.add_subplot(grid[0, 0])
    cax = fig.add_subplot(grid[1, 0])
    label_ax = fig.add_subplot(grid[0, 1])
    label_ax.set_axis_off()
    cf = ax.contourf(xx, yy, zz, levels=24, cmap="viridis")
    contour = ax.contour(xx, yy, zz, levels=10, colors="white", linewidths=.5, alpha=.6)
    ax.clabel(contour, inline=True, fontsize=7, fmt="%.2f")
    ax.scatter(x_real, y_real, s=44, facecolors="none", edgecolors="black", linewidths=1.1, label="13 validated catalyst data", zorder=4)
    order = np.argsort(y_real)
    stack = np.linspace(.04, .96, len(df))
    for pos, idx in zip(stack, order):
        label_ax.text(.08, pos, df.catalyst.iloc[idx], transform=label_ax.transAxes, ha="left", va="center",
                      fontsize=8, bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "#bbbbbb", "lw": .35, "alpha": .92})
        connector = ConnectionPatch(xyA=(x_real[idx], y_real[idx]), coordsA=ax.transData,
                                    xyB=(.02, pos), coordsB=label_ax.transAxes,
                                    axesA=ax, axesB=label_ax, color="#777777", lw=.55, zorder=2)
        fig.add_artist(connector)
    ax.set_xlabel(f"{PRETTY.get(features[0],features[0])} (eV)")
    ax.set_ylabel(f"{PRETTY.get(features[1],features[1])} (eV)")
    ax.set_title(f"Empirical linear response surface for {target.replace('_',' ')}\n{features[0]} + {features[1]} | n=13 | no convex-hull mask\nNot CatMAP continuous kinetics; rectangle is a linear-fit display")
    cb = fig.colorbar(cf, cax=cax, orientation="horizontal"); cb.set_label(f"Predicted {target.replace('_',' ')} (eV)")
    ax.legend(frameon=True, loc="best")
    fig.subplots_adjust(left=.08, right=.98, bottom=.08, top=.84)
    for ext in ["png", "pdf"]:
        fig.savefig(OUT / f"figures/{stem}_{nx}x{ny}.{ext}", dpi=220 if ext == "png" else None)
    plt.close(fig)
    return {"features": list(features), "target": target, "alignment": alignment, "n_x": nx, "n_y": ny,
            "bounds": {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax}, "coefficients": {features[0]: b1, features[1]: b2, "intercept": intercept},
            "fit_R2": float(m.iloc[0].R2_in_sample), "fit_RMSE_eV": float(m.iloc[0].RMSE_eV),
            "export_xyz": str(OUT / f"origin_exports/{stem}_{nx}x{ny}_xyz.csv"),
            "export_matrix": str(OUT / f"origin_exports/{stem}_{nx}x{ny}_matrix.csv"),
            "points_file": str(OUT / f"origin_exports/{stem}_catalyst_points.csv")}


def main() -> None:
    mkdirs()
    data, strict, source_meta = load_data()
    source_meta["workbook_crosswalk"] = export_workbook_crosswalk(data)
    data.to_csv(OUT / "data/matched_carbon_O_H_TS1_data.csv", index=False, encoding="utf-8-sig")
    source_meta["source_hashes"][str(OUT / "data/matched_carbon_O_H_TS1_data.csv")] = sha256(OUT / "data/matched_carbon_O_H_TS1_data.csv")

    # Full, explicit inventory keeps candidate transition-state structures out of the primary feature search.
    inventory = [
        ("C6H12_star_legacy", "C6H12*", "single-molecule alkene adsorbate formation energy", "historical DFT_formation_catmap.txt; DFT_Gibbs_single", "primary; pre-reaction adsorbate"),
        ("C6H12_adsorption_energy_legacy", "C6H12* adsorption", "Eform(C6H12*) - Eform(C6H12_g)", "derived from historical gas and adsorbate rows", "same site variation as C6H12* formation; constant gas shift; not a new independent feature"),
        ("C6H12O_star_legacy", "C6H12O*", "epoxide product adsorbate formation energy", "historical DFT_formation_catmap.txt; DFT_Gibbs_single", "supplementary; product-related, not a pre-reaction descriptor"),
        ("C6H12O_adsorption_energy_legacy", "C6H12O* adsorption", "Eform(C6H12O*) - Eform(C6H12O_g)", "derived from historical gas and adsorbate rows", "same site variation as C6H12O* formation; constant gas shift"),
        ("C6H12O_minus_C6H12_legacy", "C6H12O* - C6H12*", "same-surface state-energy difference", "derived; common slab cancels algebraically", "product-vs-reactant state contrast; exploratory"),
        ("O_star_legacy", "O*", "historical atomic oxygen adsorbate formation energy", "historical DFT_formation_catmap.txt; DFT_Gibbs_single", "legacy reference; distinct from H features"),
        ("Ha_star_legacy", "Ha*", "historical Ha adsorbate formation energy", "historical DFT_formation_catmap.txt; DFT_Gibbs_single", "legacy reference; not merged with Hb"),
        ("Hb_star_legacy", "Hb*", "historical Hb adsorbate formation energy", "historical DFT_formation_catmap.txt; DFT_Gibbs_single", "legacy reference; not merged with Ha"),
        ("OOH_star_legacy", "OOH*", "historical OOH adsorbate formation energy", "historical DFT_formation_catmap.txt; DFT_Gibbs_single", "historical single state; thermal/reference comparability conditional"),
        ("OH_star_legacy", "OH*", "historical OH adsorbate formation energy", "historical DFT_formation_catmap.txt; DFT_Gibbs_single", "historical single state; thermal/reference comparability conditional"),
        ("OOH_C6H12_candidate_legacy", "OOH-C6H12*", "historical candidate-state formation energy", "DFT_Gibbs_TS_candidate", "excluded from primary ranking: possible TS1 information overlap; structure identity unverified"),
        ("Hb_O_candidate_legacy", "Hb-O*", "historical candidate-state formation energy", "DFT_Gibbs_TS_candidate", "excluded: candidate transition-state-like state; not an independently verified stable adsorbate descriptor"),
        ("Rpair_strict", "Rpair*", "current six-step model initial coadsorbed-state formation energy", "current official energies.txt", "network-derived; barrier target includes subtraction of this same state"),
        ("TS1_strict", "TS1*", "current six-step model transition-state formation energy", "current official energies.txt", "target only; never a predictor"),
        ("NEB1_forward_energy_difference_TS1_minus_Rpair", "TS1 - Rpair", "literal signed NEB1 energy difference", "derived from current strict energies.txt", "separate target; no clipping of negative values"),
    ]
    inv = pd.DataFrame(inventory, columns=["feature_key", "species_or_contrast", "definition", "source", "status_and_interpretation"])
    inv.to_csv(OUT / "data/descriptor_source_inventory.csv", index=False, encoding="utf-8-sig")
    avail = []
    for _, row in inv.iterrows():
        key = row.feature_key
        col = key if key in data.columns else None
        if key.endswith("_legacy") and key.replace("_legacy", "_legacy") in data.columns: col = key
        if key == "C6H12_adsorption_energy_legacy": col = key
        if key == "C6H12O_adsorption_energy_legacy": col = key
        if key == "TS1_strict": col = "TS1_strict"
        if col:
            series = data[col]
            avail.append({"feature_key": key, "n_nonmissing": int(series.notna().sum()), "n_unique": int(series.nunique()), "min_eV": float(series.min()), "max_eV": float(series.max())})
    pd.DataFrame(avail).to_csv(OUT / "data/descriptor_availability.csv", index=False, encoding="utf-8-sig")

    metrics_df, coeff_df, pred_df = build_candidate_fits(data)
    focus_rows = []
    for target in ["TS1_formation_energy", "NEB1_forward_energy_difference_TS1_minus_Rpair"]:
        for alignment, suffix in [("legacy_reference", "legacy"), ("slabshift_conditional", "slabshift_conditional")]:
            for other in ["O_star", "Ha_star", "Hb_star", "OOH_star", "OH_star"]:
                names = f"C6H12_star_{suffix} + {other}_{suffix}"
                row = metrics_df[(metrics_df.target == target) & (metrics_df.alignment == alignment) &
                                 (metrics_df.model == "two_feature_linear") & (metrics_df.features == names)]
                if not row.empty:
                    focus_rows.append(row.iloc[0].to_dict())
    focus_df = pd.DataFrame(focus_rows)
    focus_df.to_csv(OUT / "data/C6H12star_plus_O_H_pair_comparison.csv", index=False, encoding="utf-8-sig")
    plot_performance(metrics_df)
    # Fit actual-vs-predicted panels for leading carbon-containing pairs.
    def carbon_pattern(s: str) -> bool:
        # Match the same prereaction C6H12* + O/H feature family in either
        # reference track, without mixing legacy and conditional columns.
        parts = s.split(" + ")
        if len(parts) != 2:
            return False
        carbon, partner = parts
        if carbon.endswith("_legacy"):
            suffix = "_legacy"
        elif carbon.endswith("_slabshift_conditional"):
            suffix = "_slabshift_conditional"
        else:
            return False
        allowed_partners = {f"{name}{suffix}" for name in
                            ("O_star", "Ha_star", "Hb_star", "OOH_star", "OH_star")}
        return carbon == f"C6H12_star{suffix}" and partner in allowed_partners
    carbon_any_pattern = lambda s: (any(k in s for k in ["C6H12_star", "C6H12O_star", "C6H12O_minus_C6H12"]) and "Rpair" not in s)
    ts_paths = plot_actual_predicted(metrics_df, pred_df, "TS1_formation_energy", "TS1_actual_vs_predicted_carbon_pair",
                                     "TS1 formation energy: leading carbon + O/H descriptor pair", carbon_pattern)
    ts_paths += plot_actual_predicted(metrics_df, pred_df, "TS1_formation_energy", "TS1_actual_vs_predicted_best_carbon_pair",
                                     "TS1 formation energy: best carbon-containing pair overall", carbon_any_pattern)
    barrier_paths = plot_actual_predicted(metrics_df, pred_df, "NEB1_forward_energy_difference_TS1_minus_Rpair", "NEB1_barrier_actual_vs_predicted_carbon_pair",
                                          "Signed NEB1 TS1−Rpair energy difference: leading carbon + O/H descriptor pair", carbon_pattern)

    # Select by sample RMSE among interpretable carbon-containing pairs; keep pre-reaction and product-state models separate.
    map_summaries = []
    map_requests = [
        ("TS1_formation_energy", "TS1_best_prereaction_C6H12_O_H_empirical_surface", carbon_pattern),
        ("TS1_formation_energy", "TS1_best_carbon_overall_empirical_surface", carbon_any_pattern),
        ("NEB1_forward_energy_difference_TS1_minus_Rpair", "NEB1_barrier_best_prereaction_C6H12_O_H_empirical_surface", carbon_pattern),
    ]
    for target, stem, pair_filter in map_requests:
        pool = metrics_df[(metrics_df.target == target) & (metrics_df.alignment == "legacy_reference") &
                          (metrics_df.model == "two_feature_linear") & (metrics_df.interpretation_flag == "primary_candidate")]
        pool = pool[pool.features.map(pair_filter)].sort_values(["RMSE_eV", "LOOCV_RMSE_eV"])
        if not pool.empty:
            features = tuple(pool.iloc[0].features.split(" + "))
            for nx, ny in [(200, 200), (400, 400)]:
                summary = export_response_map(data, metrics_df, coeff_df, target, "legacy_reference", features, stem, nx, ny)
                if summary: map_summaries.append(summary)

    # Reproducibility: verify baseline O+Ha from the current input records and compare it to the archived reported value.
    baseline = metrics_df[(metrics_df.target == "TS1_formation_energy") & (metrics_df.alignment == "legacy_reference") &
                         (metrics_df.model == "two_feature_linear") & (metrics_df.features == "O_star_legacy + Ha_star_legacy")]
    if baseline.empty:
        raise ValueError("could not reproduce the O+Ha comparison model")
    baseline_row = baseline.iloc[0].to_dict()

    # Two clean rank lists: primary non-Rpair and all mathematical candidates (Rpair is explicitly flagged).
    ts1_rank = metrics_df[(metrics_df.target == "TS1_formation_energy") & (metrics_df.alignment == "legacy_reference") &
                          (metrics_df.model == "two_feature_linear")].sort_values("RMSE_eV")
    ts1_conditional_rank = metrics_df[(metrics_df.target == "TS1_formation_energy") & (metrics_df.alignment == "slabshift_conditional") &
                                      (metrics_df.model == "two_feature_linear") & (metrics_df.interpretation_flag == "primary_candidate")].sort_values("RMSE_eV")
    barrier_rank = metrics_df[(metrics_df.target == "NEB1_forward_energy_difference_TS1_minus_Rpair") & (metrics_df.alignment == "legacy_reference") &
                              (metrics_df.model == "two_feature_linear")].sort_values("RMSE_eV")
    ts1_rank.to_csv(OUT / "data/TS1_two_descriptor_ranking_legacy.csv", index=False, encoding="utf-8-sig")
    barrier_rank.to_csv(OUT / "data/NEB1_signed_energy_difference_two_descriptor_ranking_legacy.csv", index=False, encoding="utf-8-sig")
    metrics_df[metrics_df.target == "TS1_formation_energy"].to_csv(OUT / "data/TS1_complete_descriptor_ranking_both_reference_tracks.csv", index=False, encoding="utf-8-sig")

    summary = {
        "cohort": SURFACES, "catalysts": [DISPLAY[s] for s in SURFACES], "n": len(data),
        "temperature_K": 333.15, "tof_recalculated": False, "catmap_recalculated": False,
        "baseline_O_plus_Ha": baseline_row,
        "best_TS1_pair_primary": ts1_rank[ts1_rank.interpretation_flag == "primary_candidate"].iloc[0].to_dict(),
        "best_TS1_pair_conditional_slabshift": ts1_conditional_rank.iloc[0].to_dict(),
        "best_TS1_pair_conditional_carbon_prereaction": ts1_conditional_rank[ts1_conditional_rank.features.map(carbon_pattern)].iloc[0].to_dict(),
        "best_TS1_pair_carbon_prereaction": ts1_rank[ts1_rank.features.map(carbon_pattern)].iloc[0].to_dict(),
        "best_NEB1_signed_difference_pair_primary": barrier_rank[barrier_rank.interpretation_flag == "primary_candidate"].iloc[0].to_dict(),
        "best_NEB1_pair_conditional_slabshift": metrics_df[(metrics_df.target == "NEB1_forward_energy_difference_TS1_minus_Rpair") &
            (metrics_df.alignment == "slabshift_conditional") & (metrics_df.model == "two_feature_linear") &
            (metrics_df.interpretation_flag == "primary_candidate")].sort_values("RMSE_eV").iloc[0].to_dict(),
        "best_NEB1_signed_difference_pair_carbon_prereaction": barrier_rank[barrier_rank.features.map(carbon_pattern)].iloc[0].to_dict(),
        "best_TS1_Rpair_coupled_diagnostic": ts1_rank[ts1_rank.interpretation_flag == "diagnostic_target_coupling"].sort_values("RMSE_eV").iloc[0].to_dict(),
        "best_barrier_Rpair_coupled_diagnostic": barrier_rank[barrier_rank.interpretation_flag == "diagnostic_target_coupling"].sort_values("RMSE_eV").iloc[0].to_dict(),
        "BEP_fit": metrics_df[(metrics_df.target == "NEB1_forward_energy_difference_TS1_minus_Rpair") & (metrics_df.model == "BEP_linear")].iloc[0].to_dict(),
        "C6H12_star_pair_focus": focus_df[focus_df.target == "TS1_formation_energy"].to_dict("records"),
        "best_TS1_single_primary": metrics_df[(metrics_df.target == "TS1_formation_energy") & (metrics_df.alignment == "legacy_reference") & (metrics_df.model == "single_linear") & (metrics_df.interpretation_flag == "primary_candidate")].sort_values("RMSE_eV").iloc[0].to_dict(),
        "best_TS1_single_Rpair_diagnostic": metrics_df[(metrics_df.target == "TS1_formation_energy") & (metrics_df.alignment == "legacy_reference") & (metrics_df.model == "single_linear") & (metrics_df.interpretation_flag == "diagnostic_target_coupling")].sort_values("RMSE_eV").iloc[0].to_dict(),
        "map_summaries": map_summaries,
        "actual_predicted_figures": ts_paths + barrier_paths,
        "source_metadata": source_meta,
    }
    (OUT / "data/analysis_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    (OUT / "reports/ANALYSIS_REPORT.md").write_text(render_report(summary, data), encoding="utf-8")
    # Manifest covers every produced file except itself; source hashes are also in analysis_summary.json.
    files = sorted(p for p in OUT.rglob("*") if p.is_file() and p.name != "SHA256SUMS.csv")
    pd.DataFrame([{"path": str(p.relative_to(OUT)), "sha256": sha256(p), "bytes": p.stat().st_size} for p in files]).to_csv(
        OUT / "reports/SHA256SUMS.csv", index=False, encoding="utf-8-sig")
    print(json.dumps({"output": str(OUT), "baseline": baseline_row,
                      "best_TS1_pair": summary["best_TS1_pair_primary"],
                      "best_TS1_carbon_pre_reaction": summary["best_TS1_pair_carbon_prereaction"],
                      "best_barrier_pair": summary["best_NEB1_signed_difference_pair_primary"],
                      "best_barrier_carbon_pre_reaction": summary["best_NEB1_signed_difference_pair_carbon_prereaction"],
                      "outputs": len(files)}, indent=2, ensure_ascii=True, default=str))


def render_report(summary: dict, df: pd.DataFrame) -> str:
    def fmt(rec):
        return f"`{rec['features']}` (N={rec['N']}, R²={rec['R2_in_sample']:.4f}, RMSE={rec['RMSE_eV']:.4f} eV, MAE={rec['MAE_eV']:.4f} eV, LOOCV RMSE={rec['LOOCV_RMSE_eV']:.4f} eV)"
    ts1 = summary["best_TS1_pair_primary"]
    ts1c = summary["best_TS1_pair_carbon_prereaction"]
    barr = summary["best_NEB1_signed_difference_pair_primary"]
    barrc = summary["best_NEB1_signed_difference_pair_carbon_prereaction"]
    baseline = summary["baseline_O_plus_Ha"]
    ts1_gain = float(ts1c["R2_in_sample"] - baseline["R2_in_sample"])
    focus_rows = summary["C6H12_star_pair_focus"]
    focus_names = [("O", "O_star"), ("Ha", "Ha_star"), ("Hb", "Hb_star"), ("OOH", "OOH_star"), ("OH", "OH_star")]
    focus_table = []
    for label, feature in focus_names:
        raw_key = f"C6H12_star_legacy + {feature}_legacy"
        shift_key = f"C6H12_star_slabshift_conditional + {feature}_slabshift_conditional"
        raw = next(r for r in focus_rows if r["alignment"] == "legacy_reference" and r["features"] == raw_key)
        shift = next(r for r in focus_rows if r["alignment"] == "slabshift_conditional" and r["features"] == shift_key)
        focus_table.append(f"| C6H12* + {label}* | {raw['R2_in_sample']:.4f} | {raw['RMSE_eV']:.4f} | {raw['LOOCV_RMSE_eV']:.4f} | {shift['R2_in_sample']:.4f} | {shift['RMSE_eV']:.4f} | {shift['LOOCV_RMSE_eV']:.4f} |")
    return f"""# 含碳吸附物种与 O/H 特征对 TS1 / NEB1 的探索性标度分析

## 范围与来源

- 13个位点：{', '.join(summary['catalysts'])}；Ti–Zr排除。
- 使用归档的原始形成能表、当前严格模型 `energies.txt`、已有13位点收敛TOF表和slab对照审计；所有连接均按 `surface_name` 一对一完成。
- 原始形成能表与归档副本SHA256一致：`{summary['source_metadata']['source_hashes'][str(RAW_LEGACY)]}`；原始转换脚本SHA256：`{summary['source_metadata']['source_hashes'][str(CONVERTER)]}`。Excel源单元格和公式见 `data/source_workbook_crosswalk.csv`。
- 没有重新运行CatMAP、稳态、条件搜索或DFT；所有图为普通线性回归的经验响应，不是CatMAP描述符网格。
- 单物种 `C6H12*` 与 `C*` 不同。本项目既有数据是完整1-己烯吸附态 `C6H12*`，并非单原子碳。

## 主要结果

- O+Ha历史参考基准复核：R²={baseline['R2_in_sample']:.4f}, RMSE={baseline['RMSE_eV']:.4f} eV (N={baseline['N']});与此前0.1382基准相符。
- TS1形成能的数值最优双特征（含网络态/所有主候选）：{fmt(summary['best_TS1_pair_primary'])}。若组合含 `Rpair*`，因目标 `TS1` 与反应初态高度相关，必须视为网络能量重建诊断，不自动视作简化独立描述符。
- TS1形成能的最佳“烯烃吸附态 + 反应前O/H稳定态”组合：{fmt(ts1c)}；相对O+Ha样本内R²变化 {ts1_gain:+.4f}。
- 仅作参考敏感性诊断的 slab 平移轨迹中，TS1双特征最佳为：{fmt(summary['best_TS1_pair_conditional_slabshift'])}；其中前驱C6H12*与O/H组合为 {fmt(summary['best_TS1_pair_conditional_carbon_prereaction'])}。这类提升不能在参考与热修正闭合前当作已验证规律。
- NEB1原始有符号差值目标 `TS1−Rpair` 的最佳双描述符：{fmt(summary['best_NEB1_signed_difference_pair_primary'])}。
- 同一条件slab平移敏感性轨迹的NEB1差值最佳双特征为：{fmt(summary['best_NEB1_pair_conditional_slabshift'])}；仅为参考敏感性比较。
- 仅在不把 `Rpair` 直接用于构造差值目标的含碳前驱候选中，NEB1差值最佳组合：{fmt(barrc)}。负差值原样保留，没有把能垒截断为零。
- NEB1反应能对有符号 `TS1−Rpair` 的单变量BEP式拟合：R²={summary['BEP_fit']['R2_in_sample']:.4f}、RMSE={summary['BEP_fit']['RMSE_eV']:.4f} eV、LOOCV RMSE={summary['BEP_fit']['LOOCV_RMSE_eV']:.4f} eV；这是同一网络能量内部的关系，不等于独立验证。
- 若将 `Rpair*` 纳入预测变量，数值最佳网络态耦合诊断为TS1：{fmt(summary['best_TS1_Rpair_coupled_diagnostic'])}；NEB1差值：{fmt(summary['best_barrier_Rpair_coupled_diagnostic'])}。这些不列为简化的独立描述符结论，因为目标本身与该网络态直接相关。
- 最佳单一历史吸附态特征TS1模型（排除Rpair网络态）：{fmt(summary['best_TS1_single_primary'])}。`Rpair*`单特征单列为网络耦合诊断：{fmt(summary['best_TS1_single_Rpair_diagnostic'])}。

### C6H12* + O/H 对TS1形成能的逐组合结果（N=13）

| 组合 | 历史参考 R² | 历史 RMSE (eV) | 历史 LOOCV RMSE (eV) | 条件slab平移 R² | 条件slab平移 RMSE (eV) | 条件slab平移 LOOCV RMSE (eV) |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(focus_table)}

条件平移列仅用于显示参考选择敏感性；不是已验证的统一物理参考。

## 物种可用性与定义

`DFT_formation_catmap.txt` 中存在13个共同位点的 `C6H12*`、`C6H12O*`、`O*`、`Ha*`、`Hb*`、`OOH*`、`OH*` 和候选 `OOH-C6H12*`、`Hb-O*`。C6H12和C6H12O气相形成能也存在。`C6H12*`吸附自由能按 `Eform(C6H12*)−Eform(C6H12_g)` 得到；所用历史气相值分别为 {summary['source_metadata']['old_gas_C6H12_formation_energy_eV']:.6f} eV 和 {summary['source_metadata']['old_gas_C6H12O_formation_energy_eV']:.6f} eV。气相项为跨位点常数，因此含截距的OLS中它仅平移截距，不改变R²/残差。没有独立C*，也未将 `C6H12*` 冒充C*。

## 能量参考与解释限制

旧单吸附态和当前六步网络能量并非已完全证明为同一物理参考。旧/当前洁净slab差值在13个位点间并不统一，已知大偏移包括Ti-Fe、Ti-Co、Ti-Ti、Ti-Ta、Ti-Ni及Ti-Mo。报告同时保留：

1. `legacy_reference`：直接使用历史表原数值；
2. `slabshift_conditional`：仅按已审计的当前减历史slab能差平移旧吸附态形成能，假设对应吸附态总G及热修正在其他方面可比。

第二种按 `Eform_conditional = Eform_legacy - (Gslab_current - Gslab_legacy)` 计算，仅为条件性参考对齐诊断，不代表热修正/结构等价性已经验证。所有排名同时保存两条轨迹。形成能的绝对目标 `TS1` 来自当前严格模型输入，仍沿用其原有条件性来源标记。

`OOH-C6H12*` 和 `Hb-O*` 在旧转换脚本中被标为 `DFT_Gibbs_TS_candidate`，注释明确其并不等同于共吸附NEB TS自由能。因其与TS1的结构身份/独立性未核实，OOH-C6H12只作为单变量信息重叠诊断，未进入主组合排名；Hb-O不作为稳定态描述符。`Rpair*`是当前网络前驱态，主排名中单独标识，并避免把它作为NEB1差值候选组合的机制结论。

## 图与Origin数据

- `figures/TS1_top_descriptor_performance.png/.pdf`：TS1形成能最佳候选模型指标。
- `figures/TS1_actual_vs_predicted_carbon_pair_*.png/.pdf`：最佳含碳前驱态对照模型，带13个位点标签及y=x线。
- `figures/TS1_actual_vs_predicted_best_carbon_pair_*.png/.pdf`：全体非循环含碳组合中的高排名对照；产品吸附态模型单独标作事后描述。
- `figures/NEB1_barrier_actual_vs_predicted_carbon_pair_*.png/.pdf`：NEB1原始有符号能量差拟合图。
- `figures/*_200x200.png/.pdf`、`*_400x400.png/.pdf`：经验线性响应面，无凸包遮罩、全矩形显示、标出全部13个数据点和名称；不是CatMAP连续网格。
- `origin_exports/*_xyz.csv` 与 `*_matrix.csv`：完整规则网格；`*_catalyst_points.csv`：真实点坐标、目标值与TOF。
- `data/fit_metrics_all_candidates.csv`、系数、逐位点预测/残差CSV和两个目标的排名CSV保存全部可用模型；`data/descriptor_source_inventory.csv`记录物种定义与来源。

## 科学解释边界

这些是13个催化剂的小样本内线性重建结果。LOOCV用于稳定性提示，不作为绘图门槛。拟合相关不证明因果机制；含产品吸附态的模型属于事后描述，含 `Rpair` 的模型对反应能垒目标有代数耦合。最优关联应结合结构核查、统一热力学参考和独立位点数据再用于TS1–TS2–TOF简化模型。当前没有证据把任何经验线性响应面称为火山峰。

## 可复现性

脚本：`scripts/analyze_carbon_ts1_scaling.py`。本分析没有重算任何TOF或CatMAP稳态。原始输入SHA256列于 `data/analysis_summary.json`；输出SHA256列于 `reports/SHA256SUMS.csv`。
"""


if __name__ == "__main__":
    main()
