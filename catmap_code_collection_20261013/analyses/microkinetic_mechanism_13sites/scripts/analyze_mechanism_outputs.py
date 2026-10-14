"""Summarize existing CatMAP 0.3.1 strict-baseline outputs for 13 sites.

No steady-state solve is performed here. Eight missing DRC values were
previously obtained with CatMAP's official rate_control output in isolated
site runs. This script joins and checks exported results only.
"""
from __future__ import annotations

import csv
import ast
import hashlib
import json
import math
import pickle
import re
from decimal import Decimal, InvalidOperation
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parents[1]
STRICT = Path(r"C:\Users\Administrator\Documents\codex探索\catmap_neb_coadsorption_model\catmap_official_run\results\official_strict_tolerance")
GRID_DIR = Path(r"H:\catmap\catmap\catmap_code_collection_20261013\catmap_official_run\TS1_TS2_descriptor_grid\runs\grid_50x50")
OLD_DRC = Path(r"C:\Users\Administrator\Documents\codex探索\catmap_neb_coadsorption_model\catmap_official_run\results\activity_analysis_round1\official_drc")
NEW_DRC = HERE / "data" / "official_drc_runs_supplemental"
OUT_DATA = HERE / "data"
OUT_FIG = HERE / "figures"
OUT_ORIGIN = HERE / "origin_exports"
OUT_REPORT = HERE / "reports"

SITES = ["tife", "timn", "tihf", "tire", "timo", "tiv", "tico", "titi", "tiw", "tita", "ticr", "ti", "tini"]
NAMES = {"tife": "Ti–Fe", "timn": "Ti–Mn", "tihf": "Ti–Hf", "tire": "Ti–Re", "timo": "Ti–Mo", "tiv": "Ti–V", "tico": "Ti–Co", "titi": "Ti–Ti", "tiw": "Ti–W", "tita": "Ti–Ta", "ticr": "Ti–Cr", "ti": "Ti", "tini": "Ti–Ni"}
ADS = ["EpoxideHHO_s", "H2O_s", "HHO_s", "OHHa_s", "Rpair_s"]
ALL_COVERAGE = ADS + ["s"]
DRC_PARAMS_OFFICIAL = ["EpoxideHHO_s", "H2O_s", "HHO_s", "OHHa_s", "Rpair_s", "TS1_s", "TS2_s"]
DRC_PARAMS_REPORT = ["Rpair_s", "TS1_s", "EpoxideHHO_s", "HHO_s", "TS2_s", "OHHa_s", "H2O_s"]
GASES = ["C6H12O_g", "C6H12_g", "H2O2_g", "H2O_g"]

REACTION_LABELS = {
    1: "* + H2O2_g + C6H12_g <-> Rpair*",
    2: "Rpair* <-> TS1* -> EpoxideHHO*",
    3: "EpoxideHHO* <-> HHO* + C6H12O_g",
    4: "HHO* <-> TS2* -> OHHa*",
    5: "OHHa* <-> H2O*",
    6: "H2O* <-> H2O_g + *",
}

plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False, "font.size": 10})


def read_csv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def read_directional_export(path: Path):
    """Read CatMAP's existing direction export, whose header has 3 fields
    while its data rows have 2: `direction,reaction,rate` vs label,value.
    Preserve the raw label and numeric text without shifting/guessing columns.
    """
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        if header != ["direction", "reaction", "rate"]:
            raise ValueError(f"Unexpected directional rate header in {path}: {header}")
        rows = []
        for line_no, row in enumerate(reader, 2):
            if len(row) != 2:
                raise ValueError(f"Expected 2 data fields under 3-column header at {path}:{line_no}; got {len(row)}")
            rows.append({"label": row[0], "rate": row[1]})
    if len(rows) != 12:
        raise ValueError(f"Expected 12 directional rows in {path}, found {len(rows)}")
    return rows


def dec(v: str) -> Decimal:
    try:
        return Decimal(str(v).strip())
    except (InvalidOperation, ValueError):
        return Decimal("NaN")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(rows, pd.DataFrame):
        rows.to_csv(path, index=False, encoding="utf-8-sig", float_format="%.12g")
    else:
        rows = list(rows)
        if fields is None:
            fields = list(rows[0].keys()) if rows else []
        with path.open("w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)


def load_drc(site: str) -> tuple[pd.DataFrame, str, Path]:
    new = NEW_DRC / site / "rate_control.csv"
    if new.is_file():
        return read_csv(new), "new_official_catmap_rate_control_run", new
    old = OLD_DRC / site / "rate_control.csv"
    if old.is_file():
        return read_csv(old), "existing_official_catmap_rate_control_run", old
    raise FileNotFoundError(f"No official DRC result for {site}")


def catalyst_rows():
    summary = read_csv(STRICT / "strict_14_summary.csv")
    summary = summary[summary.surface.isin(SITES)].copy()
    if set(summary.surface) != set(SITES):
        raise ValueError(f"strict summary surface mismatch: {sorted(set(SITES)-set(summary.surface))}")
    summary["catalyst"] = summary.surface.map(NAMES)
    summary["tof"] = summary.tof_C6H12O_g.map(float)
    summary["log10_tof"] = summary.tof.map(math.log10)
    summary["residual_float"] = summary.residual.map(float)
    summary["tolerance_float"] = summary.tolerance.map(float)
    summary["coverage_sum_float"] = summary.coverage_sum.map(float)
    summary["status"] = summary.status.astype(str)
    return summary.set_index("surface").loc[SITES].reset_index()


def collect_outputs():
    cats = catalyst_rows()
    coverage_rows, balance_rows, flux_rows, barrier_rows, gas_rows = [], [], [], [], []
    drc_rows, drc_diagnostics, source_hash_rows, const_rows = [], [], [], []
    convergence_rows = []

    summary = read_csv(STRICT / "strict_14_summary.csv").set_index("surface")
    existing_drc = read_csv(OLD_DRC / "five_site_product_drc.csv")
    existing_drc = existing_drc[existing_drc.gas_species.eq("C6H12O_g")]
    old5 = {site for site in SITES if (OLD_DRC / site / "rate_control.csv").is_file()}

    for _, cat in cats.iterrows():
        site = cat.surface
        root = STRICT / site
        cov = read_csv(root / "coverage.csv")
        cov_map = dict(zip(cov.species, cov.coverage))
        if set(cov_map) != set(ALL_COVERAGE):
            raise ValueError(f"{site}: coverage species mismatch {set(cov_map)}")
        vals = {k: dec(v) for k, v in cov_map.items()}
        total = sum(vals.values(), Decimal(0))
        ads_total = sum((vals[k] for k in ADS), Decimal(0))
        for sp in ALL_COVERAGE:
            coverage_rows.append({"surface": site, "catalyst": cat.catalyst, "species": sp, "coverage": cov_map[sp], "is_stable_adsorbate": sp != "s"})
        max_ads = max(ADS, key=lambda s: vals[s])
        max_all = max(ALL_COVERAGE, key=lambda s: vals[s])
        balance_rows.append({"surface":site,"catalyst":cat.catalyst,"status":cat.status,
            "adsorbate_coverage_sum":str(ads_total),"vacancy_coverage":cov_map["s"],
            "total_site_coverage":str(total),"site_balance_error":str(total-Decimal(1)),
            "steady_state_residual":cat.residual,"tolerance":cat.tolerance,
            "dominant_adsorbate":max_ads,"dominant_adsorbate_coverage":cov_map[max_ads],
            "dominant_including_vacancy":max_all,"dominant_all_coverage":cov_map[max_all]})

        net = read_csv(root / "net_reaction_rates.csv")
        direc = read_directional_export(root / "directional_reaction_rates.csv")
        eff = read_csv(root / "effective_barriers.csv")
        gas = read_csv(root / "gas_turnover_frequencies.csv")
        direction_map = {}
        for row in direc:
            label = str(row["label"])
            reaction, direction = label.rsplit(" ", 1)
            if direction not in {"forward", "reverse"}:
                raise ValueError(f"Unexpected directional-rate label: {label!r}")
            direction_map[(reaction, direction)] = row["rate"]
        for _, row in net.iterrows():
            step = int(row.step)
            rxn = row["reaction"]
            if (rxn, "forward") not in direction_map or (rxn, "reverse") not in direction_map:
                raise ValueError(f"{site} step {row['step']}: directional output does not exactly match net reaction string {rxn!r}")
            fw, rv = direction_map[(rxn, "forward")], direction_map[(rxn, "reverse")]
            netd, fwd, rev = dec(row.net_rate), dec(fw), dec(rv)
            diff = (fwd-rev)-netd if fwd.is_finite() and rev.is_finite() else Decimal("NaN")
            scale=max(abs(fwd),abs(rev),abs(netd)) if fwd.is_finite() and rev.is_finite() else Decimal("NaN")
            scaled_error=abs(diff)/scale if scale.is_finite() and scale!=0 else Decimal("NaN")
            stages=ast.literal_eval(rxn)
            product_nu=stages[-1].count("C6H12O_g")-stages[0].count("C6H12O_g")
            flux_rows.append({"surface":site,"catalyst":cat.catalyst,"step":step,"reaction":rxn,
                "reaction_label":REACTION_LABELS[step],"forward_rate_as_exported":fw,"reverse_rate_as_exported":rv,
                "net_rate_high_precision":row.net_rate,"forward_minus_reverse_minus_net":str(diff),"scaled_directional_closure_error":str(scaled_error),
                "CatMAP_parsed_C6H12O_stoichiometric_coefficient":product_nu,
                "forward_rate_negative":fwd<0,"reverse_rate_negative":rev<0,
                "directional_output_status":"NEGATIVE_DIRECTIONAL_RATE; NUMERICALLY INVALID AS A PHYSICAL ONE-WAY FLUX" if fwd<0 or rev<0 else "NONNEGATIVE_API_OUTPUT; closure error retained for review",
                "directional_rate_precision":"CatMAP CSV float output","net_rate_precision":"high-precision exported decimal"})
        for _, row in eff.iterrows():
            step = int(row.step)
            barrier_rows.append({"surface":site,"catalyst":cat.catalyst,"step":step,
                "reaction_label":REACTION_LABELS[step],"barrier_model":row.barrier_model,
                "delta_G_eV":row.delta_G,"forward_effective_barrier_eV":row.forward_effective_barrier,
                "reverse_effective_barrier_eV":row.reverse_effective_barrier,
                "interpretation":"explicit NEB TS" if row.barrier_model=="explicit_TS" else "CatMAP endpoint approximation; not an independently calculated TS"})
        for _, row in gas.iterrows():
            gas_rows.append({"surface":site,"catalyst":cat.catalyst,"gas_species":row.gas_species,
                "turnover_frequency":row.turnover_frequency,"log10_abs_turnover_frequency":math.log10(abs(float(row.turnover_frequency))) if float(row.turnover_frequency)!=0 else "",
                "sign":int(np.sign(float(row.turnover_frequency)))})
        ddf, provenance, dpath = load_drc(site)
        ddf = ddf[ddf.gas_species.eq("C6H12O_g")]
        if len(ddf) != 7:
            raise ValueError(f"{site}: expected 7 product DRC parameters, found {len(ddf)}")
        dm = dict(zip(ddf.energy_parameter, ddf.degree_of_rate_control))
        if set(dm) != set(DRC_PARAMS_OFFICIAL):
            raise ValueError(f"{site}: DRC parameter mismatch {set(dm)}")
        nums = [float(dm[p]) for p in DRC_PARAMS_OFFICIAL]
        if not all(np.isfinite(nums)):
            raise ValueError(f"{site}: nonfinite DRC")
        largest_idx=int(np.argmax(np.abs(nums)))
        for p in DRC_PARAMS_OFFICIAL:
            drc_rows.append({"surface":site,"catalyst":cat.catalyst,"gas_species":"C6H12O_g",
                "energy_parameter":p,"degree_of_rate_control":dm[p],"source_type":provenance,
                "source_file":str(dpath),"strict_baseline_tof":cat.tof})
        drc_diagnostics.append({"surface":site,"catalyst":cat.catalyst,"drc_parameter_count":len(dm),
            "all_finite":True,"drc_sum":sum(nums),"largest_absolute_DRC_parameter":DRC_PARAMS_OFFICIAL[largest_idx],
            "largest_absolute_DRC_value":nums[largest_idx],"source_type":provenance,
            "official_rate_control_file":str(dpath),"jacobian_fail_in_run_log":check_jacobian_fail(site, provenance),
            "interpretation_caution":"finite-difference DRC; sensitivity may be ill-conditioned at extreme low TOF" if cat.tof < 1e-20 else "finite-difference DRC; interpret as local sensitivity"})

        for p in ["forward", "reverse"]:
            for step in range(1,7):
                const_rows.append({"surface":site,"catalyst":cat.catalyst,"step":step,"direction":p,
                    "rate_constant":"","status":"NOT_EXPORTED_IN_EXISTING_BASELINE_CACHE",
                    "note":"No rate-constant column in strict baseline CSV; not reconstructed from net rates."})
        convergence_rows.append({"surface":site,"catalyst":cat.catalyst,"status":cat.status,"tof_C6H12O":cat.tof,
            "log10_tof":cat.log10_tof,"residual":cat.residual,"tolerance":cat.tolerance,
            "coverage_sum":cat.coverage_sum,"flux_spread_relative":cat.flux_spread_relative})

        for rel in ["alkene_epoxidation.mkm","energies.txt","run_catmap_official.py",f"catmap_{site}_strict.log","catmap_console.log",f"catmap_{site}_strict.pkl","convergence.csv","coverage.csv","directional_reaction_rates.csv","net_reaction_rates.csv","effective_barriers.csv","gas_turnover_frequencies.csv"]:
            f = root / rel
            source_hash_rows.append({"source_file":str(f),"sha256":sha(f),"role":"strict baseline input/output"})
        source_hash_rows.append({"source_file":str(dpath),"sha256":sha(dpath),"role":"official DRC output"})

    # Hash complete official DRC evidence directories and this analysis code.
    for drc_root,role in [(OLD_DRC,"existing official DRC evidence"),(NEW_DRC,"supplemental official DRC evidence")]:
        for f in sorted(drc_root.rglob("*")):
            if f.is_file():source_hash_rows.append({"source_file":str(f),"sha256":sha(f),"role":role})
    source_hash_rows.append({"source_file":str(Path(__file__).resolve()),"sha256":sha(Path(__file__).resolve()),"role":"analysis script"})
    # Add summary-held flux and residual fields as direct cross-checks.
    for s in SITES:
        row = summary.loc[s]
        if row.status != "CONVERGED":
            raise ValueError(f"{s}: strict source not CONVERGED ({row.status})")
    coverage = pd.DataFrame(coverage_rows)
    coverage_wide = coverage.pivot(index=["surface","catalyst"], columns="species", values="coverage").reset_index()
    coverage_wide.columns.name = None
    rates = pd.DataFrame(flux_rows)
    barriers = pd.DataFrame(barrier_rows)
    gas_df = pd.DataFrame(gas_rows)
    drc_df = pd.DataFrame(drc_rows)
    drc_wide = drc_df.pivot(index=["surface","catalyst"],columns="energy_parameter",values="degree_of_rate_control").reset_index()
    drc_wide = drc_wide[["surface","catalyst"]+DRC_PARAMS_REPORT]
    tof_flux = []
    for _, cat in cats.iterrows():
        sr = rates[rates.surface.eq(cat.surface)].sort_values("step")
        step3 = dec(sr[sr.step.eq(3)].net_rate_high_precision.iloc[0])
        stoich_product=sum((dec(row.net_rate_high_precision)*int(row.CatMAP_parsed_C6H12O_stoichiometric_coefficient) for _,row in sr.iterrows()),Decimal(0))
        gas_prod = dec(gas_df[(gas_df.surface.eq(cat.surface))&(gas_df.gas_species.eq("C6H12O_g"))].turnover_frequency.iloc[0])
        tof_flux.append({"surface":cat.surface,"catalyst":cat.catalyst,"official_product_TOF":str(gas_prod),
            "step3_net_rate":str(step3),"stoichiometric_product_rate_from_all_six_steps":str(stoich_product),
            "step3_minus_product_TOF":str(step3-gas_prod),"stoichiometric_product_rate_minus_TOF":str(stoich_product-gas_prod),
            "six_step_net_flux_spread_relative":cat.flux_spread_relative,
            "steady_state_status":cat.status})
    # Raw rate-constant status: a traceability table, not an inferred kinetic parameter.
    constants = pd.DataFrame(const_rows)
    return cats, coverage, coverage_wide, pd.DataFrame(balance_rows), rates, barriers, gas_df, drc_df, drc_wide, pd.DataFrame(drc_diagnostics), pd.DataFrame(tof_flux), constants, pd.DataFrame(convergence_rows), pd.DataFrame(source_hash_rows)


def collect_signed_neb_differences(barriers):
    """Preserve raw signed TS-minus-initial-state energy differences."""
    rows=[]
    for site in SITES:
        e=pd.read_csv(STRICT/site/"energies.txt",sep="\t",dtype=str,keep_default_na=False)
        vals={r.species_name:dec(r.formation_energy) for _,r in e.iterrows() if r.surface_name==site}
        for step,ts,initial in [(2,"TS1","Rpair"),(4,"TS2","HHO")]:
            raw=vals[ts]-vals[initial]
            b=barriers[(barriers.surface.eq(site))&(barriers.step.eq(step))].iloc[0]
            rows.append({"surface":site,"catalyst":NAMES[site],"step":step,"TS_species":ts+"_s","initial_state":initial+"_s",
                "raw_signed_TS_minus_initial_eV":str(raw),"CatMAP_forward_effective_barrier_eV":b.forward_effective_barrier_eV,
                "CatMAP_reverse_effective_barrier_eV":b.reverse_effective_barrier_eV,
                "raw_TS_minus_initial_negative":raw<0,
                "note":"Raw signed difference preserved; CatMAP effective barrier is a separate model quantity."})
    return pd.DataFrame(rows)


def collect_existing_grid_coverage():
    """Read the already-computed official CatMAP 50x50 cache; never solve."""
    pkl_path=GRID_DIR/"catmap_ts1_ts2_grid_50x50.pkl"
    with pkl_path.open("rb") as f:
        cache=pickle.load(f)
    covmap=cache.get("coverage_map")
    if not covmap or len(covmap)!=2500:
        raise ValueError("Existing TS1-TS2 cache does not contain the expected 2500 coverage nodes")
    tof=pd.read_csv(GRID_DIR/"catmap_TS1_TS2_50x50_xyz.csv")
    labels=["EpoxideHHO_s","H2O_s","HHO_s","OHHa_s","Rpair_s","s"]
    rows=[]
    for item in covmap:
        xy,cvgs=item
        if len(cvgs)!=len(labels):
            raise ValueError(f"Unexpected coverage vector length {len(cvgs)}; refusing to guess labels")
        row={"TS1_s_eV":float(xy[0]),"TS2_s_eV":float(xy[1])}
        row.update({name:str(value) for name,value in zip(labels,cvgs)})
        if not math.isclose(sum(float(v) for v in cvgs),1.0,rel_tol=0,abs_tol=1e-10):
            raise ValueError(f"Grid coverage sum != 1 at {xy}")
        rows.append(row)
    cv=pd.DataFrame(rows)
    if len(tof)!=len(cv):
        raise ValueError("Grid TOF and coverage map lengths differ")
    # Join by descriptor coordinates, not positional row order.
    cv["_xkey"]=cv.TS1_s_eV.round(8);cv["_ykey"]=cv.TS2_s_eV.round(8)
    tof["_xkey"]=tof.TS1_s_eV.round(8);tof["_ykey"]=tof.TS2_s_eV.round(8)
    cv=cv.merge(tof.drop(columns=["TS1_s_eV","TS2_s_eV"]),on=["_xkey","_ykey"],how="left",validate="one_to_one",suffixes=("","_tof"))
    cv=cv.drop(columns=["_xkey","_ykey"])
    if cv.CatMAP_log10_TOF_C6H12O.isna().any():
        raise ValueError("One or more coverage nodes lack a matching TOF node")
    summary=json.loads((GRID_DIR/"grid_summary.json").read_text(encoding="utf-8"))
    peak=summary["maximum"]
    px,py=float(peak["TS1_s_eV"]),float(peak["TS2_s_eV"])
    xs=sorted(cv.TS1_s_eV.unique());ys=sorted(cv.TS2_s_eV.unique())
    ix=min(range(len(xs)),key=lambda i:abs(xs[i]-px));iy=min(range(len(ys)),key=lambda i:abs(ys[i]-py))
    selection=[("grid_peak",ix,iy),("TS1_lower_side",max(0,ix-2),iy),("TS1_higher_side",min(len(xs)-1,ix+2),iy),
        ("TS2_lower_side",ix,max(0,iy-2)),("TS2_higher_side",ix,min(len(ys)-1,iy+2))]
    sides=[]
    for label,i,j in selection:
        x,y=xs[i],ys[j]
        row=cv[(cv.TS1_s_eV.eq(x))&(cv.TS2_s_eV.eq(y))].iloc[0].to_dict()
        row["grid_location_role"]=label
        sides.append(row)
    anchors=pd.read_csv(GRID_DIR/"validated_catalyst_points.csv")
    return cv,pd.DataFrame(sides),anchors,summary,[pkl_path,GRID_DIR/"grid_summary.json",GRID_DIR/"catmap_TS1_TS2_50x50_xyz.csv",GRID_DIR/"catmap_TS1_TS2_50x50_matrix.csv",GRID_DIR/"catmap_ts1_ts2_grid_50x50.log",GRID_DIR/"validated_catalyst_points.csv"]


def plot_grid_coverages(grid, anchors):
    names=["EpoxideHHO_s","H2O_s","HHO_s","OHHa_s","Rpair_s","s"]
    fig,axes=plt.subplots(2,3,figsize=(18,11),constrained_layout=True)
    for name,ax in zip(names,axes.flat):
        pivot=grid.pivot(index="TS2_s_eV",columns="TS1_s_eV",values=name)
        xx=pivot.columns.to_numpy(dtype=float); yy=pivot.index.to_numpy(dtype=float)
        zz=pivot.to_numpy(dtype=float)
        cf=ax.contourf(xx,yy,zz,levels=18,cmap="viridis",vmin=0,vmax=1)
        ax.scatter(anchors.TS1_s_eV,anchors.TS2_s_eV,s=24,c="white",edgecolors="black",linewidths=.7,zorder=4,label="13 validated catalysts")
        for _,p in anchors.iterrows():
            ax.annotate(p.catalyst,(p.TS1_s_eV,p.TS2_s_eV),xytext=(3,3),textcoords="offset points",fontsize=6,color="black",zorder=5)
        ax.set_title(name);ax.set_xlabel("TS1_s descriptor (eV)");ax.set_ylabel("TS2_s descriptor (eV)")
        fig.colorbar(cf,ax=ax,label="CatMAP steady-state coverage")
    fig.suptitle("Existing CatMAP 0.3.1 TS1–TS2 grid coverage (50×50 computed nodes); actual catalysts overlaid")
    for ext in ["png","pdf"]: fig.savefig(OUT_FIG/f"ts1_ts2_grid_coverage_panels.{ext}",dpi=220) if ext=="png" else fig.savefig(OUT_FIG/f"ts1_ts2_grid_coverage_panels.{ext}")
    plt.close(fig)


def check_jacobian_fail(site, provenance):
    log = (NEW_DRC/site/f"catmap_{site}_mechanism_drc.log") if provenance.startswith("new_") else (OLD_DRC/site/f"catmap_{site}_strict.log")
    if not log.exists():
        log = NEW_DRC/site/"catmap_console.log" if provenance.startswith("new_") else OLD_DRC/site/"catmap_drc_console.log"
    if not log.exists():
        return "log_not_found"
    return bool(re.search(r"jacobian_fail|Jacobian\s+fail", log.read_text(encoding="utf-8", errors="replace"), re.I))


def plot_heatmap(matrix, xlabels, ylabels, title, cbar, path, cmap="viridis", center_zero=False, annotate=False):
    arr = np.asarray(matrix, dtype=float)
    fig, ax = plt.subplots(figsize=(13, 6.5), constrained_layout=True)
    cm=plt.get_cmap(cmap).copy();cm.set_bad("#bdbdbd")
    if center_zero:
        from matplotlib.colors import TwoSlopeNorm
        lim = max(abs(np.nanmin(arr)), abs(np.nanmax(arr)))
        im = ax.imshow(arr, aspect="auto", cmap=cm, norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim))
    else:
        im = ax.imshow(arr, aspect="auto", cmap=cm)
    ax.set_xticks(range(len(xlabels)), labels=xlabels, rotation=35, ha="right")
    ax.set_yticks(range(len(ylabels)), labels=ylabels)
    ax.set_title(title)
    cb = fig.colorbar(im, ax=ax); cb.set_label(cbar)
    if annotate and arr.size < 250:
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                ax.text(j, i, f"{arr[i,j]:.2g}", ha="center", va="center", fontsize=7, color="white" if abs(arr[i,j]) > np.nanmax(abs(arr))*.5 else "black")
    fig.savefig(path.with_suffix(".png"), dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def save_matrix_csv(df, rowkey, colkeys, path):
    out = df[[rowkey]+colkeys].copy()
    write_csv(path, out)


def make_figures(cats, coverage, coverage_wide, barriers, rates, drc_wide, tof_flux):
    labels = [NAMES[s] for s in SITES]
    # Coverage stacked columns with all stable adsorbates and vacancy.
    fig, ax = plt.subplots(figsize=(15, 7), constrained_layout=True)
    x = np.arange(len(SITES)); bottom = np.zeros(len(SITES))
    colors = {"EpoxideHHO_s":"#4c78a8","HHO_s":"#f58518","OHHa_s":"#54a24b","H2O_s":"#e45756","Rpair_s":"#b279a2","s":"#dddddd"}
    for species in ["EpoxideHHO_s","HHO_s","OHHa_s","H2O_s","Rpair_s","s"]:
        vals = np.array([float(coverage_wide.set_index("surface").loc[s,species]) for s in SITES])
        ax.bar(x, vals, bottom=bottom, label=species, color=colors[species], width=.78)
        bottom += vals
    ax.set_xticks(x, labels=labels, rotation=35, ha="right"); ax.set_ylabel("Steady-state coverage"); ax.set_ylim(0,1.02)
    ax.set_title("CatMAP official steady-state coverages (13 validated sites)"); ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(.5,-.18))
    for ext in ["png","pdf"]: fig.savefig(OUT_FIG/f"coverage_stacked.{ext}",dpi=220) if ext=="png" else fig.savefig(OUT_FIG/f"coverage_stacked.{ext}")
    plt.close(fig)
    covmat = np.array([[float(coverage_wide.set_index("surface").loc[s,k]) for k in ALL_COVERAGE] for s in SITES])
    plot_heatmap(covmat, ALL_COVERAGE, labels, "Steady-state coverage by catalyst", "Coverage", OUT_FIG/"coverage_heatmap")

    cat_index = cats.set_index("surface")
    logtof = np.array([cat_index.loc[s,"log10_tof"] for s in SITES])
    for sp, fname, xlabel in [("EpoxideHHO_s","epoxide_coverage_vs_tof","θ(EpoxideHHO*)"),("HHO_s","hho_coverage_vs_tof","θ(HHO*)"),("s","vacancy_coverage_vs_tof","θ(*) vacant site")]:
        xs = np.array([float(coverage_wide.set_index("surface").loc[s,sp]) for s in SITES])
        fig,ax=plt.subplots(figsize=(9,6),constrained_layout=True); ax.scatter(xs,logtof,c="#276FBF",edgecolor="black",zorder=3)
        for s,xv,yv in zip(SITES,xs,logtof): ax.annotate(NAMES[s],(xv,yv),xytext=(4,4),textcoords="offset points",fontsize=8)
        ax.set_xlabel(xlabel);ax.set_ylabel("log10(TOF C6H12O / s^-1)");ax.set_title(f"{xlabel} and product activity")
        for ext in ["png","pdf"]: fig.savefig(OUT_FIG/f"{fname}.{ext}",dpi=220) if ext=="png" else fig.savefig(OUT_FIG/f"{fname}.{ext}")
        plt.close(fig)

    # Flux heatmaps: rates are positive/near positive; log10 absolute values preserve visibility.
    for direction in ["forward","reverse"]:
        mat=[]
        for s in SITES:
            sub=rates[(rates.surface.eq(s))].sort_values("step")
            vals=[float(sub[f"{direction}_rate_as_exported"].iloc[i]) for i in range(6)]
            mat.append([math.log10(v) if v>0 else np.nan for v in vals])
        plot_heatmap(mat,[f"Step {i}" for i in range(1,7)],labels,f"CatMAP {direction} API rates; negative entries masked as invalid","log10(positive API output); gray = negative/invalid",OUT_FIG/f"{direction}_rate_heatmap")
    netmat=[]
    for s in SITES:
        sub=rates[rates.surface.eq(s)].sort_values("step")
        vals=np.array([float(v) for v in sub.net_rate_high_precision])
        # signed log transform, explicitly reported in axis label.
        netmat.append(np.sign(vals)*np.log10(1+np.abs(vals)))
    plot_heatmap(netmat,[f"Step {i}" for i in range(1,7)],labels,"CatMAP net reaction rates by site","sign(rate) × log10(1 + |rate|)",OUT_FIG/"net_rate_heatmap",cmap="coolwarm",center_zero=True)

    for direction,col in [("forward","forward_effective_barrier_eV"),("reverse","reverse_effective_barrier_eV")]:
        mat=[]
        for s in SITES:
            sub=barriers[barriers.surface.eq(s)].sort_values("step");mat.append([float(v) for v in sub[col]])
        plot_heatmap(mat,[f"Step {i}" for i in range(1,7)],labels,f"Effective {direction} barriers (eV)","Barrier (eV)",OUT_FIG/f"{direction}_barrier_heatmap",cmap="magma")

    # Product TOF against each step's directional rates.
    fig,axes=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
    for step,ax in enumerate(axes.flat,1):
        sub=rates[(rates.step.eq(step))].set_index("surface").loc[SITES]
        vals=np.array([float(v) for v in sub.net_rate_high_precision])
        ax.scatter(np.log10(np.maximum(np.abs(vals),1e-300)),logtof,c="#138A72",edgecolor="black")
        for s,xv,yv in zip(SITES,np.log10(np.maximum(np.abs(vals),1e-300)),logtof): ax.annotate(NAMES[s],(xv,yv),xytext=(3,3),textcoords="offset points",fontsize=7)
        ax.set_xlabel(f"Step {step} log10(|net rate|)");ax.set_ylabel("log10(product TOF)");ax.set_title(REACTION_LABELS[step])
    fig.suptitle("Step flux and product TOF (same strict steady-state output)")
    for ext in ["png","pdf"]: fig.savefig(OUT_FIG/f"step_net_flux_vs_tof.{ext}",dpi=220) if ext=="png" else fig.savefig(OUT_FIG/f"step_net_flux_vs_tof.{ext}")
    plt.close(fig)

    dmat=np.array([[float(drc_wide.set_index("surface").loc[s,p]) for p in DRC_PARAMS_REPORT] for s in SITES])
    plot_heatmap(dmat,DRC_PARAMS_REPORT,labels,"Official CatMAP degree of rate control for C6H12O","Degree of rate control",OUT_FIG/"product_drc_heatmap",cmap="coolwarm",center_zero=True)
    # Per-parameter DRC vs activity: compact small multiples.
    fig,axes=plt.subplots(3,3,figsize=(15,12),constrained_layout=True)
    for i,p in enumerate(DRC_PARAMS_REPORT):
        ax=axes.flat[i];xs=np.array([float(drc_wide.set_index("surface").loc[s,p]) for s in SITES]);ax.scatter(xs,logtof,c="#8e5ea2",edgecolor="black")
        for s,xv,yv in zip(SITES,xs,logtof):ax.annotate(NAMES[s],(xv,yv),xytext=(3,3),textcoords="offset points",fontsize=7)
        ax.axvline(0,color="grey",lw=.8);ax.set_xlabel(p+" DRC");ax.set_ylabel("log10(TOF)")
    for ax in axes.flat[len(DRC_PARAMS_REPORT):]: ax.axis("off")
    fig.suptitle("Official local DRC and product activity; catalyst points are not grid points")
    for ext in ["png","pdf"]: fig.savefig(OUT_FIG/f"drc_vs_tof_small_multiples.{ext}",dpi=220) if ext=="png" else fig.savefig(OUT_FIG/f"drc_vs_tof_small_multiples.{ext}")
    plt.close(fig)

    # Highlight mechanism-sensitive states for representative sites, without asserting a universal RDS.
    focus=["ti","tico","tini","tire","tife"]
    fig,ax=plt.subplots(figsize=(11,6),constrained_layout=True)
    sel=drc_wide.set_index("surface").loc[focus,DRC_PARAMS_REPORT]
    width=.14; xpos=np.arange(len(DRC_PARAMS_REPORT))
    for j,s in enumerate(focus):ax.bar(xpos+(j-2)*width,[float(v) for v in sel.loc[s]],width,label=NAMES[s])
    ax.axhline(0,color="black",lw=.8);ax.set_xticks(xpos,labels=DRC_PARAMS_REPORT,rotation=30,ha="right");ax.set_ylabel("Official DRC for C6H12O");ax.set_title("Official DRC across five representative catalysts");ax.legend(ncol=3)
    for ext in ["png","pdf"]:fig.savefig(OUT_FIG/f"focus_site_drc.{ext}",dpi=220) if ext=="png" else fig.savefig(OUT_FIG/f"focus_site_drc.{ext}")
    plt.close(fig)


def main():
    for p in [OUT_DATA,OUT_FIG,OUT_ORIGIN,OUT_REPORT]:p.mkdir(parents=True,exist_ok=True)
    (cats,cov,covwide,balance,rates,barriers,gas,drc,drcwide,drcdiag,tofflux,constants,conv,hashes)=collect_outputs()
    signed_neb=collect_signed_neb_differences(barriers)
    reaction_energy=barriers[["surface","catalyst","step","reaction_label","barrier_model","delta_G_eV"]].copy()
    grid_cov,grid_sides,grid_anchors,grid_summary,grid_sources=collect_existing_grid_coverage()
    for f in grid_sources:hashes.loc[len(hashes)]={"source_file":str(f),"sha256":sha(f),"role":"existing CatMAP TS1-TS2 grid cache; read-only"}
    directional_diag=rates.groupby("step").agg(negative_forward_rows=("forward_rate_negative","sum"),negative_reverse_rows=("reverse_rate_negative","sum"),
        maximum_scaled_closure_error=("scaled_directional_closure_error",lambda x:max(float(v) for v in x))).reset_index()
    products={"coverage_13sites_long.csv":cov,"coverage_13sites.csv":covwide,"coverage_site_balance.csv":balance,
        "coverage_dominant_species.csv":balance[["surface","catalyst","dominant_adsorbate","dominant_adsorbate_coverage","dominant_including_vacancy","dominant_all_coverage"]],
        "directional_reaction_rates_13sites.csv":rates,"effective_barriers_13sites.csv":barriers,
        "gas_turnover_frequencies_13sites.csv":gas,"rate_control_13sites.csv":drc,"rate_control_matrix.csv":drcwide,
        "rate_control_diagnostics.csv":drcdiag,"TOF_flux_consistency.csv":tofflux,"directional_rate_api_diagnostics.csv":directional_diag,
        "reaction_energy_13sites.csv":reaction_energy,"explicit_TS_signed_energy_differences.csv":signed_neb,
        "ts1_ts2_grid_coverage_2500nodes.csv":grid_cov,"ts1_ts2_grid_peak_and_sides.csv":grid_sides,
        "reaction_rate_constants_13sites.csv":constants,"convergence_13sites.csv":conv,"source_SHA256.csv":hashes}
    for filename,df in products.items():write_csv(OUT_DATA/filename,df)
    # Origin-ready wide matrices and point tables.
    save_matrix_csv(covwide,"catalyst",ALL_COVERAGE,OUT_ORIGIN/"coverage_matrix.csv")
    save_matrix_csv(drcwide,"catalyst",DRC_PARAMS_REPORT,OUT_ORIGIN/"product_DRC_matrix.csv")
    write_csv(OUT_ORIGIN/"directional_and_net_flux_long.csv",rates)
    write_csv(OUT_ORIGIN/"effective_barriers_matrix_forward.csv",barriers.pivot(index="catalyst",columns="step",values="forward_effective_barrier_eV").reset_index())
    write_csv(OUT_ORIGIN/"effective_barriers_matrix_reverse.csv",barriers.pivot(index="catalyst",columns="step",values="reverse_effective_barrier_eV").reset_index())
    write_csv(OUT_ORIGIN/"activity_and_coverage_points.csv",cats[["surface","catalyst","tof","log10_tof","status"]].merge(covwide,on=["surface","catalyst"]))
    write_csv(OUT_ORIGIN/"ts1_ts2_grid_coverage_2500nodes_xyz.csv",grid_cov)
    for species in ["EpoxideHHO_s","H2O_s","HHO_s","OHHa_s","Rpair_s","s"]:
        m=grid_cov.pivot(index="TS2_s_eV",columns="TS1_s_eV",values=species).reset_index()
        write_csv(OUT_ORIGIN/f"ts1_ts2_{species}_coverage_matrix_50x50.csv",m)
    make_figures(cats,cov,covwide,barriers,rates,drcwide,tofflux)
    plot_grid_coverages(grid_cov,grid_anchors)
    # Compute simple deterministic evidence counts for concise report.
    balance_err=max(abs(float(v)) for v in balance.site_balance_error)
    tof_err=max(abs(float(v)) for v in tofflux.step3_minus_product_TOF)
    flux_scaled_errors=[abs(float(v)) for v in rates.scaled_directional_closure_error if str(v).lower()!="nan"]
    max_directional_scaled=max(flux_scaled_errors)
    negative_forward_count=sum(bool(v) for v in rates.forward_rate_negative)
    negative_reverse_count=sum(bool(v) for v in rates.reverse_rate_negative)
    drc_sum_min=min(float(v) for v in drcdiag.drc_sum);drc_sum_max=max(float(v) for v in drcdiag.drc_sum)
    conv_ok=int((cats.status=="CONVERGED").sum())
    ep_cov=np.array([float(covwide.set_index("surface").loc[s,"EpoxideHHO_s"]) for s in SITES])
    hho_cov=np.array([float(covwide.set_index("surface").loc[s,"HHO_s"]) for s in SITES])
    ohha_cov=np.array([float(covwide.set_index("surface").loc[s,"OHHa_s"]) for s in SITES])
    top_ep=[NAMES[SITES[i]] for i in np.argsort(ep_cov)[-3:][::-1]]
    top_hho=[NAMES[SITES[i]] for i in np.argsort(hho_cov)[-3:][::-1]]
    top_ohha=[NAMES[SITES[i]] for i in np.argsort(ohha_cov)[-3:][::-1]]
    grid_peak_row=grid_sides[grid_sides.grid_location_role.eq("grid_peak")].iloc[0]
    grid_peak_log=float(grid_peak_row.CatMAP_log10_TOF_C6H12O)
    real_distance=np.sqrt((grid_anchors.TS1_s_eV.astype(float)-float(grid_peak_row.TS1_s_eV))**2+(grid_anchors.TS2_s_eV.astype(float)-float(grid_peak_row.TS2_s_eV))**2)
    nearest_idx=int(np.argmin(real_distance));nearest_real=grid_anchors.iloc[nearest_idx]
    peak_occ={k:float(grid_peak_row[k]) for k in ["EpoxideHHO_s","H2O_s","HHO_s","OHHa_s","Rpair_s","s"]}
    report=f"""# 13个Ti–M位点 CatMAP 微观动力学机制结果

## 范围与来源

- 计算：本报告复用已完成的 CatMAP 0.3.1 严格容差六步稳态输出；未重跑稳态基线。
- 位点：13个已验证位点；Ti–Zr未纳入。
- 条件：333.15 K；固定基线储库输入 H2O2=0.185、C6H12=1.0、H2O=0.815、C6H12O=1e-20；模型中的压力描述符=1.0；其单位/物理量遵循基线配置，不在此处重新解释。
- 机理：原六步反应网络不变；第2/4步保留显式NEB TS1/TS2，第1/3/5/6步为 CatMAP 端点势垒近似。
- DRC：5个位点复用既有官方 `rate_control` 输出；另8个位点通过 CatMAP 0.3.1 官方 `ReactionModel.run(recalculate=True)`、仅请求 `rate_control` 的隔离补充运行取得。该补充没有改动能量、机理、温度或储库。
- 主要输入与输出路径及SHA256见 `../data/source_SHA256.csv`；补充DRC运行日志、配置、能量表位于 `../data/official_drc_runs_supplemental/`。

## 数值验收摘要

- 严格稳态：{conv_ok}/13 为 `CONVERGED`。
- 最大覆盖度位点守恒误差：{balance_err:.3e}；覆盖度检查同时保留高精度十进制原值。
- 第3步净速率与 C6H12O 气相净TOF最大差：{tof_err:.3e}（以导出十进制值比较）。还用CatMAP CSV中每步反应两端的实际物种计量系数求和核算了六步网络产物速率；结果见 `TOF_flux_consistency.csv`。
- 原方向速率CSV表头有3列而数据行实际为2列；脚本按原始数据行结构解析。方向API输出含负值：forward {negative_forward_count} 行、reverse {negative_reverse_count} 行；最大缩放闭合差为 {max_directional_scaled:.3e}。本机CatMAP 0.3.1 `mean_field_solver.py` 中 `get_empty_site_cvgs()` 将覆盖度转成Python `float` 后计算空位，`get_elem_ec()` 又将速率常数、覆盖度和结果转成`float`；当空位接近零时会产生负空位/负方向输出。因此第1、6步方向API数值在这些体系中不能作为可信单向物理通量；原始值保留并在诊断CSV中标记。此问题不推翻高精度稳态净速率、覆盖度或官方TOF。
- DRC输出覆盖：13/13位点 × 7个表面态/过渡态参数；DRC总和范围={drc_sum_min:.6g} 至 {drc_sum_max:.6g}。总和偏离1时不强制归一化，保留CatMAP原值。
- 已有TS1–TS2 50×50 CatMAP缓存的pickle中确实含2500个`coverage_map`节点，且与既有TOF节点按描述符坐标一对一匹配；本报告复用该缓存，没有重跑网格。覆盖度向量按CatMAP日志的5个`adsorbate_names`顺序加末尾空位`s`核验，逐节点总和为1。峰值及TS1/TS2轴向两侧代表节点详见`ts1_ts2_grid_peak_and_sides.csv`；完整覆盖度长表和Origin矩阵一并导出。虚拟网格数据与13个真实催化剂数据分开保存。
- 速率常数：基线CSV没有直接导出每步正/逆速率常数，因此 `reaction_rate_constants_13sites.csv` 将其记为 `NOT_EXPORTED_IN_EXISTING_BASELINE_CACHE`；未由净速率反推或填造。

## 机制结果摘要

- 覆盖度最大稳定吸附态按位点列于 `coverage_dominant_species.csv`（见覆盖平衡表）；代表性最高覆盖的前三个位点：EpoxideHHO*—{', '.join(top_ep)}；HHO*—{', '.join(top_hho)}；OHHa*—{', '.join(top_ohha)}。请结合逐位点覆盖度，不把最大覆盖物种直接等同于速率控制步骤。
- 已有虚拟网格离散最大点在 TS1={float(grid_peak_row.TS1_s_eV):.6f} eV、TS2={float(grid_peak_row.TS2_s_eV):.6f} eV，log10TOF={grid_peak_log:.4f}；峰处θ(EpoxideHHO,H2O,HHO,OHHa,Rpair,*)=({peak_occ['EpoxideHHO_s']:.3g},{peak_occ['H2O_s']:.3g},{peak_occ['HHO_s']:.3g},{peak_occ['OHHa_s']:.3g},{peak_occ['Rpair_s']:.3g},{peak_occ['s']:.3g})。沿TS1和TS2两轴各±2格时，H2O/Rpair/EpoxideHHO覆盖重新分配（见5点及2500节点明细）；最近真实催化剂是{nearest_real.catalyst}，描述符平面欧氏距离={float(real_distance.iloc[nearest_idx]):.3f} eV，真实log10TOF={float(nearest_real.actual_log10_TOF_C6H12O):.3f}。真实位点并未达到虚拟网格峰附近的预测活性，不能把该离散峰等同于实际催化剂的可实现性能。
- 全部13个位点的产物TOF、六步高精度净通量、方向API原始输出、有效正逆势垒和官方DRC分别见数据表与图；方向速率图把负API值遮罩为灰色，不解释为真实通量。
- 有效势垒图中第1/3/5/6步明确标为 `CatMAP_endpoint`，不能解释为DFT独立过渡态能垒；第2/4步为输入的显式NEB TS势垒。
- TS1−Rpair与TS2−HHO的原始有符号差值逐位点见 `explicit_TS_signed_energy_differences.csv`；负值原样保留，和CatMAP有效正/逆势垒分列，不把负差强制改成零。
- DRC是CatMAP官方局部能量敏感度，不等同于全局限速步骤判定。极低TOF体系中的有限差分DRC可能受数值病态影响；诊断表提供来源/日志检查及低TOF提醒。尚无导出的状态Jacobian特征值或多初值吸引域证据，故不把这些稳态结果宣称为动力学可达性证明。
- 机理读法：Ti–Fe 的OHHa覆盖约0.910且OHHa DRC约−0.906，指向对后段OHHa能量/积累敏感；Ti的TS1 DRC约+0.999，是本组中最明显的显式NEB势垒敏感例。多数异核体系（含Ti–W、Ti–Ta）EpoxideHHO接近饱和，且HHO DRC约+1、EpoxideHHO DRC约−1，显示活性响应主要受吸附态能量平衡和中间体释放/转化共同影响；这不是“最高覆盖=唯一限速步”的证明。Ti–Hf、Ti–Ta的TOF极低，DRC仅作警示性局部敏感度，不作强因果结论。

## 13个位点汇总

| 位点 | CatMAP状态 | TOF (s⁻¹) | log10(TOF) | 主导稳定吸附态 | 覆盖度 | 最大|净通量|相对离散 |
|---|---:|---:|---:|---|---:|---:|
"""
    bal=balance.set_index("surface");
    for _,row in cats.iterrows():
        b=bal.loc[row.surface]
        report+=f"| {row.catalyst} | {row.status} | {row.tof:.8g} | {row.log10_tof:.5f} | {b.dominant_adsorbate} | {float(b.dominant_adsorbate_coverage):.8g} | {float(row.flux_spread_relative):.3e} |\n"
    report+="\n## 各位点最大绝对DRC参数（仅作局部敏感度摘要）\n\n| 位点 | 参数 | DRC | TOF (s^-1) |\n|---|---|---:|---:|\n"
    ddiag=drcdiag.set_index("surface")
    for _,row in cats.iterrows():
        z=ddiag.loc[row.surface]
        report+=f"| {row.catalyst} | {z.largest_absolute_DRC_parameter} | {float(z.largest_absolute_DRC_value):.6g} | {row.tof:.6g} |\n"
    report+="""

## DRC解释边界

官方D(R)C保留原始符号：正值表示提高该能量在局部会降低产物速率；负值表示相反的局部响应。不同表面态之间存在耦合，单项DRC不等于唯一决速步。对于接近数值下限的TOF，不宜把很大的、微小或符号敏感的有限差分值作强机理断言。具体逐态数据见 `rate_control_matrix.csv` 和 `rate_control_diagnostics.csv`。

## 可复现性与图表

- 复现分析：运行 `scripts/analyze_mechanism_outputs.py`；该脚本只读取CSV、组表和绘图，不调用稳态求解器。
- 补充官方DRC：`scripts/run_missing_official_drc.py`；运行副本与日志保存在本分析目录的 `data/official_drc_runs_supplemental/`。
- H盘目录是归档副本；脚本仍记录/引用本次C盘工作区、CatMAP安装路径及既有网格缓存路径。没有做路径迁移或从H盘独立重放，因此归档副本用于审计保存，不宣称可脱离当前主机目录独立运行。
- PNG/PDF图位于 `../figures/`；Origin矩阵/长表位于 `../origin_exports/`。
- 方向速率由基线CatMAP接口以双精度导出；净通量、覆盖度和TOF保留高精度十进制文本。除格式精度差异外，本数据中第1/6步还出现负方向速率及显著通量闭合异常，根因与官方方向速率路径对近零空位的双精度计算有关；这些方向值不作物理单向通量解释。高精度稳态净通量、覆盖度及TOF的验证独立保留。
"""
    (OUT_REPORT/"MECHANISM_ANALYSIS_REPORT.md").write_text(report,encoding="utf-8")
    # Include a deterministic output SHA256 manifest, excluding itself.
    manifest=[]
    for base in [OUT_DATA,OUT_FIG,OUT_ORIGIN,OUT_REPORT,HERE/"scripts"]:
        for f in sorted(base.rglob("*")):
            if f.is_file() and f.name!="OUTPUT_SHA256.csv":manifest.append({"path":str(f.relative_to(HERE)),"size_bytes":f.stat().st_size,"sha256":sha(f)})
    write_csv(OUT_DATA/"OUTPUT_SHA256.csv",manifest)
    print(json.dumps({"sites":len(cats),"converged":conv_ok,"coverage_balance_max_abs":balance_err,"step3_product_tof_max_abs_diff":tof_err,"drc_rows":len(drc),"drc_sum_range":[drc_sum_min,drc_sum_max],"figures":len(list(OUT_FIG.glob("*.png"))),"files":len(manifest)+1,"report":str(OUT_REPORT/"MECHANISM_ANALYSIS_REPORT.md")},ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
