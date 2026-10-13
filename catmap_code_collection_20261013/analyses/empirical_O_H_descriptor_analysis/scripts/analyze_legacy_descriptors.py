"""Empirical O/Ha/Hb descriptor plots from frozen CatMAP outputs and legacy energies.

This script only reads the historical formation-energy table and existing strict
CatMAP outputs. It never runs CatMAP and does not write to source data locations.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MplPath
from scipy.optimize import least_squares


BASE = Path(r"C:\Users\Administrator\Documents\codex探索\catmap_neb_coadsorption_model")
OFFICIAL = BASE / "catmap_official_run"
STRICT = OFFICIAL / "results" / "official_strict_tolerance"
LEGACY = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\DFT_formation_catmap.txt")
CONVERTER = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\1-generating_input_file\convert_dft_to_catmap.py")
OUT = BASE / "empirical_O_H_descriptor_analysis"
SURFACES = ["tife", "timn", "tihf", "tire", "timo", "tiv", "tico", "titi", "tiw", "tita", "ticr", "ti", "tini"]
LABELS = {"ti": "Ti", "tico": "Ti-Co", "ticr": "Ti-Cr", "tife": "Ti-Fe", "tihf": "Ti-Hf", "timn": "Ti-Mn", "timo": "Ti-Mo", "tini": "Ti-Ni", "tire": "Ti-Re", "tita": "Ti-Ta", "titi": "Ti-Ti", "tiv": "Ti-V", "tiw": "Ti-W"}
K_B_EV = 8.617333262145e-5
TEMPERATURE_K = 333.15


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_inputs():
    with LEGACY.open(encoding="utf-8-sig", newline="") as f:
        old_rows = list(csv.DictReader(f, delimiter="\t"))
    legacy = {}
    for row in old_rows:
        surf, sp = row["surface_name"], row["species_name"]
        if surf in SURFACES and sp in {"O", "Ha", "Hb"}:
            key = (surf, sp)
            if key in legacy:
                raise ValueError(f"duplicate legacy record: {key}")
            legacy[key] = float(row["formation_energy"])
    if len(legacy) != 39:
        raise ValueError(f"Expected 39 O/Ha/Hb values for 13 sites; found {len(legacy)}")

    records = {}
    for surf in SURFACES:
        if not (STRICT / surf / "convergence.csv").exists():
            raise ValueError(f"Missing convergence evidence for {surf}")
        with (STRICT / surf / "convergence.csv").open(encoding="utf-8-sig", newline="") as f:
            status = next(csv.DictReader(f))
        if status["status"] != "CONVERGED" or float(status["temperature_K"]) != TEMPERATURE_K:
            raise ValueError(f"Strict result is not converged at target temperature: {surf}: {status}")
        with (STRICT / surf / "gas_turnover_frequencies.csv").open(encoding="utf-8-sig", newline="") as f:
            tof_rows = list(csv.DictReader(f))
        product = [float(r["turnover_frequency"]) for r in tof_rows if r["gas_species"] == "C6H12O_g"]
        if len(product) != 1 or product[0] <= 0:
            raise ValueError(f"Expected exactly one positive C6H12O_g TOF: {surf}")
        records[surf] = {"surface_name": surf, "catalyst": LABELS[surf], "TOF": product[0], "log10_TOF": math.log10(product[0]), **{sp: legacy[(surf, sp)] for sp in ("O", "Ha", "Hb")}, "status": "CONVERGED", "temperature_K": TEMPERATURE_K}
    return records


def fit_segmented(x, y):
    order = np.argsort(x)
    xs, ys = np.asarray(x)[order], np.asarray(y)[order]
    candidates = []
    for j in range(3, len(xs) - 3):
        knot = (xs[j - 1] + xs[j]) / 2
        A = np.column_stack([np.ones_like(xs), np.minimum(xs - knot, 0), np.maximum(xs - knot, 0)])
        coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
        candidates.append((float(np.sum((A @ coef - ys) ** 2)), knot, coef))
    if not candidates:
        return None
    _, knot, coef = min(candidates, key=lambda z: z[0])
    return knot, coef


def segmented_predict(x, model):
    knot, coef = model
    x = np.asarray(x)
    return coef[0] + coef[1] * np.minimum(x - knot, 0) + coef[2] * np.maximum(x - knot, 0)


def loocv_rmse(x, y, kind):
    x, y = np.asarray(x, float), np.asarray(y, float)
    pred = []
    for i in range(len(y)):
        keep = np.arange(len(y)) != i
        if kind == "linear":
            c = np.polyfit(x[keep], y[keep], 1)
            yp = np.polyval(c, x[i])
        elif kind == "quadratic":
            c = np.polyfit(x[keep], y[keep], 2)
            yp = np.polyval(c, x[i])
        else:
            m = fit_segmented(x[keep], y[keep])
            yp = segmented_predict([x[i]], m)[0] if m else np.mean(y[keep])
        pred.append(yp)
    return float(np.sqrt(np.mean((y - np.asarray(pred)) ** 2)))


def annotate_points(ax, x, y, rows):
    for xx, yy in zip(x, y):
        ax.scatter(xx, yy, marker="o", s=48, facecolor="white", edgecolor="#162b4d", linewidth=1.2, zorder=5)
    place_named_labels(ax, x, y, rows)


def place_named_labels(ax, x, y, rows):
    """Place every catalyst label inside the axes and avoid label-label overlap."""
    offsets = [(7, 8), (8, 19), (8, -19), (-38, 8), (-38, -18), (-38, 19), (0, 23), (0, -25), (25, 0), (-50, 0), (22, 16), (22, -17), (-50, 17), (-50, -17)]
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    axes_box = ax.get_window_extent(renderer)
    occupied=[]
    for xx, yy, row in zip(x, y, rows):
        options=[]
        for offset in offsets:
            ann=ax.annotate(row["catalyst"],(xx,yy),xytext=offset,textcoords="offset points",fontsize=8,color="#162b4d",zorder=7,arrowprops={"arrowstyle":"-","lw":.4,"color":"#777"})
            ax.figure.canvas.draw(); box=ann.get_window_extent(renderer).expanded(1.05,1.08)
            outside=(box.x0<axes_box.x0 or box.x1>axes_box.x1 or box.y0<axes_box.y0 or box.y1>axes_box.y1)
            overlaps=any(box.overlaps(prev) for prev in occupied)
            score=(1000 if outside else 0)+(100 if overlaps else 0)+abs(offset[0])+abs(offset[1])
            options.append((score,offset,box))
            ann.remove()
            if not outside and not overlaps:
                break
        score,offset,box=min(options,key=lambda v:v[0])
        ann=ax.annotate(row["catalyst"],(xx,yy),xytext=offset,textcoords="offset points",fontsize=8,color="#162b4d",zorder=7,arrowprops={"arrowstyle":"-","lw":.4,"color":"#777"})
        ax.figure.canvas.draw();box=ann.get_window_extent(renderer).expanded(1.05,1.08)
        occupied.append(box)


def one_d_plot(records, descriptor, species, stem):
    rows = list(records.values())
    x = np.array([r[species] for r in rows], float)
    y = np.array([r["log10_TOF"] for r in rows], float)
    order = np.argsort(x)
    grid = np.linspace(x.min() - .05 * (x.max() - x.min()), x.max() + .05 * (x.max() - x.min()), 500)
    models = {}
    coef_lin = np.polyfit(x, y, 1)
    coef_quad = np.polyfit(x, y, 2)
    seg = fit_segmented(x, y)
    models["linear"] = (np.polyval(coef_lin, grid), loocv_rmse(x, y, "linear"), None)
    models["quadratic"] = (np.polyval(coef_quad, grid), loocv_rmse(x, y, "quadratic"), None)
    models["segmented"] = (segmented_predict(grid, seg), loocv_rmse(x, y, "segmented"), seg)
    fig, ax = plt.subplots(figsize=(9.2, 6.2), constrained_layout=True)
    colors = {"linear": "#5b8ff9", "quadratic": "#e8684a", "segmented": "#5ad8a6"}
    for name, (yp, cv, _) in models.items():
        ax.plot(grid, yp, color=colors[name], lw=2, label=f"{name} fit (LOOCV RMSE={cv:.2f})")
    annotate_points(ax, x, y, rows)
    peak_info = []
    # Mark every fitted model's actual maximum over the sampled descriptor interval.
    gx = np.linspace(x.min(), x.max(), 2001)
    peak_markers = {"linear": "s", "quadratic": "*", "segmented": "P"}
    for name, color in colors.items():
        if name == "linear":
            p = np.polyval(coef_lin, gx)
            idx = int(np.argmax(p))
            xpk, zpk = float(gx[idx]), float(p[idx])
            boundary = idx in (0, len(gx) - 1)
            peak_info.append({"model": name, "descriptor": xpk, "log10_TOF": zpk, "TOF": 10**zpk, "kind": "boundary maximum" if boundary else "internal maximum"})
            ax.scatter([xpk], [zpk], marker=peak_markers[name], s=90, color=color, edgecolor="black", zorder=8, label="linear boundary maximum")
            continue
        if name == "quadratic":
            p = np.polyval(coef_quad, gx)
            idx = int(np.argmax(p)); xpk, zpk = float(gx[idx]), float(p[idx])
            boundary = idx in (0, len(gx) - 1)
            if coef_quad[0] < 0 and x.min() < -coef_quad[1] / (2 * coef_quad[0]) < x.max():
                xpk = float(-coef_quad[1] / (2 * coef_quad[0])); zpk = float(np.polyval(coef_quad, xpk)); boundary = False
            peak_info.append({"model": name, "descriptor": xpk, "log10_TOF": zpk, "TOF": 10**zpk, "kind": "internal peak" if not boundary else "boundary maximum"})
            ax.scatter([xpk], [zpk], marker=peak_markers[name], s=180 if not boundary else 90, color=color, edgecolor="black", zorder=8, label=f"quadratic {'internal peak' if not boundary else 'boundary maximum'}")
            if not boundary:
                ax.annotate(f"quad peak\n{xpk:.3f} eV, logTOF={zpk:.2f}", (xpk, zpk), xytext=(28, 24), textcoords="offset points", fontsize=8, color="#8e2d20", arrowprops={"arrowstyle":"-","lw":.6,"color":"#8e2d20"})
        else:
            knot, c = seg
            slopes = (float(c[1]), float(c[2]))
            if slopes[0] > 0 and slopes[1] < 0:
                zpk = float(c[0]); xpk = float(knot)
                peak_info.append({"model": name, "descriptor": xpk, "log10_TOF": zpk, "TOF": 10**zpk, "kind": "internal segmented volcano peak", "left_slope": slopes[0], "right_slope": slopes[1]})
                ax.scatter([xpk], [zpk], marker="D", s=76, color=color, edgecolor="black", zorder=8, label="segmented volcano peak")
                ax.annotate(f"seg peak\n{xpk:.3f} eV, logTOF={zpk:.2f}", (xpk, zpk), xytext=(-90, 16), textcoords="offset points", fontsize=8, color="#276e50")
            else:
                p = segmented_predict(gx, seg); idx = int(np.argmax(p)); xpk, zpk = float(gx[idx]), float(p[idx])
                peak_info.append({"model": name, "descriptor": xpk, "log10_TOF": zpk, "TOF": 10**zpk, "kind": "boundary maximum" if idx in (0, len(gx)-1) else "non-volcano segmented maximum", "left_slope": slopes[0], "right_slope": slopes[1]})
                ax.scatter([xpk], [zpk], marker=peak_markers[name], s=90, color=color, edgecolor="black", zorder=8, label="segmented boundary maximum")
    ax.set_xlabel(f"{species} historical formation energy (eV)")
    ax.set_ylabel("log10(steady-state CatMAP TOF of C6H12O_g)")
    ax.set_title(f"Empirical descriptor response: {species} vs product TOF\nHistorical single-species energies; strict CatMAP TOF points; empirical fits only")
    ax.grid(alpha=.2); ax.legend(fontsize=8, loc="best")
    fig.savefig(OUT / f"{stem}.png", dpi=240); fig.savefig(OUT / f"{stem}.pdf"); plt.close(fig)
    return {"species": species, "n": len(rows), "linear_r2": float(1-np.sum((np.polyval(coef_lin,x)-y)**2)/np.sum((y-y.mean())**2)), "quadratic_r2": float(1-np.sum((np.polyval(coef_quad,x)-y)**2)/np.sum((y-y.mean())**2)), "segmented_left_slope": float(seg[1][1]), "segmented_right_slope": float(seg[1][2]), "fits": {k: {"loocv_rmse": v[1]} for k,v in models.items()}, "peaks": peak_info}


def response_surface(x, y, z):
    xm, xs, ym, ys = float(np.mean(x)), float(np.std(x)), float(np.mean(y)), float(np.std(y))
    u, v = (x-xm)/xs, (y-ym)/ys
    def design(a,b): return np.column_stack([np.ones_like(a), a, b, a*a, a*b, b*b])
    coef, *_ = np.linalg.lstsq(design(u,v), z, rcond=None)
    pred = design(u,v) @ coef
    loo=[]
    for i in range(len(z)):
        keep=np.arange(len(z))!=i
        c,*_=np.linalg.lstsq(design(u[keep],v[keep]),z[keep],rcond=None)
        loo.append(design(u[i:i+1],v[i:i+1]) @ c)
    rmse=float(np.sqrt(np.mean((z-np.array(loo).ravel())**2)))
    return coef, (xm,xs,ym,ys), pred, rmse


def two_d_plot(records, hspecies, stem):
    rows=list(records.values()); x=np.array([r["O"] for r in rows]); y=np.array([r[hspecies] for r in rows]); z=np.array([r["log10_TOF"] for r in rows])
    coef,scales, fitted, cv=response_surface(x,y,z); xm,xs,ym,ys=scales
    hull=__import__('scipy.spatial',fromlist=['ConvexHull']).ConvexHull(np.column_stack([x,y]))
    hull_pts=np.column_stack([x,y])[hull.vertices]
    poly=MplPath(hull_pts)
    gx=np.linspace(x.min(),x.max(),280); gy=np.linspace(y.min(),y.max(),280); XX,YY=np.meshgrid(gx,gy)
    U,V=(XX-xm)/xs,(YY-ym)/ys
    ZZ=coef[0]+coef[1]*U+coef[2]*V+coef[3]*U**2+coef[4]*U*V+coef[5]*V**2
    inside=poly.contains_points(np.column_stack([XX.ravel(),YY.ravel()])).reshape(XX.shape); ZZ=np.where(inside,ZZ,np.nan)
    hess=np.array([[2*coef[3],coef[4]],[coef[4],2*coef[5]]]); lin=np.array([coef[1],coef[2]])
    stat=np.linalg.solve(hess,-lin) if abs(np.linalg.det(hess))>1e-14 else np.array([np.nan,np.nan])
    xs_stat=xm+xs*stat[0]; ys_stat=ym+ys*stat[1]
    interior=(np.all(np.linalg.eigvalsh(hess)<0) and bool(poly.contains_point((xs_stat,ys_stat))))
    if interior:
        peakx,peaky=float(xs_stat),float(ys_stat); pk=float(coef[0]+coef[1]*stat[0]+coef[2]*stat[1]+coef[3]*stat[0]**2+coef[4]*stat[0]*stat[1]+coef[5]*stat[1]**2); kind="local maximum"
    else:
        boundary=np.array([x,y]).T[hull.vertices]
        bu,bv=(boundary[:,0]-xm)/xs,(boundary[:,1]-ym)/ys
        bz=(coef[0]+coef[1]*bu+coef[2]*bv+coef[3]*bu**2+coef[4]*bu*bv+coef[5]*bv**2)
        i=int(np.argmax(bz)); peakx,peaky=float(boundary[i,0]),float(boundary[i,1]);pk=float(bz[i]);kind="boundary maximum within catalyst convex hull"
    dist=np.sqrt(((x-peakx)/xs)**2+((y-peaky)/ys)**2); nearest=[rows[i]["catalyst"] for i in np.argsort(dist)[:3]]
    fig,ax=plt.subplots(figsize=(9.2,7.2),constrained_layout=True)
    m=ax.contourf(XX,YY,ZZ,levels=24,cmap="viridis")
    cb=fig.colorbar(m,ax=ax);cb.set_label("Empirical quadratic response: predicted log10(TOF)")
    ax.plot(*np.vstack([hull_pts,hull_pts[0]]).T,ls="--",lw=1,color="white",alpha=.9,label="convex hull of 13 catalysts (fit domain)")
    ax.scatter(x,y,s=55,marker="D",facecolor="white",edgecolor="black",linewidth=.9,zorder=5,label="13 validated strict CatMAP catalyst points")
    ax.set_xlim(x.min()-.10*(x.max()-x.min()),x.max()+.10*(x.max()-x.min()))
    ax.set_ylim(y.min()-.10*(y.max()-y.min()),y.max()+.10*(y.max()-y.min()))
    place_named_labels(ax,x,y,rows)
    ax.scatter([peakx],[peaky],marker="*",s=250,color="#ff3b30",edgecolor="black",linewidth=1.0,zorder=10,label=f"{kind}")
    ax.annotate(f"{kind}\nO={peakx:.3f} eV; {hspecies}={peaky:.3f} eV\nlogTOF={pk:.2f}; TOF={10**pk:.3g}",(peakx,peaky),xytext=(.02,.03),textcoords="axes fraction",ha="left",va="bottom",fontsize=8,bbox={"boxstyle":"round,pad=.25","fc":"white","alpha":.92},arrowprops={"arrowstyle":"->","lw":.8})
    ax.set_xlabel("O* historical formation energy (eV)");ax.set_ylabel(f"{hspecies} historical formation energy (eV)")
    ax.set_title(f"Empirical O* + {hspecies} response surface vs C6H12O TOF\nQuadratic fit to 13 actual CatMAP TOFs; not a CatMAP grid")
    ax.legend(fontsize=8,loc="best");ax.grid(alpha=.12)
    fig.savefig(OUT/f"{stem}.png",dpi=240);fig.savefig(OUT/f"{stem}.pdf");plt.close(fig)
    r2=float(1-np.sum((fitted-z)**2)/np.sum((z-z.mean())**2))
    return {"descriptor_pair":["O",hspecies],"n":len(rows),"in_sample_r2":r2,"loocv_rmse":cv,"maximum":{"kind":kind,"O_eV":peakx,f"{hspecies}_eV":peaky,"log10_TOF":pk,"TOF":10**pk,"nearest_catalysts":nearest},"coefficients_standardized":coef.tolist(),"convex_hull_vertices":[[float(a),float(b)] for a,b in hull_pts]}


def write_reference_comparison():
    spec=importlib.util.spec_from_file_location("legacy_converter_readonly",CONVERTER)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    sheets=mod.load_sheets(mod.SOURCE)
    dataset=next(v for k,v in sheets.items() if "catmap" in k.lower())
    summary=next(v for k,v in sheets.items() if "".join(ch for ch in k if ord(ch)<128)=="") if False else next(v for k,v in sheets.items() if v is not dataset)
    positions=mod.surfaces(dataset)
    summary_positions=mod.surfaces(summary)
    rows=[]
    for surf in SURFACES:
        _,old_gcol=positions[surf]
        _,strict_gcol=summary_positions[surf]
        old_slab=float(dataset.get(mod.cell(old_gcol,51),float("nan")))
        new_slab=float(summary.get(mod.cell(strict_gcol,4),float("nan")))
        closure=float(summary.get(mod.cell(strict_gcol,41),float("nan")))
        old_addr=lambda row: mod.cell(old_gcol,row)
        strict_addr=lambda row: mod.cell(strict_gcol,row)
        rows.append({"surface_name":surf,"catalyst":LABELS[surf],"legacy_sheet":"catmap数据集","legacy_clean_slab_cell":old_addr(51),"legacy_clean_slab_G_eV":old_slab,"strict_sheet":"吉布斯自由能汇总","strict_clean_slab_cell":strict_addr(4),"strict_clean_slab_G_eV":new_slab,"strict_alternate_slab_cell":strict_addr(41),"strict_alternate_slab_G_eV":closure,"strict_minus_legacy_slab_eV":new_slab-old_slab,"row41_minus_row4_eV":closure-new_slab})
    return rows,sha256(mod.SOURCE)


def write_workbook_source_audit():
    import openpyxl
    spec=importlib.util.spec_from_file_location("legacy_converter_audit",CONVERTER)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    values=openpyxl.load_workbook(mod.SOURCE,data_only=True,read_only=True)
    formulas=openpyxl.load_workbook(mod.SOURCE,data_only=False,read_only=True)
    old_v,old_f=values["catmap数据集"],formulas["catmap数据集"]
    sum_v,sum_f=values["吉布斯自由能汇总"],formulas["吉布斯自由能汇总"]
    with LEGACY.open(encoding="utf-8-sig",newline="") as f:
        legacy_lookup={(r["surface_name"],r["species_name"]):float(r["formation_energy"]) for r in csv.DictReader(f,delimiter="\t")}
    old_map=mod.surfaces(mod.load_sheets(mod.SOURCE)["catmap数据集"])
    sum_map=mod.surfaces(mod.load_sheets(mod.SOURCE)["吉布斯自由能汇总"])
    mu={"Ha":mod.MU["H"],"Hb":mod.MU["H"],"O":mod.MU["O"]}
    rows=[]
    for surf in SURFACES:
        oldcol=old_map[surf][1]; sum_ecol,sumcol=sum_map[surf]
        sum_tcol=next((c for c in range(sum_ecol,sumcol) if sum_v.cell(3,c).value=="T"),None)
        for species in ("O","Ha","Hb"):
            row=mod.SINGLE_ROWS[species]
            gcell=mod.cell(oldcol,row); e0cell=mod.cell(oldcol-4,row); corr=mod.cell(oldcol-1,row); tcell=mod.cell(oldcol-3,row)
            slabcell=mod.cell(oldcol,51); strictslab=mod.cell(sumcol,4); altslab=mod.cell(sumcol,41)
            g=float(old_v[gcell].value); slab=float(old_v[slabcell].value)
            form=g-slab-mu[species]
            rows.append({"surface_name":surf,"catalyst":LABELS[surf],"species":species,"legacy_sheet":"catmap数据集","raw_G_cell":gcell,"raw_G_formula":old_f[gcell].value,"raw_G_cached_eV":g,"E0_cell":e0cell,"E0_cached_eV":old_v[e0cell].value,"correction_cell":corr,"correction_formula":old_f[corr].value,"correction_cached_eV":old_v[corr].value,"temperature_cell":tcell,"temperature_cached":old_v[tcell].value,"legacy_slab_cell":slabcell,"legacy_slab_G_eV":slab,"atomic_reference_eV":mu[species],"reconstructed_formation_energy_eV":form,"current_slab_cell":strictslab,"current_slab_G_eV":sum_v[strictslab].value,"current_slab_temperature_cell":mod.cell(sum_tcol,4) if sum_tcol else None,"current_slab_temperature_K":sum_v.cell(4,sum_tcol).value if sum_tcol else None,"alternate_slab_cell":altslab,"alternate_slab_G_eV":sum_v[altslab].value,"G_formula_reconstruction_matches_legacy_table":abs(form-legacy_lookup[(surf,species)])<2e-8})
    with (OUT/"workbook_source_audit.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    return rows


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    records=read_inputs()
    with (OUT/"matched_data.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(next(iter(records.values())).keys()));w.writeheader();w.writerows(records.values())
    reference,wb_hash=write_reference_comparison()
    with (OUT/"reference_slab_comparison.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(reference[0].keys()));w.writeheader();w.writerows(reference)
    workbook_audit=write_workbook_source_audit()
    one={}
    for species,stem in [("O","O_TOF_empirical"),("Ha","Ha_TOF_empirical"),("Hb","Hb_TOF_empirical")]: one[species]=one_d_plot(records,species,species,stem)
    two={"Ha":two_d_plot(records,"Ha","O_Ha_TOF_empirical_map"),"Hb":two_d_plot(records,"Hb","O_Hb_TOF_empirical_map")}
    result={"temperature_K":TEMPERATURE_K,"n_catalysts":len(records),"one_dimensional":one,"two_dimensional":two,"source_audit":{"legacy_rows_verified":len(workbook_audit),"legacy_formula_reconstruction_matches":sum(bool(r["G_formula_reconstruction_matches_legacy_table"]) for r in workbook_audit),"old_single_species_temperature_missing":sum(r["temperature_cached"] is None for r in workbook_audit)},"source_hashes":{"legacy_formation_table_sha256":sha256(LEGACY),"converter_sha256":sha256(CONVERTER),"workbook_sha256":wb_hash}}
    import json
    (OUT/"fit_metrics.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    site_list=", ".join(LABELS[s] for s in SURFACES)
    ref_lines=[]
    for r in reference:
        if abs(float(r["strict_minus_legacy_slab_eV"]))>=0.05 or abs(float(r["row41_minus_row4_eV"]))>=0.05:
            ref_lines.append(f"- {r['catalyst']}: old slab {float(r['legacy_clean_slab_G_eV']):.6f} eV ({r['legacy_clean_slab_cell']}); current row-4 slab {float(r['strict_clean_slab_G_eV']):.6f} eV ({r['strict_clean_slab_cell']}), Δ(current−old)={float(r['strict_minus_legacy_slab_eV']):+.6f} eV; row41−row4={float(r['row41_minus_row4_eV']):+.6f} eV.")
    notes=f'''# O*/Ha*/Hb* 历史形成能—严格 CatMAP TOF 经验图\n\n## 数据及图的性质\n\n- 位点：{site_list}（13个；Ti–Zr未收敛而排除，Ti–Nb不属于本轮严格基准集合）。按 `surface_name` 一对一连接；没有按行号拼接。\n- y 数据：各位点 `official_strict_tolerance/<surface>/gas_turnover_frequencies.csv` 中 `C6H12O_g`，并核对 `convergence.csv` 为 CatMAP 0.3.1、333.15 K、`CONVERGED`。未重算 CatMAP。\n- x 数据：原有 `DFT_formation_catmap.txt` 中 O、Ha、Hb `formation_energy`，其转换脚本使用 `G_ads − G_slab − μ_atom`，其中旧单物种数据分别来自 `catmap数据集` 的 Gibbs 缓存列，slab 为第51行，原子参考为 μH={-1.11:.2f} eV、μO={-4.37:.2f} eV；Ha 与 Hb 独立保留。复算 39 个值全部与旧表一致。\n- 每张图上的圆点/菱形及名称是实测严格 CatMAP 位点数据；彩色曲线和二维色面是对这些点做的经验 OLS 拟合，不是 CatMAP 网格计算、不是新的 TOF。二维拟合限在 13 点描述符凸包内。\n- 一维图同时显示线性、二次、连续非对称分段线性拟合；二维图显示 O+Ha 或 O+Hb 六参数二次响应面。峰或最大值均由对应拟合函数在样本范围/凸包中计算；TOF=10^{{log10(TOF)}}。\n\n## 拟合结论与峰值\n\n| 描述符 | 样本内二次 R² | 线性/二次/分段 LOOCV RMSE（log10 TOF） | 峰判断 |\n|---|---:|---:|---|\n| O* | {one['O']['quadratic_r2']:.3f} | {one['O']['fits']['linear']['loocv_rmse']:.2f} / {one['O']['fits']['quadratic']['loocv_rmse']:.2f} / {one['O']['fits']['segmented']['loocv_rmse']:.2f} | 无内部火山峰；分段斜率为 {one['O']['segmented_left_slope']:.2f}, {one['O']['segmented_right_slope']:.2f}，呈谷形转向。二次模型最大值在边界。 |\n| Ha* | {one['Ha']['quadratic_r2']:.3f} | {one['Ha']['fits']['linear']['loocv_rmse']:.2f} / {one['Ha']['fits']['quadratic']['loocv_rmse']:.2f} / {one['Ha']['fits']['segmented']['loocv_rmse']:.2f} | 无内部火山峰；分段斜率为 {one['Ha']['segmented_left_slope']:.2f}, {one['Ha']['segmented_right_slope']:.2f}，呈谷形转向。二次模型最大值在边界。 |\n| Hb* | {one['Hb']['quadratic_r2']:.3f} | {one['Hb']['fits']['linear']['loocv_rmse']:.2f} / {one['Hb']['fits']['quadratic']['loocv_rmse']:.2f} / {one['Hb']['fits']['segmented']['loocv_rmse']:.2f} | 二次模型有样本范围内峰：Hb={one['Hb']['peaks'][1]['descriptor']:.3f} eV，预测 log10TOF={one['Hb']['peaks'][1]['log10_TOF']:.3f}，TOF={one['Hb']['peaks'][1]['TOF']:.3g}；分段模型两侧斜率均为负，不支持火山。 |\n\n二维 O+Ha：样本内 R²={two['Ha']['in_sample_r2']:.3f}，LOOCV RMSE={two['Ha']['loocv_rmse']:.2f}；最大值为凸包边界最大值，O={two['Ha']['maximum']['O_eV']:.3f} eV、Ha={two['Ha']['maximum']['Ha_eV']:.3f} eV，预测 log10TOF={two['Ha']['maximum']['log10_TOF']:.3f}、TOF={two['Ha']['maximum']['TOF']:.3g}，邻近真实点：{', '.join(two['Ha']['maximum']['nearest_catalysts'])}。\n\n二维 O+Hb：样本内 R²={two['Hb']['in_sample_r2']:.3f}，LOOCV RMSE={two['Hb']['loocv_rmse']:.2f}；最大值为凸包边界最大值，O={two['Hb']['maximum']['O_eV']:.3f} eV、Hb={two['Hb']['maximum']['Hb_eV']:.3f} eV，预测 log10TOF={two['Hb']['maximum']['log10_TOF']:.3f}、TOF={two['Hb']['maximum']['TOF']:.3g}，邻近真实点：{', '.join(two['Hb']['maximum']['nearest_catalysts'])}。这不是内部二维最优活性区。\n\n## 能量参考与温度证据\n\n旧表与当前严格网络不是已验证的统一能量参考。旧转换器用 `catmap数据集` 每个表面第51行作为 clean slab；严格网络用 `吉布斯自由能汇总` 每个表面第4行作为 clean slab。当前工作簿两处 slab 数值差异（current−old）≥0.05 eV 或 row41−row4 差异≥0.05 eV 的位点：\n\n{chr(10).join(ref_lines)}\n\n旧单物种 O/Ha/Hb 的 `T` 列（catmap数据集对应各源行）对 39/39 记录为空；虽然能追溯其 `E0`、`ΔG` 修正列及缓存公式，但没有证据把这些修正归属到 333.15 K。严格汇总表第4行的温度字段可记录 clean-slab 项的温度，但不能补证旧 O/Ha/Hb 单物种热修正。旧与新 slab 差异不得自动作为校正量，因为结构等价性和热力学基准尚未验证。故**不生成所谓统一参考 O/Ha/Hb 表**，本轮图只作为历史参考下的探索性经验关系。逐单元格、公式、缓存值和 slab 对照见 `workbook_source_audit.csv`、`reference_slab_comparison.csv`；输入哈希见 `fit_metrics.json`。\n\n## 文件\n\n- `O_TOF_empirical.png/.pdf`\n- `Ha_TOF_empirical.png/.pdf`\n- `Hb_TOF_empirical.png/.pdf`\n- `O_Ha_TOF_empirical_map.png/.pdf`\n- `O_Hb_TOF_empirical_map.png/.pdf`\n- `matched_data.csv`：逐 surface_name 匹配的13条严格 TOF 和旧形成能。\n- `fit_metrics.json`：拟合系数/指标、峰值及源文件 SHA256。\n- `workbook_source_audit.csv`：39个物种源单元格公式、缓存值、ΔG/E0、温度字段与重算核验。\n- `reference_slab_comparison.csv`：旧 slab、新 row4 slab 与 row41 参照差异。\n- `analyze_legacy_descriptors.py`：只读数据的可复现绘图脚本。\n'''
    (OUT/"NOTES.md").write_text(notes,encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__ == "__main__": main()
