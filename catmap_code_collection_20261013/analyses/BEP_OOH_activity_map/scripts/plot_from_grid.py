"""Plot committed, independently solved Ti-Co-anchored BEP grid; no CatMAP run."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
GRID = BASE/"data"/"grid_BEP_OOH_dG2_TiCo_conditional.csv"
OUT = BASE/"figures"
ROOT = Path(__file__).resolve().parents[3]
SITES = ROOT/"analyses"/"BEP_energy_correction_optimization"/"data"/"source_carbon_O_H_TS1.csv"

def main():
    with GRID.open(encoding="utf-8-sig", newline="") as f: rows=list(csv.DictReader(f))
    xs=sorted({float(r["deltaG2_eV"]) for r in rows})
    zs=sorted({float(r["OOH_legacy_formation_eV"]) for r in rows})
    assert len(rows)==len(xs)*len(zs)
    vals=np.full((len(zs),len(xs)),np.nan);cut=np.ones_like(vals,dtype=bool)
    ix={v:i for i,v in enumerate(xs)};iz={v:i for i,v in enumerate(zs)}
    for r in rows:
        i=iz[float(r["OOH_legacy_formation_eV"])]
        j=ix[float(r["deltaG2_eV"])]
        vals[i,j]=float(r["log10_TOF_analytic"])
        cut[i,j]=bool(int(r["TS1_endpoint_clipped"]))
    fig,ax=plt.subplots(figsize=(9.2,5.2))
    im=ax.contourf(xs,zs,np.ma.masked_where(cut,vals),levels=np.arange(-38,5,2),cmap="viridis",extend="both")
    ax.contourf(xs,zs,cut.astype(float),levels=[.5,1.5],colors=["#eeeeee"],alpha=.94)
    ax.contour(xs,zs,cut.astype(float),levels=[.5],colors=["#c5604c"],linestyles=":",linewidths=.8)
    with SITES.open(encoding="utf-8-sig",newline="") as f:
        sites=list(csv.DictReader(f))
    for r in sites:
        x=float(r["NEB1_reaction_energy_EpoxideHHO_minus_Rpair"])
        z=float(r["OOH_star_legacy"])
        if r["surface_name"]=="tico":
            ax.scatter(x,z,s=160,marker="*",c="#d34548",edgecolor="white",zorder=5)
            ax.annotate("Ti-Co anchor",(x,z),xytext=(-78,-21),textcoords="offset points",color="#d34548",fontsize=8)
        else:
            ax.scatter(x,z,s=21,facecolor="none",edgecolor="#ced7df")
            ax.annotate(r["catalyst"],(x,z),xytext=(4,3),textcoords="offset points",fontsize=6.5,color="#eeeeee")
    c=fig.colorbar(im,ax=ax,pad=.02)
    c.set_label("Conditional log10 TOF (s-1)")
    ax.set(xlabel="Delta G2 = G(EpoxideHHO*) - G(Rpair*) (eV)",
           ylabel="Legacy E(OOH*) formation energy (eV)",
           title="Extended-BEP activity map (Ti-Co anchored; independent six-step solution)")
    fig.text(.5,.012,"Grey: TS1 below initial state; other energies fixed at Ti-Co. Not a new official CatMAP grid.",ha="center",fontsize=8)
    fig.subplots_adjust(left=.12,right=.89,bottom=.17,top=.9)
    OUT.mkdir(parents=True,exist_ok=True)
    for typ in ("png","pdf","svg"):
        fig.savefig(OUT/f"01_BEP_OOH_dG2_TiCo_conditional_activity.{typ}",dpi=300,facecolor="white")
    plt.close(fig)

if __name__=="__main__":main()
