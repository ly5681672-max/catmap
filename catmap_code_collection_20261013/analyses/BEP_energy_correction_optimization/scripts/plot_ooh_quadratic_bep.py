"""Reproduce 4-coefficient quadratic OOH*-corrected NEB1 BEP figure.
Archival post-processing ONLY: no DFT, NEB, CatMAP steady state or TOF runs.
Run from anywhere: python plot_ooh_quadratic_bep.py
Requires numpy, pandas, matplotlib (all present in typical Anaconda installs).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data"/"source_carbon_O_H_TS1.csv"
OUT=ROOT/"figures"
X="NEB1_reaction_energy_EpoxideHHO_minus_Rpair"
Z="OOH_star_legacy"
Y="NEB1_forward_energy_difference_TS1_minus_Rpair"

def main():
    df=pd.read_csv(SOURCE,encoding="utf-8-sig")
    assert len(df)==13 and df.surface_name.nunique()==13,"Expected 13 unique sites"
    x=df[X].to_numpy(float);z=df[Z].to_numpy(float);y=df[Y].to_numpy(float)
    xc=x-x.mean();zc=z-z.mean()
    A=np.column_stack([np.ones(len(y)),xc,zc,zc**2])
    coeff=np.linalg.lstsq(A,y,rcond=None)[0]
    fit=A@coeff
    loo=np.array([A[i]@np.linalg.lstsq(A[np.arange(len(y))!=i],y[np.arange(len(y))!=i],rcond=None)[0] for i in range(len(y))])
    basic=np.column_stack([np.ones(len(y)),xc])
    baseline=basic@np.linalg.lstsq(basic,y,rcond=None)[0]
    r2=1-np.sum((y-fit)**2)/np.sum((y-y.mean())**2)
    rmse=np.sqrt(np.mean((y-fit)**2))
    cv=np.sqrt(np.mean((y-loo)**2))
    print(f"R2={r2:.7f} RMSE={rmse:.6f} eV LOOCV RMSE={cv:.6f} eV")
    print("centered coefficients",coeff.tolist(),"mean DG2",x.mean(),"mean OOH",z.mean())
    OUT.mkdir(parents=True,exist_ok=True)
    tab=df[["surface_name","catalyst",X,Z,Y]].copy()
    tab["fitted_Y_eV"]=fit;tab["LOOCV_Y_eV"]=loo;tab["residual_eV"]=y-fit
    tab.to_csv(OUT/"08_OOH_quadratic_BEP_predictions.csv",index=False)
    plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,5.7),gridspec_kw={"width_ratios":[1,1.15]})
    lim=(-1.35,1.05)
    ax.plot(lim,lim,"--",lw=1.2,color="#8794a7")
    offsets={"Ti-Fe":(8,-12),"Ti-Mn":(8,14),"Ti-Hf":(-42,-24),"Ti-Re":(9,-12),"Ti-Mo":(-44,13),"Ti-V":(-42,-14),"Ti-Co":(8,7),"Ti-Ti":(-43,18),"Ti-W":(8,-13),"Ti-Ta":(8,7),"Ti-Cr":(8,-15),"Ti":(-21,9),"Ti-Ni":(8,10)}
    for i,r in df.iterrows():
        site=r["catalyst"];sp=site in {"Ti-W","Ti-Ta"}
        ax.scatter(y[i],fit[i],s=72 if sp else 47,marker="^" if sp else "o",c="#be4a50" if sp else "#24599b",edgecolor="white",linewidth=.7,zorder=3)
        ax.annotate(site,(y[i],fit[i]),xytext=offsets[site],textcoords="offset points",fontsize=8)
    ax.set(xlim=lim,ylim=lim,xlabel="Calculated signed TS1 − Rpair (eV)",ylabel="Fitted signed TS1 − Rpair (eV)",title="(a) Quadratic OOH*-corrected BEP")
    ax.set_aspect("equal",adjustable="box")
    ax.grid(alpha=.15)
    ax.text(.03,.97,f"R²={r2:.4f} | RMSE={rmse:.4f} eV\nLOOCV RMSE={cv:.4f} eV",transform=ax.transAxes,va="top",fontsize=9,bbox=dict(boxstyle="round",fc="white",ec="#d7e0ee"))
    pos=np.arange(len(y))
    bx.barh(pos-.18,y-baseline,height=.32,color="#a7aeb9",label="Traditional BEP")
    bx.barh(pos+.18,y-fit,height=.32,color="#245e9d",label="Quadratic OOH* BEP")
    bx.axvline(0,c="#667085",lw=1)
    bx.set(yticks=pos,yticklabels=df.catalyst,xlabel="Observed − fitted (eV)",title="(b) Residuals across all 13 sites")
    bx.invert_yaxis();bx.grid(axis="x",alpha=.16);bx.legend(frameon=False,fontsize=8,loc="upper right")
    fig.suptitle("Quadratic OOH* correction to NEB1 BEP",fontsize=15,weight="bold")
    fig.text(.5,.012, f"Y = {coeff[0]:.4f} {coeff[1]:+.4f} u {coeff[2]:+.4f} z {coeff[3]:+.4f} z²;  u = ΔG₂ + {-x.mean():.4f}, z = E(OOH*) + {-z.mean():.4f}",ha="center",fontsize=9)
    fig.subplots_adjust(left=.08,right=.97,top=.86,bottom=.16,wspace=.33)
    for ext in ("png","pdf","svg"):
        fig.savefig(OUT/f"08_OOH_quadratic_BEP_fit.{ext}",dpi=300,facecolor="white",bbox_inches="tight")
    plt.close(fig)

if __name__=="__main__":
    main()
