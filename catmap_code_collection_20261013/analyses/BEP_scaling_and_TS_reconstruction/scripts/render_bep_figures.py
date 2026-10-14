from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
D,O,F,R=ROOT/"data",ROOT/"origin_exports",ROOT/"figures",ROOT/"reports"
for p in (O,F,R): p.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "figure.dpi":150,"savefig.dpi":300,"axes.grid":True,"grid.alpha":0.2})
B1=pd.read_csv(O/"BEP1_13sites.csv")
B2=pd.read_csv(O/"BEP2_13sites.csv")
MET=pd.read_csv(O/"official_BEP_comparison.csv")
EXT=pd.read_csv(O/"BEP_extended_model_ranking.csv")
RECON=pd.read_csv(O/"BEP_TS1_TS2_reconstruction.csv")
CHAIN=pd.read_csv(D/"source_prior_chain_predictions.csv")
GRID=pd.read_csv(D/"source_catmap_VectorMap_400x400.csv")

def save(fig,stem):
    fig.savefig(F/(stem+".png"),bbox_inches="tight")
    fig.savefig(F/(stem+".pdf"),bbox_inches="tight")
    plt.close(fig)

def label_points(ax,x,y,names,offset=5):
    ax.scatter(x,y,s=36,marker="o",facecolor="#0072B2",edgecolor="white",linewidth=.6,zorder=3,label="13 validated catalyst data")
    manual={"Ti-Ti":(-22,17),"Ti-Hf":(-34,-9),"Ti":(-14,12),"Ti-Fe":(-27,12),"Ti-V":(10,-16),
            "Ti-Re":(8,-13),"Ti-Ni":(8,8),"Ti-Mn":(8,9),"Ti-Co":(8,8),"Ti-Mo":(9,-12),
            "Ti-Cr":(8,9),"Ti-W":(8,8),"Ti-Ta":(8,9)}
    for xx,yy,name in zip(x,y,names):
        ax.annotate(name,(xx,yy),xytext=manual.get(name,(4,offset)),textcoords="offset points",fontsize=7,
                    bbox={"boxstyle":"round,pad=0.08","fc":"white","ec":"none","alpha":0.75})

for data, num, title in [(B1,"NEB1","NEB1: Rpair → TS1 → EpoxideHHO"),(B2,"NEB2","NEB2: HHO → TS2 → OHHa")]:
    fig,ax=plt.subplots(figsize=(7.6,5.1))
    x=data.reaction_energy_eV.to_numpy(); y=data.signed_TS_minus_IS_eV.to_numpy()
    label_points(ax,x,y,data.catalyst)
    grid=np.linspace(x.min()-.08*(x.max()-x.min()),x.max()+.08*(x.max()-x.min()),250)
    row=MET.loc[MET.path.eq("NEB1_path2" if num=="NEB1" else "NEB2_path4")].iloc[0]
    ax.plot(grid,row.numpy_alpha*grid+row.numpy_beta,color="#D55E00",lw=2,label=f"OLS BEP: α={row.numpy_alpha:.4f}, β={row.numpy_beta:.4f} eV")
    neg=data.loc[data.signed_TS_minus_IS_eV<0]
    if len(neg): ax.scatter(neg.reaction_energy_eV,neg.signed_TS_minus_IS_eV,marker="X",s=64,color="#CC3311",zorder=4,label="negative signed TS−IS retained")
    ax.axhline(0,color="black",lw=.7)
    ax.set(xlabel="Reaction energy ΔG_rxn = G_FS − G_IS (eV)",ylabel="Signed TS − IS energy difference (eV)",title=title+" — signed values retained")
    ax.legend(frameon=False,loc="best"); save(fig,f"{num}_path_BEP")

fig,ax=plt.subplots(figsize=(8.2,5.2))
for dat,label,color,marker in [(B1,"NEB1 / path 2","#0072B2","o"),(B2,"NEB2 / path 4","#D55E00","s")]:
    ax.scatter(dat.catalyst,dat.BEP_residual_eV,label=label,color=color,marker=marker,s=42)
ax.axhline(0,color="black",lw=.8); ax.set_ylabel("In-sample BEP residual (eV)"); ax.set_xlabel("Catalyst")
ax.set_title("Signed BEP residuals by catalyst"); ax.tick_params(axis="x",rotation=45); ax.legend(frameon=False)
save(fig,"BEP_residuals_by_catalyst")

fig,axs=plt.subplots(1,2,figsize=(12,5.2))
for ax,metric,title in [(axs[0],"R2","In-sample R²"),(axs[1],"RMSE_eV","In-sample RMSE (eV)")]:
    for path,color,mark in [("NEB1_path2","#0072B2","o"),("NEB2_path4","#D55E00","s")]:
        sub=EXT.loc[EXT.path.eq(path)].sort_values(metric,ascending=(metric!="R2")).head(5)
        base=MET.loc[MET.path.eq(path)].iloc[0]
        baseval=base.R2 if metric=="R2" else base.RMSE_eV
        ax.scatter([f"{path}: ΔG only"],[baseval],color=color,marker="X",s=58,label=f"{path}: baseline")
        ax.scatter([f"{path}: +{sub.iloc[0].added_descriptor}"],[sub.iloc[0][metric]],color=color,marker=mark,s=48,label=f"{path}: best added term")
        # All tested one-covariate corrections are retained in the CSV; show rank labels compactly.
        for rank,(_,v) in enumerate(sub.iterrows(),1):
            ax.scatter([f"{path}: +{v.added_descriptor}"],[v[metric]],color=color,marker=mark,s=30,alpha=.55)
    ax.set_title(title); ax.tick_params(axis="x",rotation=75,labelsize=7); ax.legend(frameon=False,fontsize=7)
fig.suptitle("One-extra-descriptor BEP diagnostics (exploratory; historic references flagged)")
fig.tight_layout(); save(fig,"BEP_extended_model_comparison")

fig,axs=plt.subplots(1,2,figsize=(11,5.2))
for ax,path,actual,col in [(axs[0],"NEB1_path2","TS1","#0072B2"),(axs[1],"NEB2_path4","TS2","#D55E00")]:
    sub=RECON.loc[RECON.path.eq(path)]
    label_points(ax,sub.actual_TS_eV,sub.BEP_in_sample_TS_reconstruction_eV,sub.catalyst)
    lo=min(sub.actual_TS_eV.min(),sub.BEP_in_sample_TS_reconstruction_eV.min()); hi=max(sub.actual_TS_eV.max(),sub.BEP_in_sample_TS_reconstruction_eV.max())
    ax.plot([lo,hi],[lo,hi],"k--",lw=1,label="identity")
    ax.set_xlabel(f"Actual {actual} formation energy (eV)"); ax.set_ylabel("BEP in-sample reconstructed TS energy (eV)"); ax.set_title(path); ax.legend(frameon=False,fontsize=7)
save(fig,"BEP_TS_reconstruction_parity")

# Display actual catalyst points and the prior O/H descriptor-chain predictions in separate panels.
maps=pd.read_csv(O/"BEP_TOF_mapping_summary.csv")
fig,axs=plt.subplots(1,2,figsize=(13.5,6.0),sharex=True,sharey=True)
xcol=next(c for c in GRID.columns if "TS1" in c and "eV" in c); ycol=next(c for c in GRID.columns if "TS2" in c and "eV" in c); zcol=next(c for c in GRID.columns if "log10_TOF" in c)
piv=GRID.pivot(index=ycol,columns=xcol,values=zcol).sort_index()
levels=np.linspace(np.nanpercentile(piv.values,5),np.nanpercentile(piv.values,95),18)
cf=None
points=pd.read_csv(D/"source_validated_catalyst_points.csv")
px=next(c for c in points.columns if "TS1" in c and "eV" in c); py=next(c for c in points.columns if "TS2" in c and "eV" in c)
for ax in axs:
    cf=ax.contourf(piv.columns,piv.index,piv.values,levels=levels,cmap="viridis",extend="both")
inside=maps.loc[maps.existing_grid_bounds_supported]
axs[0].scatter(points[px],points[py],s=42,facecolor="white",edgecolor="black",label="Validated catalyst coordinates",zorder=4)
for _,v in points.iterrows():
    axs[0].annotate(v.catalyst,(v[px],v[py]),xytext=(4,4),textcoords="offset points",fontsize=7,color="black",
                    bbox={"boxstyle":"round,pad=0.12","fc":"white","ec":"none","alpha":0.72})
axs[1].scatter(inside.TS1_LOOCV_prediction_eV,inside.TS2_LOOCV_prediction_eV,marker="D",s=34,facecolor="#E69F00",edgecolor="black",label="O/H chain LOOCV-predicted coordinates",zorder=5)
offsets={"Ti-Ta":(-12,12),"Ti-Hf":(9,13),"Ti-W":(9,-15),"Ti-Re":(8,13),"Ti-V":(8,-14),
         "Ti-Cr":(-23,10),"Ti-Mo":(9,7),"Ti-Ti":(9,-13),"Ti-Co":(9,12),"Ti-Fe":(9,-8),
         "Ti-Mn":(-19,-12),"Ti-Ni":(8,10),"Ti":(9,4)}
for _,v in inside.iterrows():
    axs[1].annotate(v.catalyst,(v.TS1_LOOCV_prediction_eV,v.TS2_LOOCV_prediction_eV),xytext=offsets.get(v.catalyst,(4,4)),textcoords="offset points",fontsize=7,color="black")
axs[0].set_title("Actual validated catalyst coordinates")
axs[1].set_title("Descriptor-chain LOOCV predictions")
for ax in axs:
    ax.set_xlabel("TS1_s coordinate (eV)")
    dx=(piv.columns.max()-piv.columns.min())*.025; dy=(piv.index.max()-piv.index.min())*.025
    ax.set_xlim(piv.columns.min()-dx,piv.columns.max()+dx)
    ax.set_ylim(piv.index.min()-dy,piv.index.max()+dy)
    ax.legend(frameon=True,fontsize=7,loc="best")
axs[0].set_ylabel("TS2_s coordinate (eV)")
fig.colorbar(cf,ax=axs.ravel().tolist(),label="Saved CatMAP VectorMap-interpolated log10(TOF C6H12O)",shrink=.84)
fig.suptitle("Predicted TS coordinates mapped to saved CatMAP VectorMap\nBackground: prior interpolation only; no new CatMAP solves")
save(fig,"BEP_chain_coordinates_on_existing_CatMAP_surface")

# Compact Origin-ready combined tables and hash manifest for inputs/outputs.
pd.concat([B1.assign(path="NEB1_path2"),B2.assign(path="NEB2_path4")]).to_csv(O/"BEP_all_signed_points_origin.csv",index=False)
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()
rows=[]
for p in sorted(ROOT.rglob("*")):
    if p.is_file() and p.name!="SHA256SUMS.csv": rows.append({"relative_path":str(p.relative_to(ROOT)),"sha256":sha(p),"bytes":p.stat().st_size})
pd.DataFrame(rows).to_csv(R/"SHA256SUMS.csv",index=False)
print(f"Rendered {len(list(F.glob('*.png')))} PNG and {len(list(F.glob('*.pdf')))} PDF files")
