"""Plot archived Rpair* coverage--TOF regression and the official product DRC heatmap.
No new CatMAP calculation; regenerate PNG/PDF/SVG from joined data.
"""
from pathlib import Path
import csv
import math
import numpy as np
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parents[1]
SOURCE=BASE/"data"/"rpair_13sites_joined.csv"
OUT=BASE/"figures"
SITES={"ti":"#17529b","tico":"#df9e13","tini":"#0a977e","tire":"#ae75ab","tife":"#ce6538"}
def main():
    with SOURCE.open(encoding="utf-8-sig",newline="") as f:
        data=list(csv.DictReader(f))
    assert len(data)==13 and len(set(r["surface"] for r in data))==13
    OUT.mkdir(exist_ok=True,parents=True)
    xx=np.array([float(r["log_theta_Rpair"]) for r in data])
    yy=np.array([float(r["log_TOF"]) for r in data])
    coeff=np.polyfit(xx,yy,1)
    yp=np.polyval(coeff,xx)
    r2=1-np.sum((yy-yp)**2)/np.sum((yy-yy.mean())**2)
    plt.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"svg.fonttype":"none"})
    fig,ax=plt.subplots(figsize=(9.5,5.3),layout="constrained")
    xxgrid=np.linspace(xx.min(),xx.max(),200)
    ax.plot(xxgrid,np.polyval(coeff,xxgrid),"--",color="#516580",label="OLS (descriptive only)")
    for i,row in enumerate(data):
        ax.scatter(xx[i],yy[i],s=60,color=SITES.get(row["surface"],"#8b93a4"),edgecolor="white",zorder=3)
        ax.annotate(row["catalyst"],(xx[i],yy[i]),xytext=(5,7),textcoords="offset points",fontsize=8)
    ax.text(.02,.98,f"R²={r2:.4f}   slope={coeff[0]:.3f}",transform=ax.transAxes,va="top",bbox=dict(facecolor="white",edgecolor="#dde4ee"))
    ax.set(xlabel="log10 θ(Rpair*)",ylabel="log10 TOF(epoxide, s−1)",title="Rpair* coverage and product TOF (13 sites)")
    ax.grid(alpha=.17);ax.legend(frameon=False)
    for fmt in ("png","pdf","svg"):
        fig.savefig(OUT/f"01_Rpair_coverage_vs_TOF.{fmt}",dpi=300)
    plt.close(fig)
    keys=["DRC_Rpair","DRC_TS1","DRC_Epoxide","DRC_HHO","DRC_OHHa"]
    z=np.array([[float(r[k]) for k in keys] for r in data])
    fig,ax=plt.subplots(figsize=(9.6,5.6),layout="constrained")
    im=ax.imshow(z,vmin=-1.05,vmax=1.05,cmap="coolwarm",aspect="auto")
    ax.set_yticks(np.arange(len(data)),[r["catalyst"] for r in data])
    ax.set_xticks(np.arange(5),["Rpair*","TS1*","EpoxideHHO*","HHO*","OHHa*"])
    for i in range(len(data)):
        for j in range(5):
            ax.text(j,i,f"{z[i,j]:.3f}" if abs(z[i,j])>1e-4 else "~0",ha="center",va="center",fontsize=7.4)
    fig.colorbar(im,ax=ax,label="Product TOF degree of rate control")
    ax.set_title("Official CatMAP DRC (Rpair* near zero)")
    for fmt in ("png","pdf","svg"):
        fig.savefig(OUT/f"02_Rpair_DRC_heatmap.{fmt}",dpi=300)
    plt.close(fig)
    print("Rpair coverage--TOF fit R2",r2)

if __name__=="__main__":
    main()
