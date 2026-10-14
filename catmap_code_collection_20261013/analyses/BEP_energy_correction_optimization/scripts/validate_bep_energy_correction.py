from pathlib import Path
import json
import pandas as pd

PKG=Path(__file__).resolve().parents[1]
DATA=PKG/"data"; FIG=PKG/"figures"; EXP=PKG/"origin_exports"
expected={
    ("13site_refit","traditional_BEP"):(.20542753321490248,.41769472237702737,.7134306722655978),
    ("13site_refit","model_A_OH_Ha_OOH"):(.6305746995970117,.2848102070492015,.7135941784373467),
    ("13site_refit","model_B_OH_Ha_OHHa"):(.6756588733318574,.2668659977374489,.6710370996105984),
    ("13site_refit","variable_slope_OH"):(.7563994313961024,.23127643312808743,.41967788966321284),
    ("11site_refit","model_A_OH_Ha_OOH"):(.726991148088036,.0803895330784755,.24524902087147954),
    ("11site_refit","model_B_OH_Ha_OHHa"):(.704215871469061,.08367554989520552,.1420619388540508),
}
m=pd.read_csv(DATA/"BEP_MODEL_COMPARISON.csv")
for key,vals in expected.items():
    row=m[(m.subset==key[0])&(m.model==key[1])].iloc[0]
    got=(row.R2,row.RMSE_eV,row.LOOCV_RMSE_eV)
    if any(abs(float(a)-b)>1e-9 for a,b in zip(got,vals)):
        raise AssertionError(f"metric reproduction failed {key}: {got} vs {vals}")
points=pd.read_csv(EXP/"BEP_source_points_13catalysts.csv")
assert len(points)==13 and points.surface_name.nunique()==13
assert not any("tizr" in str(x).lower() for x in points.surface_name)
comp=pd.read_csv(DATA/"BEP_13SITE_11SITE_COMPARISON.csv")
assert set(comp[comp.subset=="11site_refit"].n)=={11}
assert set(comp[comp.subset=="13site_refit"].n)=={13}
pred=pd.read_csv(DATA/"BEP_TS1_RECONSTRUCTION.csv")
assert len(pred[pred.subset=="13site_refit"])==52
assert abs((pred[pred.subset=="13site_refit"].predicted_TS1_formation_energy_eV-(pred[pred.subset=="13site_refit"].Rpair_strict_eV+pred[pred.subset=="13site_refit"].predicted_Y_eV)).abs().max())<1e-10
figs=list(FIG.glob("*.png")); pdfs=list(FIG.glob("*.pdf"))
assert len(figs)==7 and len(pdfs)==7,(len(figs),len(pdfs))
for f in figs+pdfs:
    assert f.stat().st_size>5000,f
for f in EXP.glob("*.csv"):
    d=pd.read_csv(f)
    assert len(d)>0,f
print(json.dumps({"status":"PASS","catalysts":13,"subset_11_excludes":["tiw","tita"],"figures_png":len(figs),"figures_pdf":len(pdfs),"origin_csv_count":len(list(EXP.glob('*.csv'))),"metrics_checked":len(expected)},indent=2))
