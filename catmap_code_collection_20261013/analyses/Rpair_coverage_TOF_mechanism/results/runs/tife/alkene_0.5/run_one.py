import json
from catmap import ReactionModel, __version__
if __version__ != "0.3.1":
    raise RuntimeError("CatMAP 0.3.1 is required; found "+str(__version__))
model=ReactionModel(setup_file="alkene_epoxidation.mkm")
model.run(recalculate=True)
tof_map=model.turnover_frequency_map
cover_map=model.coverage_map
assert len(tof_map)==1 and len(cover_map)==1
gas=list(model.gas_names)
names=list(model.adsorbate_names)+[s for s in model.site_names if s!="g"]
p=tof_map[0][1]
c=cover_map[0][1]
assert "C6H12O_g" in gas and "Rpair_s" in names
tof=p[gas.index("C6H12O_g")]
theta=c[names.index("Rpair_s")]
assert tof>0 and theta>=0
res=model.solver.get_residual(c)
assert res<=model.solver.tolerance
with open("result.json","w",encoding="utf-8") as f:
    json.dump(dict(tof=str(tof),theta_Rpair=str(theta),residual=str(res),
                   version=__version__,surface=model.surface_names[0]),f,indent=2)
