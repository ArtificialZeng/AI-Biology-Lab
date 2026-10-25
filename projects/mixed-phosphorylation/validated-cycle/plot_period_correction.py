"""Scientific illustration: floating finite-K data versus certified asymptotic coefficient."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
DATA=ROOT.parent/"evidence"
T=float(json.loads((DATA/"period_response.json").read_text())["T_limit"])
rows=[]
for name in ["baseline_K100.0_LSODA.json","baseline_K1000.0_LSODA.json","baseline_K10000.0_LSODA.json"]:
    d=json.loads((DATA/name).read_text());k=float(d["K"])
    rows.append((k,k*(d["period_mean"]-T)))
k,y=np.array(rows).T
fig,ax=plt.subplots(figsize=(8.4,5.0),layout="constrained")
ax.set_xscale("log")
ax.axhspan(7.65,7.69,color="#3174B5",alpha=.25,label="Interval-certified asymptotic coefficient: [7.65, 7.69]")
ax.plot(k,y,"o--",color="#686868",markersize=7,label="Finite-K simulation (floating point)")
for xx,yy in rows:ax.annotate(f"{yy:.4f}",(xx,yy),xytext=(5,8),textcoords="offset points",fontsize=9)
ax.set(xlabel="Total kinase K (model units)",ylabel=r"Scaled period excess $K[T(K)-T_\infty]$",ylim=(7.53,8.9),xlim=(70,16000))
ax.set_title("The period approaches its limit with a nonzero correction",loc="left",pad=13,fontweight="bold")
ax.grid(alpha=.15,which="both")
ax.spines[["top","right"]].set_visible(False)
ax.legend(loc="upper right",frameon=False,fontsize=9)
fig.supxlabel("Fixed published rates, P=8 and S=40. The blue band bounds T1, not finite-K errors.",fontsize=9)
output=ROOT/"figures";output.mkdir(exist_ok=True)
fig.savefig(output/"period_correction.png",dpi=180)
print(output/"period_correction.png")
