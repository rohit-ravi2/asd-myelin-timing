#!/usr/bin/env python3
"""PART 3 (partial): dimensional test. HBN-POD2 has dki_awf imaging for 1,879
children but NO diagnosis (that is what the CMI DUA gates). The EEG subcohort
however carries dimensional psychopathology factor scores, and 179 subjects have
both. This tests the FORM of the spectrum claim -- does myelin-sensitive white
matter scale with symptom load -- without being ASD-specific."""
import json
import numpy as np, pandas as pd
from scipy import stats
from pathlib import Path
H=Path("data/asd-corpus/imaging/hbn_pod2"); P=Path(".")
eeg=pd.read_csv(H/"hbn_eeg_participants.tsv",sep="\t")
eeg["pid"]=eeg.participant_id.astype(str).str.replace("^sub-","",regex=True)
keep=set(eeg.pid)
print("loading tract profiles (543 MB)...", flush=True)
tp=pd.read_csv(H/"combined_tract_profiles.csv",
               usecols=["tractID","nodeID","dki_fa","dki_md","dki_mk","dki_awf","subjectID"])
tp["subjectID"]=tp.subjectID.astype(str)
tp=tp[tp.subjectID.isin(keep)]
print(f"rows for EEG-cohort subjects: {len(tp):,}; subjects: {tp.subjectID.nunique()}")
# drop the 10 end-nodes each side (endpoint artefacts), then mean per subject x tract
tp=tp[(tp.nodeID>=10)&(tp.nodeID<=89)]
agg=tp.groupby(["subjectID","tractID"])[["dki_fa","dki_md","dki_mk","dki_awf"]].mean().reset_index()
wide=agg.pivot(index="subjectID",columns="tractID",values="dki_awf")
glob=agg.groupby("subjectID")[["dki_fa","dki_md","dki_mk","dki_awf"]].mean()
ph=eeg.set_index(eeg.pid)[["age","sex","p_factor","attention","internalizing","externalizing"]]
D=glob.join(ph,how="inner").dropna(subset=["p_factor","age"])
print(f"\nanalysis sample n={len(D)}  age {D.age.min():.1f}-{D.age.max():.1f}")
print(f"\n{'metric':<10}{'factor':<16}{'partial r (age,sex adj)':>24}{'p':>10}")
res=[]
for met in ["dki_awf","dki_fa","dki_mk","dki_md"]:
    for fac in ["p_factor","attention","internalizing","externalizing"]:
        sub=D.dropna(subset=[met,fac])
        X=np.column_stack([np.ones(len(sub)),sub.age.values,(sub.sex.astype(str)=="M").astype(float)])
        ry=sub[met].values-X@np.linalg.lstsq(X,sub[met].values,rcond=None)[0]
        rx=sub[fac].values-X@np.linalg.lstsq(X,sub[fac].values,rcond=None)[0]
        r,p=stats.pearsonr(rx,ry)
        star="**" if p<0.01 else "*" if p<0.05 else ""
        print(f"{met:<10}{fac:<16}{r:>+24.3f}{p:>10.4f} {star}")
        res.append({"metric":met,"factor":fac,"n":len(sub),"partial_r":round(float(r),3),"p":round(float(p),4)})
# tract-level, dki_awf vs p_factor
print(f"\nTract-level: dki_awf vs p_factor (age/sex adjusted), top 8 of {wide.shape[1]} tracts")
W=wide.join(ph,how="inner").dropna(subset=["p_factor","age"])
rows=[]
for t in wide.columns:
    sub=W.dropna(subset=[t])
    if len(sub)<60: continue
    X=np.column_stack([np.ones(len(sub)),sub.age.values,(sub.sex.astype(str)=="M").astype(float)])
    ry=sub[t].values-X@np.linalg.lstsq(X,sub[t].values,rcond=None)[0]
    rx=sub.p_factor.values-X@np.linalg.lstsq(X,sub.p_factor.values,rcond=None)[0]
    r,p=stats.pearsonr(rx,ry); rows.append((p,t,r,len(sub)))
rows.sort()
m=len(rows)
for i,(p,t,r,nn) in enumerate(rows[:8]):
    q=min(1,p*m/(i+1))
    print(f"  {t:<22} n={nn:>3}  r={r:+.3f}  p={p:.4f}  q={q:.3f}")
print(f"  ({m} tracts tested; min q = {min(min(1,p*m/(i+1)) for i,(p,_,_,_) in enumerate(rows)):.3f})")
json.dump(res,open(P/"outputs/apoptosis/hbn_dimensional.json","w"),indent=1)
