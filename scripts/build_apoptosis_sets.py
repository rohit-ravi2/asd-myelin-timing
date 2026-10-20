#!/usr/bin/env python3
"""Gene sets to test the apoptosis x oligodendrocyte convergence hypothesis.

Hypothesis: failed apoptosis -> hyperconnectivity; hypomyelination -> E/I
imbalance; and the genes for the two overlap. The decisive object is therefore
the INTERSECTION, tested against both parent sets and against size-matched
random sets drawn from the same universe.
"""
import json, random
from pathlib import Path
C = Path("data/asd-corpus"); P = Path(__file__).resolve().parent.parent

sym2e, universe = {}, set()
for line in open(C/"tools/magma_geneloc/NCBI37.3.gene.loc"):
    f = line.split()
    if len(f) >= 6:
        sym2e[f[5]] = f[0]; universe.add(f[0])
print(f"MAGMA universe: {len(universe)} genes")

gmt = {}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f = line.rstrip("\n").split("\t"); gmt[f[0]] = f[2:]
def go(term):
    return {sym2e[s] for s in gmt.get(term,[]) if s in sym2e}

# ---- existing sets from the prior pipeline (reuse verbatim) -----------------
prior = {}
for line in open(P/"annot/asd_genesets.magma.set"):
    f = line.split(); prior[f[0]] = set(f[1:])

APO = {
 "APOP_PROCESS":        go("GOBP_APOPTOTIC_PROCESS"),
 "APOP_NEURON":         go("GOBP_NEURON_APOPTOTIC_PROCESS"),
 "APOP_NEURON_NEGREG":  go("GOBP_NEGATIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"),
 "APOP_INTRINSIC":      go("GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY"),
 "APOP_IN_DEVELOPMENT": go("GOBP_APOPTOTIC_PROCESS_INVOLVED_IN_DEVELOPMENT"),
 "APOP_GLIAL":          go("GOBP_GLIAL_CELL_APOPTOTIC_PROCESS"),
}
PRUNE = {
 "MACROAUTOPHAGY":      go("GOBP_MACROAUTOPHAGY"),
 "TOR_SIGNALING":       go("GOBP_TOR_SIGNALING"),
 "SYNAPSE_PRUNING":     go("GOBP_SYNAPSE_PRUNING") | go("GOBP_REGULATION_OF_SYNAPSE_PRUNING"),
 "OLIGO_DIFFERENTIATION": go("GOBP_OLIGODENDROCYTE_DIFFERENTIATION"),
 "ENSHEATHMENT":        go("GOBP_ENSHEATHMENT_OF_NEURONS"),
}
OL_REF = prior["OLIGODENDROCYTE_VELM"]          # cell-type marker set, data-driven
MY_REF = prior["MYELIN"]

APOP_UNION = set().union(*APO.values())
sets = {}
sets.update({k: v for k, v in APO.items()})
sets.update({k: v for k, v in PRUNE.items()})
sets["APOP_UNION"] = APOP_UNION

# ---- the decisive objects: intersections and their complements -------------
sets["APOPxOLIGO"]      = APOP_UNION & OL_REF
sets["APOPxMYELIN"]     = APOP_UNION & MY_REF
sets["APOP_not_OLIGO"]  = APOP_UNION - OL_REF - MY_REF
sets["OLIGO_not_APOP"]  = OL_REF - APOP_UNION
sets["MYELIN_not_APOP"] = MY_REF - APOP_UNION

# ---- size-matched random controls for every intersection set ---------------
rng = random.Random(11)
pool = sorted(universe)
for tgt in ["APOPxOLIGO","APOPxMYELIN"]:
    n = len(sets[tgt])
    for i in range(5):
        sets[f"RANDMATCH_{tgt}_{i}"] = set(rng.sample(pool, n))
# random controls drawn from the apoptosis union (tests whether the OL slice is special)
for i in range(5):
    n = len(sets["APOPxOLIGO"])
    sets[f"RANDWITHIN_APOP_{i}"] = set(rng.sample(sorted(APOP_UNION), min(n, len(APOP_UNION))))

sets = {k: (v & universe) for k, v in sets.items()}
sets = {k: v for k, v in sets.items() if len(v) >= 5}

out = P/"annot/apoptosis_sets.magma.set"
with open(out,"w") as f:
    for k, v in sets.items():
        f.write(k + "\t" + "\t".join(sorted(v)) + "\n")
print(f"\nwrote {out} with {len(sets)} sets")
for k, v in sets.items():
    if not k.startswith(("RANDMATCH","RANDWITHIN")): print(f"  {k:<24} n={len(v)}")
print(f"\n  overlap check: |APOP_UNION|={len(APOP_UNION)}  |OLIGO|={len(OL_REF)}  "
      f"|intersection|={len(sets.get('APOPxOLIGO',[]))}")
exp = len(APOP_UNION)*len(OL_REF)/len(universe)
print(f"  expected overlap by chance = {exp:.1f}; observed = {len(sets.get('APOPxOLIGO',[]))}"
      f"  -> ratio {len(sets.get('APOPxOLIGO',[]))/exp:.2f}x")
json.dump({k: sorted(v) for k,v in sets.items()}, open(P/"outputs/apoptosis/sets.json","w"))
