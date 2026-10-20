#!/usr/bin/env python3
"""Which top ASD genes do NOT fall into the clusters we've established, and
could any of them route into the PV/oligodendrocyte pipeline?"""
import csv
from pathlib import Path
C=Path("data/asd-corpus"); S=Path("data/dvad-corpus/asd_scope/rare_variant")
gmt={}
for fn in ["c5.go.bp.v2024.1.Hs.symbols.gmt","c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    for line in open(C/"pathways"/fn):
        f=line.rstrip("\n").split("\t"); gmt[f[0]]=set(f[2:])
g=lambda t:gmt.get(t,set())
sat={r["hugoGene"] for r in csv.DictReader(open(S/"Satterstrom2020_102genes.csv")) if r.get("hugoGene")}
fu={r["HGNC"] for r in csv.DictReader(open(S/"Fu2022_72_studywide_genes.csv")) if r.get("HGNC")}
GOLD=sorted(sat|fu)
CL={
 "SYNAPTIC":   g("GOCC_SYNAPSE")|g("GOBP_SYNAPSE_ORGANIZATION")|g("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING"),
 "CHROMATIN":  g("GOBP_CHROMATIN_REMODELING")|g("GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION")|g("GOCC_HISTONE_DEACETYLASE_COMPLEX")|g("GOBP_HISTONE_MODIFICATION"),
 "OL_MYELIN":  g("GOBP_ENSHEATHMENT_OF_NEURONS")|g("GOBP_OLIGODENDROCYTE_DIFFERENTIATION")|g("GOCC_MYELIN_SHEATH"),
 "BCAT_TCF":   g("GOCC_BETA_CATENIN_TCF_COMPLEX"),
 "mTOR_RAS":   {"MTOR","PTEN","TSC1","TSC2","NF1","HRAS","AKT3","PIK3CA","RHEB","DEPDC5","BRAF","KRAS","SOS1","SYNGAP1","RASA1"},
 "ION_CHANNEL":g("GOCC_ION_CHANNEL_COMPLEX")|g("GOBP_REGULATION_OF_MEMBRANE_POTENTIAL"),
 "TRANSCRIPTION": g("GOBP_REGULATION_OF_DNA_TEMPLATED_TRANSCRIPTION"),
 "RNA_PROCESS": g("GOBP_MRNA_PROCESSING")|g("GOBP_RNA_SPLICING"),
 "CYTOSKELETON": g("GOBP_CYTOSKELETON_ORGANIZATION"),
}
assign={}
for gene in GOLD:
    hits=[k for k,v in CL.items() if gene in v]
    assign[gene]=hits
n_un=[gene for gene,h in assign.items() if not h]
core=[gene for gene,h in assign.items() if set(h)&{"SYNAPTIC","CHROMATIN"}]
print(f"GOLD standard genes: {len(GOLD)}")
for k in CL:
    n=sum(1 for gene in GOLD if k in assign[gene])
    print(f"  {k:<15} {n:>3}")
print(f"\nIn SYNAPTIC or CHROMATIN (the two dominant clusters): {len(core)}")
print(f"UNASSIGNED to any cluster: {len(n_un)}")
print(f"  {sorted(n_un)}")
# which unassigned genes touch the pipeline via less obvious routes
ROUTES={
 "OPC_marker_Velmeshev":None,
 "GLIAL_any": g("GOBP_GLIAL_CELL_DIFFERENTIATION")|g("GOBP_GLIOGENESIS")|g("GOBP_ASTROCYTE_DIFFERENTIATION"),
 "AXON_GUIDANCE_myelin_adjacent": g("GOBP_AXON_ENSHEATHMENT")|g("GOCC_NODE_OF_RANVIER")|g("GOCC_PARANODE_REGION_OF_AXON")|g("GOCC_AXON_INITIAL_SEGMENT"),
 "LIPID_metabolism": g("GOBP_LIPID_BIOSYNTHETIC_PROCESS")|g("GOBP_STEROL_BIOSYNTHETIC_PROCESS"),
 "WNT_broad": g("GOBP_CANONICAL_WNT_SIGNALING_PATHWAY")|g("GOBP_WNT_SIGNALING_PATHWAY"),
 "NOTCH": g("GOBP_NOTCH_SIGNALING_PATHWAY"),
 "CELL_CYCLE_proliferation": g("GOBP_REGULATION_OF_CELL_CYCLE")|g("GOBP_CELL_POPULATION_PROLIFERATION"),
 "GABAergic": g("GOBP_GAMMA_AMINOBUTYRIC_ACID_SIGNALING_PATHWAY")|g("GOBP_INHIBITORY_SYNAPSE_ASSEMBLY"),
}
import gzip,glob,os
V="data/dvad-corpus/asd_scope/velmeshev2019/markers"
cm={}
for f in glob.glob(V+"/*.tsv.gz"):
    nm=os.path.basename(f).replace(".tsv.gz","")
    s=set()
    with gzip.open(f,"rt") as fh:
        fh.readline()
        for line in fh:
            p=line.split("\t")
            if len(p)>1: s.add(p[1])
    cm[nm]=s
ROUTES["OPC_marker_Velmeshev"]=cm.get("OPC",set())
ROUTES["PV_marker_Velmeshev"]=cm.get("INPV",set())
print(f"\n--- UNASSIGNED genes: do they route into the pipeline by a less obvious path? ---")
routed=0
for gene in sorted(n_un):
    hits=[k for k,v in ROUTES.items() if v and gene in v]
    if hits:
        routed+=1; print(f"  {gene:<10} -> {', '.join(hits)}")
print(f"\n  {routed}/{len(n_un)} unassigned genes route in; {len(n_un)-routed} remain unexplained")
print(f"  truly unexplained: {sorted(gene for gene in n_un if not [k for k,v in ROUTES.items() if v and gene in v])}")
print(f"\n--- ALL gold genes that are OPC or PV markers (potential pipeline entry points) ---")
print(f"  OPC markers: {sorted(set(GOLD)&ROUTES['OPC_marker_Velmeshev'])}")
print(f"  PV  markers: {sorted(set(GOLD)&ROUTES['PV_marker_Velmeshev'])}")
