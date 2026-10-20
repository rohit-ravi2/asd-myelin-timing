#!/usr/bin/env python3
"""
Build the gene-set panel for the ASD lipid/myelin analyses.

Design discipline: every set is declared here BEFORE any association result is
seen, and each is tagged with its role -- PRIMARY (hypothesis), POSCTRL (must be
enriched or the instrument is broken), NEGCTRL (must NOT be enriched), or
RANDOM (size-matched null calibration).

Outputs a MAGMA-format .set file (one line per set: NAME gene1 gene2 ...) using
ENTREZ ids, plus a symbol-keyed JSON manifest for reporting.
"""
import csv, gzip, json, random, re, sys
from pathlib import Path
from collections import defaultdict

CORPUS = Path("data/asd-corpus")
SCOPE  = Path("data/dvad-corpus/asd_scope")
OUT    = Path("./annot")
OUT.mkdir(parents=True, exist_ok=True)
random.seed(20260916)

# ---------------------------------------------------------------- gene universe
# MAGMA gene.loc is the authoritative universe: only genes MAGMA can test.
sym2entrez, entrez2sym, gene_chr = {}, {}, {}
with open(CORPUS / "tools/magma_geneloc/NCBI37.3.gene.loc") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            p = line.split()
        eid, chrom, sym = p[0], p[1], p[5]
        sym2entrez[sym] = eid
        entrez2sym[eid] = sym
        gene_chr[eid] = chrom
UNIVERSE = set(entrez2sym)
print(f"[universe] {len(UNIVERSE)} genes in MAGMA NCBI37.3 gene.loc", file=sys.stderr)


def to_entrez(symbols):
    """Map symbols -> entrez, keeping only genes in the MAGMA universe."""
    out = set()
    for s in symbols:
        s = s.strip().upper()
        if s in sym2entrez:
            out.add(sym2entrez[s])
    return out


# ---------------------------------------------------------------- GMT loading
def load_gmt(path):
    sets = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            sets[p[0]] = [g for g in p[2:] if g]
    return sets


gmts = {}
for name in ["c2.cp.v2024.1.Hs.symbols.gmt",
             "c5.go.bp.v2024.1.Hs.symbols.gmt",
             "c5.go.cc.v2024.1.Hs.symbols.gmt"]:
    p = CORPUS / "pathways" / name
    if p.exists():
        gmts[name] = load_gmt(p)
        print(f"[gmt] {name}: {len(gmts[name])} sets", file=sys.stderr)
ALLGMT = {}
for g in gmts.values():
    ALLGMT.update(g)


def gmt_union(patterns, exclude=()):
    """Union of all MSigDB sets whose NAME matches any regex in `patterns`."""
    genes, used = set(), []
    for name, members in ALLGMT.items():
        if any(re.search(p, name, re.I) for p in patterns) and \
           not any(re.search(e, name, re.I) for e in exclude):
            genes |= set(members)
            used.append(name)
    return genes, sorted(used)


# ---------------------------------------------------------------- PRIMARY sets
PANEL = {}   # name -> dict(role=, genes=set(entrez), source=, note=)


def add(name, role, symbols, source, note="", srcsets=None):
    e = to_entrez(symbols)
    PANEL[name] = dict(role=role, genes=e, source=source, note=note,
                       n_input=len(set(symbols)), n_mapped=len(e),
                       srcsets=srcsets or [])
    print(f"[set] {name:38s} {role:8s} {len(e):5d} genes "
          f"(of {len(set(symbols))} input)", file=sys.stderr)


# --- 1. cholesterol biosynthesis (the master hypothesis set)
chol_syn, s1 = gmt_union([
    r"CHOLESTEROL_BIOSYNTH", r"STEROL_BIOSYNTH",
    r"REACTOME_CHOLESTEROL_BIOSYNTHESIS", r"KEGG.*STEROID_BIOSYNTHESIS",
    r"GOBP_CHOLESTEROL_BIOSYNTHETIC_PROCESS",
    r"GOBP_STEROL_BIOSYNTHETIC_PROCESS",
])
add("CHOL_BIOSYNTHESIS", "PRIMARY", chol_syn, "MSigDB C2CP+C5BP union", srcsets=s1)

# --- 2. cholesterol/sterol transport + homeostasis
chol_tr, s2 = gmt_union([
    r"CHOLESTEROL_TRANSPORT", r"STEROL_TRANSPORT", r"CHOLESTEROL_EFFLUX",
    r"CHOLESTEROL_HOMEOSTASIS", r"STEROL_HOMEOSTASIS",
    r"LIPOPROTEIN_PARTICLE", r"PLASMA_LIPOPROTEIN",
])
add("CHOL_TRANSPORT", "PRIMARY", chol_tr, "MSigDB union", srcsets=s2)

# --- 3. broader lipid transport (the original "lipid-transport proteins" idea)
lip_tr, s3 = gmt_union([r"LIPID_TRANSPORT", r"LIPID_LOCALIZATION",
                        r"PHOSPHOLIPID_TRANSPORT", r"FATTY_ACID_TRANSPORT"])
add("LIPID_TRANSPORT", "PRIMARY", lip_tr, "MSigDB union", srcsets=s3)

# --- 4. myelin structural / myelination
mye, s4 = gmt_union([r"MYELIN", r"ENSHEATHMENT_OF_NEURONS",
                     r"AXON_ENSHEATHMENT"])
add("MYELIN", "PRIMARY", mye, "MSigDB union", srcsets=s4)

# --- 5. canonical myelin structural genes (hand-curated, high confidence)
MYELIN_CORE = """MBP PLP1 MOG MAG CNP MOBP MAL CLDN11 GJB1 GJC2 PMP22 MPZ
SOX10 OLIG1 OLIG2 MYRF ERMN ASPA UGT8 GAL3ST1 FA2H CERS2 ASPA TF CA2 QKI
SIRT2 NKX6-2 ZNF488 PLLP OPALIN""".split()
add("MYELIN_CORE_CURATED", "PRIMARY", MYELIN_CORE, "hand-curated canonical")

# --- 6. cholesterol synthesis enzymes, hand-curated (the pathway proper)
CHOL_ENZ = """HMGCR HMGCS1 MVK PMVK MVD IDI1 IDI2 FDPS FDFT1 SQLE LSS
CYP51A1 TM7SF2 MSMO1 NSDHL HSD17B7 EBP SC5D DHCR7 DHCR24 SREBF1 SREBF2
SCAP MBTPS1 MBTPS2 INSIG1 INSIG2 ACAT1 ACAT2 SOAT1 CYP46A1 ABCA1 ABCG1
ABCA2 APOE APOA1 LDLR LRP1 LRP8 VLDLR NPC1 NPC2 STARD3 OSBP OSBPL1A
SCARB1 CETP PLTP SORT1""".split()
add("CHOL_ENZ_CURATED", "PRIMARY", CHOL_ENZ, "hand-curated synth+transport")

# --- 7/8. oligodendrocyte + OPC markers (Velmeshev 2019, human ASD cortex)
def velmeshev(fname):
    p = SCOPE / "velmeshev2019/markers" / fname
    op = gzip.open if fname.endswith(".gz") else open
    with op(p, "rt") as f:
        rd = csv.DictReader(f, delimiter="\t")
        col = "symbol" if "symbol" in (rd.fieldnames or []) else rd.fieldnames[1]
        return [r[col] for r in rd if r.get(col)]

add("OLIGODENDROCYTE_VELM", "PRIMARY", velmeshev("OL.tsv.gz"),
    "Velmeshev2019 OL markers", note="human ASD cortex snRNA-seq")
add("OPC_VELM", "PRIMARY", velmeshev("OPC.tsv.gz"),
    "Velmeshev2019 OPC markers")

# ---------------------------------------------------------------- POSITIVE controls
# ASD rare-variant architecture converges on chromatin (~57%) and synaptic (~24%).
chrom, s5 = gmt_union([r"CHROMATIN_REMODEL", r"CHROMATIN_ORGANIZATION",
                       r"HISTONE_MODIFICATION", r"REGULATION_OF_TRANSCRIPTION_BY_RNA_POLYMERASE_II"])
add("CHROMATIN_REG", "POSCTRL", chrom, "MSigDB union", srcsets=s5)

syn, s6 = gmt_union([r"SYNAP", r"POSTSYNAPTIC_DENSITY", r"NEUROTRANSMITTER_SECRETION"])
add("SYNAPTIC", "POSCTRL", syn, "MSigDB union", srcsets=s6)

# SFARI high-confidence
sfari_p = SCOPE / "sfari/SFARI-Gene_genes_2026-09-10.csv"
sf_all, sf_s1 = [], []
with open(sfari_p) as f:
    for r in csv.DictReader(f):
        g = r.get("gene-symbol") or r.get("gene_symbol") or r.get("Gene Symbol")
        sc = (r.get("gene-score") or r.get("gene_score") or "").strip()
        if g:
            sf_all.append(g)
            if sc == "1":
                sf_s1.append(g)
add("SFARI_ALL", "POSCTRL", sf_all, "SFARI Gene 2026-09-10")
add("SFARI_SCORE1", "POSCTRL", sf_s1, "SFARI score-1 (high confidence)")

# excitatory neuron markers (Velmeshev) -- ASD's strongest cell type
exc = []
for f in ["L2_3.tsv.gz", "L4.tsv.gz", "L5.tsv.gz", "L6.tsv.gz"]:
    if (SCOPE / "velmeshev2019/markers" / f).exists():
        exc += velmeshev(f)
add("EXC_NEURON_VELM", "POSCTRL", exc, "Velmeshev2019 L2/3+L4+L5+L6")

# ---------------------------------------------------------------- NEGATIVE controls
add("ENDOTHELIAL_VELM", "NEGCTRL", velmeshev("END.tsv.gz"), "Velmeshev2019 endothelial")
add("MICROGLIA_VELM", "NEGCTRL", velmeshev("MG.tsv.gz"), "Velmeshev2019 microglia")

olf, s7 = gmt_union([r"OLFACTORY"])
add("OLFACTORY_RECEPTOR", "NEGCTRL", olf, "MSigDB union", srcsets=s7)
ker, s8 = gmt_union([r"KERATIN", r"CORNIFIED"])
add("KERATINIZATION", "NEGCTRL", ker, "MSigDB union", srcsets=s8)

# ---------------------------------------------------------------- APOE contrast
add("APOE_ONLY", "CONTRAST", ["APOE"], "single gene")
add("APOE_CLUSTER", "CONTRAST", "APOE APOC1 APOC2 APOC4 TOMM40 NECTIN2 BCL3 EXOC3L2".split(),
    "chr19q13.32 APOE locus genes")
add("LIPOPROTEIN_NONAPOE", "CONTRAST",
    [g for g in CHOL_ENZ if g != "APOE"], "cholesterol set minus APOE")

# ---------------------------------------------------------------- RANDOM null sets
# size-matched to the primary sets, to calibrate what "enrichment" looks like by chance
prim_sizes = sorted({len(PANEL[k]["genes"]) for k in PANEL
                     if PANEL[k]["role"] == "PRIMARY"})
univ = sorted(UNIVERSE)
nrand_per_size = 20
for sz in prim_sizes:
    if sz < 5:
        continue
    for i in range(nrand_per_size):
        PANEL[f"RANDOM_{sz}_{i:02d}"] = dict(
            role="RANDOM", genes=set(random.sample(univ, sz)),
            source="size-matched random", note="", n_input=sz, n_mapped=sz, srcsets=[])
print(f"[set] {nrand_per_size} random sets x {len(prim_sizes)} sizes", file=sys.stderr)

# ---------------------------------------------------------------- write
setfile = OUT / "asd_genesets.magma.set"
with open(setfile, "w") as f:
    for name, d in PANEL.items():
        if len(d["genes"]) < 3:
            print(f"[skip] {name}: only {len(d['genes'])} genes", file=sys.stderr)
            continue
        f.write(name + "\t" + "\t".join(sorted(d["genes"])) + "\n")

manifest = {name: dict(role=d["role"], n_genes=len(d["genes"]),
                       n_input=d["n_input"], source=d["source"], note=d["note"],
                       symbols=sorted(entrez2sym[g] for g in d["genes"])
                       if d["role"] != "RANDOM" else [],
                       member_sets=d["srcsets"])
            for name, d in PANEL.items()}
with open(OUT / "geneset_manifest.json", "w") as f:
    json.dump(manifest, f, indent=1)

named = [k for k in PANEL if PANEL[k]["role"] != "RANDOM"]
print(f"\n[done] {len(PANEL)} sets written ({len(named)} named + "
      f"{len(PANEL)-len(named)} random) -> {setfile}", file=sys.stderr)
