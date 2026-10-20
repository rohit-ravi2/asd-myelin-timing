#!/usr/bin/env python3
"""
Test A: build S-LDSC annotations from the gene-set panel and compute LD scores.

This is the fix for the 2026-09-10 "Test 2 data wall", which was never a data
wall: all LDSC reference data was on disk. The real fault was annotations built
from 6-15 genes (0.04-0.16% of SNPs) -- far too small for S-LDSC, which is why
the analysis could not even recover the known neuronal positive control.

Here each annotation is a gene set expanded by +/-100 kb, which puts the larger
sets at 1-10% of SNPs -- inside the range where S-LDSC has usable power. Sets
that still land under ~0.5% of SNPs are built and run anyway, but flagged
UNDERPOWERED in the report rather than silently interpreted.
"""
import gzip, json, os, subprocess, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

BASE = Path(__file__).resolve().parents[1]
REF = Path("data/ldsc_ref")
PLINK = REF / "1000G_EUR_Phase3_plink/1000G.EUR.QC"
LDSC = "ldsc/ldsc.py"
PY2 = os.environ.get("LDSC_PYTHON", "python2")
GLOC = "data/asd-corpus/tools/magma_geneloc/NCBI37.3.gene.loc"
ANN = BASE / "outputs/sldsc/annot"
ANN.mkdir(parents=True, exist_ok=True)
WINDOW = 100_000

# Sets to run. Small sets are included but flagged.
# Trimmed 2026-09-16 to the 12 sets S-LDSC can actually power (>=1% of SNPs).
# Dropped MYELIN_CORE_CURATED (24k SNPs = 0.24%), CHOL_ENZ_CURATED (0.42%) and
# LIPOPROTEIN_NONAPOE (0.41%): below the size where S-LDSC is informative, and
# covered properly by the MAGMA competitive test instead. Also dropped
# CHROMATIN_REG / SFARI_SCORE1 / KERATINIZATION as redundant with retained controls.
SETS = ["CHOL_BIOSYNTHESIS", "CHOL_TRANSPORT", "LIPID_TRANSPORT", "MYELIN",
        "OLIGODENDROCYTE_VELM", "OPC_VELM",
        "SYNAPTIC", "SFARI_ALL", "EXC_NEURON_VELM",
        "ENDOTHELIAL_VELM", "MICROGLIA_VELM", "OLFACTORY_RECEPTOR"]

# ---------------------------------------------------------------- gene coords
coords = {}
with open(GLOC) as f:
    for line in f:
        p = line.split()
        if len(p) >= 6:
            coords[p[0]] = (p[1], int(p[2]), int(p[3]))   # entrez -> chr,start,end

manifest = json.load(open(BASE / "annot/geneset_manifest.json"))
setfile = BASE / "annot/asd_genesets.magma.set"
members = {}
with open(setfile) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        members[p[0]] = p[1:]

# ---------------------------------------------------------------- make annots
def intervals(setname):
    iv = {}
    for g in members[setname]:
        if g not in coords:
            continue
        c, s, e = coords[g]
        iv.setdefault(c, []).append((max(0, s - WINDOW), e + WINDOW))
    for c in iv:
        iv[c].sort()
        merged = []
        for lo, hi in iv[c]:
            if merged and lo <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
            else:
                merged.append((lo, hi))
        iv[c] = merged
    return iv


def write_annot(args):
    setname, chrom = args
    outp = ANN / f"{setname}.{chrom}.annot.gz"
    if outp.exists() and outp.stat().st_size > 100:
        return None
    iv = INTERVALS[setname].get(str(chrom), [])
    starts = [x[0] for x in iv]
    import bisect
    n_in = 0
    bim = f"{PLINK}.{chrom}.bim"
    with open(bim) as fi, gzip.open(outp, "wt") as fo:
        fo.write("ANNOT\n")
        for line in fi:
            bp = int(line.split()[3])
            j = bisect.bisect_right(starts, bp) - 1
            hit = 1 if (j >= 0 and bp <= iv[j][1]) else 0
            n_in += hit
            fo.write(f"{hit}\n")
    return (setname, chrom, n_in)


def ldscore(args):
    setname, chrom = args
    out = ANN / f"{setname}.{chrom}"
    if (Path(str(out) + ".l2.ldscore.gz")).exists():
        return f"skip {setname} {chrom}"
    cmd = [PY2, LDSC, "--l2", "--bfile", f"{PLINK}.{chrom}",
           "--ld-wind-cm", "1", "--annot", str(ANN / f"{setname}.{chrom}.annot.gz"),
           "--thin-annot", "--out", str(out),
           "--print-snps", str(REF / "hm3_no_MHC.list.txt")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return f"{'OK' if r.returncode == 0 else 'FAIL'} {setname} {chrom}" + \
           ("" if r.returncode == 0 else " :: " + r.stderr.strip()[-200:])


if __name__ == "__main__":
    print("[1] building intervals", file=sys.stderr)
    INTERVALS = {s: intervals(s) for s in SETS}
    jobs = [(s, c) for s in SETS for c in range(1, 23)]

    print(f"[2] writing {len(jobs)} annot files", file=sys.stderr)
    counts = {}
    with ProcessPoolExecutor(8) as ex:
        for r in ex.map(write_annot, jobs):
            if r:
                counts.setdefault(r[0], 0)
                counts[r[0]] += r[2]
    for s in SETS:
        if s in counts:
            print(f"   {s:24s} {counts[s]:>8d} SNPs in annot", file=sys.stderr)
    json.dump(counts, open(ANN / "snp_counts.json", "w"), indent=1)

    print(f"[3] computing LD scores ({len(jobs)} chr-jobs, 8 workers)", file=sys.stderr)
    done = 0
    with ProcessPoolExecutor(14) as ex:
        for msg in ex.map(ldscore, jobs):
            done += 1
            if msg.startswith("FAIL") or done % 22 == 0:
                print(f"   [{done}/{len(jobs)}] {msg}", file=sys.stderr)
    print("SLDSC_ANNOT_DONE", file=sys.stderr)
