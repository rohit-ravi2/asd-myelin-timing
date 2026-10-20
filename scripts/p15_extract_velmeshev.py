#!/usr/bin/env python3
"""
p15 STEP 1 -- extract the gene rows needed to apply the atlas-2 OL maturation clock
to the Velmeshev 2019 ASD cortex snRNA-seq matrix.

INPUT  (on disk, data/dvad-corpus/asd_scope/velmeshev2019/):
  raw/matrix.mtx    65,217 genes x 104,559 nuclei, 238,268,142 nonzeros, raw 10x UMI counts
  raw/genes.tsv     ENSG<TAB>symbol
  raw/barcodes.tsv  cell barcode, same order as meta.tsv rows
  meta.tsv          cell / cluster / sample / individual / region / age / sex / diagnosis /
                    Capbatch / Seqbatch / PMI / RIN / genes / UMIs / mito% / ribo%

WHY A SEPARATE STEP: the .mtx is 3.3 GB of text. We never load it into Python. A single
mawk pass (a) accumulates per-gene total counts + nnz and per-cell total counts, and
(b) dumps ONLY the rows we need: the union of the atlas-2 stage markers plus a
seed-fixed random pool of 6,000 genes used later to build expression-matched NULL
maturation scores. Everything is written to outputs/phase15/extract/ so the analysis
step is re-runnable without touching the 3.3 GB file again.

SEED 20260921 fixed before any gene was looked at.
"""
import subprocess, sys, time
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]
VEL = Path("data/dvad-corpus/asd_scope/velmeshev2019")
RAW = VEL / "raw"
OUT = BASE / "outputs/phase15/extract"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20260921
N_RANDOM_POOL = 6000

# ---------------------------------------------------------------- stage marker panels
# VERBATIM from scripts/p6_atlas2_replication.py (the script whose age-monotonicity is the
# clock's internal validation). p10's slightly trimmed panels are carried as a sensitivity.
STAGE_P6 = {
    "OPC":  "PDGFRA CSPG4 PTPRZ1 LHFPL3 SOX6 OLIG1 DSCAM GPR17".split(),
    "COP":  "GPR17 BCAS1 ENPP6 TMEM235 NKX2-2 SIRT2 FYN".split(),
    "NFOL": "BCAS1 ENPP6 ITPR2 CNP TCF7L2 NKX6-2".split(),
    "MFOL": "MAL MOG OPALIN CTPS1 CNP PLP1".split(),
    "MOL":  "PLP1 MBP MOBP CLDN11 APOD KLK6 ASPA CARNS1 SLC5A11".split(),
}
STAGE_P10 = {
    "OPC":  "PDGFRA CSPG4 PTPRZ1 LHFPL3 SOX6 OLIG1 DSCAM".split(),
    "COP":  "GPR17 BCAS1 ENPP6 TMEM235 NKX2-2 SIRT2 FYN".split(),
    "NFOL": "BCAS1 ENPP6 ITPR2 CNP TCF7L2 NKX6-2".split(),
    "MFOL": "MAL MOG OPALIN CTPS1 CNP PLP1".split(),
    "MOL":  "PLP1 MBP MOBP CLDN11 APOD KLK6 ASPA CARNS1".split(),
}
# analogous astrocyte maturation panel (immature -> mature), used as a lineage-matched
# "analogous score" control. Committed before any result was seen.
ASTRO_PANEL = {
    "A0": "VIM NES SOX2 HES1 FABP7 PTPRZ1 TNC".split(),          # immature / radial-glia-like
    "A1": "SPARCL1 ATP1A2 SLC1A3 GJA1 CLU".split(),               # intermediate
    "A2": "AQP4 GJB6 SLC1A2 GLUL MT1X MT2A CD44 ALDH1L1".split(), # mature protoplasmic/fibrous
}
# analogous excitatory-neuron maturation panel
NEURON_PANEL = {
    "N0": "SOX11 SOX4 DCX NEUROD2 NEUROD6 TBR1 BCL11B".split(),
    "N1": "MEF2C SATB2 RORB NEFL".split(),
    "N2": "SNAP25 SYT1 NEFM NEFH CAMK2A SYN1 GRIN1 CUX2".split(),
}

ALL_MARKERS = sorted(set(
    [g for v in STAGE_P6.values() for g in v] +
    [g for v in STAGE_P10.values() for g in v] +
    [g for v in ASTRO_PANEL.values() for g in v] +
    [g for v in NEURON_PANEL.values() for g in v]
))


def main():
    t0 = time.time()
    # ---------------- genes
    sym, ens = [], []
    with open(RAW / "genes.tsv") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            ens.append(p[0])
            sym.append(p[1] if len(p) > 1 else p[0])
    sym = np.array(sym)
    n_genes = len(sym)
    print(f"[genes] {n_genes}", file=sys.stderr)

    # symbol -> row index (1-based for mtx). Duplicated symbols: keep the row with the
    # highest total count, decided AFTER the stats pass; for now keep all candidate rows.
    marker_rows = np.flatnonzero(np.isin(sym, ALL_MARKERS)) + 1
    missing = sorted(set(ALL_MARKERS) - set(sym.tolist()))
    print(f"[markers] {len(ALL_MARKERS)} requested, {len(marker_rows)} matrix rows matched, "
          f"{len(missing)} symbols absent: {missing}", file=sys.stderr)

    rng = np.random.default_rng(SEED)
    rand_rows = rng.choice(np.arange(1, n_genes + 1), size=N_RANDOM_POOL, replace=False)
    keep = np.union1d(marker_rows, rand_rows)
    print(f"[keep] {len(keep)} matrix rows dumped", file=sys.stderr)

    krows = OUT / "keep_rows.txt"
    krows.write_text("\n".join(map(str, keep)) + "\n")

    # ---------------- one mawk pass
    ent = OUT / "entries.tsv"
    gstat = OUT / "gene_stats.tsv"
    cstat = OUT / "cell_totals.tsv"
    if ent.exists() and gstat.exists() and cstat.exists():
        print("[awk] outputs present, skipping pass", file=sys.stderr)
    else:
        prog = r'''
FNR==NR { K[$1]=1; next }
FNR<=3  { next }
{ gs[$1]+=$3; gn[$1]++; ct[$2]+=$3; cn[$2]++
  if ($1 in K) print $1"\t"$2"\t"$3 > ENT }
END {
  for (g=1; g<=NG; g++) print g"\t"(gs[g]+0)"\t"(gn[g]+0) > GST
  for (c=1; c<=NC; c++) print c"\t"(ct[c]+0)"\t"(cn[c]+0) > CST
}
'''
        cmd = ["mawk", "-v", f"ENT={ent}", "-v", f"GST={gstat}", "-v", f"CST={cstat}",
               "-v", f"NG={n_genes}", "-v", "NC=104559", prog,
               str(krows), str(RAW / "matrix.mtx")]
        print("[awk] streaming 238M entries ...", file=sys.stderr)
        subprocess.run(cmd, check=True)
    print(f"[awk] done in {time.time()-t0:.0f}s", file=sys.stderr)

    # ---------------- pack to npz
    gs = np.loadtxt(gstat, dtype=np.int64)
    cs = np.loadtxt(cstat, dtype=np.int64)
    e = np.loadtxt(ent, dtype=np.int64)
    np.savez_compressed(
        OUT / "velmeshev_subset.npz",
        ent_gene=e[:, 0].astype(np.int32), ent_cell=e[:, 1].astype(np.int32),
        ent_val=e[:, 2].astype(np.int32),
        gene_total=gs[:, 1], gene_nnz=gs[:, 2],
        cell_total=cs[:, 1], cell_nnz=cs[:, 2],
        symbols=sym, ensembl=np.array(ens), keep_rows=keep.astype(np.int32),
    )
    print(f"[pack] {e.shape[0]} entries kept -> velmeshev_subset.npz "
          f"({(OUT/'velmeshev_subset.npz').stat().st_size/1e6:.0f} MB) "
          f"total {time.time()-t0:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
