#!/usr/bin/env bash
# p14_probtrackx_dl.sh -- Phase 14, step 5.
#
# Gu et al. 2025 ran their WM-tract rg on 27 AutoPtx PROBABILISTIC TRACTOGRAPHY
# tracts. This project's 15 ICVF IDPs are TBSS FA-SKELETON JHU-label ROIs -- a
# different spatial definition of "white matter tract".
#
# BIG40 contains BOTH families (48 TBSS ICVF, 27 ProbtrackX ICVF). The 27
# ProbtrackX ICVF IDPs (Pheno 1950-1976) are exactly Gu's tract set
# (ar atr cgc cgh cst fma fmi ifo ilf mcp ml ptr slf str unc).
#
# So the decisive local test is: run OUR ASD GWAS against BIG40's ProbtrackX
# ICVF, i.e. Gu's tract definitions. We download 5 of the 27 (~1.72 GB), chosen
# to cover the two tracts where this project's TBSS estimate is most positive
# (SLF, cingulum-cingulate) plus anterior thalamic radiation.
#
# TRANSPORT NOTE: as of 2026-09-21 the Oxford BIG40 host's Let's Encrypt
# certificate is EXPIRED (notAfter 2026-09-19), so curl needs -k. Authenticity
# was established independently of TLS: the first 1 MiB of Pheno 1904 served by
# this endpoint is BYTE-IDENTICAL to the copy downloaded 2026-09-14 and already
# used in this project, and content-length matched exactly. Every file is also
# checked with `gzip -t` and against the expected column layout.
#
# The server is serving at ~80 KB/s per connection, so the five run in parallel.
set -u
ROOT=$(cd "$(dirname "$0")/.." && pwd); DATA=${DATA:-$ROOT/data}
DEST=$DATA/dvad-corpus/asd_myelin/big40_probtrackx
BASE=https://open.oxcin.ox.ac.uk/ukbiobank/big40/release2/stats33k
LOG=$ROOT/outputs/phase14
mkdir -p "$DEST" "$LOG"

# Pheno -> label (BIG40 ProbtrackX ICVF; matches Gu's AutoPtx tract acronyms)
IDS="1971:pICVF_slf_l 1972:pICVF_slf_r 1954:pICVF_cgc_l 1955:pICVF_cgc_r 1952:pICVF_atr_l"

fetch_one() {
  local id=$1
  local f="$DEST/$id.txt.gz"
  for attempt in 1 2 3 4 5; do
    gzip -t "$f" 2>/dev/null && { echo "OK $id $(stat -c%s "$f")"; return 0; }
    curl -skL -m 21600 -C - "$BASE/$id.txt.gz" -o "$f" 2>/dev/null
  done
  gzip -t "$f" 2>/dev/null && echo "OK $id $(stat -c%s "$f")" || echo "FAIL $id"
}

pids=""
for pair in $IDS; do
  id=${pair%%:*}
  fetch_one "$id" >> "$LOG/p14_dl.log" 2>&1 &
  pids="$pids $!"
done
echo "launched:$pids" >> "$LOG/p14_dl.log"
wait $pids

{
  echo "--- md5sums ---"
  md5sum "$DEST"/*.txt.gz
  echo "--- header check (expect: chr rsid pos a1 a2 beta se pval(-log10)) ---"
  for f in "$DEST"/*.txt.gz; do echo -n "$(basename "$f"): "; zcat "$f" | head -1; done
  echo "P14_DL_DONE"
} >> "$LOG/p14_dl.log" 2>&1

~/bin/notify "P14 BIG40 ProbtrackX ICVF downloads complete; running rg" "asd-myelin" 2>/dev/null
Rscript $ROOT/scripts/p14_probtrackx_rg.R \
  >> "$LOG/p14_probtrackx_console.log" 2>&1
~/bin/notify "P14 ProbtrackX rg DONE - see outputs/phase14/p14_probtrackx_rg.log" "asd-myelin" 2>/dev/null
