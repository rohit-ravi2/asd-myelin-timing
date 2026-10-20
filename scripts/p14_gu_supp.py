#!/usr/bin/env python3
"""
p14_gu_supp.py -- Phase 14, step 4.

Read the supplementary workbook of Gu Y et al., Mol Psychiatry 2025;30(8):3393-3403
(PMID 39994426, PMC12240833) and extract:
  * the Index, so each supplementary table is identified by its stated contents;
  * every table reporting genetic correlations (rg) between autism and ICVF;
  * in particular the WHITE MATTER TRACT rg table, whose SIGN the main text does
    not report -- this is the number that decides whether this project's positive
    WM-tract ICVF rg conflicts with Gu or merely sits in their null region.

Dependency-free xlsx reader (no openpyxl on this machine).
Outputs into outputs/phase14/: gu2025_index.tsv, gu2025_rg_tables/*.tsv, p14_gu_supp.log
"""
import zipfile, re, os, sys, html

XLSX = "./outputs/phase14/gu_MOESM2_ESM.xlsx"
OUT = "./outputs/phase14"
os.makedirs(f"{OUT}/gu2025_rg_tables", exist_ok=True)

z = zipfile.ZipFile(XLSX)
wb = z.read("xl/workbook.xml").decode("utf8")
sheets = re.findall(r'<sheet[^>]*name="([^"]+)"[^>]*r:id="(rId\d+)"', wb)
rels = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"',
                       z.read("xl/_rels/workbook.xml.rels").decode("utf8")))

# shared strings
sst = []
if "xl/sharedStrings.xml" in z.namelist():
    ss = z.read("xl/sharedStrings.xml").decode("utf8")
    for si in re.findall(r"<si>(.*?)</si>", ss, re.S):
        sst.append(html.unescape("".join(re.findall(r"<t[^>]*>(.*?)</t>", si, re.S))))


def col_idx(ref):
    m = re.match(r"([A-Z]+)", ref)
    n = 0
    for ch in m.group(1):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_sheet(target):
    path = "xl/" + target.lstrip("/").replace("worksheets/", "worksheets/")
    if path not in z.namelist():
        path = target.lstrip("/")
        if not path.startswith("xl/"):
            path = "xl/" + path
    if path not in z.namelist():
        return []
    xml = z.read(path).decode("utf8")
    rows = []
    for rx in re.findall(r"<row[^>]*>(.*?)</row>", xml, re.S):
        cells = {}
        for c in re.finditer(r'<c r="([A-Z]+\d+)"([^>]*)>(.*?)</c>', rx, re.S):
            ref, attrs, body = c.group(1), c.group(2), c.group(3)
            t = re.search(r't="(\w+)"', attrs)
            t = t.group(1) if t else "n"
            v = re.search(r"<v>(.*?)</v>", body, re.S)
            if t == "s" and v:
                val = sst[int(v.group(1))] if int(v.group(1)) < len(sst) else ""
            elif t == "inlineStr":
                val = html.unescape("".join(re.findall(r"<t[^>]*>(.*?)</t>", body, re.S)))
            elif v:
                val = v.group(1)
            else:
                val = ""
            cells[col_idx(ref)] = val
        if cells:
            w = max(cells) + 1
            rows.append([cells.get(i, "") for i in range(w)])
    return rows


log = open(f"{OUT}/p14_gu_supp.log", "w")
def out(*a):
    print(*a); print(*a, file=log)

out("=== Gu et al. 2025 (PMID 39994426) supplementary tables ===")
out(f"workbook: {XLSX}\nsheets: {len(sheets)}\n")

sheet_rows = {}
for name, rid in sheets:
    sheet_rows[name] = read_sheet(rels.get(rid, ""))

# ---- index
out("--- Index sheet ---")
with open(f"{OUT}/gu2025_index.tsv", "w") as f:
    for r in sheet_rows.get("Index", []):
        line = "\t".join(str(x) for x in r).strip()
        if line:
            out("  " + line[:200])
            f.write(line + "\n")

# ---- find rg tables
out("")
out("--- sheets whose header mentions genetic correlation / rg ---")
cands = []
for name, rows in sheet_rows.items():
    blob = " ".join(" ".join(str(c) for c in r) for r in rows[:6]).lower()
    if re.search(r"genetic correlation|\brg\b|gencor", blob):
        cands.append(name)
        out(f"  {name}: {blob[:220]}")

# ---- white matter tract sheets
out("")
out("--- sheets mentioning white matter tract acronyms (AutoPtx) ---")
TRACTS = ["mcp", "ml", "cst", "ar", "atr", "str", "ptr", "slf", "ilf", "ifo",
          "unc", "cgc", "cgh", "fmi", "fma"]
wm = []
for name, rows in sheet_rows.items():
    blob = " ".join(" ".join(str(c) for c in r) for r in rows).lower()
    hits = sum(1 for t in TRACTS if re.search(rf"\b{t}(_[lr])?\b", blob))
    if hits >= 8:
        wm.append((name, hits))
        out(f"  {name}: {hits}/15 tract acronyms, {len(rows)} rows")

# ---- dump every candidate
for name in sorted(set(cands) | {n for n, _ in wm}):
    rows = sheet_rows[name]
    with open(f"{OUT}/gu2025_rg_tables/{name}.tsv", "w") as f:
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")

out("")
out(f"dumped {len(set(cands) | {n for n, _ in wm})} sheets to {OUT}/gu2025_rg_tables/")
log.close()
