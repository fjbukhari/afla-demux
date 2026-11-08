#!/usr/bin/env python3
"""One-page CFTR report for one sample (education and research use only, not for clinical use).

  afla_cftr.py --sample S1 --vcf S1.annotated.vcf.gz --polyt S1.polyt.json --cnv cnv.json \\
      --thresholds S1.targets.thresholds.bed.gz --regions S1.targets.regions.bed.gz \\
      --fasta hg38.fa --cftr2 CFTR2_30January2026.xlsx --out S1.cftr.html --json S1.cftr.json

Variants are matched to the CFTR2 list by GRCh38 position, reference and alternative base (taken from its
"Genomic coordinates" sheet and left-aligned like a normalised VCF), so naming differences do not matter.
cDNA names are written from the genome on NM_000492.4 with the HGVS 3' rule. ClinVar is read from the VEP
annotation when it is present. A statement that nothing was found is always limited to the regions that were
adequately covered. Standard library only (the Excel file is read directly).
"""
import argparse
import gzip
import html
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import afla_cftr_hgvs as hg  # noqa: E402

GENE_START, GENE_END = 117470000, 117670000  # generous CFTR window, GRCh38
CF_CAUSING, VCC, NON_CF, NO_INTERP = ("CF-causing", "Varying clinical consequence", "Non CF-causing",
                                      "No interpretation available")
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


# ---------------------------------------------------------------- CFTR2 list
def read_xlsx(path):
    """{sheet name: [rows]}; each row is a list of cell values (text or number, None when empty)."""
    with zipfile.ZipFile(path) as z:
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
                strings.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
        book = ET.fromstring(z.read("xl/workbook.xml"))
        rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
        out = {}
        for sheet in book.find("m:sheets", NS):
            target = rels[sheet.get(f"{{{NS['r']}}}id")].lstrip("/")
            name = target if target.startswith("xl/") else "xl/" + target
            rows = []
            for row in ET.fromstring(z.read(name)).find("m:sheetData", NS):
                vals = []
                for c in row:
                    col = 0
                    for ch in re.match(r"[A-Z]+", c.get("r")).group(0):
                        col = col * 26 + ord(ch) - 64
                    vals.extend([None] * (col - len(vals)))
                    t, v = c.get("t"), c.find("m:v", NS)
                    if t == "inlineStr":
                        val = "".join(x.text or "" for x in c.iter(f"{{{NS['m']}}}t"))
                    elif v is None or v.text is None:
                        val = None
                    elif t == "s":
                        val = strings[int(v.text)]
                    elif t in ("str", "b", "e"):
                        val = v.text
                    else:
                        val = float(v.text)
                        val = int(val) if val.is_integer() else val
                    vals[col - 1] = val
                rows.append(vals)
            out[sheet.get("name")] = rows
    return out


def _header(rows, first):
    for i, r in enumerate(rows):
        if r and isinstance(r[0], str) and r[0].startswith(first):
            return i, [str(c or "") for c in r]
    raise ValueError(f"header starting with {first!r} not found in the CFTR2 file")


def _cell(row, idx):
    return row[idx] if idx < len(row) else None


def load_cftr2(path, fasta):
    """Returns ({(pos, ref, alt): record}, {legacy name: record}, date text)."""
    book = read_xlsx(path)
    main = next(r for n, r in book.items() if "legacy" in n.lower())
    coords = next(r for n, r in book.items() if "coord" in n.lower())
    date = next((str(r[0]) for r in main[:5] if r and str(r[0] or "").startswith("Date")), "")
    h, cols = _header(main, "Variant legacy name")
    det = next(i for i, c in enumerate(cols) if "final determination" in c and "current" in c)
    cd, lg, pr, al = cols.index("Variant cDNA name"), cols.index("Variant legacy name"), cols.index(
        "Variant protein name"), next(i for i, c in enumerate(cols) if c.startswith("# alleles"))
    by_c, by_legacy = {}, {}
    for r in main[h + 1:]:
        cdna = _cell(r, cd)
        if not cdna or not _cell(r, det):
            continue
        rec = {"cdna": str(cdna).strip(), "legacy": str(_cell(r, lg) or "").strip(),
               "protein": str(_cell(r, pr) or "").strip(), "alleles": _cell(r, al),
               "determination": str(_cell(r, det)).strip()}
        by_c[rec["cdna"]] = rec
        by_legacy[rec["legacy"]] = rec
    h, cols = _header(coords, "Variant cDNA name")
    pos_i, ref_i, alt_i = (cols.index("grch38_" + k) for k in ("pos", "ref", "alt"))
    by_pos = {}
    for r in coords[h + 1:]:
        rec = by_c.get(str(_cell(r, 0) or "").strip())
        pos, ref, alt = _cell(r, pos_i), _cell(r, ref_i), _cell(r, alt_i)
        if rec and pos and ref and alt:
            by_pos[hg.left_align(fasta, int(pos), str(ref).upper(), str(alt).upper())] = rec
    return by_pos, by_legacy, date.replace("Date:", "").strip()


# ---------------------------------------------------------------- variant names
def hgvs_name(tx, fasta, pos, ref, alt):
    """cDNA name on NM_000492.4 (+ strand) following the HGVS 3' rule."""
    c = tx.g_to_c
    if len(ref) == 1 and len(alt) == 1:
        return f"c.{c(pos)}{ref}>{alt}"
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    base = lambda p: fasta.fetch(hg.CHROM, p, p)  # noqa: E731
    if not alt:  # deletion: move as far 3' as the sequence allows
        n = len(ref)
        while base(pos + n) == base(pos):
            pos += 1
        return f"c.{c(pos)}del" if n == 1 else f"c.{c(pos)}_{c(pos + n - 1)}del"
    if not ref:  # insertion before pos
        while base(pos) == alt[0]:
            alt, pos = alt[1:] + alt[0], pos + 1
        n = len(alt)
        if fasta.fetch(hg.CHROM, pos - n, pos - 1) == alt:
            return f"c.{c(pos - 1)}dup" if n == 1 else f"c.{c(pos - n)}_{c(pos - 1)}dup"
        return f"c.{c(pos - 1)}_{c(pos)}ins{alt}"
    return f"c.{c(pos)}delins{alt}" if len(ref) == 1 else f"c.{c(pos)}_{c(pos + len(ref) - 1)}delins{alt}"


# ---------------------------------------------------------------- inputs
def _open(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def read_vcf(path, fasta):
    """CFTR calls of the (single) sample in a normalised VCF, annotation fields included when present."""
    fields, out = [], []
    if Path(path).name.startswith("NO_"):
        return out
    with _open(path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = re.search(r"Format: ([^\">]+)", line).group(1).strip().split("|")
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            pos = int(p[1])
            if p[0] != hg.CHROM or not GENE_START <= pos <= GENE_END or len(p) < 10 or p[4] == "*":
                continue
            fmt = dict(zip(p[8].split(":"), p[9].split(":")))
            alleles = re.split(r"[/|]", fmt.get("GT", "."))
            if "." in alleles or set(alleles) == {"0"}:
                continue
            ad = [int(x) for x in fmt.get("AD", "").split(",") if x.isdigit()]
            dp = int(fmt["DP"]) if fmt.get("DP", "").isdigit() else (sum(ad) or None)
            if len(ad) >= 2 and sum(ad):
                vaf = ad[1] / sum(ad)
            else:
                try:
                    vaf = float(fmt.get("VAF", "").split(",")[0])
                except ValueError:
                    vaf = None
            csq = {}
            info = dict(x.split("=", 1) for x in p[7].split(";") if "=" in x)
            if fields and "CSQ" in info:
                csq = dict(zip(fields, info["CSQ"].split(",")[0].split("|")))
            gt = fmt.get("GT", "")
            out.append({"pos": pos, "ref": p[3], "alt": p[4], "filter": p[6], "zygosity": "het" if len(set(alleles)) > 1
                        else "hom", "vaf": vaf, "dp": dp, "csq": csq, "gt": gt, "ps": fmt.get("PS", "") if "|" in gt else "",
                        "key": hg.left_align(fasta, pos, p[3].upper(), p[4].upper())})
    return out


def read_coverage(thresholds, regions, min_depth):
    """Per-target share of bases at >= min_depth, and mean depth. Thresholds come from mosdepth."""
    means = {}
    if regions and not Path(regions).name.startswith("NO_"):
        with _open(regions) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 5:
                    means[p[3]] = float(p[4])
    rows = []
    with _open(thresholds) as f:
        head = f.readline().rstrip("\n").split("\t")
        cols = {int(c[:-1]): i for i, c in enumerate(head) if re.fullmatch(r"\d+X", c)}
        usable = [d for d in cols if d <= min_depth]
        if not usable:
            raise ValueError(f"no depth column at or below {min_depth}x in {thresholds}")
        depth = max(usable)
        for line in f:
            p = line.rstrip("\n").split("\t")
            width = int(p[2]) - int(p[1])
            pct = 100.0 * int(p[cols[depth]]) / width if width else 0.0
            rows.append({"target": p[3], "start": int(p[1]) + 1, "end": int(p[2]), "pct": round(pct, 1),
                         "mean": round(means.get(p[3], 0.0), 0), "gap": pct < 100.0})
    return rows, depth


# ---------------------------------------------------------------- content
def classify(var, tx, fasta, by_pos):
    rec = by_pos.get(var["key"])
    csq = var["csq"]
    clinvar = (csq.get("ClinVar_CLNSIG") or csq.get("CLIN_SIG") or "").replace("_", " ").replace("&", ", ")
    det = rec["determination"] if rec else ""
    pathogenic = "pathogenic" in clinvar.lower() and "conflicting" not in clinvar.lower()
    benign = "benign" in clinvar.lower() and "pathogenic" not in clinvar.lower()
    tier = "other"
    if det in (CF_CAUSING, VCC) or (det == NO_INTERP and not benign) or (not rec and pathogenic):
        tier = "report"
    protein = rec["protein"] if rec and rec["protein"] else (csq.get("HGVSp", "").split(":")[-1] if csq.get("HGVSp") else "")
    return {"hgvs_c": hg_name_safe(tx, fasta, var), "protein": protein, "legacy": rec["legacy"] if rec else "",
            "region": tx.region_label(var["pos"]), "zygosity": var["zygosity"],
            "vaf": None if var["vaf"] is None else round(var["vaf"], 3), "depth": var["dp"], "cftr2": det or "not listed",
            "clinvar": clinvar or "-", "consequence": (csq.get("Consequence") or "").replace("_", " ").replace("&", ", "),
            "filter": var["filter"], "pos": var["pos"], "ref": var["ref"], "alt": var["alt"], "tier": tier,
            "gt": var.get("gt", ""), "ps": var.get("ps", "")}


def hg_name_safe(tx, fasta, var):
    try:
        return hgvs_name(tx, fasta, var["pos"], var["ref"].upper(), var["alt"].upper())
    except (ValueError, IndexError):
        return f"g.{var['pos']}{var['ref']}>{var['alt']}"


def conclusions(rows, polyt, cnv, gaps, cftr2_given, primers_clipped, platform):
    """Plain-language notes; they describe findings and limits and never give a diagnosis."""
    notes = []
    rep = [r for r in rows if r["tier"] == "report"]
    strong = [r for r in rep if r["cftr2"] in (CF_CAUSING, VCC) or (r["cftr2"] == "not listed")]
    het = [r for r in strong if r["zygosity"] == "het"]
    unphased = False
    for i, a in enumerate(het):
        for b in het[i + 1:]:
            if a.get("ps") and a.get("ps") == b.get("ps"):
                where = "the same chromosome (cis)" if a["gt"] == b["gt"] else "opposite chromosomes (trans)"
                notes.append(f"{a['hgvs_c']} and {b['hgvs_c']}: reads place them on {where}. This is read-based "
                             "phasing; confirm with parental samples or another method.")
            else:
                unphased = True
    if unphased:
        notes.append("Two or more heterozygous variants not joined by the same reads: whether they are on the same "
                     "chromosome (cis) or on opposite chromosomes (trans) is not known from this test. Parental "
                     "samples are needed.")
    for r in strong:
        if r["zygosity"] == "hom":
            notes.append(f"{r['hgvs_c']} looks homozygous: a large deletion on the other chromosome, or a primer-site "
                         "variant causing allele drop-out, can look the same. Confirm with parents or another method.")
    if any(r["hgvs_c"] == "c.350G>A" for r in rows):
        notes.append("R117H (c.350G>A) present: its effect depends on the poly-T tract on the same chromosome. "
                     "Amplicons cannot join exon 4 and intron 9, so cis/trans needs parental samples.")
    if polyt.get("notes"):
        notes.extend(polyt["notes"])
    if not rep and not (cnv or {}).get("calls"):
        where = "in the regions adequately covered" if gaps else "in the targeted regions"
        notes.append(f"No CF-causing, varying-consequence or uninterpreted CFTR2 variant was detected {where}. "
                     "This lowers but does not remove the chance of carrying a CFTR variant.")
    if gaps:
        notes.append(f"{len(gaps)} target(s) below full coverage (listed in red): a variant there could be missed.")
    if not cftr2_given:
        notes.append("CFTR2 list not supplied: variants could not be classified by CFTR2.")
    if not primers_clipped:
        notes.append("Primer positions were not supplied, so primer bases were not removed: a real variant under a "
                     "primer can be hidden or under-counted.")
    if platform == "ont":
        notes.append("Nanopore data: indels in homopolymers and repeat counts are less reliable.")
    return notes


# ---------------------------------------------------------------- HTML
CSS = """
:root{--ink:#1c2733;--muted:#5d6b78;--line:#d9e0e6;--accent:#0b6e75;--soft:#eef5f6;--bad:#b3261e;--badbg:#fbeceb;--warn:#8a5a00;--warnbg:#fff5e0}
*{box-sizing:border-box}body{font:12px/1.35 system-ui,Segoe UI,Arial,sans-serif;color:var(--ink);background:#fff;margin:0 auto;max-width:880px;padding:16px}
h1{font-size:18px;margin:0;color:var(--accent)}h2{font-size:12.5px;margin:12px 0 4px;color:var(--accent);border-bottom:1px solid var(--line);padding-bottom:2px}
.banner{background:var(--badbg);color:var(--bad);font-weight:600;padding:5px 8px;border-radius:4px;margin:6px 0}
.meta{color:var(--muted);margin-top:2px}.wrap{overflow-x:auto}table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:3px 5px;border-bottom:1px solid var(--line);vertical-align:top}
th{background:var(--soft);font-weight:600}.CF{color:var(--bad);font-weight:700}.VCC,.NI{color:var(--warn);font-weight:700}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(84px,1fr));gap:3px}.cell{border:1px solid var(--line);border-radius:3px;padding:2px 4px;font-size:10.5px}
.cell.gap{background:var(--badbg);border-color:var(--bad);color:var(--bad);font-weight:600}.cell b{display:block}
ul{margin:3px 0 0;padding-left:16px}.note{background:var(--warnbg);padding:5px 8px;border-radius:4px;margin:0}details{margin-top:8px;color:var(--muted)}.foot{color:var(--muted);font-size:10px;margin-top:10px}
@media print{@page{size:A4;margin:9mm}body{padding:0;max-width:none}details:not([open]){display:none}}
"""


def _e(x):
    return html.escape("" if x is None else str(x))


def _cls(det):
    return {CF_CAUSING: "CF", VCC: "VCC", NO_INTERP: "NI"}.get(det, "")


def _var_table(rows):
    if not rows:
        return "<p>None.</p>"
    out = ["<div class=wrap><table><tr><th>cDNA (NM_000492.4)</th><th>Protein</th><th>Legacy</th><th>Location</th>"
           "<th>Zygosity</th><th>VAF</th><th>Depth</th><th>CFTR2</th><th>ClinVar</th></tr>"]
    for r in rows:
        vaf = "" if r["vaf"] is None else f"{r['vaf'] * 100:.0f}%"
        flag = "" if r["filter"] in ("PASS", ".") else f" <small>({_e(r['filter'])})</small>"
        out.append(f"<tr><td>{_e(r['hgvs_c'])}{flag}</td><td>{_e(r['protein'])}</td><td>{_e(r['legacy'])}</td><td>{_e(r['region'])}</td>"
                   f"<td>{_e(r['zygosity'])}</td><td>{vaf}</td><td>{_e(r['depth'])}</td>"
                   f"<td class='{_cls(r['cftr2'])}'>{_e(r['cftr2'])}</td><td>{_e(r['clinvar'])}</td></tr>")
    out.append("</table></div>")
    return "".join(out)


def _label(target):
    m = re.match(r"CFTR_exon(\d+)", target)
    if m:
        n = int(m.group(1))
        leg = hg.LEGACY_EXON[n]
        return f"Exon {n}" + (f" ({leg})" if leg != str(n) else "")
    return target.replace("CFTR_", "")


def render(d):
    rep = [r for r in d["variants"] if r["tier"] == "report"]
    other = [r for r in d["variants"] if r["tier"] != "report"]
    p = d["polyt"]
    pt = (f"<b>{_e(p.get('genotype'))}</b> ({_e(p.get('hgvs'))}); {_e(p.get('spanning_reads'))} spanning reads"
          if p.get("genotype") else _e(p.get("status") or "not assessed"))
    if p.get("confidence"):
        pt += f"; confidence: {_e(p['confidence'])}"
    if d.get("polyt_cftr2"):
        pt += "<br>CFTR2: " + "; ".join(f"{_e(k)} = {_e(v)}" for k, v in d["polyt_cftr2"].items())
    c = d["cnv"] or {}
    if c.get("status") == "assessed":
        cn = ("<ul>" + "".join(
            f"<li><b>{_e(x['zygosity'])} {_e(x['kind'])}, {_e(x['exons'])}</b> ({_e(x['legacy_exons'])}), ratio {_e(x['mean_ratio'])}"
            f"{' = ' + _e(x['known']) if x['known'] else ''}. {_e(x['confirm'])}.</li>" for x in c["calls"]) + "</ul>"
              if c["calls"] else f"No exon deletion or duplication suggested ({_e(c['exons_assessed'])} of {_e(c['exons_total'])} exons assessed).")
    else:
        cn = _e(c.get("status") or "not assessed")
    cov = "".join(f"<div class='cell{' gap' if t['gap'] else ''}'><b>{_e(_label(t['target']))}</b>{t['pct']:g}% &middot; {t['mean']:.0f}x</div>"
                  for t in d["coverage"])
    notes = "".join(f"<li>{_e(n)}</li>" for n in d["notes"])
    return f"""<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport content="width=device-width,initial-scale=1">
<title>AFLA CFTR report</title><style>{CSS}</style></head><body>
<h1>AFLA CFTR report</h1>
<div class=meta>Sample <b>{_e(d['sample'])}</b> &middot; {_e(d['platform_label'])} &middot; {_e(d['date'])} &middot; GRCh38, NM_000492.4 &middot; workflow {_e(d['version'])}</div>
<div class=banner>Education and research use only. Not for clinical use.</div>
<h2>Variants of interest</h2>{_var_table(rep)}
<h2>Poly-T / TG tract (c.1210-34TG[n]T[n], legacy IVS8)</h2><p style="margin:0">{pt}</p>
<h2>Exon deletions and duplications</h2><div>{cn}</div>
<h2>Interpretation notes</h2><div class=note><ul>{notes}</ul></div>
<h2>Coverage: share of bases at {_e(d['depth_used'])}x or more, mean depth</h2><div class=grid>{cov}</div>
<details><summary>Other variants ({len(other)}; benign, non CF-causing or not listed)</summary>{_var_table(other)}</details>
<div class=foot>Sources: CFTR2 ({_e(d['cftr2_date'] or 'not supplied')}), ClinVar via Ensembl VEP when annotation is on. Legacy names follow CFTR2. Variant calls: {_e(d['caller'])}. Analysis limited to the targeted regions.</div>
</body></html>
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--polyt", required=True)
    ap.add_argument("--cnv", required=True, help="JSON from afla_cftr_cnv.py (all samples)")
    ap.add_argument("--thresholds", required=True, help="mosdepth <prefix>.thresholds.bed.gz of the named targets")
    ap.add_argument("--regions", help="mosdepth <prefix>.regions.bed.gz of the named targets")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--cftr2", help="CFTR2 variant list (.xlsx)")
    ap.add_argument("--exons", default=str(hg.DEFAULT_EXONS))
    ap.add_argument("--platform", choices=["illumina", "ont"], default="illumina")
    ap.add_argument("--min-depth", type=int, default=20)
    ap.add_argument("--caller", default="")
    ap.add_argument("--version", default="")
    ap.add_argument("--date", default="")
    ap.add_argument("--primers-clipped", action="store_true")
    ap.add_argument("--out", required=True)
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    tx = hg.Transcript.load(a.exons)
    fasta = hg.Fasta(a.fasta)
    has_cftr2 = bool(a.cftr2) and not Path(a.cftr2).name.startswith("NO_")
    by_pos, by_legacy, cftr2_date = load_cftr2(a.cftr2, fasta) if has_cftr2 else ({}, {}, "")
    variants = [classify(v, tx, fasta, by_pos) for v in read_vcf(a.vcf, fasta)]
    order = {CF_CAUSING: 0, VCC: 1, NO_INTERP: 2}
    variants.sort(key=lambda r: (r["tier"] != "report", order.get(r["cftr2"], 3), r["pos"]))
    polyt = json.loads(Path(a.polyt).read_text())
    cnv = json.loads(Path(a.cnv).read_text()).get(a.sample)
    polyt_cftr2 = {}
    for al in polyt.get("alleles", []):  # most specific CFTR2 entry first: "5T;TG12", then "5T"
        for name in (f"{al['t']}T;TG{al['tg']}", f"{al['t']}T"):
            if name in by_legacy:
                polyt_cftr2[name] = by_legacy[name]["determination"]
                break
    coverage, depth = read_coverage(a.thresholds, a.regions, a.min_depth)
    gaps = [t for t in coverage if t["gap"]]
    data = {"sample": a.sample, "platform": a.platform, "platform_label": "Illumina" if a.platform == "illumina" else "Nanopore",
            "date": a.date, "version": a.version, "caller": a.caller, "cftr2_date": cftr2_date, "variants": variants, "polyt": polyt,
            "polyt_cftr2": polyt_cftr2, "cnv": cnv, "coverage": coverage, "depth_used": depth,
            "notes": conclusions(variants, polyt, cnv, gaps, has_cftr2, a.primers_clipped, a.platform)}
    Path(a.out).write_text(render(data))
    if a.json:
        Path(a.json).write_text(json.dumps(data, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
