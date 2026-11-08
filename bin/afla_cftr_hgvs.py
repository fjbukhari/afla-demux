#!/usr/bin/env python3
"""CFTR coordinates for AFLA (education and research use only, not for clinical use).

Converts between CFTR cDNA names on the MANE transcript NM_000492.4 (ENST00000003084) and GRCh38 positions,
and gives exon/intron names in both the current and the legacy numbering (legacy is still used in CF reports:
"IVS8 poly-T" is intron 9 today, "exon 10" of F508del is exon 11).

  afla_cftr_hgvs.py exons --ncbi data_report.jsonl --out cftr/resources/CFTR_NM_000492.4_exons.tsv
  afla_cftr_hgvs.py c2g c.1521_1523del c.1210-34 c.3718-2477C>T
  afla_cftr_hgvs.py g2c 117559591 117548795
  afla_cftr_hgvs.py check-acmg "ACMG CFTR 23 Variants list.xlsx" --fasta hg38.fa
  afla_cftr_hgvs.py targets --pad 25 --out cftr_targets.bed      (exons, deep intronic sites, poly-T/TG)

The exon table in cftr/resources/ was made with the first command from the NCBI Datasets gene report for CFTR
(GeneID 1080, NCBI annotation release 109, GRCh38.p13): public-domain data. CFTR is on the + strand of chr7.
Standard library only (openpyxl just for check-acmg on .xlsx).
"""
import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESOURCES = HERE.parent / "cftr" / "resources"
DEFAULT_EXONS = RESOURCES / "CFTR_NM_000492.4_exons.tsv"
DEFAULT_DEEP = RESOURCES / "CFTR_deep_intronic.tsv"
POLYT_C = ("1210-34", "1210-6")  # (TG)n(T)n tract, legacy IVS8
CHROM = "chr7"
TRANSCRIPT = "NM_000492.4"

# current exon number -> legacy exon name (CFTR2 / Cystic Fibrosis Mutation Database numbering)
LEGACY_EXON = {1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "6a", 7: "6b", 8: "7", 9: "8", 10: "9", 11: "10",
               12: "11", 13: "12", 14: "13", 15: "14a", 16: "14b", 17: "15", 18: "16", 19: "17a", 20: "17b",
               21: "18", 22: "19", 23: "20", 24: "21", 25: "22", 26: "23", 27: "24"}
# legacy intron names follow the legacy exon before them: intron 9 (between exons 9 and 10) is legacy IVS8
LEGACY_INTRON = {n: LEGACY_EXON[n] for n in range(1, 27)}


class Transcript:
    """Exons of a + strand transcript with its CDS start/end in transcript coordinates."""

    def __init__(self, exons, cds_start_t, cds_end_t):
        self.exons = sorted(exons)  # [(number, g_start, g_end)], 1-based inclusive
        self.cds_start_t, self.cds_end_t = cds_start_t, cds_end_t
        self.t_start = []  # transcript coordinate of each exon's first base
        t = 1
        for _, a, b in self.exons:
            self.t_start.append(t)
            t += b - a + 1
        self.length = t - 1

    @classmethod
    def load(cls, path=DEFAULT_EXONS):
        exons, cds = [], {}
        for line in open(path).read().splitlines(keepends=True):
            if line.startswith("#"):
                m = re.match(r"#\s*(cds_start_t|cds_end_t)\s*=\s*(\d+)", line)
                if m:
                    cds[m.group(1)] = int(m.group(2))
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "exon":
                continue
            exons.append((int(f[0]), int(f[2]), int(f[3])))
        return cls(exons, cds["cds_start_t"], cds["cds_end_t"])

    # transcript coordinate <-> c. position (no offset)
    def _t_to_c(self, t):
        if t < self.cds_start_t:
            return f"-{self.cds_start_t - t}"
        if t > self.cds_end_t:
            return f"*{t - self.cds_end_t}"
        return str(t - self.cds_start_t + 1)

    def _c_to_t(self, c):
        if c.startswith("-"):
            return self.cds_start_t - int(c[1:])
        if c.startswith("*"):
            return self.cds_end_t + int(c[1:])
        return self.cds_start_t + int(c) - 1

    def _t_to_g(self, t):
        for (n, a, b), ts in zip(self.exons, self.t_start):
            if ts <= t <= ts + (b - a):
                return a + (t - ts)
        raise ValueError(f"transcript position {t} is outside {TRANSCRIPT} (1-{self.length})")

    def c_pos_to_g(self, cpos):
        """'1521', '1210-34', '-8', '*12', '3718-2477', '489+1' -> GRCh38 position."""
        m = re.fullmatch(r"([-*]?\d+)([+-]\d+)?", cpos.strip())
        if not m:
            raise ValueError(f"cannot read the c. position {cpos!r}")
        g = self._t_to_g(self._c_to_t(m.group(1)))
        if m.group(2):
            off = int(m.group(2))
            # an intronic offset must start from an exon edge in the right direction
            g += off
        return g

    def g_to_c(self, g):
        """GRCh38 position -> c. position string (e.g. '1210-34', '489+1', '*1251')."""
        for i, ((n, a, b), ts) in enumerate(zip(self.exons, self.t_start)):
            if a <= g <= b:
                return self._t_to_c(ts + g - a)
            if i + 1 < len(self.exons) and b < g < self.exons[i + 1][1]:
                na = self.exons[i + 1][1]
                up, down = g - b, na - g  # distance from the exon before / after
                if up <= down:
                    return f"{self._t_to_c(ts + b - a)}+{up}"
                return f"{self._t_to_c(self.t_start[i + 1])}-{down}"
        if g < self.exons[0][1]:
            return f"-{self.cds_start_t - 1 + self.exons[0][1] - g}"
        return f"*{self.length - self.cds_end_t + g - self.exons[-1][2]}"

    def region(self, g):
        """('exon', 11, '10') or ('intron', 9, '8') with current and legacy numbers, or ('outside', None, None)."""
        for i, (n, a, b) in enumerate(self.exons):
            if a <= g <= b:
                return "exon", n, LEGACY_EXON[n]
            if i + 1 < len(self.exons) and b < g < self.exons[i + 1][1]:
                return "intron", n, LEGACY_INTRON[n]
        return "outside", None, None

    def region_label(self, g):
        kind, n, leg = self.region(g)
        if kind == "exon":
            return f"exon {n}" + (f" (legacy {leg})" if leg != str(n) else "")
        if kind == "intron":
            return f"intron {n} (legacy IVS{leg})"
        return "outside the transcript"


HGVS_RE = re.compile(r"^c\.(?P<a>[-*]?\d+(?:[+-]\d+)?)(?:_(?P<b>[-*]?\d+(?:[+-]\d+)?))?(?P<rest>.*)$")


def c_to_g(tx, name):
    """'c.1521_1523del' -> (start, end) GRCh38 of the named bases; 'c.350G>A' -> (pos, pos)."""
    m = HGVS_RE.match(name.strip().split(":")[-1])
    if not m:
        raise ValueError(f"cannot read {name!r} (expected c. notation, e.g. c.1521_1523del)")
    a = tx.c_pos_to_g(m.group("a"))
    b = tx.c_pos_to_g(m.group("b")) if m.group("b") else a
    return min(a, b), max(a, b)


def read_tsv(path):
    """Rows of a resource TSV as dicts (lines starting with # are comments; first other line is the header)."""
    header, rows = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
            else:
                rows.append(dict(zip(header, f + [""] * (len(header) - len(f)))))
    return rows


def targets(tx, pad=25, deep=DEFAULT_DEEP, cds_only=False):
    """Named BED regions (0-based start): every exon +/- pad (coding part only with cds_only), each deep
    intronic site +/- pad, and the poly-T/TG tract. Overlapping regions are kept separate (they are named)."""
    out = []
    cds_a, cds_b = tx.c_pos_to_g("1"), tx.c_pos_to_g(f"{tx.cds_end_t - tx.cds_start_t + 1}")
    for n, a, b in tx.exons:
        if cds_only:
            a, b = max(a, cds_a), min(b, cds_b)
            if a > b:
                continue
        leg = LEGACY_EXON[n]
        out.append((CHROM, a - pad - 1, b + pad, f"CFTR_exon{n:02d}" + (f"_legacy{leg}" if leg != str(n) else "")))
    if deep and Path(deep).exists():
        for r in read_tsv(deep):
            a, b = c_to_g(tx, r["cdna"])
            out.append((CHROM, a - pad - 1, b + pad, "CFTR_" + re.sub(r"[^A-Za-z0-9+._-]", "", r["cdna"].replace(">", "_"))))
    a, b = tx.c_pos_to_g(POLYT_C[0]), tx.c_pos_to_g(POLYT_C[1])
    out.append((CHROM, a - pad - 1, b + pad, "CFTR_polyTG_polyT_IVS8"))
    return sorted(out, key=lambda r: (r[1], r[2]))


def make_exon_table(ncbi_jsonl, out, transcript=TRANSCRIPT):
    for line in open(ncbi_jsonl):
        d = json.loads(line)
        for t in d.get("transcripts", []):
            if t["accessionVersion"] != transcript:
                continue
            cds = t["cds"]["range"][0]
            exons = sorted((int(e["order"]), int(e["begin"]), int(e["end"])) for e in t["exons"]["range"])
            if t["exons"]["range"][0].get("orientation", "plus") != "plus":
                sys.exit("only + strand transcripts are handled")
            ann = d.get("annotations", [{}])[0]
            with open(out, "w") as f:
                f.write(f"# {transcript} ({t.get('ensemblTranscript', '')}) exons on GRCh38 {CHROM}, 1-based inclusive\n")
                f.write(f"# source: NCBI Datasets gene report, {ann.get('releaseName', '')} "
                        f"({ann.get('releaseDate', '')}); public-domain data\n")
                f.write(f"# cds_start_t = {cds['begin']}\n# cds_end_t = {cds['end']}\n")
                f.write("exon\tlegacy_exon\tstart\tend\tlength\n")
                for n, a, b in exons:
                    f.write(f"{n}\t{LEGACY_EXON.get(n, '')}\t{a}\t{b}\t{b - a + 1}\n")
            return len(exons)
    sys.exit(f"{transcript} not found in {ncbi_jsonl}")


# ---------------------------------------------------------------- reference sequence (no samtools needed)
class Fasta:
    def __init__(self, path):
        self.path = path
        self.idx = {}
        for line in open(str(path) + ".fai"):
            name, length, offset, lb, lw = line.split("\t")[:5]
            self.idx[name] = (int(length), int(offset), int(lb), int(lw))

    def fetch(self, chrom, start, end):
        """1-based inclusive."""
        if chrom not in self.idx and chrom.replace("chr", "") in self.idx:
            chrom = chrom.replace("chr", "")
        length, offset, lb, lw = self.idx[chrom]
        start, end = max(1, start), min(length, end)
        out = []
        with open(self.path, "rb") as f:
            p = start - 1
            while p < end:
                line_i, col = divmod(p, lb)
                f.seek(offset + line_i * lw + col)
                n = min(lb - col, end - p)
                out.append(f.read(n).decode())
                p += n
        return "".join(out).upper()


def left_align(fasta, pos, ref, alt):
    """Left-align and trim a VCF-style indel (as bcftools norm does)."""
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1] and (len(ref) > 1 and len(alt) > 1 or pos > 1):
            # shift left: drop the common last base, add the base before at the front
            if len(ref) == 1 or len(alt) == 1:
                prev = fasta.fetch(CHROM, pos - 1, pos - 1)
                ref, alt, pos = prev + ref[:-1], prev + alt[:-1], pos - 1
            else:
                ref, alt = ref[:-1], alt[:-1]
            changed = True
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0] and not changed:
            ref, alt, pos = ref[1:], alt[1:], pos + 1
        if not changed:
            return pos, ref, alt


HGVS_CHANGE_RE = re.compile(r"^(?:(?P<snv_ref>[ACGT])>(?P<snv_alt>[ACGT])|(?P<kind>delins|del|dup|ins)(?P<seq>[ACGT]*))$")


def hgvs_to_vcf(tx, fasta, name):
    """'c.1521_1523del' -> ('chr7', 117559590, 'ATCT', 'A'), left-aligned like a normalised VCF.
    Handles substitutions, del, dup, ins and delins with explicit or implied sequence."""
    m = HGVS_RE.match(name.strip().split(":")[-1])
    if not m:
        raise ValueError(f"cannot read {name!r}")
    a, b = c_to_g(tx, name)
    ch = HGVS_CHANGE_RE.match(m.group("rest").strip())
    if not ch:
        raise ValueError(f"{name}: only substitutions, del, dup, ins and delins are converted")
    if ch.group("snv_ref"):
        if fasta.fetch(CHROM, a, a) != ch.group("snv_ref"):
            raise ValueError(f"{name}: reference base at {CHROM}:{a} is {fasta.fetch(CHROM, a, a)}")
        return CHROM, a, ch.group("snv_ref"), ch.group("snv_alt")
    kind, seq = ch.group("kind"), ch.group("seq")
    if kind == "del":
        prev = fasta.fetch(CHROM, a - 1, a - 1)
        pos, ref, alt = a - 1, prev + fasta.fetch(CHROM, a, b), prev
    elif kind == "dup":
        dup = fasta.fetch(CHROM, a, b)
        pos, ref, alt = b, fasta.fetch(CHROM, b, b), fasta.fetch(CHROM, b, b) + dup
    elif kind == "ins":
        if not seq:
            raise ValueError(f"{name}: inserted sequence not given")
        pos, ref, alt = a, fasta.fetch(CHROM, a, a), fasta.fetch(CHROM, a, a) + seq
    else:  # delins
        if not seq:
            raise ValueError(f"{name}: inserted sequence not given")
        pos, ref, alt = a, fasta.fetch(CHROM, a, b), seq
        if len(ref) == len(alt) or ref[0] != alt[0] and ref[-1] != alt[-1]:
            return CHROM, pos, ref, alt
    pos, ref, alt = left_align(fasta, pos, ref, alt)
    return CHROM, pos, ref, alt


def check_acmg(path, tx, fasta=None):
    """Compare our c. -> GRCh38 conversion with the VCF positions given in the ACMG 2023 CFTR table."""
    rows = []
    if str(path).lower().endswith(".xlsx"):
        import openpyxl
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        for ws in wb.worksheets:
            for r in ws.iter_rows(values_only=True):
                cells = [str(c) for c in r if c is not None]
                cn = next((c for c in cells if c.startswith("c.")), None)
                vcf = next((c for c in cells if re.fullmatch(r"7:\d+:[ACGT]+:[ACGT]+", c)), None)
                if cn and vcf:
                    rows.append((cn, vcf))
            if rows:
                break
    else:
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0].startswith("c."):
                rows.append((f[0], f[1]))
    ok, bad = 0, []
    for cn, vcf in rows:
        _, pos, ref, alt = vcf.split(":")
        pos = int(pos)
        try:
            a, b = c_to_g(tx, cn)
        except ValueError as e:
            bad.append((cn, vcf, str(e)))
            continue
        if len(ref) == len(alt) == 1:
            good = a == pos
            if good and fasta:
                good = fasta.fetch(CHROM, pos, pos) == ref
        else:
            # VCF: anchor base before the change, left-aligned; HGVS: shifted to the 3' side. The change must lie
            # within the stretch that the left shift can move across.
            # (delins and ins have no extra anchor: the VCF position is then the first changed or flanking base)
            span = max(len(ref), len(alt)) + 1
            good = pos <= a and a - pos <= 60 + span
            if good and fasta:
                good = fasta.fetch(CHROM, pos, pos + len(ref) - 1) == ref
        if good and fasta:
            try:
                key = hgvs_to_vcf(tx, fasta, cn)
            except ValueError as e:
                key = str(e)
            if key != (CHROM, pos, ref, alt):
                good = False
                b = f"{b}; VCF from name {key}"
        if good:
            ok += 1
        else:
            bad.append((cn, vcf, f"ours {a}-{b}"))
    return ok, bad


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("exons")
    p.add_argument("--ncbi", required=True)
    p.add_argument("--out", default=str(DEFAULT_EXONS))
    p = sub.add_parser("c2g")
    p.add_argument("names", nargs="+")
    p.add_argument("--exons", default=str(DEFAULT_EXONS))
    p = sub.add_parser("g2c")
    p.add_argument("positions", nargs="+", type=int)
    p.add_argument("--exons", default=str(DEFAULT_EXONS))
    p = sub.add_parser("targets")
    p.add_argument("--pad", type=int, default=25)
    p.add_argument("--deep", default=str(DEFAULT_DEEP))
    p.add_argument("--cds-only", action="store_true", help="coding part of the exons only (no UTRs)")
    p.add_argument("--out", required=True)
    p.add_argument("--exons", default=str(DEFAULT_EXONS))
    p = sub.add_parser("check-acmg")
    p.add_argument("table")
    p.add_argument("--fasta")
    p.add_argument("--exons", default=str(DEFAULT_EXONS))
    a = ap.parse_args(argv)
    if a.cmd == "exons":
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        print(f"{make_exon_table(a.ncbi, a.out)} exons written to {a.out}")
        return 0
    tx = Transcript.load(a.exons)
    if a.cmd == "c2g":
        for n in a.names:
            s, e = c_to_g(tx, n)
            print(f"{n}\t{CHROM}:{s}" + (f"-{e}" if e != s else "") + f"\t{tx.region_label(s)}")
    elif a.cmd == "g2c":
        for g in a.positions:
            print(f"{CHROM}:{g}\tc.{tx.g_to_c(g)}\t{tx.region_label(g)}")
    elif a.cmd == "targets":
        rows = targets(tx, a.pad, a.deep, a.cds_only)
        with open(a.out, "w") as f:
            for r in rows:
                f.write("\t".join(map(str, r)) + "\n")
        print(f"{len(rows)} regions, {sum(r[2] - r[1] for r in rows):,} bases, written to {a.out}")
    elif a.cmd == "check-acmg":
        ok, bad = check_acmg(a.table, tx, Fasta(a.fasta) if a.fasta else None)
        print(f"{ok} of {ok + len(bad)} ACMG table variants agree with our conversion")
        for b in bad:
            print("  differs:", *b)
        return 1 if bad else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
