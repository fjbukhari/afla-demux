#!/usr/bin/env python3
"""AFLA demultiplexer for barcoded Nanopore runs (education and research use only, not for clinical use).

Sorts reads into one FASTQ per sample by finding the sample's barcode sequence at the start and/or the end of
each read. It reads the barcodes from your sample sheet, so it works for any kit or home-made barcodes whose reads
look like  [adapter][START sequence][insert][END sequence][adapter]. For standard ONT barcoding kits, MinKNOW/Dorado
already does this; use this tool for custom or third-party barcodes.

  # what is in this run? (no sample sheet needed; writes a draft sample sheet you can edit)
  afla_demux.py --discover --input fastq_pass/ --out demux/

  # sort the reads
  afla_demux.py --sample-sheet sample_sheet.csv --input fastq_pass/ --out demux/ --threads 8

  # live: keep sorting new files while MinKNOW is still writing them (Ctrl+C to stop)
  afla_demux.py --sample-sheet sample_sheet.csv --input /path/to/run/fastq_pass --out demux/ --watch

  # only check that the sample sheet is understood
  afla_demux.py --sample-sheet sample_sheet.xlsx --check-sheet

Inputs: FASTQ files (.fastq, .fq, optionally .gz), folders of them (searched recursively), or .zip archives of
them (read in place, nothing is extracted). Sample sheets: a CSV or XLSX with columns sample,start[,end] (headers
such as "Start sequence bases", "End sequence bases" and "SAMPLE_NAME" are also understood; only rows with a sample
name are used). If END is missing it is taken as the reverse complement of START (symmetric kits).

Outputs in --out: <sample>/<sample>.fastq.gz, unclassified/unclassified.fastq.gz, demux_summary.tsv,
demux_report.html and (with --read-table) demux_reads.tsv.gz.

Needs edlib (pip install edlib).
"""
import argparse
import csv
import gzip
import html
import io
import json
import math
import os
import re
import signal
import sys
import time
import zipfile
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

try:
    import edlib
except ImportError:  # checked in main() so --help still works
    edlib = None

VERSION = "0.1.0"
FASTQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
SAFE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Generic, public adapter/flank sequences used by --discover to recognise a layout. They are the standard ONT
# native-barcode flanks and Illumina/Ion Torrent adapter parts, not kit barcodes (barcodes are never stored here).
LAYOUTS = [
    # name, left flank, right flank, typical barcode length
    ("ONT native-barcode style", "AAGGTTAA", "CAGCACCT", 24),
    ("Illumina i7 style", "CAAGCAGAAGACGGCATACGAGAT", "GTCTCGTGGGCTCGG", 8),
    ("Illumina i5 style", "AATGATACGGCGACCACCGAGATCTACAC", "TCGTCGGCAGCGTC", 8),
    ("Ion Torrent style", "CCATCTCATCCCTGCGTGTCTCCGACTCAG", "GATAGATGTGTATAAGAGACAG", 11),
]
ELEMENTS = [
    ("ONT ligation adapter", "CCTGTACTTCGTTCAGTTACGTATTGCT"),
    ("Nextera mosaic end (ME)", "AGATGTGTATAAGAGACAG"),
]


def rc(s):
    return s.translate(COMP)[::-1]


def die(msg):
    sys.exit(f"afla_demux: {msg}")


# ---------------------------------------------------------------- sample sheet
class Barcode:
    __slots__ = ("sample", "name", "well", "start", "end", "kmax_start", "kmax_end", "core_start", "core_end",
                 "kcore_start", "kcore_end")

    def __init__(self, sample, name, well, start, end):
        self.sample, self.name, self.well = sample, name, well
        self.start, self.end = start, end
        self.kmax_start = self.kmax_end = self.kcore_start = self.kcore_end = 0
        self.core_start = (0, len(start))
        self.core_end = (0, len(end))


def _clean_seq(s):
    s = re.sub(r"\s+", "", str(s or "")).upper()
    if s and not re.fullmatch(r"[ACGTN]+", s):
        return None
    return s


def _norm_header(h):
    return re.sub(r"[^a-z0-9]", "", str(h or "").lower())


# header aliases (normalised) -> field
HEADER_ALIASES = {
    "startsequencebases": "start", "start": "start", "startseq": "start", "barcodesequence": "start",
    "sequence": "start", "barcodeseq": "start",
    "endsequencebases": "end", "end": "end", "endseq": "end",
    "samplename": "sample", "sample": "sample", "alias": "sample", "sampleid": "sample",
    "startsequencename": "name", "barcode": "name", "barcodename": "name",
    "index": "well", "indice": "well", "well": "well",
}


def _rows_from_table(rows, source):
    """rows: list of lists (first non-empty row containing a start/sample column is the header)."""
    header_i = None
    for i, r in enumerate(rows[:20]):
        fields = {HEADER_ALIASES.get(_norm_header(c)) for c in r}
        if "start" in fields:
            header_i = i
            break
    if header_i is None:
        die(f"{source}: no column with the barcode sequence (expected 'Start sequence bases', 'start' or 'sequence')")
    col = {}
    for j, c in enumerate(rows[header_i]):
        f = HEADER_ALIASES.get(_norm_header(c))
        if f and f not in col:
            col[f] = j
    out = []
    for r in rows[header_i + 1:]:
        get = lambda f: (str(r[col[f]]).strip() if f in col and col[f] < len(r) and r[col[f]] is not None else "")
        rec = {f: get(f) for f in ("start", "end", "sample", "name", "well")}
        if not any(rec.values()):
            continue
        out.append(rec)
    return out


def _read_xlsx(path):
    from afla_cftr import read_xlsx  # standard-library reader, so the container needs no extra package
    sheets = list(read_xlsx(path).values())
    # main sheet = first sheet that has a sample-name column; lookup table = a sheet whose first column is wells
    main = None
    for rows in sheets:
        for r in rows[:20]:
            fields = {HEADER_ALIASES.get(_norm_header(c)) for c in r}
            if "start" in fields and "sample" in fields:
                main = rows
                break
        if main:
            break
    if main is None:
        main = next((rows for rows in sheets if any(
            "start" in {HEADER_ALIASES.get(_norm_header(c)) for c in r} for r in rows[:20])), None)
    if main is None:
        die(f"{path}: no sheet with a barcode sequence column")
    lookup = {}
    for rows in sheets:
        for r in rows:
            if len(r) >= 3 and isinstance(r[0], str) and re.fullmatch(r"[A-Ha-h](1[0-2]|[1-9])", r[0].strip()):
                s, e = _clean_seq(r[1]), _clean_seq(r[2])
                if s:
                    lookup[r[0].strip().upper()] = (s, e or "")
    recs = _rows_from_table(main, path)
    for rec in recs:  # formulas whose results were not saved: resolve from the well lookup table
        if (not _clean_seq(rec["start"]) or str(rec["start"]).startswith("=")) and rec["well"].upper() in lookup:
            rec["start"], rec["end"] = lookup[rec["well"].upper()]
    return recs


def _read_csv(path):
    raw = Path(path).read_bytes()
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    first = next((ln for ln in text.splitlines() if ln.strip()), "")
    delim = max([",", ";", "\t"], key=first.count)
    return _rows_from_table(list(csv.reader(io.StringIO(text), delimiter=delim)), path)


def read_sample_sheet(path):
    """Returns (barcodes, notes). Only rows with a sample name are used, unless no row has one."""
    p = Path(path)
    if not p.exists():
        die(f"sample sheet not found: {path}")
    recs = _read_xlsx(p) if p.suffix.lower() in (".xlsx", ".xlsm") else _read_csv(p)
    notes = []
    named = [r for r in recs if r["sample"]]
    if not named:
        named = [r for r in recs if _clean_seq(r["start"])]
        for r in named:
            r["sample"] = re.sub(r"[^A-Za-z0-9_.-]", "_", r["name"] or r["well"] or "sample")
        if named:
            notes.append("no sample names in the sheet: every barcode row is used, named after its barcode")
    bcs, seen_names, seen_seq = [], {}, {}
    for r in named:
        start = _clean_seq(r["start"])
        if not start:
            die(f"sample {r['sample']!r}: missing or invalid start sequence ({str(r['start'])[:20]!r}); "
                "if the sheet uses formulas, open it in Excel and save it again, or save as CSV")
        end = _clean_seq(r["end"]) if r["end"] else ""
        if end is None:
            die(f"sample {r['sample']!r}: invalid end sequence")
        if not end:
            end = rc(start)
        elif end != rc(start):
            notes.append(f"{r['sample']}: end sequence is not the reverse complement of the start (asymmetric layout)")
        name = r["sample"].strip()
        if not SAFE_NAME_RE.match(name):
            die(f"sample name {name!r}: use only letters, digits, '-', '_' or '.' (no spaces)")
        if name.lower() in seen_names:
            die(f"sample name {name!r} is used twice (rows {seen_names[name.lower()]} and {r['name'] or r['well']})")
        if start in seen_seq:
            die(f"samples {seen_seq[start]!r} and {name!r} have the same barcode")
        if name.lower() in ("unclassified",):
            die("'unclassified' is reserved; choose another sample name")
        seen_names[name.lower()] = r["name"] or r["well"]
        seen_seq[start] = name
        bcs.append(Barcode(name, r["name"], r["well"], start, end))
    if not bcs:
        die(f"{path}: no usable rows (fill in SAMPLE_NAME for the wells you used)")
    return bcs, notes


def _lcp(a, b):
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def barcode_core(seq, others, min_core=6):
    """(a, b): the part of seq that is the barcode proper, i.e. without the flanks it shares with the other
    barcodes of the sheet or with a known public layout. Errors are counted separately in this part, so that a
    barcode that is not in the sheet is not given to the nearest one that is."""
    # flank = what seq shares with every member of its family (barcodes sharing >= 6 bases at that end); the
    # minimum over the family avoids counting bases that two barcodes share by chance
    fam_p = [x for x in (_lcp(seq, o) for o in others) if x >= 6]
    fam_s = [x for x in (_lcp(seq[::-1], o[::-1]) for o in others) if x >= 6]
    pre = min(fam_p) if fam_p else 0
    suf = min(fam_s) if fam_s else 0
    for _, left, right, _ in LAYOUTS:
        for l_, r_ in ((left, right), (rc(right), rc(left))):
            i = seq.find(l_)
            if i >= 0:
                pre = max(pre, i + len(l_))
            j = seq.rfind(r_)
            if j >= 0:
                suf = max(suf, len(seq) - j)
    if len(seq) - pre - suf < min_core:
        return 0, len(seq)
    return pre, len(seq) - suf


def set_thresholds(bcs, max_err, max_core_err=0.2):
    for b in bcs:
        b.kmax_start = max(1, int(math.floor(len(b.start) * max_err)))
        b.kmax_end = max(1, int(math.floor(len(b.end) * max_err)))
        b.core_start = barcode_core(b.start, [o.start for o in bcs if o is not b])
        if b.end == rc(b.start):  # symmetric: the barcode part of the end mirrors the start
            b.core_end = (len(b.end) - b.core_start[1], len(b.end) - b.core_start[0])
        else:
            b.core_end = barcode_core(b.end, [o.end for o in bcs if o is not b])
        b.kcore_start = max(1, int(math.floor((b.core_start[1] - b.core_start[0]) * max_core_err)))
        b.kcore_end = max(1, int(math.floor((b.core_end[1] - b.core_end[0]) * max_core_err)))


def close_pairs(bcs, limit):
    """Pairs of samples whose barcode sequences differ by fewer than `limit` edits, closest first."""
    out = []
    for i in range(len(bcs)):
        for j in range(i + 1, len(bcs)):
            d = edlib.align(bcs[i].start, bcs[j].start, mode="NW")["editDistance"]
            if d < limit:
                out.append((d, bcs[i].sample, bcs[j].sample))
    return sorted(out)


# ---------------------------------------------------------------- reading FASTQ
def win_path(p):
    """C:\\Users\\me\\run -> /mnt/c/Users/me/run when running in WSL."""
    m = re.match(r"^([A-Za-z]):[\\/](.*)$", str(p))
    if m and Path("/mnt").is_dir():
        return Path(f"/mnt/{m.group(1).lower()}/" + m.group(2).replace("\\", "/"))
    return Path(p)


def list_inputs(paths):
    files = []
    for p in paths:
        p = win_path(p)
        if p.is_dir():
            files += sorted(f for f in p.rglob("*") if f.is_file() and (FASTQ_RE.search(f.name) or
                                                                         f.suffix.lower() == ".zip"))
        elif p.is_file():
            files.append(p)
        else:
            die(f"input not found: {p}")
    if not files:
        die("no FASTQ (.fastq/.fq[.gz]) or .zip files found in the input")
    bad = [f for f in files if f.suffix.lower() in (".bam", ".pod5", ".fast5")]
    if bad:
        die(f"{bad[0]}: only FASTQ is read; in MinKNOW choose FASTQ output, or convert BAM with 'samtools fastq'")
    return files


def _fastq_records(handle, source):
    while True:
        h = handle.readline()
        if not h:
            return
        if not h.strip():
            continue
        s = handle.readline().rstrip("\r\n")
        plus = handle.readline()
        q = handle.readline().rstrip("\r\n")
        if not h.startswith("@") or not plus.startswith("+") or len(s) != len(q):
            raise ValueError(f"{source}: not a valid FASTQ record near {h[:60]!r} (file may be truncated)")
        yield h.rstrip("\r\n"), s, q


def iter_file(path):
    """Yields (header, seq, qual) from a FASTQ, FASTQ.gz or a .zip of them."""
    path = Path(path)
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as z:
            for m in sorted(z.namelist()):
                if not FASTQ_RE.search(m):
                    continue
                with z.open(m) as raw:
                    fh = gzip.GzipFile(fileobj=raw) if m.lower().endswith(".gz") else raw
                    yield from _fastq_records(io.TextIOWrapper(fh, encoding="ascii", errors="replace"), f"{path}:{m}")
        return
    op = gzip.open if path.name.lower().endswith(".gz") else open
    with op(path, "rt", encoding="ascii", errors="replace") as fh:
        yield from _fastq_records(fh, path)


def chunks(it, n):
    buf = []
    for x in it:
        buf.append(x)
        if len(buf) >= n:
            yield buf
            buf = []
    if buf:
        yield buf


# ---------------------------------------------------------------- matching
_W = {}  # worker globals


def _init_worker(bcs, opts):
    signal.signal(signal.SIGINT, signal.SIG_IGN)  # the main process handles Ctrl+C
    _W["bcs"], _W["o"] = bcs, opts
    _W["asym"] = any(b.end != rc(b.start) for b in bcs)


def _best_two(window, seqs_k, margin):
    """Edit distance of every query inside window. Returns (best_i, best_d, second_i, second_d); a second
    barcode further than best_d + margin is not looked for (second_d is then 10**6)."""
    best_i, best_d, second_i, second_d = -1, 10 ** 6, -1, 10 ** 6
    for i, (q, kmax) in enumerate(seqs_k):
        k = min(kmax, best_d + margin) if best_i >= 0 else kmax
        d = edlib.align(q, window, mode="HW", k=k)["editDistance"]
        if d < 0:
            continue
        if d < best_d:
            best_i, best_d, second_i, second_d = i, d, best_i, best_d
        elif d < second_d:
            second_i, second_d = i, d
    return best_i, best_d, second_i, second_d


CIGAR_RE = re.compile(r"(\d+)([=XID])")


def core_edits(cigar, a, b):
    """Edits that fall in query positions a..b-1 along an edlib extended CIGAR (I = extra base in the query,
    D = extra base in the read)."""
    q = n = 0
    for length, op in CIGAR_RE.findall(cigar):
        length = int(length)
        if op == "D":
            if a < q < b:
                n += length
            continue
        lo, hi = max(q, a), min(q + length, b)
        if op in "XI" and hi > lo:
            n += hi - lo
        q += length
    return n


def _locate(query, core, kmax, kcore, window):
    """Where query sits in window and whether its barcode part has at most kcore errors.
    Returns (ok, start, end) in window coordinates."""
    r = edlib.align(query, window, mode="HW", task="path", k=kmax)
    if r["editDistance"] < 0 or not r["locations"]:
        return False, 0, 0
    ls, le = r["locations"][0]
    a, b = core
    if (a, b) != (0, len(query)) and core_edits(r["cigar"], a, b) > kcore:
        return False, ls, le
    return True, ls, le


def classify(seq):
    """Returns dict: idx (or -1), reason, ends ('both'/'start'/'end'), d_start, d_end, trim (a, b).

    Each end: the closest barcode, and how many edits closer it is than the runner-up (its margin, capped at
    --min-margin). A read is assigned when the margins of the ends that agree add up to --min-margin, and the
    barcode part itself (without the shared flanks) has few enough errors at every end used."""
    bcs, o = _W["bcs"], _W["o"]
    W, M = o["window"], o["min_margin"]
    head, tail = seq[:W], seq[-W:]
    tail_off = max(0, len(seq) - W)
    si, sd, s2i, s2 = _best_two(head, [(b.start, b.kmax_start) for b in bcs], M)
    ei, ed, e2i, e2 = _best_two(tail, [(b.end, b.kmax_end) for b in bcs], M)
    sm, em = min(s2 - sd, M), min(e2 - ed, M)
    res = {"idx": -1, "reason": "", "ends": "", "d_start": sd if si >= 0 else "", "d_end": ed if ei >= 0 else ""}
    s_ok = e_ok = False
    core_fail = False
    if si >= 0:
        b = bcs[si]
        s_ok, s_ls, s_le = _locate(b.start, b.core_start, b.kmax_start, b.kcore_start, head)
        core_fail |= not s_ok
    if ei >= 0:
        b = bcs[ei]
        e_ok, e_ls, e_le = _locate(b.end, b.core_end, b.kmax_end, b.kcore_end, tail)
        core_fail |= not e_ok

    def ambiguous(i, j):
        names = sorted({bcs[i].sample, bcs[j].sample}) if j >= 0 else [bcs[i].sample]
        res["reason"] = "barcode ambiguous: " + " or ".join(names)
        return res

    if s_ok and e_ok:
        if si == ei:
            if sm + em < M:
                return ambiguous(si, s2i if sm <= em else e2i)
            idx, ends = si, "both"
        elif sm >= M and em >= M:
            res["reason"] = "different barcodes at the two ends"
            return res
        elif sm >= M and e2i == si:  # the end is torn between its best and the start's barcode
            idx, ends = si, "start"
        elif em >= M and s2i == ei:
            idx, ends = ei, "end"
        elif sm < M and em < M:
            return ambiguous(si, ei)
        else:
            res["reason"] = "different barcodes at the two ends"
            return res
    elif s_ok or e_ok:
        idx, m, i2 = (si, sm, s2i) if s_ok else (ei, em, e2i)
        if m < M:
            return ambiguous(idx, i2)
        ends = "start" if s_ok else "end"
    else:
        if core_fail:
            res["reason"] = "barcode part too different from every sample's barcode (not in the sheet?)"
        else:
            res["reason"] = "no barcode found"
        return res
    if o["require_both"] and ends != "both":
        res["reason"] = "barcode found at one end only (--require-both)"
        return res
    b = bcs[idx]
    a, z = 0, len(seq)
    if ends in ("both", "start"):
        a = s_le + 1
    if ends in ("both", "end"):
        z = tail_off + e_ls
    if o["check_middle"] and z - a > 2 * W:
        mid = seq[a + W // 2: z - W // 2]
        for q, k in ((b.start, b.kmax_start), (b.end, b.kmax_end)):
            if edlib.align(q, mid, mode="HW", k=k)["editDistance"] >= 0:
                res["reason"] = "barcode also inside the read (joined/chimeric read)"
                return res
    if z - a < o["min_length"]:
        res["reason"] = f"shorter than {o['min_length']} bases after trimming"
        return res
    if o["max_length"] and z - a > o["max_length"]:
        res["reason"] = f"longer than {o['max_length']} bases"
        return res
    res.update(idx=idx, ends=ends, trim=(a, z))
    return res


def _work(batch):
    o, bcs, out = _W["o"], _W["bcs"], []
    for h, s, q in batch:
        r = classify(s)
        if _W["asym"] and r["idx"] < 0:  # asymmetric layouts: the reverse strand reads as rc(end)...rc(start)
            r2 = classify(rc(s))
            if r2["idx"] >= 0:
                r, s, q = r2, rc(s), q[::-1]
        if r["idx"] >= 0 and o["trim"]:
            a, z = r["trim"]
            s, q = s[a:z], q[a:z]
        out.append((h, s, q, r["idx"], r["reason"], r["ends"], r["d_start"], r["d_end"]))
    return out


# ---------------------------------------------------------------- statistics and outputs
class Stats:
    def __init__(self, n):
        self.reads = [0] * n
        self.bases = [0] * n
        self.both = [0] * n
        self.lengths = [Counter() for _ in range(n)]  # length in 10-base bins
        self.unclassified = Counter()
        self.total = 0
        self.files_done = []

    def to_json(self):
        return {"reads": self.reads, "bases": self.bases, "both": self.both,
                "lengths": [dict(c) for c in self.lengths], "unclassified": dict(self.unclassified),
                "total": self.total, "files_done": self.files_done}

    @classmethod
    def from_json(cls, d, n):
        st = cls(n)
        if len(d.get("reads", [])) != n:
            die("the saved state belongs to a different sample sheet; use a new --out folder")
        st.reads, st.bases, st.both = d["reads"], d["bases"], d["both"]
        st.lengths = [Counter({int(k): v for k, v in c.items()}) for c in d["lengths"]]
        st.unclassified, st.total, st.files_done = Counter(d["unclassified"]), d["total"], d["files_done"]
        return st


def _median_from_bins(c):
    n = sum(c.values())
    if not n:
        return 0
    acc = 0
    for k in sorted(c):
        acc += c[k]
        if acc * 2 >= n:
            return k * 10 + 5
    return 0


class Writers:
    def __init__(self, out, bcs, keep_unclassified, append, level):
        self.files = {}
        self.out, self.bcs, self.level = Path(out), bcs, level
        self.mode = "ab" if append else "wb"
        self.keep_unclassified = keep_unclassified

    def _get(self, key):
        if key not in self.files:
            name = self.bcs[key].sample if key >= 0 else "unclassified"
            d = self.out / name
            d.mkdir(parents=True, exist_ok=True)
            self.files[key] = gzip.open(d / f"{name}.fastq.gz", self.mode, compresslevel=self.level)
        return self.files[key]

    def write(self, key, h, s, q):
        if key < 0 and not self.keep_unclassified:
            return
        self._get(key).write(f"{h}\n{s}\n+\n{q}\n".encode())

    def close(self):
        for f in self.files.values():
            f.close()
        self.files = {}


def write_summary(out, bcs, st, opts, notes, final):
    out = Path(out)
    assigned = sum(st.reads)
    rows = []
    with open(out / "demux_summary.tsv", "w") as f:
        f.write("sample\tbarcode\twell\treads\tpercent_of_all\tbases\tboth_ends_percent\tmedian_length\n")
        for i, b in enumerate(bcs):
            pct = 100 * st.reads[i] / st.total if st.total else 0
            both = 100 * st.both[i] / st.reads[i] if st.reads[i] else 0
            med = _median_from_bins(st.lengths[i])
            rows.append((b, st.reads[i], pct, st.bases[i], both, med))
            f.write(f"{b.sample}\t{b.name}\t{b.well}\t{st.reads[i]}\t{pct:.2f}\t{st.bases[i]}\t{both:.1f}\t{med}\n")
        un = st.total - assigned
        f.write(f"unclassified\t\t\t{un}\t{100 * un / st.total if st.total else 0:.2f}\t\t\t\n")
    with open(out / "demux_unclassified_reasons.tsv", "w") as f:
        f.write("reason\treads\n")
        for k, v in st.unclassified.most_common():
            f.write(f"{k}\t{v}\n")
    # HTML
    maxr = max([r[1] for r in rows] + [1])
    expected = assigned / len(bcs) if bcs else 0
    trs = []
    for b, n, pct, bases, both, med in rows:
        flag = ""
        if final and st.total and n < 0.1 * expected:
            flag = '<span class="flag">low</span>'
        trs.append(
            f"<tr><td>{html.escape(b.sample)}</td><td>{html.escape(b.name or '')}</td><td>{html.escape(b.well or '')}</td>"
            f"<td class=n>{n:,}</td><td class=n>{pct:.1f}%</td><td><div class=bar style=\"width:{100 * n / maxr:.1f}%\"></div></td>"
            f"<td class=n>{both:.0f}%</td><td class=n>{med}</td><td>{flag}</td></tr>")
    un = st.total - assigned
    reasons = "".join(f"<tr><td>{html.escape(k)}</td><td class=n>{v:,}</td><td class=n>"
                      f"{100 * v / st.total if st.total else 0:.1f}%</td></tr>" for k, v in st.unclassified.most_common())
    note_html = "".join(f"<li>{html.escape(n)}</li>" for n in notes)
    status = "finished" if final else "running (this page is refreshed as files are processed)"
    refresh = "" if final else '<meta http-equiv="refresh" content="30">'
    settings = ", ".join(f"{k}={v}" for k, v in sorted(opts.items()))
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">{refresh}
<meta name="viewport" content="width=device-width,initial-scale=1"><title>AFLA demultiplexing</title>
<style>
:root{{--bg:#fff;--fg:#1d2733;--mut:#5b6773;--line:#dde3ea;--acc:#1f6f8b;--warn:#b54708}}
@media (prefers-color-scheme:dark){{:root{{--bg:#12171d;--fg:#e6edf3;--mut:#9aa7b4;--line:#2c3540;--acc:#5ab0cf;--warn:#f0a35e}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.45 system-ui,sans-serif;margin:0 auto;max-width:1000px;padding:16px}}
h1{{font-size:1.4em;margin:.2em 0}} .mut{{color:var(--mut)}} table{{border-collapse:collapse;width:100%;margin:12px 0}}
td,th{{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.bar{{background:var(--acc);height:10px;border-radius:2px;min-width:1px}} .flag{{color:var(--warn);font-weight:600}}
.box{{border:1px solid var(--line);border-radius:6px;padding:8px 12px;margin:10px 0}} .wrap{{overflow-x:auto}}
</style></head><body>
<p class="mut">AFLA demultiplexer {VERSION}: education and research use only, not for clinical use.</p>
<h1>Demultiplexing summary</h1>
<p>Status: <b>{status}</b>. Reads read: <b>{st.total:,}</b>; assigned to a sample: <b>{assigned:,}</b>
({100 * assigned / st.total if st.total else 0:.1f}%); unclassified: <b>{un:,}</b>. Files processed: {len(st.files_done)}.</p>
<div class="wrap"><table><tr><th>Sample</th><th>Barcode</th><th>Well</th><th>Reads</th><th>% of all</th><th></th>
<th>Both ends</th><th>Median length</th><th></th></tr>{''.join(trs)}</table></div>
<p class="mut">"Both ends": reads whose barcode was found at the start and at the end (the rest had it at one end,
usually because the read is truncated). "low": fewer than 10% of the average reads per sample.</p>
<h2>Why reads were not assigned</h2><div class="wrap"><table><tr><th>Reason</th><th>Reads</th><th>% of all</th></tr>{reasons}</table></div>
{f'<div class="box"><b>Notes</b><ul>{note_html}</ul></div>' if note_html else ''}
<p class="mut">Settings: {html.escape(settings)}</p></body></html>"""
    (out / "demux_report.html").write_text(page)


def process_files(files, bcs, opts, out, st, writers, read_table, pool, state_path=None, notes=(), stop=None):
    """stop: dict set by Ctrl+C. Without state_path (a one-off run) it stops between batches; with state_path
    (--watch, resumable) only between files, so that no file is ever half counted."""
    for fp in files:
        if stop and stop["now"]:
            return
        t0, n0 = time.time(), st.total
        for batch in pool.imap(_work, chunks(iter_file(fp), 2000)):
            if stop and stop["now"] and not state_path:
                print("  interrupted: the outputs hold only the reads processed so far", flush=True)
                return
            for h, s, q, idx, reason, ends, ds, de in batch:
                st.total += 1
                if idx >= 0:
                    st.reads[idx] += 1
                    st.bases[idx] += len(s)
                    st.lengths[idx][len(s) // 10] += 1
                    if ends == "both":
                        st.both[idx] += 1
                else:
                    st.unclassified[reason] += 1
                writers.write(idx, h, s, q)
                if read_table:
                    rid = h[1:].split()[0]
                    read_table.write(f"{rid}\t{bcs[idx].sample if idx >= 0 else 'unclassified'}\t{ends}\t"
                                     f"{ds}\t{de}\t{len(s)}\t{reason}\n")
        st.files_done.append(str(fp))
        dt = time.time() - t0
        print(f"  {fp.name}: {st.total - n0:,} reads in {dt:.0f} s; total {st.total:,}, "
              f"assigned {100 * sum(st.reads) / st.total if st.total else 0:.1f}%", flush=True)
        if state_path:
            for f in writers.files.values():
                f.flush()
            state_path.write_text(json.dumps(st.to_json()))
            write_summary(out, bcs, st, opts, notes, final=False)


# ---------------------------------------------------------------- discovery
def discover(files, out, n_reads, window):
    """Survey a run without a sample sheet: which adapter/barcode layout, which barcodes, draft sample sheet."""
    reads = []
    for fp in files:
        for h, s, q in iter_file(fp):
            reads.append(s)
            if len(reads) >= n_reads:
                break
        if len(reads) >= n_reads:
            break
    n = len(reads)
    if not n:
        die("no reads found")
    lens = sorted(len(s) for s in reads)
    print(f"Read {n:,} reads from {len(files)} file(s); median length {lens[n // 2]}, "
          f"10-90% range {lens[n // 10]}-{lens[9 * n // 10]}.\n")
    print("Known adapter parts (% of reads where found near the start / near the end):")
    for name, e in ELEMENTS:
        k = max(1, len(e) // 6)
        hs = sum(edlib.align(e, s[:window], mode="HW", k=k)["editDistance"] >= 0 for s in reads)
        ts = sum(edlib.align(e, rc(s[-window:]), mode="HW", k=k)["editDistance"] >= 0 for s in reads)
        print(f"  {name:28s} {100 * hs / n:5.1f}% / {100 * ts / n:5.1f}%")
    found = []
    for name, left, right, bl in LAYOUTS:
        pat = re.compile(f"{left}([ACGT]{{{max(4, bl - 3)},{bl + 3}}}){right}")
        c, tail_c = Counter(), Counter()
        for s in reads:
            m = pat.search(s[:window])
            if m:
                c[m.group(1)] += 1
            m = pat.search(rc(s[-window:]))
            if m:
                tail_c[m.group(1)] += 1
        tot = sum(c.values())
        if tot < max(20, 0.01 * n):
            continue
        # merge 1-2 edit variants into the most common exact sequences
        clusters = []
        for seq, cnt in c.most_common():
            for cl in clusters:
                if edlib.align(seq, cl[0], mode="NW", k=2)["editDistance"] >= 0:
                    cl[1] += cnt
                    break
            else:
                clusters.append([seq, cnt])
        main = [cl for cl in clusters if cl[1] >= max(10, 0.01 * tot) and len(cl[0]) == bl]
        print(f"\nLayout: {name} ({left}...barcode...{right}); exact flanks found at the start of "
              f"{100 * tot / n:.1f}% of reads.")
        print(f"  Barcodes carrying >=1% of these reads: {len(main)}")
        for seq, cnt in main:
            print(f"   {seq}  {cnt:6,} reads at start ({100 * cnt / n:.1f}% of all), {tail_c.get(seq, 0):,} exact at end")
        found += [(name, left + seq + right, cnt) for seq, cnt in main]
    if not found:
        print("\nNo known barcode layout was recognised. If you have the kit's sample sheet, run with --sample-sheet "
              "and look at the unclassified reasons.")
        return
    Path(out).mkdir(parents=True, exist_ok=True)
    draft = Path(out) / "discovered_sample_sheet.csv"
    with open(draft, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "barcode", "start", "reads_in_survey"])
        for i, (name, seq, cnt) in enumerate(found, 1):
            w.writerow([f"S{i:02d}", f"found{i:02d}", seq, cnt])
    print(f"\nDraft sample sheet: {draft}\n  Replace S01, S02... with your sample names (no spaces), delete rows you do "
          "not want, then run with --sample-sheet. The start sequence includes the flanks; the end is taken as its "
          "reverse complement.\n  Barcodes found in a run are not a substitute for the kit's sheet: check them against "
          "the wells you loaded.")


# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-sheet", help="sample sheet (.csv or .xlsx) with sample,start[,end]")
    ap.add_argument("--input", nargs="+", help="FASTQ files, folders (e.g. MinKNOW fastq_pass) or .zip archives")
    ap.add_argument("--out", help="output folder")
    ap.add_argument("--threads", type=int, default=max(1, min(8, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--max-error", type=float, default=0.2,
                    help="largest edit distance accepted, as a fraction of the barcode sequence length (default 0.2)")
    ap.add_argument("--max-barcode-error", type=float, default=0.2,
                    help="largest edit distance in the barcode part itself (without the shared flanks), as a "
                         "fraction of its length (default 0.2: 1 error in 8 bases, 4 in 24)")
    ap.add_argument("--min-margin", type=int, default=3,
                    help="the best barcode must be at least this many edits closer than the second best (default 3)")
    ap.add_argument("--window", type=int, default=150, help="bases searched at each end of the read (default 150)")
    ap.add_argument("--require-both", action="store_true", help="assign only reads with the barcode at both ends")
    ap.add_argument("--no-trim", action="store_true", help="keep the barcode/adapter sequence in the output reads")
    ap.add_argument("--no-chimera-check", action="store_true", help="do not look for the barcode inside the read")
    ap.add_argument("--min-length", type=int, default=50, help="shortest read kept after trimming (default 50)")
    ap.add_argument("--max-length", type=int, default=0, help="longest read kept (default: no limit)")
    ap.add_argument("--no-unclassified", action="store_true", help="do not write the unclassified reads")
    ap.add_argument("--read-table", action="store_true", help="also write demux_reads.tsv.gz (one line per read)")
    ap.add_argument("--compress-level", type=int, default=3, choices=range(1, 10), metavar="1-9")
    ap.add_argument("--force", action="store_true", help="overwrite an existing output folder")
    ap.add_argument("--watch", action="store_true", help="keep watching the input folder(s) for new files")
    ap.add_argument("--poll", type=int, default=30, help="--watch: seconds between checks (default 30)")
    ap.add_argument("--settle", type=int, default=60,
                    help="--watch: a file is read once unchanged for this many seconds (default 60)")
    ap.add_argument("--stop-when-idle", type=float, default=0,
                    help="--watch: stop after this many minutes without new files (default: run until Ctrl+C)")
    ap.add_argument("--discover", action="store_true", help="survey the run (no sample sheet needed)")
    ap.add_argument("--discover-reads", type=int, default=20000)
    ap.add_argument("--check-sheet", action="store_true", help="read the sample sheet, report, and stop")
    ap.add_argument("--version", action="version", version=VERSION)
    a = ap.parse_args(argv)

    if edlib is None:
        die("the edlib module is missing: pip install edlib (or run inside the AFLA container)")

    if a.discover:
        if not a.input:
            die("--discover needs --input")
        discover(list_inputs(a.input), a.out or ".", a.discover_reads, a.window)
        return 0

    if not a.sample_sheet:
        die("--sample-sheet is required (or use --discover to see what barcodes the run contains)")
    bcs, notes = read_sample_sheet(win_path(a.sample_sheet))
    set_thresholds(bcs, a.max_error, a.max_barcode_error)
    print(f"Sample sheet: {len(bcs)} samples; barcode sequence lengths "
          f"{sorted({len(b.start) for b in bcs})}; edit distance accepted up to {a.max_error:.0%} of the length.")
    close = close_pairs(bcs, 2 * a.min_margin)
    for d, x, y in close[:10]:
        notes.append(f"barcodes of {x} and {y} differ by only {d} edit(s): many of their reads will be called "
                     "ambiguous, mostly those with the barcode at one end only. If you can, do not use these "
                     "two wells in the same run.")
    if len(close) > 10:
        notes.append(f"... and {len(close) - 10} more close pairs")
    for n_ in notes:
        print(f"  Note: {n_}")
    if a.check_sheet:
        for b in bcs:
            print(f"  {b.sample:20s} {b.name or '':10s} {b.well or '':4s} start {len(b.start)} bases "
                  f"(barcode part {b.core_start[1] - b.core_start[0]} bases at {b.core_start[0] + 1}-{b.core_start[1]}), "
                  f"end {len(b.end)} bases")
        return 0
    if not a.input or not a.out:
        die("--input and --out are required")

    out = win_path(a.out)
    state_path = out / ".demux_state.json"
    resume = a.watch and state_path.exists() and not a.force
    if out.exists() and any(out.iterdir()) and not (a.force or resume):
        die(f"{out} already has files; use --force to overwrite or choose another --out")
    out.mkdir(parents=True, exist_ok=True)
    if a.force:
        for f in list(out.glob("*/*.fastq.gz")) + [state_path, out / "demux_reads.tsv.gz"]:
            if f.exists():
                f.unlink()

    opts = {"window": a.window, "min_margin": a.min_margin, "require_both": a.require_both, "trim": not a.no_trim,
            "check_middle": not a.no_chimera_check, "min_length": a.min_length, "max_length": a.max_length,
            "max_error": a.max_error, "max_barcode_error": a.max_barcode_error}
    st = Stats.from_json(json.loads(state_path.read_text()), len(bcs)) if resume else Stats(len(bcs))
    writers = Writers(out, bcs, not a.no_unclassified, append=resume, level=a.compress_level)
    read_table = None
    if a.read_table:
        read_table = gzip.open(out / "demux_reads.tsv.gz", "at" if resume else "wt", compresslevel=a.compress_level)
        if not resume:
            read_table.write("read_id\tsample\tbarcode_ends\tdist_start\tdist_end\tlength\tunclassified_reason\n")

    stop = {"now": False}

    def _stop(signum, frame):
        stop["now"] = True
        print("\nStopping after the current file...", flush=True)

    pool = Pool(a.threads, initializer=_init_worker, initargs=(bcs, opts))
    signal.signal(signal.SIGINT, _stop)
    signal.signal(signal.SIGTERM, _stop)
    t0 = time.time()
    try:
        if not a.watch:
            files = list_inputs(a.input)
            print(f"Reading {len(files)} file(s) with {a.threads} worker(s)...", flush=True)
            process_files(files, bcs, opts, out, st, writers, read_table, pool, notes=notes, stop=stop)
        else:
            dirs = [win_path(p) for p in a.input]
            if not all(d.is_dir() for d in dirs):
                die("--watch needs folder(s) as --input")
            done = set(st.files_done)
            last_new = time.time()
            print(f"Watching {', '.join(map(str, dirs))} (every {a.poll} s). Ctrl+C to stop.", flush=True)
            while not stop["now"]:
                now = time.time()
                ready = []
                for d in dirs:
                    for f in sorted(d.rglob("*")):
                        if (f.is_file() and FASTQ_RE.search(f.name) and str(f) not in done
                                and now - f.stat().st_mtime >= a.settle):
                            ready.append(f)
                if ready:
                    process_files(ready, bcs, opts, out, st, writers, read_table, pool, state_path, notes, stop)
                    done.update(st.files_done)
                    last_new = time.time()
                elif a.stop_when_idle and time.time() - last_new > a.stop_when_idle * 60:
                    print(f"No new files for {a.stop_when_idle:g} min: stopping.")
                    break
                for _ in range(a.poll):
                    if stop["now"]:
                        break
                    time.sleep(1)
    finally:
        pool.terminate()
        writers.close()
        if read_table:
            read_table.close()
        if a.watch:
            state_path.write_text(json.dumps(st.to_json()))
    write_summary(out, bcs, st, opts, notes, final=True)
    assigned = sum(st.reads)
    print(f"\nDone in {time.time() - t0:.0f} s: {st.total:,} reads, {assigned:,} assigned "
          f"({100 * assigned / st.total if st.total else 0:.1f}%).")
    for b, n_ in sorted(zip(bcs, st.reads), key=lambda x: -x[1]):
        print(f"  {b.sample:20s} {n_:>10,}")
    print(f"  {'unclassified':20s} {st.total - assigned:>10,}")
    for k, v in st.unclassified.most_common():
        print(f"      {v:>10,}  {k}")
    print(f"Summary: {out / 'demux_report.html'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
