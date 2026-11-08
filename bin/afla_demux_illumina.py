#!/usr/bin/env python3
"""AFLA demultiplexer for Illumina FASTQ files (education and research use only, not for clinical use).

Sorts reads into one set of FASTQ files per sample using the index in each read's name, or in separate I1/I2 files.
It starts from FASTQ (for example Undetermined_S0_R1_001.fastq.gz); converting raw BCL files needs Illumina's own software.

  afla_demux_illumina.py --sample-sheet sheet.csv --input fastq_folder/ --out demux/

Sample sheets: a plain CSV with columns sample,index[,index2], or an Illumina SampleSheet.csv ([Data] or [BCLConvert_Data]
section with Sample_ID/Sample_Name, index, index2); .xlsx also works. Reads are matched with up to --mismatches differences
(default 1). The sheet is refused if two samples could be confused at that setting. The direction of index 2 differs
between instruments, so it is detected from the reads unless --i5 is given. Reads whose index 1 and index 2 belong to
different samples are reported as possible index hopping and never assigned.

Outputs in --out: <sample>/<sample>_R1.fastq.gz (and _R2), Undetermined/, demux_summary.tsv, demux_report.html.
Standard library only.
"""
import argparse
import csv
import gzip
import html
import io
import re
import sys
import time
from collections import Counter
from pathlib import Path

VERSION = "0.1.0"
FASTQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
SAFE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
READ_TAG_RE = re.compile(r"(_R?)1((?:_\d{3})?\.(?:fastq|fq)(?:\.gz)?)$", re.I)
HEADER_INDEX_RE = re.compile(r"\s\d:[YN]:\d+:([ACGTN]+)(?:[+_-]([ACGTN]+))?\s*$")
COMP = str.maketrans("ACGTN", "TGCAN")
RESERVED = {"undetermined", "unclassified"}
AUTO_READS = 100000
HOPPED, UNKNOWN, AMBIGUOUS = ("index 1 and index 2 belong to different samples (possible index hopping)",
                              "index not in the sample sheet", "index matches more than one sample")


def die(msg):
    sys.exit(f"afla_demux_illumina: {msg}")


def rc(s):
    return s.translate(COMP)[::-1]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


# ---------------------------------------------------------------- sample sheet
def _rows(path):
    p = Path(path)
    if not p.exists():
        die(f"sample sheet not found: {path}")
    if p.suffix.lower() in (".xlsx", ".xlsm"):
        from afla_cftr import read_xlsx
        book = read_xlsx(p)
        return [[("" if c is None else str(c)) for c in r] for r in next(iter(book.values()))]
    raw = p.read_bytes()
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    first = next((ln for ln in text.splitlines() if ln.strip()), "")
    delim = max([",", ";", "\t"], key=first.count)
    return list(csv.reader(io.StringIO(text), delimiter=delim))


def _key(c):
    return re.sub(r"[^a-z0-9]", "", str(c or "").lower())


def read_sheet(path):
    """Returns [(sample, index, index2)] with index2 '' for single-index runs."""
    rows = [r for r in _rows(path) if any(str(c).strip() for c in r)]
    start = 0
    sections = [i for i, r in enumerate(rows) if str(r[0]).strip().startswith("[")]
    data = [i for i in sections if _key(rows[i][0]) in ("data", "bclconvertdata")]
    if sections and not data:
        die(f"{path}: Illumina sheet without a [Data] or [BCLConvert_Data] section")
    end = len(rows)
    if data:
        start = data[0] + 1
        end = next((i for i in sections if i > data[0]), len(rows))
    head = [_key(c) for c in rows[start]]
    if "index" not in head:
        die(f"{path}: no 'index' column (expected sample,index[,index2] or an Illumina SampleSheet)")
    names = [head.index(k) for k in ("samplename", "sample", "sampleid", "alias") if k in head]
    if not names:
        die(f"{path}: no sample column (expected sample, Sample_Name or Sample_ID)")
    i7, i5 = head.index("index"), (head.index("index2") if "index2" in head else None)
    out, seen = [], set()
    for r in rows[start + 1: end]:
        cell = lambda j: str(r[j]).strip() if j is not None and j < len(r) else ""  # noqa: E731
        name = next((cell(j) for j in names if cell(j)), "")
        a, b = cell(i7).upper(), cell(i5).upper()
        if not name and not a:
            continue
        if not name or not a:
            die(f"{path}: a row has a sample without an index, or an index without a sample ({name or a})")
        if not SAFE_NAME_RE.match(name):
            die(f"sample name {name!r}: use only letters, digits, '-', '_' or '.' (no spaces)")
        if name.lower() in RESERVED:
            die(f"{name!r} is reserved; choose another sample name")
        if name.lower() in seen:
            die(f"sample name {name!r} is used twice")
        if not re.fullmatch("[ACGT]+", a) or (b and not re.fullmatch("[ACGT]+", b)):
            die(f"sample {name}: the index may contain only A, C, G, T")
        seen.add(name.lower())
        out.append((name, a, b))
    if not out:
        die(f"{path}: no samples found")
    if len({len(x[1]) for x in out}) > 1 or len({len(x[2]) for x in out}) > 1:
        die("all samples in one run must have index sequences of the same length, and all or none need index 2")
    return out


def check_collisions(samples, mm):
    """Two samples are confusable when both of their indexes could match the same read at this mismatch setting."""
    bad = []
    for i, (s1, a1, b1) in enumerate(samples):
        for s2, a2, b2 in samples[i + 1:]:
            if a1 == a2 and b1 == b2:
                die(f"samples {s1} and {s2} have the same index")
            if hamming(a1, a2) <= 2 * mm and (not b1 or hamming(b1, b2) <= 2 * mm):
                bad.append((s1, s2))
    return bad


def variants(seq, mm):
    """Every sequence within mm mismatches of seq (N counts as a mismatch)."""
    found = {seq}
    frontier = {seq}
    for _ in range(mm):
        nxt = set()
        for s in frontier:
            for i in range(len(s)):
                for c in "ACGTN":
                    if c != s[i]:
                        nxt.add(s[:i] + c + s[i + 1:])
        found |= nxt
        frontier = nxt
    return found


class Matcher:
    def __init__(self, samples, mm, i5_rc):
        self.dual = bool(samples[0][2])
        self.len7, self.len5 = len(samples[0][1]), len(samples[0][2])
        self.m7, self.m5 = {}, {}
        for n, (_, a, b) in enumerate(samples):
            for v in variants(a, mm):
                self.m7.setdefault(v, set()).add(n)
            if self.dual:
                for v in variants(rc(b) if i5_rc else b, mm):
                    self.m5.setdefault(v, set()).add(n)
        self.cache = {}

    def match(self, i7, i5):
        """Returns (sample number or -1, reason)."""
        key = (i7, i5)
        hit = self.cache.get(key)
        if hit is not None:
            return hit
        a = self.m7.get(i7[:self.len7], set())
        res = (-1, UNKNOWN)
        if not self.dual:
            if len(a) == 1:
                res = (next(iter(a)), "")
            elif a:
                res = (-1, AMBIGUOUS)
        else:
            b = self.m5.get(i5[:self.len5], set())
            both = a & b
            if len(both) == 1:
                res = (next(iter(both)), "")
            elif both:
                res = (-1, AMBIGUOUS)
            elif a and b:
                res = (-1, HOPPED)
        if len(self.cache) < 2000000:
            self.cache[key] = res
        return res


# ---------------------------------------------------------------- reading
def find_runs(paths):
    """[(R1, R2 or None, I1 or None, I2 or None)] from files and folders."""
    files = []
    for p in map(Path, paths):
        if p.is_dir():
            files += sorted(f for f in p.rglob("*") if f.is_file() and FASTQ_RE.search(f.name))
        elif p.is_file():
            files.append(p)
        else:
            die(f"input not found: {p}")
    r1s = [f for f in files if READ_TAG_RE.search(f.name)]
    if not r1s:
        die("no read 1 FASTQ found (expected names like Undetermined_S0_L001_R1_001.fastq.gz or x_R1.fastq.gz)")
    runs = []
    have = set(files)

    def partner(f, tag):
        q = f.with_name(READ_TAG_RE.sub(lambda m: tag + m.group(2), f.name))
        return q if q in have else None

    for f in r1s:
        prefix = READ_TAG_RE.search(f.name).group(1)
        has_index_files = prefix.upper() == "_R"    # only the _R1 naming style has _I1/_I2 companions
        runs.append((f, partner(f, f"{prefix}2"), partner(f, "_I1") if has_index_files else None,
                     partner(f, "_I2") if has_index_files else None))
    return runs


def _open(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="ascii", errors="replace") if gz else open(path, encoding="ascii", errors="replace")


def records(path):
    if path is None:
        while True:
            yield None
    with _open(path) as f:
        while True:
            h = f.readline()
            if not h:
                return
            s, p, q = f.readline(), f.readline(), f.readline()
            if not h.startswith("@") or not p.startswith("+") or len(s) != len(q):
                die(f"{path}: not a valid FASTQ record near {h[:60]!r} (file may be truncated)")
            yield h, s, p, q


def read_id(h):
    return h[1:].split(None, 1)[0].rsplit("/", 1)[0]


def index_of(rec, i1, i2, source):
    if i1 is not None:
        return i1[1].strip().upper(), (i2[1].strip().upper() if i2 is not None else "")
    m = HEADER_INDEX_RE.search(rec[0])
    if not m:
        die(f"{source}: no index sequence in the read name ({rec[0].strip()[:80]!r}). These reads may already be "
            "demultiplexed; otherwise supply I1/I2 FASTQ files next to them.")
    return m.group(1), m.group(2) or ""


def detect_i5(runs, samples, n=AUTO_READS):
    """True when the reads carry the reverse complement of the sheet's index 2."""
    fwd = {s[2] for s in samples}
    rev = {rc(s[2]) for s in samples}
    a = b = seen = 0
    for r1, _, i1f, i2f in runs:
        for rec, i1, i2 in zip(records(r1), records(i1f), records(i2f)):
            x = index_of(rec, i1, i2, r1)[1][:len(samples[0][2])]
            a += x in fwd
            b += x in rev
            seen += 1
            if seen >= n:
                return b > a, a, b
    return b > a, a, b


class Writers:
    def __init__(self, out, names, keep_undetermined, level):
        self.out, self.names, self.keep, self.level, self.h = Path(out), names, keep_undetermined, level, {}

    def _get(self, key):
        if key not in self.h:
            name = self.names[key[0]] if key[0] >= 0 else "Undetermined"
            d = self.out / name
            d.mkdir(parents=True, exist_ok=True)
            self.h[key] = gzip.open(d / f"{name}_R{key[1]}.fastq.gz", "wb", compresslevel=self.level)
        return self.h[key]

    def write(self, n, r, rec):
        if n < 0 and not self.keep:
            return
        self._get((n, r)).write("".join(rec).encode())

    def close(self):
        for f in self.h.values():
            f.close()


def demultiplex(runs, samples, matcher, out, keep_undetermined=True, level=3):
    w = Writers(out, [s[0] for s in samples], keep_undetermined, level)
    reads, reasons, unknown, total = [0] * len(samples), Counter(), Counter(), 0
    try:
        for r1, r2, i1f, i2f in runs:
            for rec, mate, i1, i2 in zip(records(r1), records(r2), records(i1f), records(i2f)):
                if mate is not None and read_id(rec[0]) != read_id(mate[0]):
                    die(f"{r1.name} and {r2.name} are not in the same order at {read_id(rec[0])}")
                a, b = index_of(rec, i1, i2, r1)
                n, why = matcher.match(a, b)
                total += 1
                if n >= 0:
                    reads[n] += 1
                else:
                    reasons[why] += 1
                    unknown[f"{a}+{b}" if b else a] += 1
                w.write(n, 1, rec)
                if mate is not None:
                    w.write(n, 2, mate)
    finally:
        w.close()
    return reads, reasons, unknown, total


# ---------------------------------------------------------------- outputs
def write_outputs(out, samples, reads, reasons, unknown, total, notes):
    out = Path(out)
    assigned = sum(reads)
    un = total - assigned
    pc = lambda n: 100 * n / total if total else 0  # noqa: E731
    with open(out / "demux_summary.tsv", "w") as f:
        f.write("sample\tindex\tindex2\treads\tpercent_of_all\n")
        for (s, a, b), n in zip(samples, reads):
            f.write(f"{s}\t{a}\t{b}\t{n}\t{pc(n):.2f}\n")
        f.write(f"Undetermined\t\t\t{un}\t{pc(un):.2f}\n")
    top = max(reads + [1])
    e = html.escape
    rows = "".join(f"<tr><td>{e(s)}</td><td><code>{a}{'+' + b if b else ''}</code></td><td class=n>{n:,}</td>"
                   f"<td class=n>{pc(n):.1f}%</td><td><div class=bar style='width:{100 * n / top:.1f}%'></div></td></tr>"
                   for (s, a, b), n in zip(samples, reads))
    why = "".join(f"<tr><td>{e(k)}</td><td class=n>{v:,}</td></tr>" for k, v in reasons.most_common())
    unk = "".join(f"<tr><td><code>{e(k)}</code></td><td class=n>{v:,}</td></tr>" for k, v in unknown.most_common(20))
    note = "".join(f"<li>{e(n)}</li>" for n in notes)
    page = f"""<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport content="width=device-width,initial-scale=1">
<title>AFLA demultiplexing</title><style>
:root{{--bg:#fff;--fg:#1d2733;--mut:#5b6773;--line:#dde3ea;--acc:#1f6f8b}}
@media (prefers-color-scheme:dark){{:root{{--bg:#12171d;--fg:#e6edf3;--mut:#9aa7b4;--line:#2c3540;--acc:#5ab0cf}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.45 system-ui,sans-serif;margin:0 auto;max-width:1000px;padding:16px}}
h1{{font-size:1.4em;margin:.2em 0}}.mut{{color:var(--mut)}}.wrap{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;margin:12px 0}}
td,th{{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left}}td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.bar{{background:var(--acc);height:10px;border-radius:2px;min-width:1px}}</style></head><body>
<p class=mut>AFLA demultiplexer (Illumina) {VERSION}: education and research use only, not for clinical use.</p>
<h1>Demultiplexing summary</h1>
<p>Read pairs: <b>{total:,}</b>; assigned: <b>{assigned:,}</b> ({pc(assigned):.1f}%); undetermined: <b>{un:,}</b> ({pc(un):.1f}%).</p>
<div class=wrap><table><tr><th>Sample</th><th>Index</th><th>Reads</th><th>% of all</th><th></th></tr>{rows}</table></div>
<h2>Why reads were not assigned</h2><div class=wrap><table><tr><th>Reason</th><th>Reads</th></tr>{why}</table></div>
<h2>Most common unassigned indexes</h2><div class=wrap><table><tr><th>Index</th><th>Reads</th></tr>{unk}</table></div>
{'<h2>Notes</h2><ul>' + note + '</ul>' if note else ''}</body></html>"""
    (out / "demux_report.html").write_text(page)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-sheet", required=True)
    ap.add_argument("--input", nargs="+", required=True, help="FASTQ files or folders")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mismatches", type=int, default=1, choices=[0, 1, 2], help="allowed per index (default 1)")
    ap.add_argument("--i5", choices=["auto", "forward", "rc"], default="auto",
                    help="direction of index 2 in the reads relative to the sheet (default: detect from the reads)")
    ap.add_argument("--no-undetermined", action="store_true", help="do not write the unassigned reads")
    ap.add_argument("--compress-level", type=int, default=3, choices=range(1, 10), metavar="1-9")
    ap.add_argument("--force", action="store_true", help="write into an output folder that already has files")
    ap.add_argument("--version", action="version", version=VERSION)
    a = ap.parse_args(argv)
    samples = read_sheet(a.sample_sheet)
    bad = check_collisions(samples, a.mismatches)
    if bad:
        x, y = bad[0]
        die(f"samples {x} and {y} (and {len(bad) - 1} other pair(s)) could be confused with {a.mismatches} mismatch(es): "
            "use a lower --mismatches or check the sheet")
    runs = find_runs(a.input)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()) and not a.force:
        die(f"{out} already has files; use --force or choose another --out")
    out.mkdir(parents=True, exist_ok=True)
    notes = []
    i5_rc = a.i5 == "rc"
    if samples[0][2] and a.i5 == "auto":
        i5_rc, fwd, rev = detect_i5(runs, samples)
        notes.append(f"Index 2 direction detected from the reads: {'reverse complement' if i5_rc else 'as in the sheet'} "
                     f"({fwd:,} exact matches as written, {rev:,} as reverse complement).")
        if not fwd and not rev:
            notes.append("No read's index 2 exactly matched the sheet in either direction: check the sample sheet.")
    t0 = time.time()
    reads, reasons, unknown, total = demultiplex(runs, samples, Matcher(samples, a.mismatches, i5_rc), out,
                                                 not a.no_undetermined, a.compress_level)
    if not total:
        die("no reads found in the input")
    write_outputs(out, samples, reads, reasons, unknown, total, notes)
    print(f"Done in {time.time() - t0:.0f} s: {total:,} reads, {sum(reads):,} assigned ({100 * sum(reads) / total:.1f}%).")
    for n in notes:
        print(f"  {n}")
    for (s, _, _), c in sorted(zip(samples, reads), key=lambda x: -x[1]):
        print(f"  {s:24s} {c:>12,}")
    print(f"Summary: {out / 'demux_report.html'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
