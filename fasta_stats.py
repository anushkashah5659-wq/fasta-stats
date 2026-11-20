#!/usr/bin/env python3
"""Report basic properties (length, GC content, base composition) of every
sequence in a FASTA file. Standard library only. Supports .gz input."""

import argparse
import gzip
import sys
from collections import Counter


def open_text(path):
    """Open a plain or gzip-compressed text file (or '-' for stdin)."""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_fasta(handle):
    """Yield (record_id, description, sequence) one record at a time.

    Handles multi-line sequences, blank lines and lowercase bases.
    Streams the file, so large files do not need to fit in memory.
    """
    header, chunks = None, []
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield (*split_header(header), "".join(chunks))
            header, chunks = line[1:], []
        else:
            if header is None:
                raise ValueError("Sequence data found before first '>' header; not valid FASTA.")
            chunks.append(line)
    if header is not None:
        yield (*split_header(header), "".join(chunks))


def split_header(header):
    parts = header.split(None, 1)
    return (parts[0] if parts else "", parts[1] if len(parts) > 1 else "")


def analyse(seq):
    """Return length, GC% and per-base counts for one sequence."""
    counts = Counter(seq.upper())
    a, c, g, t = (counts.get(b, 0) for b in "ACGT")
    n = counts.get("N", 0)
    other = len(seq) - (a + c + g + t + n)
    acgt = a + c + g + t
    # GC% is computed over unambiguous A/C/G/T bases only (N and IUPAC codes excluded).
    gc = 100 * (g + c) / acgt if acgt else float("nan")
    return {"length": len(seq), "gc": gc, "A": a, "C": c, "G": g, "T": t, "N": n, "other": other}


COLUMNS = ["id", "length", "gc_percent", "A", "C", "G", "T", "N", "other"]


def row(name, r):
    gc = "NA" if r["gc"] != r["gc"] else f"{r['gc']:.2f}"
    return [name, r["length"], gc, r["A"], r["C"], r["G"], r["T"], r["N"], r["other"]]


def main():
    p = argparse.ArgumentParser(description="Length, GC content and composition for each FASTA sequence.")
    p.add_argument("fasta", help="FASTA file (.fa/.fasta, optionally .gz), or '-' for stdin")
    p.add_argument("-o", "--output", help="write TSV report here instead of stdout")
    args = p.parse_args()

    totals = Counter()
    rows = []
    with open_text(args.fasta) as fh:
        for rec_id, _desc, seq in read_fasta(fh):
            r = analyse(seq)
            rows.append(row(rec_id, r))
            for k in ("length", "A", "C", "G", "T", "N", "other"):
                totals[k] += r[k]

    if not rows:
        sys.exit("No FASTA records found in input.")

    acgt = sum(totals[b] for b in "ACGT")
    overall_gc = 100 * (totals["G"] + totals["C"]) / acgt if acgt else float("nan")
    total_r = {**totals, "gc": overall_gc}
    rows.append(row("TOTAL", total_r))

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        out.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
        out.write(f"# {len(rows) - 1} sequence(s), {totals['length']} bases total\n")
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    main()
