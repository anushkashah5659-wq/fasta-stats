# fasta-stats

A small Python script that reads a FASTA file and reports, for each sequence:
**length**, **GC content (%)** and **base composition** (A, C, G, T, N, other),
plus a TOTAL row for the whole file.

Standard library only (Python 3.7+). No install needed.

## Run it

```bash
python3 fasta_stats.py path/to/file.fasta                 # prints a TSV table
python3 fasta_stats.py path/to/file.fasta.gz -o report.tsv  # gzip works, save to file
```

Works on any FASTA file: single or multi-record, wrapped or unwrapped lines,
upper or lower case.

## Example on a public dataset

The bacteriophage phiX174 genome (NCBI RefSeq NC_001422.1), downloaded from NCBI:

```bash
curl -o data/phix174.fasta \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta"
python3 fasta_stats.py data/phix174.fasta -o output/phix174_report.tsv
```

The output is in `output/phix174_report.tsv`. It is generated automatically by GitHub Actions (`.github/workflows/run.yml`) each time the repo is pushed, so the report always comes from the real script run on the real data.
A tiny synthetic file (`data/example.fasta`) is also included so the script can be tested offline.

## Output columns

`id, length, gc_percent, A, C, G, T, N, other`

## Decisions

- **GC% is calculated over A/C/G/T only.** N and other ambiguity codes are
  excluded from the denominator so gaps do not dilute the result. They are
  still counted in the `N` / `other` columns.
- **Streaming parser**: records are read one at a time, so large files are fine.
- **Bad input**: a file with sequence before any `>` header, or with no records,
  stops with a clear error instead of giving silent wrong output.
- I wrote the FASTA parser myself (about 15 lines) rather than using Biopython,
  to keep it dependency-free.

## Limitations

DNA/RNA focus: `U` in RNA would be counted under `other`. Protein FASTA files
will run but the composition columns are not meaningful for them.
