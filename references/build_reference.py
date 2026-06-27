#!/usr/bin/env python3
"""Build the alignment panel from per-virus FASTAs in references/raw/.

fetch_references.sh downloads each manifest row into references/raw/<taxon>.fasta
(a multi-FASTA for segmented/multi-strain viruses). This script groups every
contig in a file under that file's <taxon>, so flu's 8 segments stay separate
contigs (per-segment coverage) but share one taxon label (per-virus calls).

Writes:
  references/panel.fasta        contigs renamed <taxon>__<accession>
  references/panel.lengths.tsv  contig_label  taxon  length
"""
from __future__ import annotations
import glob, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
PANEL = os.path.join(HERE, "panel.fasta")
LENGTHS = os.path.join(HERE, "panel.lengths.tsv")


def iter_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def accession_of(header: str) -> str:
    # first whitespace token, strip version and any db prefix
    tok = header.split()[0]
    return tok.split("|")[-1].split(".")[0]


def main():
    raws = []
    for ext in ("*.fasta", "*.fa", "*.fna"):
        raws += glob.glob(os.path.join(RAW, ext))
    raws = sorted(set(raws))
    if not raws:
        sys.exit(f"No FASTAs in {RAW}. Run fetch_references.sh (or download manually) first.")

    seen, n_taxa = set(), 0
    with open(PANEL, "w") as out, open(LENGTHS, "w") as lf:
        lf.write("contig_label\ttaxon\tlength\n")
        for fa in raws:
            taxon = os.path.splitext(os.path.basename(fa))[0]   # filename = taxon
            n_taxa += 1
            n_contigs = 0
            for header, seq in iter_fasta(fa):
                if not seq:
                    continue
                acc = accession_of(header)
                label = f"{taxon}__{acc}"
                if label in seen:
                    label = f"{taxon}__{acc}_{n_contigs}"   # disambiguate dups
                seen.add(label)
                out.write(f">{label}\n")
                for i in range(0, len(seq), 70):
                    out.write(seq[i:i + 70] + "\n")
                lf.write(f"{label}\t{taxon}\t{len(seq)}\n")
                n_contigs += 1
            print(f"  {taxon}: {n_contigs} contig(s)")
    print(f"Wrote {PANEL} ({n_taxa} taxa) and {LENGTHS}")


if __name__ == "__main__":
    main()
