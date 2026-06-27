"""Align reads to the curated reference panel with mappy (minimap2, map-ont).

mappy is minimap2's official Python binding: identical alignment engine and
the same `map-ont` preset you'd get on the command line. We keep only the
primary alignment per read and record the fields needed for confidence
scoring (identity, mapq, query/ref spans) and for coverage (cigar).
"""
from __future__ import annotations

from dataclasses import dataclass
import mappy as mp


@dataclass
class Hit:
    read_id: str
    read_len: int
    ctg: str          # reference (taxon) name aligned to
    r_st: int         # 0-based reference start
    r_en: int
    q_st: int
    q_en: int
    strand: int
    mapq: int
    mlen: int         # matching bases
    blen: int         # alignment block length
    cigar: list       # list of (length, op)

    @property
    def identity(self) -> float:
        return self.mlen / self.blen if self.blen else 0.0

    @property
    def query_cov(self) -> float:
        return (self.q_en - self.q_st) / self.read_len if self.read_len else 0.0

    def as_cov_hit(self) -> dict:
        return {"ctg": self.ctg, "r_st": self.r_st, "cigar": self.cigar}


def build_aligner(reference_fa: str, preset: str = "map-ont", n_threads: int = 3) -> mp.Aligner:
    aln = mp.Aligner(reference_fa, preset=preset, n_threads=n_threads)
    if not aln:
        raise RuntimeError(f"Failed to load/index reference: {reference_fa}")
    return aln


def align_reads(aligner: mp.Aligner, fasta_path: str, min_mapq: int = 0):
    """Yield the best (primary) Hit per read that aligns to the panel.

    Reads with no alignment are skipped (they remain 'unclassified' for this
    confirmation step; discovery via kraken2/sourmash still sees them).
    """
    for name, seq, _qual in mp.fastx_read(fasta_path, read_comment=False):
        best = None
        for h in aligner.map(seq):  # cs/MD off by default; cigar uses M/I/D
            if not h.is_primary:
                continue
            if best is None or h.mlen > best.mlen:
                best = h
        if best is None or best.mapq < min_mapq:
            continue
        yield Hit(
            read_id=name, read_len=len(seq), ctg=best.ctg,
            r_st=best.r_st, r_en=best.r_en, q_st=best.q_st, q_en=best.q_en,
            strand=best.strand, mapq=best.mapq, mlen=best.mlen, blen=best.blen,
            cigar=list(best.cigar),
        )
