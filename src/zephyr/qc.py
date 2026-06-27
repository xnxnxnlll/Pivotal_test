"""Lightweight FASTA QC (read counts, length distribution, N50).

Pure-Python so it runs anywhere; for richer ONT QC you can additionally run
NanoPlot on your server, but this is enough to characterise a pool and to
decide whether a yield is high enough to expect genome-wide coverage.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np
import mappy as mp


@dataclass
class PoolQC:
    pool: str
    n_reads: int
    total_bp: int
    min_len: int
    max_len: int
    median_len: int
    mean_len: float
    n50: int

    def to_row(self) -> dict:
        return asdict(self)


def _n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    half = sum(s) / 2
    acc = 0
    for L in s:
        acc += L
        if acc >= half:
            return L
    return s[-1]


def qc_fasta(pool: str, fasta_path: str) -> PoolQC:
    lengths = [len(seq) for _n, seq, _q in mp.fastx_read(fasta_path)]
    if not lengths:
        return PoolQC(pool, 0, 0, 0, 0, 0, 0.0, 0)
    arr = np.array(lengths)
    return PoolQC(
        pool=pool, n_reads=len(lengths), total_bp=int(arr.sum()),
        min_len=int(arr.min()), max_len=int(arr.max()),
        median_len=int(np.median(arr)), mean_len=round(float(arr.mean()), 1),
        n50=_n50(lengths),
    )
