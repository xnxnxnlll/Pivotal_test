"""Per-reference coverage from alignment intervals.

Coverage is computed by walking each primary alignment's CIGAR over the
reference, marking only reference-consuming, base-present operations
(M/=/X) as covered and skipping deletions/ref-skips (D/N). This matches
`samtools depth` semantics (depth = number of reads with a base at a
position) closely enough for ONT respiratory-virus work, and is exact and
unit-testable without writing a BAM.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

# CIGAR op codes (BAM spec, as returned by mappy): 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X
_COVER_OPS = {0, 7, 8}   # consume reference AND carry a base
_SKIP_OPS = {2, 3}       # consume reference, no base (deletion / ref-skip)


@dataclass
class RefCoverage:
    name: str
    length: int
    depth: np.ndarray = field(repr=False)

    @property
    def covered_bp(self) -> int:
        return int((self.depth > 0).sum())

    @property
    def breadth(self) -> float:
        """Fraction of reference positions with depth >= 1."""
        return self.covered_bp / self.length if self.length else 0.0

    @property
    def mean_depth(self) -> float:
        return float(self.depth.mean()) if self.length else 0.0

    @property
    def mean_depth_covered(self) -> float:
        cov = self.depth[self.depth > 0]
        return float(cov.mean()) if cov.size else 0.0

    @property
    def max_depth(self) -> int:
        return int(self.depth.max()) if self.length else 0

    @property
    def evenness_cv(self) -> float:
        """Coefficient of variation of depth over COVERED positions.
        Lower = more even (genuine genome-wide signal); high = piled on a
        few loci (often cross-mapping or a conserved region only)."""
        cov = self.depth[self.depth > 0]
        if cov.size < 2 or cov.mean() == 0:
            return 0.0
        return float(cov.std() / cov.mean())

    @property
    def largest_gap(self) -> int:
        """Longest run of zero-depth positions."""
        zero = self.depth == 0
        if not zero.any():
            return 0
        # run-length over the boolean mask
        best = run = 0
        for z in zero:
            run = run + 1 if z else 0
            if run > best:
                best = run
        return best


def _covered_ranges(r_start: int, cigar):
    """Yield (start, end) reference intervals carrying a base for one alignment."""
    pos = r_start
    for length, op in cigar:
        if op in _COVER_OPS:
            yield pos, pos + length
            pos += length
        elif op in _SKIP_OPS:
            pos += length
        # I/S/H/P consume query/clip only -> no reference advance


def coverage_from_hits(ref_lengths: dict[str, int], hits) -> dict[str, RefCoverage]:
    """Build RefCoverage per reference from an iterable of hit dicts.

    Each hit must provide: ctg (ref name), r_st (0-based start), cigar
    (list of (length, op)). Only hits whose ctg is in ref_lengths are used.
    """
    depths = {name: np.zeros(L, dtype=np.int32) for name, L in ref_lengths.items()}
    for h in hits:
        ctg = h["ctg"]
        if ctg not in depths:
            continue
        arr = depths[ctg]
        for s, e in _covered_ranges(h["r_st"], h["cigar"]):
            s = max(0, s)
            e = min(arr.size, e)
            if e > s:
                arr[s:e] += 1
    return {name: RefCoverage(name, ref_lengths[name], depths[name]) for name in ref_lengths}
