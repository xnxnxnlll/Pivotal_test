"""Correctness tests for the alignment + coverage core (no external data)."""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from zephyr.align import build_aligner, align_reads
from zephyr.coverage import coverage_from_hits
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "make_synthetic", os.path.join(os.path.dirname(__file__), "make_synthetic.py"))
ms = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ms)


def _setup():
    d = tempfile.mkdtemp()
    reads = os.path.join(d, "reads.fasta")
    facts = ms.make_case(d, reads)
    aligner = build_aligner(os.path.join(d, "panel.fasta"))
    hits = list(align_reads(aligner, reads))
    cov = coverage_from_hits(facts["ref_lengths"], (h.as_cov_hit() for h in hits))
    return facts, hits, cov


def test_reads_align():
    _facts, hits, _cov = _setup()
    assert len(hits) == 3
    assert all(h.ctg == "virusX" for h in hits)
    assert all(h.identity > 0.85 for h in hits)  # ONT-like 5% error


def test_breadth_and_depth():
    facts, _hits, cov = _setup()
    vx = cov["virusX"]
    assert abs(vx.covered_bp - facts["virusX_covered_bp"]) <= 4
    assert vx.max_depth == 2
    assert 0.30 < vx.breadth < 0.34


def test_offpanel_virus_has_zero_coverage():
    _facts, _hits, cov = _setup()
    assert cov["virusY"].covered_bp == 0
    assert cov["virusY"].breadth == 0.0


if __name__ == "__main__":
    test_reads_align()
    test_breadth_and_depth()
    test_offpanel_virus_has_zero_coverage()
    print("all tests passed")
