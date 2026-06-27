#!/usr/bin/env python3
"""Process every pool FASTA in a directory and aggregate results.

Usage:
  python scripts/run_batch.py --input data/ --config config.yaml

Input files may be .fasta or .fasta.gz (mappy reads both). Pool names and a
run key (date+site, for cross-contamination flagging) are parsed from the
SecureBio filename convention YYMMDD-Site-NAS[-P#].respiratory.fasta[.gz].
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import yaml
import pandas as pd

from zephyr import confidence, plots
from zephyr.align import build_aligner
from zephyr.pipeline import load_panel, process_pool


def pool_name(path: str) -> str:
    b = os.path.basename(path)
    b = re.sub(r"\.respiratory", "", b)
    b = re.sub(r"\.(fasta|fa|fna)(\.gz)?$", "", b)
    return b


def run_key(pool: str) -> str:
    """date+site, e.g. 250107-BoDT-NAS-P1 -> 250107-BoDT (shared sequencing run proxy)."""
    parts = pool.split("-")
    return "-".join(parts[:2]) if len(parts) >= 2 else pool


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="directory of pool FASTAs")
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    os.makedirs(cfg["out_dir"], exist_ok=True)

    files = []
    for ext in ("*.fasta", "*.fa", "*.fna", "*.fasta.gz", "*.fa.gz", "*.fna.gz"):
        files += glob.glob(os.path.join(args.input, ext))
    files = sorted(set(files))
    if not files:
        sys.exit(f"No FASTAs found in {args.input}")
    print(f"Found {len(files)} pool(s)")

    ref_lengths, taxon_map = load_panel(cfg["panel_lengths"])
    print(f"Loading panel index ({len(ref_lengths)} contigs)...")
    aligner = build_aligner(cfg["panel_fasta"], n_threads=cfg.get("threads", 4))

    qc_rows, aln_frames, kr_frames, sm_frames = [], [], [], []
    cov_store = {}
    for f in files:
        pool = pool_name(f)
        print(f"  - {pool}")
        res = process_pool(pool, f, aligner, ref_lengths, taxon_map, cfg)
        qc_rows.append(res["qc_row"])
        cov_store[pool] = res["cov"]
        if not res["aln_df"].empty:
            res["aln_df"]["run"] = run_key(pool)
            aln_frames.append(res["aln_df"])
        if not res["kraken_df"].empty:
            kr_frames.append(res["kraken_df"])
        if not res["sourmash_df"].empty:
            sm_frames.append(res["sourmash_df"])

    out = cfg["out_dir"]
    pd.DataFrame(qc_rows).to_csv(os.path.join(out, "qc.tsv"), sep="\t", index=False)
    if kr_frames:
        pd.concat(kr_frames).to_csv(os.path.join(out, "kraken2.tsv"), sep="\t", index=False)
    if sm_frames:
        pd.concat(sm_frames).to_csv(os.path.join(out, "sourmash.tsv"), sep="\t", index=False)

    if not aln_frames:
        print("No alignment-confirmed detections across pools "
              "(check discovery tables / panel coverage).")
        return

    master = pd.concat(aln_frames, ignore_index=True)
    master = confidence.flag_cross_contamination(master, run_key="run",
                                                 ratio=cfg.get("contam_ratio", 0.01))
    master = confidence.assign_confidence(
        master,
        min_high_reads=cfg.get("min_high_reads", 10),
        min_high_breadth=cfg.get("min_high_breadth", 0.10),
        min_identity=cfg.get("min_identity", 0.85))

    # taxon-level master (group segments -> virus) for the headline table
    master.to_csv(os.path.join(out, "detections_by_contig.tsv"), sep="\t", index=False)
    taxon_tbl = (master.groupby(["pool", "taxon_group"])
                 .agg(aln_reads=("aln_reads", "sum"),
                      mean_identity=("mean_identity", "mean"),
                      best_breadth=("breadth", "max"),
                      mean_depth=("mean_depth", "mean"),
                      confidence=("confidence",
                                  lambda s: _best_tier(s)))
                 .reset_index()
                 .rename(columns={"taxon_group": "taxon"}))
    taxon_tbl.to_csv(os.path.join(out, "detections_by_taxon.tsv"), sep="\t", index=False)

    # figures
    fig_dir = os.path.join(out, "figures")
    for pool, cov in cov_store.items():
        for contig, c in cov.items():
            if c.covered_bp > 0:
                plots.plot_coverage(c, pool, os.path.join(fig_dir, "coverage"))
    plots.plot_detection_heatmap(taxon_tbl.assign(taxon=taxon_tbl["taxon"]),
                                 os.path.join(fig_dir, "detection_heatmap.png"))

    _print_answers(taxon_tbl)


def _best_tier(series):
    order = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
    return max(series, key=lambda t: order.get(t, 0))


def _print_answers(taxon_tbl):
    print("\n================ SUMMARY ================")
    print("Q1: What viruses are present, and how confident?\n")
    for tier in ("HIGH", "MEDIUM", "LOW"):
        sub = taxon_tbl[taxon_tbl["confidence"] == tier]
        if sub.empty:
            continue
        print(f"  [{tier}]")
        for _, r in sub.sort_values("aln_reads", ascending=False).iterrows():
            print(f"    {r['taxon']:<16} pool={r['pool']:<18} "
                  f"reads={int(r['aln_reads']):>6}  id={r['mean_identity']:.2f}  "
                  f"breadth={r['best_breadth']:.0%}")
    print("\nQ2: Genome coverage for common respiratory viruses")
    print("    see results/figures/coverage/*.png and detections_by_contig.tsv")
    print("    (breadth/mean_depth/evenness_cv per reference; per-segment for flu)")
    print("=========================================")


if __name__ == "__main__":
    main()
