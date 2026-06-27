from __future__ import annotations
import os
import pandas as pd

from . import qc, discovery, confidence
from .align import build_aligner, align_reads
from .coverage import coverage_from_hits


def load_panel(lengths_tsv: str):
    """Return (ref_lengths {contig: len}, taxon_map {contig: taxon})."""
    df = pd.read_csv(lengths_tsv, sep="\t")
    ref_lengths = dict(zip(df["contig_label"], df["length"].astype(int)))
    taxon_map = dict(zip(df["contig_label"], df["taxon"]))
    return ref_lengths, taxon_map


def process_pool(pool, fasta_path, aligner, ref_lengths, taxon_map, cfg):
    """Run QC + discovery + alignment + coverage for one pool.

    Returns dict with: qc_row, aln_df (per-contig, with coverage + concordance),
    cov (per-contig RefCoverage for plotting), kraken_df, sourmash_df.
    """
    out = cfg["out_dir"]
    qc_row = qc.qc_fasta(pool, fasta_path).to_row()

    kraken_df = discovery.run_kraken2(
        fasta_path, cfg.get("kraken_db", ""), os.path.join(out, "discovery"),
        pool, confidence=cfg.get("kraken_confidence", 0.10),
        threads=cfg.get("threads", 4))
    sourmash_df = discovery.run_sourmash_gather(
        fasta_path, cfg.get("sourmash_db", ""), os.path.join(out, "discovery"),
        pool, ksize=cfg.get("sourmash_k", 31), scaled=cfg.get("sourmash_scaled", 1000))

    hits = list(align_reads(aligner, fasta_path, min_mapq=cfg.get("min_mapq", 0)))
    cov = coverage_from_hits(ref_lengths, (h.as_cov_hit() for h in hits))

    aln_df = confidence.summarise_alignment(pool, hits)
    if not aln_df.empty:
        aln_df["taxon_group"] = aln_df["taxon"].map(taxon_map).fillna(aln_df["taxon"])
        aln_df = confidence.attach_coverage(aln_df, {pool: cov})
        aln_df = confidence.add_concordance(aln_df, kraken_df, sourmash_df,
                                            name_map=taxon_map)
    return {"qc_row": qc_row, "aln_df": aln_df, "cov": cov,
            "kraken_df": kraken_df, "sourmash_df": sourmash_df}
