"""Discovery classification: kraken2 (+bracken) and sourmash gather.

These two have different failure modes, so agreement between them (and with
the minimap2 confirmation step) is itself a confidence signal. Both are
invoked as subprocesses because neither has a maintained pip API; install
them via conda/bioconda on your server (see environment.yml). Each wrapper
degrades gracefully: if the tool or DB is missing it returns an empty result
and the pipeline continues with alignment-only evidence.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import pandas as pd


def _have(tool: str) -> bool:
    return shutil.which(tool) is not None


# ---------------------------------------------------------------- kraken2 ---
def run_kraken2(fasta_path: str, db: str, out_dir: str, pool: str,
                confidence: float = 0.10, threads: int = 4) -> pd.DataFrame:
    """Run kraken2 and parse the report into a tidy per-taxon table.

    `confidence` is raised above the default 0 because noisy long reads
    over-assign at species level; 0.05-0.15 is a sensible ONT range.
    Returns columns: pool, taxid, rank, name, kraken_reads (clade-level).
    """
    if not _have("kraken2") or not db or not os.path.exists(db):
        return pd.DataFrame(columns=["pool", "taxid", "rank", "name", "kraken_reads"])

    os.makedirs(out_dir, exist_ok=True)
    report = os.path.join(out_dir, f"{pool}.kraken2.report.tsv")
    out = os.path.join(out_dir, f"{pool}.kraken2.out")
    subprocess.run(
        ["kraken2", "--db", db, "--threads", str(threads),
         "--confidence", str(confidence), "--report", report,
         "--output", out, fasta_path],
        check=True, capture_output=True, text=True,
    )
    rows = []
    with open(report) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            clade_reads, _direct, rank, taxid, name = int(f[1]), f[2], f[3], f[4], f[5].strip()
            if clade_reads > 0 and rank.startswith(("S", "G", "F")):  # species/genus/family
                rows.append({"pool": pool, "taxid": taxid, "rank": rank,
                             "name": name, "kraken_reads": clade_reads})
    return pd.DataFrame(rows)


# --------------------------------------------------------------- sourmash ---
def run_sourmash_gather(fasta_path: str, sig_db: str, out_dir: str, pool: str,
                        ksize: int = 31, scaled: int = 1000) -> pd.DataFrame:
    """Sketch the reads and run `sourmash gather` against a viral signature DB.

    `gather` reports the minimum set of references that explain the k-mer
    content, with a containment-based match — a robust, error-tolerant
    second opinion. Returns: pool, name, containment, sm_bp.
    """
    if not _have("sourmash") or not sig_db or not os.path.exists(sig_db):
        return pd.DataFrame(columns=["pool", "name", "containment", "sm_bp"])

    os.makedirs(out_dir, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        sig = os.path.join(td, f"{pool}.sig")
        subprocess.run(
            ["sourmash", "sketch", "dna", "-p", f"k={ksize},scaled={scaled}",
             "-o", sig, fasta_path],
            check=True, capture_output=True, text=True,
        )
        csv = os.path.join(out_dir, f"{pool}.gather.csv")
        proc = subprocess.run(
            ["sourmash", "gather", sig, sig_db, "-o", csv, "--threshold-bp", "0"],
            capture_output=True, text=True,
        )
        if proc.returncode != 0 or not os.path.exists(csv):
            return pd.DataFrame(columns=["pool", "name", "containment", "sm_bp"])

    g = pd.read_csv(csv)
    if g.empty:
        return pd.DataFrame(columns=["pool", "name", "containment", "sm_bp"])
    g = g.rename(columns={"name": "name",
                          "f_match": "containment",
                          "unique_intersect_bp": "sm_bp"})
    g["pool"] = pool
    return g[["pool", "name", "containment", "sm_bp"]]
