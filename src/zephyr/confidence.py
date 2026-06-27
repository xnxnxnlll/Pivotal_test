from __future__ import annotations
import pandas as pd


def summarise_alignment(pool: str, hits) -> pd.DataFrame:
    """Per-taxon aggregate of alignment hits for one pool."""
    rows = []
    by_taxon: dict[str, list] = {}
    for h in hits:
        by_taxon.setdefault(h.ctg, []).append(h)
    for taxon, hs in by_taxon.items():
        n = len(hs)
        rows.append({
            "pool": pool,
            "taxon": taxon,
            "aln_reads": n,
            "mean_identity": round(sum(h.identity for h in hs) / n, 4),
            "mean_mapq": round(sum(h.mapq for h in hs) / n, 1),
            "aligned_bp": sum(h.r_en - h.r_st for h in hs),
        })
    return pd.DataFrame(rows)


def attach_coverage(aln_df: pd.DataFrame, cov_by_pool: dict) -> pd.DataFrame:
    """Add breadth / depth / evenness columns from coverage objects.

    cov_by_pool: {pool: {taxon: RefCoverage}}
    """
    out = aln_df.copy()
    for col in ["ref_len", "breadth", "mean_depth", "mean_depth_covered",
                "max_depth", "evenness_cv", "largest_gap"]:
        out[col] = 0.0
    for i, r in out.iterrows():
        cov = cov_by_pool.get(r["pool"], {}).get(r["taxon"])
        if cov is None:
            continue
        out.at[i, "ref_len"] = cov.length
        out.at[i, "breadth"] = round(cov.breadth, 4)
        out.at[i, "mean_depth"] = round(cov.mean_depth, 3)
        out.at[i, "mean_depth_covered"] = round(cov.mean_depth_covered, 3)
        out.at[i, "max_depth"] = cov.max_depth
        out.at[i, "evenness_cv"] = round(cov.evenness_cv, 3)
        out.at[i, "largest_gap"] = cov.largest_gap
    return out


def flag_cross_contamination(df: pd.DataFrame, run_key: str = "run",
                             ratio: float = 0.01) -> pd.DataFrame:
    """Flag a detection as possible index-hop if, within the same run, the
    same taxon is also seen at >> higher read count elsewhere.

    A detection is flagged when its aln_reads is below `ratio` * (max
    aln_reads for that taxon in the same run) AND that max is substantially
    larger (>=100 reads). `run_key` defaults to a 'run' column you supply
    (e.g. date+site); if absent, no flags are raised.
    """
    out = df.copy()
    out["contam_flag"] = False
    if run_key not in out.columns:
        return out
    for (run, taxon), grp in out.groupby([run_key, "taxon"]):
        mx = grp["aln_reads"].max()
        if mx < 100:
            continue
        small = grp["aln_reads"] < ratio * mx
        out.loc[grp.index[small], "contam_flag"] = True
    return out


def assign_confidence(df: pd.DataFrame,
                      min_high_reads: int = 10,
                      min_high_breadth: float = 0.10,
                      min_identity: float = 0.85,
                      concordance_col: str | None = "concordant") -> pd.DataFrame:
    """Assign HIGH / MEDIUM / LOW per row from the evidence axes.

    Defaults are deliberately conservative and documented in the README;
    treat them as a starting point, not a standard.
    """
    out = df.copy()

    def tier(r) -> str:
        if r.get("contam_flag", False):
            return "LOW"
        good_aln = r["mean_identity"] >= min_identity
        concordant = bool(r.get(concordance_col, False)) if concordance_col else False
        if (r["aln_reads"] >= min_high_reads and r["breadth"] >= min_high_breadth
                and good_aln and (concordant or r["aln_reads"] >= 3 * min_high_reads)):
            return "HIGH"
        if r["aln_reads"] >= 3 and good_aln and (r["breadth"] >= 0.02 or concordant):
            return "MEDIUM"
        return "LOW"

    out["confidence"] = out.apply(tier, axis=1)
    return out


def add_concordance(aln_df: pd.DataFrame, kraken_df: pd.DataFrame,
                    sourmash_df: pd.DataFrame, name_map=None) -> pd.DataFrame:
    """Mark whether an alignment-confirmed taxon was also seen by kraken2 or
    sourmash. `name_map` optionally maps panel reference names -> the taxon
    strings used by the discovery DBs (genus/species). Without a map we do a
    case-insensitive substring match, which is conservative.
    """
    out = aln_df.copy()
    out["in_kraken"] = False
    out["in_sourmash"] = False

    def norm(s: str) -> str:
        return str(s).lower()

    for i, r in out.iterrows():
        key = norm(name_map.get(r["taxon"], r["taxon"]) if name_map else r["taxon"])
        if not kraken_df.empty:
            k = kraken_df[kraken_df["pool"] == r["pool"]]
            out.at[i, "in_kraken"] = any(key in norm(n) or norm(n) in key
                                         for n in k["name"])
        if not sourmash_df.empty:
            s = sourmash_df[sourmash_df["pool"] == r["pool"]]
            out.at[i, "in_sourmash"] = any(key in norm(n) or norm(n) in key
                                           for n in s["name"])
    out["concordant"] = out["in_kraken"] | out["in_sourmash"]
    return out
