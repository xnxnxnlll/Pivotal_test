from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_coverage(cov, pool: str, out_dir: str, window: int = 50):
    """Depth-along-genome for one RefCoverage. Returns the saved path."""
    os.makedirs(out_dir, exist_ok=True)
    d = cov.depth
    if d.size == 0:
        return None
    # smooth with a non-overlapping mean window for readability
    if d.size >= window:
        n = (d.size // window) * window
        sm = d[:n].reshape(-1, window).mean(axis=1)
        x = np.arange(sm.size) * window
    else:
        sm, x = d, np.arange(d.size)

    fig, ax = plt.subplots(figsize=(9, 2.6))
    ax.fill_between(x, sm, step="mid", alpha=0.7)
    ax.set_title(f"{pool} — {cov.name}  "
                 f"(breadth {cov.breadth:.0%}, mean {cov.mean_depth:.1f}x, "
                 f"max {cov.max_depth}x)", fontsize=9)
    ax.set_xlabel("genome position (bp)", fontsize=8)
    ax.set_ylabel("depth", fontsize=8)
    ax.margins(x=0)
    fig.tight_layout()
    path = os.path.join(out_dir, f"{pool}__{cov.name}.coverage.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def plot_detection_heatmap(master: "pd.DataFrame", out_path: str,
                           value: str = "aln_reads"):
    """viruses (rows) x pools (cols) heatmap of read support, annotated with
    confidence tier initials. `master` must have pool, taxon, value,
    confidence columns.
    """
    import pandas as pd
    if master.empty:
        return None
    mat = master.pivot_table(index="taxon", columns="pool", values=value,
                             aggfunc="sum", fill_value=0)
    tiers = master.pivot_table(index="taxon", columns="pool", values="confidence",
                               aggfunc="first")
    fig, ax = plt.subplots(figsize=(1.1 * mat.shape[1] + 3, 0.5 * mat.shape[0] + 2))
    logmat = np.log10(mat.replace(0, np.nan))
    im = ax.imshow(logmat, aspect="auto", cmap="viridis")
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=7)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            if v > 0:
                t = tiers.iloc[i, j]
                lab = {"HIGH": "H", "MEDIUM": "M", "LOW": "L"}.get(t, "")
                ax.text(j, i, f"{int(v)}\n{lab}", ha="center", va="center",
                        fontsize=6, color="white")
    cb = fig.colorbar(im, ax=ax, fraction=0.025)
    cb.set_label("log10 aligned reads", fontsize=8)
    ax.set_title("Respiratory virus detections (read support + confidence)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path
