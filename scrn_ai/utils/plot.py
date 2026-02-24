#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

Visualization utilities — QC, expression, pseudotime, and UMAP plots.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import logging
from pathlib import Path
from typing import List, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)


# ── UMAP helpers ─────────────────────────────────────────────────────────

def umap_pseudotime(
    adata: AnnData,
    key: str = "pseudotime",
    save: Optional[str] = None,
) -> None:
    """Plot UMAP coloured by pseudotime.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with an ``X_umap`` embedding.
    key : str
        Column in ``adata.obs`` containing pseudotime values.
    save : str, optional
        Path to save the figure.  If *None*, returns without saving.
    """
    if "X_umap" not in adata.obsm:
        logger.info("X_umap not found — computing UMAP")
        sc.tl.umap(adata)
    fig = sc.pl.umap(adata, color=key, cmap="viridis", show=False, return_fig=True)
    if save:
        _save_fig(fig, save)
    plt.close(fig)


# ── QC plots ─────────────────────────────────────────────────────────────

def qc_violin(
    adata: AnnData,
    keys: Sequence[str] = ("n_genes_by_counts", "total_counts", "pct_counts_mt"),
    groupby: Optional[str] = None,
    save: Optional[str] = None,
) -> None:
    """Violin plots of QC metrics.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with QC columns in ``.obs``.
    keys : sequence of str
        Observation keys to plot.
    groupby : str, optional
        Observation key to group violins by (e.g. sample or batch).
    save : str, optional
        Path to save the figure.
    """
    present = [k for k in keys if k in adata.obs.columns]
    if not present:
        logger.warning("No QC metrics found in adata.obs — run sc.pp.calculate_qc_metrics first")
        return
    fig, axes = plt.subplots(1, len(present), figsize=(4 * len(present), 4))
    if len(present) == 1:
        axes = [axes]
    for ax, key in zip(axes, present):
        sc.pl.violin(adata, key, groupby=groupby, ax=ax, show=False)
        ax.set_title(key)
    fig.tight_layout()
    if save:
        _save_fig(fig, save)
    plt.close(fig)


def qc_scatter(
    adata: AnnData,
    x: str = "total_counts",
    y: str = "n_genes_by_counts",
    color: str = "pct_counts_mt",
    save: Optional[str] = None,
) -> None:
    """Scatter plot of two QC metrics coloured by a third.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    x, y : str
        Keys in ``adata.obs`` for the x and y axes.
    color : str
        Key in ``adata.obs`` for point colouring.
    save : str, optional
        Path to save the figure.
    """
    needed = {x, y, color}
    missing = needed - set(adata.obs.columns)
    if missing:
        logger.warning("Missing obs columns for scatter: %s", missing)
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(
        adata.obs[x], adata.obs[y],
        c=adata.obs[color], cmap="RdYlGn_r", s=3, alpha=0.6,
    )
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    fig.colorbar(scatter, ax=ax, label=color)
    fig.tight_layout()
    if save:
        _save_fig(fig, save)
    plt.close(fig)


# ── Gene expression plots ───────────────────────────────────────────────

def dotplot(
    adata: AnnData,
    genes: List[str],
    groupby: str = "leiden",
    save: Optional[str] = None,
) -> None:
    """Dot plot of gene expression across clusters.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    genes : list of str
        Gene names to display.
    groupby : str
        Observation key for grouping.
    save : str, optional
        Path to save the figure.
    """
    fig = sc.pl.dotplot(adata, var_names=genes, groupby=groupby, show=False, return_fig=True)
    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
    plt.close("all")


def stacked_violin(
    adata: AnnData,
    genes: List[str],
    groupby: str = "leiden",
    save: Optional[str] = None,
) -> None:
    """Stacked violin plot of gene expression across clusters.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    genes : list of str
        Gene names to display.
    groupby : str
        Observation key for grouping.
    save : str, optional
        Path to save the figure.
    """
    fig = sc.pl.stacked_violin(
        adata, var_names=genes, groupby=groupby, show=False, return_fig=True,
    )
    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
    plt.close("all")


# ── Pseudotime plots ────────────────────────────────────────────────────

def pseudotime_heatmap(
    adata: AnnData,
    genes: List[str],
    time_key: str = "pseudotime",
    n_bins: int = 50,
    save: Optional[str] = None,
) -> None:
    """Heatmap of gene expression ordered by pseudotime.

    Cells are binned along the pseudotime axis and mean expression is shown
    per gene per bin.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    genes : list of str
        Genes to show on the y-axis.
    time_key : str
        Observation column with pseudotime values.
    n_bins : int
        Number of pseudotime bins.
    save : str, optional
        Path to save the figure.
    """
    if time_key not in adata.obs.columns:
        logger.warning("'%s' not found in adata.obs — cannot plot heatmap", time_key)
        return

    valid_genes = [g for g in genes if g in adata.var_names]
    if not valid_genes:
        logger.warning("No requested genes found in adata.var_names")
        return

    order = np.argsort(adata.obs[time_key].values)
    bins = np.array_split(order, n_bins)

    from scipy import sparse
    X = adata[:, valid_genes].X
    if sparse.issparse(X):
        X = X.toarray()

    mean_expr = np.vstack([X[b].mean(axis=0) for b in bins]).T  # genes × bins

    fig, ax = plt.subplots(figsize=(10, max(3, 0.3 * len(valid_genes))))
    im = ax.imshow(mean_expr, aspect="auto", cmap="magma", interpolation="nearest")
    ax.set_yticks(range(len(valid_genes)))
    ax.set_yticklabels(valid_genes, fontsize=8)
    ax.set_xlabel(f"Pseudotime bins ({time_key})")
    fig.colorbar(im, ax=ax, label="Mean expression")
    fig.tight_layout()
    if save:
        _save_fig(fig, save)
    plt.close(fig)


# ── Internal helpers ─────────────────────────────────────────────────────

def _save_fig(fig: plt.Figure, path: str) -> None:
    """Save a matplotlib figure, creating parent directories as needed."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    logger.info("Figure saved to %s", out)

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
