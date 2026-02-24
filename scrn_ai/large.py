#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

Large-scale trajectory workflow using STAVIA (VIA 2.0).

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging

import anndata as ad
import pathlib as p
import scanpy as sc

logger = logging.getLogger(__name__)


def run(infile, root_cell, outfile):
    adata = ad.read_h5ad(infile)
    logger.info("Loaded %d cells × %d genes", adata.n_obs, adata.n_vars)

    try:
        from pyVIA.core import VIA
    except ImportError:
        raise ImportError(
            "pyVIA is required for large-scale trajectory analysis. "
            "Install with: pip install pyVIA"
        )

    # Ensure PCA is computed
    if 'X_pca' not in adata.obsm:
        logger.info("Computing PCA")
        sc.pp.pca(adata, n_comps=50)

    logger.info("Running VIA (root_cell=%s)", root_cell)
    via = VIA(
        data=adata.obsm['X_pca'],
        true_label=adata.obs['leiden'].values if 'leiden' in adata.obs else None,
        root_user=root_cell,
        jac_std_global=0.15,
    )
    via.run_VIA()

    # Store outputs
    adata.obs["via_pseudotime"] = via.single_cell_pt_markov
    adata.obs["via_cluster"] = via.labels
    _save(adata, outfile)
    logger.info("Saved results to %s", outfile)

def _save(adata, out):
    out = p.Path(out)
    if out.suffix == ".h5ad":
        adata.write(out)
    else:
        adata.obs.to_csv(out, sep="\t")

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
