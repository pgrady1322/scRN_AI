#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v0.1.0

UMAP embedding generation and visualization.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import logging

import scanpy as sc
import anndata as ad
import pathlib as p
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def run_umap(input_file, output_file, color_by="leiden", n_neighbors=15, min_dist=0.1, cell_types=None):
    """
    Generate UMAP visualization for dimensional reduction.
    
    Parameters
    ----------
    input_file : str
        Path to input normalized .h5ad file
    output_file : str
        Path to output image file (.png, .pdf, etc.)
    color_by : str
        Observation key to color by (default: leiden)
    n_neighbors : int
        Number of neighbors for UMAP (default: 15)
    min_dist : float
        Minimum distance for UMAP (default: 0.1)
    cell_types : str, optional
        Path to CSV file with cell type annotations
    """
    logger.info("Loading data from %s", input_file)
    adata = ad.read_h5ad(input_file)
    
    # Load cell types if provided
    if cell_types is not None:
        logger.info("Loading cell type annotations from %s", cell_types)
        import pandas as pd
        ct_df = pd.read_csv(cell_types, index_col=0)
        # Merge cell types into adata.obs
        if 'cell_type' in ct_df.columns:
            adata.obs['cell_type'] = ct_df['cell_type']
            color_by = 'cell_type'
        elif 'annotation' in ct_df.columns:
            adata.obs['cell_type'] = ct_df['annotation']
            color_by = 'cell_type'
    
    # Compute neighbors if not already done
    if 'neighbors' not in adata.uns:
        logger.info("Computing PCA...")
        if 'X_pca' not in adata.obsm:
            sc.pp.pca(adata, n_comps=50)
        
        logger.info("Computing neighbors (n_neighbors=%d)...", n_neighbors)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    
    # Compute UMAP if not already done
    if 'X_umap' not in adata.obsm:
        logger.info("Computing UMAP (min_dist=%.2f)...", min_dist)
        sc.tl.umap(adata, min_dist=min_dist)
    
    # Compute leiden clustering if coloring by leiden and not present
    if color_by == 'leiden' and 'leiden' not in adata.obs:
        logger.info("Computing Leiden clustering...")
        sc.tl.leiden(adata)
    
    # Generate plot
    logger.info("Generating plot colored by '%s'...", color_by)
    fig = sc.pl.umap(adata, color=color_by, show=False, return_fig=True)
    
    # Save figure
    output_path = p.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Saving to %s", output_file)
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    logger.info("UMAP visualization complete")

# scRN_AI v0.1.0
# Any usage is subject to this software's license.
