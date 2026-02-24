#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

Unified pseudotime trajectory analysis (DPT, diffusion, BLTSA, VIA).

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
import os

import scanpy as sc
import anndata as ad
import pathlib as p
import numpy as np

logger = logging.getLogger(__name__)


def run(input_file, output_file, method="dpt", scale="small", root_cell=None):
    """
    Perform pseudotime trajectory analysis.
    
    Parameters
    ----------
    input_file : str
        Path to input normalized .h5ad file
    output_file : str
        Path to output .h5ad file or directory
    method : str
        Pseudotime method: dpt, diffusion, bltsa, via
    scale : str
        Dataset scale: small (<50k cells) or large (>50k cells)
    root_cell : str, optional
        Root cell ID for pseudotime calculation
    """
    logger.info("Loading data from %s", input_file)
    adata = ad.read_h5ad(input_file)
    
    logger.info("Method: %s, Scale: %s", method, scale)
    logger.info("Data: %d cells × %d genes", adata.n_obs, adata.n_vars)
    
    method = method.lower()
    scale = scale.lower()
    
    # Route to appropriate analysis based on scale and method
    if scale == "small":
        _run_small_scale(adata, method, root_cell)
    elif scale == "large":
        _run_large_scale(adata, method, root_cell)
    else:
        raise ValueError(f"Unknown scale: {scale}. Use 'small' or 'large'.")
    
    # Save results
    logger.info("Saving to %s", output_file)
    _save(adata, output_file)
    logger.info("Pseudotime analysis complete")


def _run_small_scale(adata, method, root_cell=None):
    """Small-scale pseudotime analysis (<50k cells)."""
    
    # Ensure preprocessing is done
    if 'X_pca' not in adata.obsm:
        logger.info("Computing PCA...")
        sc.pp.pca(adata, n_comps=50)
    
    if 'neighbors' not in adata.uns:
        logger.info("Computing neighbors...")
        sc.pp.neighbors(adata)
    
    if method == "dpt":
        logger.info("Running Diffusion Pseudotime (DPT)...")
        
        # Compute diffusion map
        sc.tl.diffmap(adata)
        
        # Set root cell if provided
        if root_cell is not None:
            if root_cell in adata.obs_names:
                adata.uns['iroot'] = np.where(adata.obs_names == root_cell)[0][0]
            else:
                logger.warning("root_cell '%s' not found, using automatic selection", root_cell)
        
        # Compute DPT
        sc.tl.dpt(adata)
        adata.obs['pseudotime'] = adata.obs['dpt_pseudotime']
        
    elif method == "diffusion":
        logger.info("Running Diffusion Maps...")
        sc.tl.diffmap(adata)
        adata.obs['pseudotime'] = adata.obsm['X_diffmap'][:, 0]  # First diffusion component
        
    elif method == "bltsa":
        logger.info("Running BLTSA trajectory inference...")
        try:
            # BLTSA requires R and specific packages
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            from scipy import sparse
            pandas2ri.activate()
            
            # Transfer data to R
            ro.globalenv["expression_matrix"] = adata.X if not sparse.issparse(adata.X) else adata.X.toarray()
            
            bltsa_path = os.environ.get("BLTSA_PATH", "/opt/BLTSA/BLTSA.R")
            logger.info("Loading BLTSA from %s", bltsa_path)
            
            ro.r(f'''
                source("{bltsa_path}")
                result <- BLTSA_trajectory(expression_matrix)
                pseudotime <- result$pseudotime
            ''')
            
            adata.obs['pseudotime'] = np.array(ro.r['pseudotime'])
            adata.obs['bltsa_pseudotime'] = adata.obs['pseudotime']
            
        except Exception as e:
            logger.warning("BLTSA failed (%s), falling back to DPT", e)
            _run_small_scale(adata, "dpt", root_cell)
    
    elif method == "via":
        logger.warning("VIA is designed for large-scale data. Consider using scale='large'")
        _run_large_scale(adata, "via", root_cell)
    
    else:
        raise ValueError(f"Unknown method for small-scale: {method}")


def _run_large_scale(adata, method, root_cell=None):
    """Large-scale pseudotime analysis (>50k cells) using VIA/STAVIA."""
    
    logger.info("Running VIA/STAVIA for large-scale trajectory analysis...")
    
    try:
        from pyVIA.core import VIA
        
        # Ensure PCA is computed
        if 'X_pca' not in adata.obsm:
            logger.info("Computing PCA...")
            sc.pp.pca(adata, n_comps=50)
        
        # Run VIA
        logger.info("Initializing VIA (root_cell=%s)...", root_cell)
        v = VIA(
            data=adata.obsm['X_pca'],
            true_label=adata.obs['leiden'].values if 'leiden' in adata.obs else None,
            root_user=root_cell,
            jac_std_global=0.15
        )
        
        logger.info("Running VIA...")
        v.run_VIA()
        
        # Store results
        adata.obs['via_pseudotime'] = v.single_cell_pt_markov
        adata.obs['via_cluster'] = v.labels
        adata.obs['pseudotime'] = v.single_cell_pt_markov
        
        # Store additional VIA results
        adata.uns['via'] = {
            'root': root_cell,
            'terminal_clusters': v.terminal_clusters if hasattr(v, 'terminal_clusters') else []
        }
        
    except ImportError:
        logger.error("pyVIA not installed. Install with: pip install pyVIA")
        logger.warning("Falling back to DPT for large-scale analysis...")
        _run_small_scale(adata, "dpt", root_cell)


def _save(adata, output_file):
    """Save AnnData to .h5ad file."""
    output_path = p.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if output_path.suffix == ".h5ad":
        adata.write_h5ad(output_path)
    else:
        # Save as directory with multiple outputs
        output_path.mkdir(exist_ok=True)
        adata.write_h5ad(output_path / "results.h5ad")
        adata.obs[['pseudotime']].to_csv(output_path / "pseudotime.csv")
        logger.info("Results saved to directory: %s", output_path)

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
