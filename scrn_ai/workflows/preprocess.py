#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

QC preprocessing with multi-format input support.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging

import scanpy as sc
import anndata as ad
import pathlib as p
import numpy as np

logger = logging.getLogger(__name__)


def run(input_file, output_file, min_genes=200, min_cells=3, max_genes=None, max_mito_pct=None):
    """
    Preprocess raw scRNA-seq data with quality control filtering.
    
    Parameters
    ----------
    input_file : str
        Path to input file (.mtx, .h5ad, .loom, .csv)
    output_file : str
        Path to output .h5ad file
    min_genes : int
        Minimum number of genes per cell (default: 200)
    min_cells : int
        Minimum number of cells per gene (default: 3)
    max_genes : int, optional
        Maximum number of genes per cell (filter doublets)
    max_mito_pct : float, optional
        Maximum mitochondrial percentage (e.g., 20.0 for 20%)
    """
    logger.info("Loading data from %s", input_file)
    adata = _read_any(input_file)
    
    n_cells_start = adata.n_obs
    n_genes_start = adata.n_vars
    
    logger.info("Initial data: %d cells × %d genes", n_cells_start, n_genes_start)
    
    # Calculate QC metrics
    logger.info("Calculating QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # mitochondrial genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Filter cells by gene count
    logger.info("Filtering cells: min_genes=%d, max_genes=%s", min_genes, max_genes)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_genes is not None:
        adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    
    # Filter cells by mitochondrial content
    if max_mito_pct is not None:
        logger.info("Filtering cells: max_mito_pct=%.1f%%", max_mito_pct)
        adata = adata[adata.obs.pct_counts_mt < max_mito_pct, :]
    
    # Filter genes by cell count
    logger.info("Filtering genes: min_cells=%d", min_cells)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    n_cells_end = adata.n_obs
    n_genes_end = adata.n_vars
    
    logger.info("Final data: %d cells × %d genes", n_cells_end, n_genes_end)
    logger.info("Removed %d cells (%.1f%%)", n_cells_start - n_cells_end, 100*(n_cells_start - n_cells_end)/n_cells_start)
    logger.info("Removed %d genes (%.1f%%)", n_genes_start - n_genes_end, 100*(n_genes_start - n_genes_end)/n_genes_start)
    
    # Save preprocessed data
    logger.info("Saving to %s", output_file)
    _save(adata, output_file)
    logger.info("Preprocessing complete")


def _read_any(path):
    """
    Read various scRNA-seq file formats into AnnData.
    
    Supports: .h5ad, .h5 (10X), .loom, .mtx (10X directory), .csv
    """
    path = p.Path(path)
    
    if path.suffix == ".h5ad":
        return ad.read_h5ad(path)
    elif path.suffix == ".h5":
        return sc.read_10x_h5(path)
    elif path.suffix == ".loom":
        return ad.read_loom(path)
    elif path.suffix == ".csv":
        import pandas as pd
        df = pd.read_csv(path, index_col=0)
        return ad.AnnData(df)
    elif path.is_dir():
        # Assume 10X mtx directory
        return sc.read_10x_mtx(path, gex_only=True)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")


def _save(adata, output_file):
    """Save AnnData to .h5ad file."""
    output_path = p.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
