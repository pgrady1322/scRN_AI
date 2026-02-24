"""
Preprocessing workflow for scRNA-seq data.

Handles QC filtering, multiple input formats, and basic quality control.
"""
import scanpy as sc
import anndata as ad
import pathlib as p
import numpy as np


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
    print(f"[preprocess] Loading data from {input_file}")
    adata = _read_any(input_file)
    
    n_cells_start = adata.n_obs
    n_genes_start = adata.n_vars
    
    print(f"[preprocess] Initial data: {n_cells_start} cells × {n_genes_start} genes")
    
    # Calculate QC metrics
    print("[preprocess] Calculating QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # mitochondrial genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Filter cells by gene count
    print(f"[preprocess] Filtering cells: min_genes={min_genes}, max_genes={max_genes}")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_genes is not None:
        adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    
    # Filter cells by mitochondrial content
    if max_mito_pct is not None:
        print(f"[preprocess] Filtering cells: max_mito_pct={max_mito_pct}%")
        adata = adata[adata.obs.pct_counts_mt < max_mito_pct, :]
    
    # Filter genes by cell count
    print(f"[preprocess] Filtering genes: min_cells={min_cells}")
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    n_cells_end = adata.n_obs
    n_genes_end = adata.n_vars
    
    print(f"[preprocess] Final data: {n_cells_end} cells × {n_genes_end} genes")
    print(f"[preprocess] Removed {n_cells_start - n_cells_end} cells ({100*(n_cells_start - n_cells_end)/n_cells_start:.1f}%)")
    print(f"[preprocess] Removed {n_genes_start - n_genes_end} genes ({100*(n_genes_start - n_genes_end)/n_genes_start:.1f}%)")
    
    # Save preprocessed data
    print(f"[preprocess] Saving to {output_file}")
    _save(adata, output_file)
    print("[preprocess] Complete!")


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
