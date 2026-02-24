#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v0.1.0

Multi-method normalization workflow (Seurat, JMP, basic).

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import scanpy as sc
import anndata as ad
import pathlib as p
import numpy as np


def run(input_file, output_file, method="seurat", algorithm="LogNormalize", scale_factor=10000):
    """
    Normalize count data using various methods.
    
    Parameters
    ----------
    input_file : str
        Path to input .h5ad file
    output_file : str
        Path to output .h5ad file
    method : str
        Normalization method: seurat, jmp, log1p, scran, sctransform
    algorithm : str
        Specific algorithm within method:
        - seurat: LogNormalize, SCTransform
        - jmp: TMM, RLE, UpperQuartile
    scale_factor : float
        Scaling factor for normalization (default: 10000)
    """
    print(f"[normalize] Loading data from {input_file}")
    adata = ad.read_h5ad(input_file)
    
    print(f"[normalize] Method: {method}, Algorithm: {algorithm}")
    
    method = method.lower()
    algorithm_lower = algorithm.lower()
    
    # Route to appropriate normalization method
    if method == "seurat":
        _normalize_seurat(adata, algorithm_lower, scale_factor)
    elif method == "jmp":
        _normalize_jmp(adata, algorithm_lower)
    elif method == "log1p":
        _normalize_log1p(adata, scale_factor)
    elif method == "scran":
        _normalize_scran(adata)
    elif method == "sctransform":
        _normalize_sctransform(adata)
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    # Save normalized data
    print(f"[normalize] Saving to {output_file}")
    _save(adata, output_file)
    print("[normalize] Complete!")


def _normalize_seurat(adata, algorithm, scale_factor):
    """Normalize using Seurat methods (via R)."""
    print(f"[normalize] Running Seurat {algorithm} normalization...")
    
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        
        # Transfer data to R
        ro.globalenv["counts"] = adata.to_df() if adata.raw is None else adata.raw.to_adata().to_df()
        
        if algorithm == "lognormalize":
            ro.r(f'''
                suppressMessages(library(Seurat))
                seurat_obj <- CreateSeuratObject(counts = t(as.matrix(counts)))
                seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                                           scale.factor = {scale_factor})
                normalized_data <- GetAssayData(seurat_obj, slot = "data")
            ''')
        elif algorithm == "sctransform":
            ro.r('''
                suppressMessages(library(Seurat))
                seurat_obj <- CreateSeuratObject(counts = t(as.matrix(counts)))
                seurat_obj <- SCTransform(seurat_obj)
                normalized_data <- GetAssayData(seurat_obj, slot = "data")
            ''')
        else:
            raise ValueError(f"Unknown Seurat algorithm: {algorithm}")
        
        # Get normalized data back from R
        normalized = np.array(ro.r["normalized_data"]).T
        adata.X = normalized
        
    except ImportError:
        print("[normalize] Warning: rpy2 not available, falling back to scanpy methods")
        if algorithm == "lognormalize":
            _normalize_log1p(adata, scale_factor)
        elif algorithm == "sctransform":
            _normalize_sctransform(adata)


def _normalize_jmp(adata, algorithm):
    """Normalize using JMP/edgeR methods (TMM, RLE, UpperQuartile)."""
    print(f"[normalize] Running JMP {algorithm} normalization...")
    
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        
        # Transfer data to R
        ro.globalenv["counts"] = adata.to_df() if adata.raw is None else adata.raw.to_adata().to_df()
        
        if algorithm == "tmm":
            method_str = "TMM"
        elif algorithm == "rle":
            method_str = "RLE"
        elif algorithm == "upperquartile":
            method_str = "upperquartile"
        else:
            raise ValueError(f"Unknown JMP algorithm: {algorithm}")
        
        ro.r(f'''
            suppressMessages(library(edgeR))
            dge <- DGEList(counts = t(as.matrix(counts)))
            dge <- calcNormFactors(dge, method = "{method_str}")
            # Apply normalization factors
            normalized_data <- cpm(dge, log = TRUE)
        ''')
        
        # Get normalized data back from R
        normalized = np.array(ro.r["normalized_data"]).T
        adata.X = normalized
        
    except ImportError:
        print("[normalize] Warning: rpy2/edgeR not available, falling back to log1p")
        _normalize_log1p(adata, 1e6)


def _normalize_log1p(adata, scale_factor):
    """Simple log(x+1) normalization."""
    print(f"[normalize] Running log1p normalization (scale_factor={scale_factor})...")
    sc.pp.normalize_total(adata, target_sum=scale_factor)
    sc.pp.log1p(adata)


def _normalize_scran(adata):
    """Deconvolution-based normalization using scran."""
    print("[normalize] Running scran normalization...")
    
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        
        ro.globalenv["counts"] = adata.to_df() if adata.raw is None else adata.raw.to_adata().to_df()
        
        ro.r('''
            suppressMessages(library(scran))
            sce <- SingleCellExperiment(assays = list(counts = t(as.matrix(counts))))
            clusters <- quickCluster(sce)
            sce <- computeSumFactors(sce, clusters = clusters)
            sce <- logNormCounts(sce)
            normalized_data <- logcounts(sce)
        ''')
        
        normalized = np.array(ro.r["normalized_data"]).T
        adata.X = normalized
        
    except ImportError:
        print("[normalize] Warning: rpy2/scran not available, falling back to log1p")
        _normalize_log1p(adata, 1e4)


def _normalize_sctransform(adata):
    """Variance-stabilizing transformation."""
    print("[normalize] Running sctransform normalization...")
    
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        
        ro.globalenv["counts"] = adata.to_df() if adata.raw is None else adata.raw.to_adata().to_df()
        
        ro.r('''
            suppressMessages(library(sctransform))
            vst_out <- vst(as.matrix(counts))
            normalized_data <- vst_out$y
        ''')
        
        normalized = np.array(ro.r["normalized_data"])
        adata.X = normalized
        
    except ImportError:
        print("[normalize] Warning: rpy2/sctransform not available, using scanpy approximation")
        # Use Scanpy's regress_out as approximation
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.regress_out(adata, ['total_counts'])


def _save(adata, output_file):
    """Save AnnData to .h5ad file."""
    output_path = p.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

# scRN_AI v0.1.0
# Any usage is subject to this software's license.
