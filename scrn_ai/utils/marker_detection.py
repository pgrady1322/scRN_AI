#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

Marker gene detection for cell type identification.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import logging

import numpy as np
import pandas as pd
from typing import List, Dict, Optional, Union
import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)


def identify_variable_genes(
    adata: AnnData,
    n_top_genes: int = 2000,
    flavor: str = "seurat_v3",
    subset: bool = False
) -> List[str]:
    """
    Identify highly variable genes in the dataset.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_top_genes : int
        Number of top variable genes to identify.
    flavor : str
        Method for identifying variable genes (seurat, seurat_v3, cell_ranger).
    subset : bool
        Whether to subset the data to variable genes.
    
    Returns
    -------
    list of str
        Names of highly variable genes.
    """
    # Calculate highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        subset=subset
    )
    
    # Extract gene names
    if 'highly_variable' in adata.var.columns:
        variable_genes = adata.var_names[adata.var['highly_variable']].tolist()
    else:
        # Fallback: use top genes by variance
        variances = np.var(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X, axis=0)
        top_indices = np.argsort(variances)[-n_top_genes:]
        variable_genes = adata.var_names[top_indices].tolist()
    
    return variable_genes


def find_cluster_markers(
    adata: AnnData,
    cluster_key: str = "leiden",
    n_genes: int = 20,
    method: str = "wilcoxon",
    min_fold_change: float = 1.5,
    max_pval: float = 0.05,
    only_positive: bool = True
) -> Dict[str, List[str]]:
    """
    Identify marker genes for each cluster.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with clustering results.
    cluster_key : str
        Key in adata.obs containing cluster assignments.
    n_genes : int
        Number of top marker genes per cluster.
    method : str
        Statistical test to use (wilcoxon, t-test, logreg).
    min_fold_change : float
        Minimum fold change to consider a gene as marker.
    max_pval : float
        Maximum adjusted p-value to consider significant.
    only_positive : bool
        Only return upregulated markers.
    
    Returns
    -------
    dict
        Dictionary mapping cluster_id -> list of marker genes.
    """
    # Ensure cluster key exists
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
    
    # Run differential expression
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method=method,
        use_raw=False,
        n_genes=n_genes * 2  # Get more initially for filtering
    )
    
    # Extract markers for each cluster
    cluster_markers = {}
    clusters = adata.obs[cluster_key].unique()
    
    for cluster in clusters:
        cluster_id = str(cluster)
        
        # Get results for this cluster
        genes = pd.DataFrame({
            'names': adata.uns['rank_genes_groups']['names'][cluster_id],
            'scores': adata.uns['rank_genes_groups']['scores'][cluster_id],
            'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][cluster_id],
            'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][cluster_id]
        })
        
        # Filter by significance and fold change
        if only_positive:
            genes = genes[genes['logfoldchanges'] > 0]
        
        genes = genes[
            (genes['pvals_adj'] < max_pval) &
            (np.abs(genes['logfoldchanges']) > np.log2(min_fold_change))
        ]
        
        # Take top n_genes
        top_genes = genes.nlargest(n_genes, 'scores')['names'].tolist()
        
        cluster_markers[cluster_id] = top_genes
    
    return cluster_markers


def filter_marker_genes(
    marker_genes: List[str],
    remove_ribosomal: bool = True,
    remove_mitochondrial: bool = True,
    remove_hsp: bool = True,
    custom_exclude: Optional[List[str]] = None
) -> List[str]:
    """
    Filter marker genes to remove non-informative genes.
    
    Parameters
    ----------
    marker_genes : list of str
        List of marker genes to filter.
    remove_ribosomal : bool
        Remove ribosomal protein genes (RPL*, RPS*).
    remove_mitochondrial : bool
        Remove mitochondrial genes (MT-*).
    remove_hsp : bool
        Remove heat shock protein genes (HSP*).
    custom_exclude : list of str, optional
        Additional gene patterns to exclude.
    
    Returns
    -------
    list of str
        Filtered list of marker genes.
    """
    filtered_genes = []
    
    for gene in marker_genes:
        gene_upper = gene.upper()
        
        # Skip ribosomal genes
        if remove_ribosomal and (gene_upper.startswith('RPL') or gene_upper.startswith('RPS')):
            continue
        
        # Skip mitochondrial genes
        if remove_mitochondrial and gene_upper.startswith('MT-'):
            continue
        
        # Skip heat shock proteins
        if remove_hsp and gene_upper.startswith('HSP'):
            continue
        
        # Skip custom exclusions
        if custom_exclude:
            if any(pattern.upper() in gene_upper for pattern in custom_exclude):
                continue
        
        filtered_genes.append(gene)
    
    return filtered_genes


def get_top_markers_per_cluster(
    adata: AnnData,
    cluster_key: str = "leiden",
    n_markers: int = 10,
    filter_genes: bool = True,
    **kwargs
) -> Dict[str, List[str]]:
    """
    Convenience function to get filtered top marker genes for each cluster.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cluster_key : str
        Key in adata.obs containing cluster assignments.
    n_markers : int
        Number of marker genes per cluster.
    filter_genes : bool
        Whether to filter out ribosomal/mitochondrial genes.
    **kwargs
        Additional arguments passed to find_cluster_markers.
    
    Returns
    -------
    dict
        Dictionary mapping cluster_id -> filtered list of marker genes.
    """
    # Find markers (get extra to account for filtering)
    n_genes_initial = n_markers * 3 if filter_genes else n_markers
    
    all_markers = find_cluster_markers(
        adata,
        cluster_key=cluster_key,
        n_genes=n_genes_initial,
        **kwargs
    )
    
    # Filter if requested
    if filter_genes:
        filtered_markers = {}
        for cluster_id, genes in all_markers.items():
            filtered = filter_marker_genes(genes)
            # Take only n_markers after filtering
            filtered_markers[cluster_id] = filtered[:n_markers]
    else:
        filtered_markers = {
            cluster_id: genes[:n_markers]
            for cluster_id, genes in all_markers.items()
        }
    
    return filtered_markers


def validate_marker_genes(
    marker_genes: Union[List[str], Dict[str, List[str]]],
    adata: AnnData
) -> Union[List[str], Dict[str, List[str]]]:
    """
    Validate that marker genes exist in the dataset.
    
    Parameters
    ----------
    marker_genes : list or dict
        Marker genes to validate (either list or dict of cluster -> genes).
    adata : AnnData
        Annotated data matrix to validate against.
    
    Returns
    -------
    list or dict
        Validated marker genes (removes genes not in dataset).
    """
    available_genes = set(adata.var_names)
    
    if isinstance(marker_genes, dict):
        validated = {}
        for cluster_id, genes in marker_genes.items():
            valid_genes = [g for g in genes if g in available_genes]
            if len(valid_genes) < len(genes):
                missing = set(genes) - set(valid_genes)
                logger.warning("Cluster %s — %d marker genes not found in dataset", cluster_id, len(missing))
            validated[cluster_id] = valid_genes
        return validated
    else:
        valid_genes = [g for g in marker_genes if g in available_genes]
        if len(valid_genes) < len(marker_genes):
            missing = set(marker_genes) - set(valid_genes)
            logger.warning("%d marker genes not found in dataset", len(missing))
        return valid_genes


def get_marker_expression_summary(
    adata: AnnData,
    marker_genes: List[str],
    cluster_key: str = "leiden"
) -> pd.DataFrame:
    """
    Get expression summary statistics for marker genes across clusters.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    marker_genes : list of str
        Marker genes to summarize.
    cluster_key : str
        Key in adata.obs containing cluster assignments.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with mean expression per gene per cluster.
    """
    # Validate genes
    marker_genes = validate_marker_genes(marker_genes, adata)
    
    # Get expression data
    gene_indices = [adata.var_names.tolist().index(g) for g in marker_genes]
    
    if hasattr(adata.X, 'toarray'):
        expr_data = adata.X[:, gene_indices].toarray()
    else:
        expr_data = adata.X[:, gene_indices]
    
    # Calculate mean expression per cluster
    clusters = adata.obs[cluster_key].values
    unique_clusters = np.unique(clusters)
    
    summary = []
    for cluster in unique_clusters:
        cluster_mask = clusters == cluster
        cluster_expr = expr_data[cluster_mask, :]
        mean_expr = np.mean(cluster_expr, axis=0)
        
        for gene_idx, gene in enumerate(marker_genes):
            summary.append({
                'cluster': str(cluster),
                'gene': gene,
                'mean_expression': mean_expr[gene_idx],
                'n_cells': cluster_mask.sum()
            })
    
    return pd.DataFrame(summary)

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
