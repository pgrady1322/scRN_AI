#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

AI-powered cell type annotation using OpenAI GPT models.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import logging
import os
from pathlib import Path
from typing import Optional, List, Dict
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ..utils.openai_client import OpenAIClient, CellTypePrediction
from ..utils.marker_detection import (
    get_top_markers_per_cluster,
    validate_marker_genes,
    get_marker_expression_summary
)

logger = logging.getLogger(__name__)


def run(
    input_path: str,
    output_dir: str,
    timing: str = "pre_analysis",
    model: str = "gpt-4",
    confidence_threshold: float = 0.7,
    marker_genes: Optional[str] = None,
    n_markers: int = 10,
    max_clusters: int = 50,
    species: str = "human",
    tissue: Optional[str] = None,
    cluster_key: str = "leiden",
    api_key: Optional[str] = None
):
    """
    Run AI-powered cell type identification.
    
    Parameters
    ----------
    input_path : str
        Path to input .h5ad file.
    output_dir : str
        Output directory for annotations.
    timing : str
        When to run: 'pre_analysis', 'post_analysis', or 'both'.
    model : str
        OpenAI model to use (gpt-4, gpt-4-turbo, gpt-3.5-turbo).
    confidence_threshold : float
        Minimum confidence score to accept (0.0-1.0).
    marker_genes : str, optional
        Path to custom marker gene CSV file.
    n_markers : int
        Number of top marker genes per cluster.
    max_clusters : int
        Maximum number of clusters to annotate.
    species : str
        Species (human, mouse, etc.).
    tissue : str, optional
        Tissue type if known.
    cluster_key : str
        Key in adata.obs containing cluster assignments.
    api_key : str, optional
        OpenAI API key (if not in environment).
    """
    logger.info("=" * 60)
    logger.info("AI-Powered Cell Type Identification (AItyping)")
    logger.info("=" * 60)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load data
    logger.info("Loading data from: %s", input_path)
    adata = sc.read_h5ad(input_path)
    logger.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)
    
    # Check for clustering
    if cluster_key not in adata.obs.columns:
        logger.warning("No '%s' clustering found. Running Leiden clustering...", cluster_key)
        # Compute neighbors if not present
        if 'neighbors' not in adata.uns:
            logger.info("Computing neighbors...")
            sc.pp.neighbors(adata, n_neighbors=15)
        # Run clustering
        logger.info("Running Leiden clustering...")
        sc.tl.leiden(adata, key_added=cluster_key)
    
    n_clusters = adata.obs[cluster_key].nunique()
    logger.info("Found %d clusters in '%s'", n_clusters, cluster_key)
    
    if n_clusters > max_clusters:
        logger.warning("%d clusters exceeds max_clusters=%d", n_clusters, max_clusters)
        logger.info("Only processing first %d clusters", max_clusters)
        # Get top N clusters by size
        cluster_sizes = adata.obs[cluster_key].value_counts()
        clusters_to_process = cluster_sizes.head(max_clusters).index.tolist()
    else:
        clusters_to_process = adata.obs[cluster_key].unique().tolist()
    
    # Get marker genes
    logger.info("Identifying marker genes...")
    if marker_genes:
        logger.info("Loading custom markers from: %s", marker_genes)
        custom_markers = pd.read_csv(marker_genes)
        # Expect format: cluster,gene1,gene2,gene3,...
        cluster_markers = {}
        for _, row in custom_markers.iterrows():
            cluster_id = str(row['cluster'])
            genes = [g for g in row[1:] if pd.notna(g)]
            cluster_markers[cluster_id] = genes[:n_markers]
    else:
        logger.info("Auto-detecting top %d markers per cluster...", n_markers)
        cluster_markers = get_top_markers_per_cluster(
            adata,
            cluster_key=cluster_key,
            n_markers=n_markers,
            filter_genes=True
        )
    
    # Filter to clusters we're processing
    cluster_markers = {
        str(c): markers for c, markers in cluster_markers.items()
        if c in [str(x) for x in clusters_to_process]
    }
    
    # Validate markers
    cluster_markers = validate_marker_genes(cluster_markers, adata)
    
    # Print marker summary
    logger.info("Marker genes identified for %d clusters:", len(cluster_markers))
    for cluster_id, markers in sorted(cluster_markers.items())[:5]:
        logger.info("  Cluster %s: %s...", cluster_id, ', '.join(markers[:5]))
    if len(cluster_markers) > 5:
        logger.info("  ... and %d more clusters", len(cluster_markers) - 5)
    
    # Initialize OpenAI client
    logger.info("Initializing OpenAI client (model: %s)...", model)
    try:
        client = OpenAIClient(
            api_key=api_key,
            model=model,
            temperature=0.3
        )
    except Exception as e:
        logger.error("Failed to initialize OpenAI client: %s", e)
        logger.error("Ensure: (1) openai is installed, (2) OPENAI_API_KEY is set")
        raise
    
    # Run predictions
    logger.info("Running cell type predictions...")
    logger.info("Making %d API calls (rate limited)...", len(cluster_markers))
    
    context_info = f"Timing: {timing}"
    if tissue:
        context_info += f", Tissue: {tissue}"
    
    predictions = client.batch_predict_cell_types(
        clusters=cluster_markers,
        species=species,
        tissue=tissue,
        context=context_info
    )
    
    # Filter by confidence
    logger.info("Filtering predictions by confidence threshold: %.2f", confidence_threshold)
    high_confidence = {
        cid: pred for cid, pred in predictions.items()
        if pred.confidence >= confidence_threshold
    }
    low_confidence = {
        cid: pred for cid, pred in predictions.items()
        if pred.confidence < confidence_threshold
    }
    
    logger.info("High confidence: %d clusters", len(high_confidence))
    logger.info("Low confidence: %d clusters", len(low_confidence))
    
    # Save results
    _save_results(
        predictions=predictions,
        high_confidence=high_confidence,
        low_confidence=low_confidence,
        cluster_markers=cluster_markers,
        output_dir=output_path,
        timing=timing,
        adata=adata,
        cluster_key=cluster_key
    )
    
    logger.info("AItyping complete! Results saved to: %s", output_dir)
    logger.info("Output: %s_annotations.csv, %s_confidence_scores.csv, %s_marker_genes.csv, %s_reasoning.txt",
                timing, timing, timing, timing)
    if low_confidence:
        logger.info("Low confidence clusters written to: %s_low_confidence.csv", timing)


def _save_results(
    predictions: Dict[str, CellTypePrediction],
    high_confidence: Dict[str, CellTypePrediction],
    low_confidence: Dict[str, CellTypePrediction],
    cluster_markers: Dict[str, List[str]],
    output_dir: Path,
    timing: str,
    adata: AnnData,
    cluster_key: str
):
    """Save AItyping results to various output files."""
    
    # 1. Main annotations file
    annotations = []
    for cluster_id, pred in predictions.items():
        annotations.append({
            'cluster': cluster_id,
            'predicted_cell_type': pred.predicted_type,
            'confidence': pred.confidence,
            'n_cells': (adata.obs[cluster_key] == cluster_id).sum(),
            'status': 'high_confidence' if cluster_id in high_confidence else 'low_confidence'
        })
    
    annotations_df = pd.DataFrame(annotations)
    annotations_df.to_csv(output_dir / f"{timing}_annotations.csv", index=False)
    
    # 2. Confidence scores with alternatives
    confidence_data = []
    for cluster_id, pred in predictions.items():
        row = {
            'cluster': cluster_id,
            'primary_type': pred.predicted_type,
            'primary_confidence': pred.confidence
        }
        
        # Add alternative types
        if pred.alternative_types:
            for idx, (alt_type, alt_conf) in enumerate(pred.alternative_types[:3], 1):
                row[f'alternative_{idx}'] = alt_type
                row[f'alternative_{idx}_confidence'] = alt_conf
        
        confidence_data.append(row)
    
    confidence_df = pd.DataFrame(confidence_data)
    confidence_df.to_csv(output_dir / f"{timing}_confidence_scores.csv", index=False)
    
    # 3. Marker genes used
    marker_data = []
    for cluster_id, markers in cluster_markers.items():
        marker_data.append({
            'cluster': cluster_id,
            'marker_genes': ', '.join(markers),
            'n_markers': len(markers)
        })
    
    markers_df = pd.DataFrame(marker_data)
    markers_df.to_csv(output_dir / f"{timing}_marker_genes.csv", index=False)
    
    # 4. Reasoning text file
    with open(output_dir / f"{timing}_reasoning.txt", 'w') as f:
        f.write(f"AI Cell Type Identification Reasoning ({timing})\n")
        f.write("=" * 60 + "\n\n")
        
        for cluster_id in sorted(predictions.keys()):
            pred = predictions[cluster_id]
            f.write(f"Cluster {cluster_id}: {pred.predicted_type}\n")
            f.write(f"Confidence: {pred.confidence:.2f}\n")
            f.write(f"Marker genes: {', '.join(pred.marker_genes_used)}\n")
            f.write(f"\nReasoning:\n{pred.reasoning}\n")
            
            if pred.alternative_types:
                f.write(f"\nAlternative types:\n")
                for alt_type, alt_conf in pred.alternative_types:
                    f.write(f"  - {alt_type} (confidence: {alt_conf:.2f})\n")
            
            f.write("\n" + "-" * 60 + "\n\n")
    
    # 5. Low confidence clusters (if any)
    if low_confidence:
        low_conf_data = []
        for cluster_id, pred in low_confidence.items():
            low_conf_data.append({
                'cluster': cluster_id,
                'predicted_type': pred.predicted_type,
                'confidence': pred.confidence,
                'marker_genes': ', '.join(pred.marker_genes_used[:5]),
                'note': 'Review recommended - confidence below threshold'
            })
        
        low_conf_df = pd.DataFrame(low_conf_data)
        low_conf_df.to_csv(output_dir / f"{timing}_low_confidence.csv", index=False)
    
    # 6. Add annotations to AnnData object
    cluster_to_type = {
        cluster_id: pred.predicted_type
        for cluster_id, pred in predictions.items()
    }
    
    # Create cell type column
    cell_types = adata.obs[cluster_key].astype(str).map(cluster_to_type)
    adata.obs[f'aitype_{timing}'] = cell_types
    adata.obs[f'aitype_{timing}_confidence'] = adata.obs[cluster_key].astype(str).map({
        cluster_id: pred.confidence
        for cluster_id, pred in predictions.items()
    })
    
    # Save updated AnnData
    output_h5ad = output_dir / f"{timing}_annotated.h5ad"
    adata.write_h5ad(output_h5ad)
    logger.info("%s_annotated.h5ad: Updated AnnData with cell type annotations", timing)

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
