#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v0.1.0

AI-powered cell type annotation using OpenAI GPT models.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

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
    print("=" * 60)
    print("AI-Powered Cell Type Identification (AItyping)")
    print("=" * 60)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print(f"\nLoading data from: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"  Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Check for clustering
    if cluster_key not in adata.obs.columns:
        print(f"\nWarning: No '{cluster_key}' clustering found. Running Leiden clustering...")
        # Compute neighbors if not present
        if 'neighbors' not in adata.uns:
            print("  Computing neighbors...")
            sc.pp.neighbors(adata, n_neighbors=15)
        # Run clustering
        print("  Running Leiden clustering...")
        sc.tl.leiden(adata, key_added=cluster_key)
    
    n_clusters = adata.obs[cluster_key].nunique()
    print(f"\nFound {n_clusters} clusters in '{cluster_key}'")
    
    if n_clusters > max_clusters:
        print(f"Warning: {n_clusters} clusters exceeds max_clusters={max_clusters}")
        print(f"  Only processing first {max_clusters} clusters")
        # Get top N clusters by size
        cluster_sizes = adata.obs[cluster_key].value_counts()
        clusters_to_process = cluster_sizes.head(max_clusters).index.tolist()
    else:
        clusters_to_process = adata.obs[cluster_key].unique().tolist()
    
    # Get marker genes
    print(f"\nIdentifying marker genes...")
    if marker_genes:
        print(f"  Loading custom markers from: {marker_genes}")
        custom_markers = pd.read_csv(marker_genes)
        # Expect format: cluster,gene1,gene2,gene3,...
        cluster_markers = {}
        for _, row in custom_markers.iterrows():
            cluster_id = str(row['cluster'])
            genes = [g for g in row[1:] if pd.notna(g)]
            cluster_markers[cluster_id] = genes[:n_markers]
    else:
        print(f"  Auto-detecting top {n_markers} markers per cluster...")
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
    print(f"\n  Marker genes identified for {len(cluster_markers)} clusters:")
    for cluster_id, markers in sorted(cluster_markers.items())[:5]:
        print(f"    Cluster {cluster_id}: {', '.join(markers[:5])}...")
    if len(cluster_markers) > 5:
        print(f"    ... and {len(cluster_markers) - 5} more clusters")
    
    # Initialize OpenAI client
    print(f"\nInitializing OpenAI client (model: {model})...")
    try:
        client = OpenAIClient(
            api_key=api_key,
            model=model,
            temperature=0.3
        )
    except Exception as e:
        print(f"Error: Failed to initialize OpenAI client: {e}")
        print("\nPlease ensure:")
        print("  1. openai package is installed: pip install openai")
        print("  2. OPENAI_API_KEY environment variable is set")
        raise
    
    # Run predictions
    print(f"\nRunning cell type predictions...")
    print(f"  This will make {len(cluster_markers)} API calls (rate limited)...")
    
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
    print(f"\nFiltering predictions by confidence threshold: {confidence_threshold}")
    high_confidence = {
        cid: pred for cid, pred in predictions.items()
        if pred.confidence >= confidence_threshold
    }
    low_confidence = {
        cid: pred for cid, pred in predictions.items()
        if pred.confidence < confidence_threshold
    }
    
    print(f"  High confidence: {len(high_confidence)} clusters")
    print(f"  Low confidence: {len(low_confidence)} clusters")
    
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
    
    print(f"\n✅ AItyping complete! Results saved to: {output_dir}")
    print("\nOutput files:")
    print(f"  - {timing}_annotations.csv: Cell type predictions")
    print(f"  - {timing}_confidence_scores.csv: Confidence metrics")
    print(f"  - {timing}_marker_genes.csv: Marker genes used")
    print(f"  - {timing}_reasoning.txt: AI reasoning for each cluster")
    if low_confidence:
        print(f"  - {timing}_low_confidence.csv: Clusters needing review")


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
    print(f"  - {timing}_annotated.h5ad: Updated AnnData with cell type annotations")

# scRN_AI v0.1.0
# Any usage is subject to this software's license.
