#!/usr/bin/env python3
"""
scRN_AI v1.0.0

Agentic, evidence-based cell type annotation powered by CyteType.

CyteType (https://github.com/NygenAnalytics/CyteType) uses a multi-agent
AI architecture that outperforms naive GPT-based approaches by +388%.
No API keys are required for the default configuration.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd
import scanpy as sc
from anndata import AnnData

from ..utils.cytetype_client import CellTypePrediction, CyteTypeClient

logger = logging.getLogger(__name__)


def run(
    input_path: str,
    output_dir: str,
    timing: str = "pre_analysis",
    confidence_threshold: float = 0.7,
    n_top_genes: int = 100,
    max_clusters: int = 50,
    species: str = "human",
    tissue: Optional[str] = None,
    cluster_key: str = "leiden",
    study_context: Optional[str] = None,
):
    """
    Run CyteType cell type annotation.

    Parameters
    ----------
    input_path : str
        Path to input .h5ad file.
    output_dir : str
        Output directory for annotations.
    timing : str
        When to run: 'pre_analysis', 'post_analysis', or 'both'.
    confidence_threshold : float
        Minimum confidence score to accept (0.0-1.0).
    n_top_genes : int
        Number of top marker genes per cluster for CyteType (default: 100).
    max_clusters : int
        Maximum number of clusters to annotate.
    species : str
        Species (human, mouse, etc.).
    tissue : str, optional
        Tissue type if known — folded into study_context.
    cluster_key : str
        Key in adata.obs containing cluster assignments.
    study_context : str, optional
        Free-text study context passed to CyteType (e.g.,
        "Human PBMC from healthy donor").
    """
    logger.info("=" * 60)
    logger.info("CyteType Cell Type Annotation (AItyping)")
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
        if "neighbors" not in adata.uns:
            logger.info("Computing neighbors...")
            sc.pp.neighbors(adata, n_neighbors=15)
        logger.info("Running Leiden clustering...")
        sc.tl.leiden(adata, key_added=cluster_key)

    n_clusters = adata.obs[cluster_key].nunique()
    logger.info("Found %d clusters in '%s'", n_clusters, cluster_key)

    if n_clusters > max_clusters:
        logger.warning(
            "%d clusters exceeds max_clusters=%d — only the %d largest will be processed",
            n_clusters,
            max_clusters,
            max_clusters,
        )

    # Ensure rank_genes_groups has been computed for the cluster key
    rank_key = f"rank_genes_{cluster_key}"
    if rank_key not in adata.uns:
        logger.info("Computing marker genes (rank_genes_groups) for '%s'...", cluster_key)
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")
        # CyteType expects the key format "rank_genes_<group_key>"
        adata.uns[rank_key] = adata.uns["rank_genes_groups"]

    # Build study context
    ctx = study_context or ""
    if species:
        ctx = f"{species.capitalize()} sample. {ctx}".strip()
    if tissue:
        ctx = f"{ctx} Tissue: {tissue}.".strip()
    if timing:
        ctx = f"{ctx} Analysis stage: {timing}.".strip()

    # Initialise CyteType client
    logger.info("Initialising CyteType (n_top_genes=%d)...", n_top_genes)
    try:
        client = CyteTypeClient(n_top_genes=n_top_genes, study_context=ctx or None)
    except ImportError as e:
        logger.error("Failed to initialise CyteType: %s", e)
        logger.error("Install with: pip install cytetype")
        raise

    # Run annotation
    logger.info("Running CyteType annotation...")
    adata, predictions = client.annotate(
        adata,
        cluster_key=cluster_key,
        rank_key=rank_key,
        study_context=ctx or None,
    )

    # Filter by confidence
    logger.info("Filtering predictions by confidence threshold: %.2f", confidence_threshold)
    high_confidence = {
        cid: pred for cid, pred in predictions.items() if pred.confidence >= confidence_threshold
    }
    low_confidence = {
        cid: pred for cid, pred in predictions.items() if pred.confidence < confidence_threshold
    }

    logger.info("High confidence: %d clusters", len(high_confidence))
    logger.info("Low confidence: %d clusters", len(low_confidence))

    # Save results
    _save_results(
        predictions=predictions,
        high_confidence=high_confidence,
        low_confidence=low_confidence,
        output_dir=output_path,
        timing=timing,
        adata=adata,
        cluster_key=cluster_key,
    )

    logger.info("CyteType annotation complete! Results saved to: %s", output_dir)
    logger.info(
        "Output: %s_annotations.csv, %s_confidence_scores.csv, %s_reasoning.txt",
        timing,
        timing,
        timing,
    )
    if low_confidence:
        logger.info("Low confidence clusters written to: %s_low_confidence.csv", timing)


def _save_results(
    predictions: dict[str, CellTypePrediction],
    high_confidence: dict[str, CellTypePrediction],
    low_confidence: dict[str, CellTypePrediction],
    output_dir: Path,
    timing: str,
    adata: AnnData,
    cluster_key: str,
):
    """Save CyteType annotation results to various output files."""

    # 1. Main annotations file
    annotations = []
    for cluster_id, pred in predictions.items():
        annotations.append(
            {
                "cluster": cluster_id,
                "predicted_cell_type": pred.predicted_type,
                "confidence": pred.confidence,
                "cell_ontology_id": pred.cell_ontology_id or "",
                "n_cells": (adata.obs[cluster_key].astype(str) == cluster_id).sum(),
                "status": "high_confidence" if cluster_id in high_confidence else "low_confidence",
            }
        )

    annotations_df = pd.DataFrame(annotations)
    annotations_df.to_csv(output_dir / f"{timing}_annotations.csv", index=False)

    # 2. Confidence scores with alternatives
    confidence_data = []
    for cluster_id, pred in predictions.items():
        row = {
            "cluster": cluster_id,
            "primary_type": pred.predicted_type,
            "primary_confidence": pred.confidence,
        }

        if pred.alternative_types:
            for idx, (alt_type, alt_conf) in enumerate(pred.alternative_types[:3], 1):
                row[f"alternative_{idx}"] = alt_type
                row[f"alternative_{idx}_confidence"] = alt_conf

        confidence_data.append(row)

    confidence_df = pd.DataFrame(confidence_data)
    confidence_df.to_csv(output_dir / f"{timing}_confidence_scores.csv", index=False)

    # 3. Reasoning text file
    with open(output_dir / f"{timing}_reasoning.txt", "w") as f:
        f.write(f"CyteType Cell Type Annotation Reasoning ({timing})\n")
        f.write("=" * 60 + "\n\n")

        for cluster_id in sorted(predictions.keys()):
            pred = predictions[cluster_id]
            f.write(f"Cluster {cluster_id}: {pred.predicted_type}\n")
            f.write(f"Confidence: {pred.confidence:.2f}\n")
            if pred.cell_ontology_id:
                f.write(f"Cell Ontology ID: {pred.cell_ontology_id}\n")
            if pred.marker_genes_used:
                f.write(f"Marker genes: {', '.join(pred.marker_genes_used)}\n")
            f.write(f"\nReasoning:\n{pred.reasoning}\n")

            if pred.literature_refs:
                f.write("\nLiterature references:\n")
                for ref in pred.literature_refs:
                    f.write(f"  - {ref}\n")

            if pred.alternative_types:
                f.write("\nAlternative types:\n")
                for alt_type, alt_conf in pred.alternative_types:
                    f.write(f"  - {alt_type} (confidence: {alt_conf:.2f})\n")

            f.write("\n" + "-" * 60 + "\n\n")

    # 4. Low confidence clusters (if any)
    if low_confidence:
        low_conf_data = []
        for cluster_id, pred in low_confidence.items():
            low_conf_data.append(
                {
                    "cluster": cluster_id,
                    "predicted_type": pred.predicted_type,
                    "confidence": pred.confidence,
                    "cell_ontology_id": pred.cell_ontology_id or "",
                    "note": "Review recommended - confidence below threshold",
                }
            )

        low_conf_df = pd.DataFrame(low_conf_data)
        low_conf_df.to_csv(output_dir / f"{timing}_low_confidence.csv", index=False)

    # 5. Add annotations to AnnData object
    cluster_to_type = {cluster_id: pred.predicted_type for cluster_id, pred in predictions.items()}

    cell_types = adata.obs[cluster_key].astype(str).map(cluster_to_type)
    adata.obs[f"aitype_{timing}"] = cell_types
    adata.obs[f"aitype_{timing}_confidence"] = (
        adata.obs[cluster_key]
        .astype(str)
        .map({cluster_id: pred.confidence for cluster_id, pred in predictions.items()})
    )

    # Save updated AnnData
    output_h5ad = output_dir / f"{timing}_annotated.h5ad"
    adata.write_h5ad(output_h5ad)
    logger.info("%s_annotated.h5ad: Updated AnnData with cell type annotations", timing)


# scRN_AI v1.0.0
# Any usage is subject to this software's license.
