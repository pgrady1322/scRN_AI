#!/usr/bin/env python3
"""
scRN_AI v1.0.0

CyteType wrapper for agentic, evidence-based cell type annotation.

Replaces the former OpenAI GPT wrapper with CyteType
(https://github.com/NygenAnalytics/CyteType), a multi-agent AI architecture
that delivers transparent, ontology-mapped annotations with literature evidence
— no API keys required by default.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd
from anndata import AnnData

logger = logging.getLogger(__name__)

try:
    from cytetype import CyteType

    CYTETYPE_AVAILABLE = True
except ImportError:
    CYTETYPE_AVAILABLE = False
    logger.warning("cytetype package not installed. Install with: pip install cytetype")


@dataclass
class CellTypePrediction:
    """Container for cell type prediction results.

    This dataclass preserves the same interface used by the rest of
    scRN_AI so that downstream code (saving results, adding to .obs)
    continues to work without changes.
    """

    cluster_id: str
    predicted_type: str
    confidence: float
    reasoning: str
    marker_genes_used: list[str]
    alternative_types: list[tuple[str, float]] = field(default_factory=list)
    cell_ontology_id: Optional[str] = None
    literature_refs: list[str] = field(default_factory=list)


class CyteTypeClient:
    """
    Client wrapping CyteType for evidence-based cell type annotation.

    CyteType uses a multi-agent AI architecture to:
      - Analyse marker genes per cluster
      - Retrieve literature evidence
      - Map annotations to Cell Ontology (CL) IDs
      - Return confidence scores and full reasoning

    No API keys are required for the default configuration.
    See https://github.com/NygenAnalytics/CyteType/blob/master/docs/configuration.md
    for advanced LLM configuration.
    """

    def __init__(
        self,
        n_top_genes: int = 100,
        study_context: Optional[str] = None,
    ):
        """
        Initialise the CyteType client.

        Parameters
        ----------
        n_top_genes : int
            Number of top marker genes per cluster passed to CyteType
            (default: 100, as recommended by CyteType docs).
        study_context : str, optional
            Free-text context about the experiment (e.g.,
            "Human PBMC from healthy donor").  Passed to
            ``CyteType.run(study_context=...)``.
        """
        if not CYTETYPE_AVAILABLE:
            raise ImportError(
                "cytetype package not installed. Install with: pip install cytetype"
            )

        self.n_top_genes = n_top_genes
        self.study_context = study_context

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def annotate(
        self,
        adata: AnnData,
        cluster_key: str = "leiden",
        rank_key: Optional[str] = None,
        study_context: Optional[str] = None,
    ) -> tuple[AnnData, dict[str, CellTypePrediction]]:
        """
        Run CyteType annotation on an AnnData object.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix with clustering already performed.
        cluster_key : str
            Key in ``adata.obs`` that holds cluster labels.
        rank_key : str, optional
            Key in ``adata.uns`` from ``sc.tl.rank_genes_groups``.
            If *None*, defaults to ``"rank_genes_{cluster_key}"``.
        study_context : str, optional
            Override the study context set at init time.

        Returns
        -------
        adata : AnnData
            The input object with CyteType annotation columns added
            (e.g., ``cytetype_annotation_{cluster_key}``).
        predictions : dict[str, CellTypePrediction]
            Mapping of cluster id → prediction result, compatible with
            the rest of the scRN_AI pipeline.
        """
        context = study_context or self.study_context
        rk = rank_key or f"rank_genes_{cluster_key}"

        logger.info("Running CyteType annotation (group_key=%s, rank_key=%s)", cluster_key, rk)

        annotator = CyteType(
            adata,
            group_key=cluster_key,
            rank_key=rk,
            n_top_genes=self.n_top_genes,
        )

        if context:
            adata = annotator.run(study_context=context)
        else:
            adata = annotator.run()

        # ---- Build CellTypePrediction objects from CyteType output ----
        annotation_col = f"cytetype_annotation_{cluster_key}"

        predictions = self._extract_predictions(
            adata=adata,
            cluster_key=cluster_key,
            annotation_col=annotation_col,
        )

        return adata, predictions

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_predictions(
        adata: AnnData,
        cluster_key: str,
        annotation_col: str,
    ) -> dict[str, CellTypePrediction]:
        """
        Convert CyteType's output columns into ``CellTypePrediction`` objects.

        CyteType typically adds several columns to ``adata.obs``:
          - ``cytetype_annotation_<key>``   – predicted cell type label
          - ``cytetype_confidence_<key>``   – confidence score (0–1)
          - ``cytetype_cl_id_<key>``        – Cell Ontology CL ID

        If a column is missing, the helper falls back to sensible
        defaults so the rest of the pipeline never crashes.
        """
        conf_col = f"cytetype_confidence_{cluster_key}"
        cl_col = f"cytetype_cl_id_{cluster_key}"

        predictions: dict[str, CellTypePrediction] = {}

        clusters = adata.obs[cluster_key].unique()
        for cluster_id in clusters:
            cid = str(cluster_id)
            mask = adata.obs[cluster_key] == cluster_id

            # Cell type label
            if annotation_col in adata.obs.columns:
                predicted_type = adata.obs.loc[mask, annotation_col].iloc[0]
            else:
                predicted_type = "Unknown"

            # Confidence
            if conf_col in adata.obs.columns:
                confidence = float(adata.obs.loc[mask, conf_col].iloc[0])
            else:
                confidence = 0.0

            # Cell Ontology ID
            if cl_col in adata.obs.columns:
                cl_id = str(adata.obs.loc[mask, cl_col].iloc[0])
            else:
                cl_id = None

            predictions[cid] = CellTypePrediction(
                cluster_id=cid,
                predicted_type=str(predicted_type),
                confidence=confidence,
                reasoning="Annotated by CyteType multi-agent pipeline",
                marker_genes_used=[],
                cell_ontology_id=cl_id,
            )

        return predictions


# scRN_AI v1.0.0
# Any usage is subject to this software's license.
