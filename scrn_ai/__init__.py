"""
scrn_ai - Single-cell RNA-seq analysis toolkit

A comprehensive toolkit for single-cell RNA-seq data analysis including:
- Quality control and preprocessing
- Multiple normalization methods (Seurat, JMP, basic)
- Dimensionality reduction (PCA, UMAP)
- Clustering and cell type annotation
- Pseudotime trajectory analysis
- AI-powered cell typing
"""

__version__ = "0.1.0"
__author__ = "Patrick Grady"
__email__ = ""

from . import cli
from . import utils
from . import workflows

__all__ = ["cli", "utils", "workflows"]
