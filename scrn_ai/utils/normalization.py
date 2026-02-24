#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

In-place AnnData normalization (log1p, scran, sctransform, size-factor).

Thin wrapper around workflows.normalization for shared methods;
adds the lightweight size_factor option for quick CPM-style scaling.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging

import scanpy as sc
import anndata as ad
import pathlib as p

logger = logging.getLogger(__name__)

# Methods handled by the full normalization workflow
_DELEGATED_METHODS = {"log1p", "scran", "sctransform"}


def run(infile, outfile, method):
    """Normalize an AnnData file and write to disk.

    Parameters
    ----------
    infile : str
        Path to input .h5ad file.
    outfile : str
        Path to output .h5ad file.
    method : str
        One of: log1p, scran, sctransform, size_factor.
    """
    m = method.lower()

    if m in _DELEGATED_METHODS:
        # Delegate to the full-featured workflow module
        from ..workflows import normalization as wf_norm

        logger.info("Delegating '%s' to workflows.normalization", m)
        wf_norm.run(infile, outfile, method=m)
        return

    if m == "size_factor":
        logger.info("Running size-factor normalization (CPM)...")
        adata = ad.read_h5ad(infile)
        sc.pp.normalize_total(adata, target_sum=1e6)
        p.Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(outfile)
        logger.info("Saved to %s", outfile)
        return

    raise ValueError(f"Unknown method: {method}. Choose from: log1p, scran, sctransform, size_factor")

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
