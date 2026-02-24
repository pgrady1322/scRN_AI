#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v0.1.0

UMAP pseudotime plotting utility.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import scanpy as sc
def umap_pseudotime(adata, key="pseudotime"):
    if "X_umap" not in adata.obsm.keys():
        sc.tl.umap(adata)
    sc.pl.umap(adata, color=key, cmap="viridis", show=False)

# scRN_AI v0.1.0
# Any usage is subject to this software's license.
