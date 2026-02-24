#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v0.1.0

Merge multiple AnnData files into one.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import anndata as ad
import pathlib as p


def run(infiles, outfile):
    """Merge multiple AnnData files into one."""
    adatas = [ad.read_h5ad(f) for f in infiles]
    merged = ad.concat(adatas, join="outer")
    merged.obs_names_make_unique()
    p.Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(outfile)

# scRN_AI v0.1.0
# Any usage is subject to this software's license.
