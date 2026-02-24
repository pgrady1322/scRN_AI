#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

Merge multiple AnnData files into one.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import logging
from pathlib import Path
from typing import Sequence, Union

import anndata as ad

logger = logging.getLogger(__name__)


def run(infiles: Sequence[Union[str, Path]], outfile: Union[str, Path]) -> None:
    """Merge multiple .h5ad files into a single AnnData object.

    Parameters
    ----------
    infiles : sequence of str or Path
        Paths to the .h5ad files to merge.
    outfile : str or Path
        Path to the merged output .h5ad file.
    """
    logger.info("Merging %d files", len(infiles))
    adatas = [ad.read_h5ad(f) for f in infiles]
    merged = ad.concat(adatas, join="outer")
    merged.obs_names_make_unique()
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(outfile)
    logger.info("Merged result: %d cells × %d genes → %s", merged.n_obs, merged.n_vars, outfile)

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
