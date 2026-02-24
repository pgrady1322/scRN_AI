#!/usr/bin/env python3
"""
scRN_AI v1.0.0

AnnData export to loom, mtx, and csv formats.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
from pathlib import Path
from typing import Union

import anndata as ad

logger = logging.getLogger(__name__)


def run(infile: Union[str, Path], outdir: Union[str, Path], fmt: str) -> None:
    """Export an AnnData .h5ad file to *loom*, *mtx*, or *csv*.

    Parameters
    ----------
    infile : str or Path
        Path to the input .h5ad file.
    outdir : str or Path
        Directory for exported output files.
    fmt : str
        Target format — one of ``loom``, ``mtx``, ``csv``.
    """
    fmt = fmt.lower()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(infile)
    logger.info("Exporting %s (%d cells × %d genes) to %s", infile, adata.n_obs, adata.n_vars, fmt)

    if fmt == "loom":
        loom_path = outdir / "adata.loom"
        adata.write_loom(loom_path)
    elif fmt == "mtx":
        from scipy import io, sparse

        io.mmwrite(outdir / "matrix.mtx", sparse.csr_matrix(adata.X))
        adata.obs_names.to_series().to_csv(outdir / "barcodes.tsv", index=False, header=False)
        adata.var_names.to_series().to_csv(outdir / "features.tsv", index=False, header=False)
    elif fmt == "csv":
        adata.to_df().to_csv(outdir / "matrix.csv")
        adata.obs.to_csv(outdir / "obs.csv")
        adata.var.to_csv(outdir / "var.csv")
    else:
        raise ValueError(f"Unsupported format '{fmt}'. Choose from: loom, mtx, csv")

    logger.info("Export complete → %s/", outdir)


# scRN_AI v1.0.0
# Any usage is subject to this software's license.
