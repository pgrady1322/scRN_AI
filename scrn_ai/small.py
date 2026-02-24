#!/usr/bin/env python3
"""
scRN_AI v1.0.0

Small-scale pseudotime workflow (DPT, DTFLOW, BLTSA).

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import logging
import os
import pathlib as p

import anndata as ad
import scanpy as sc

logger = logging.getLogger(__name__)


def run(infile, species, method, run_bltsa, outfile):
    adata = _read_any(infile)
    adata.uns["species"] = species

    # Basic QC & normalisation
    logger.info("Preprocessing %d cells × %d genes", adata.n_obs, adata.n_vars)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
    sc.pp.scale(adata)

    # Neighbours
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata)

    if method.lower() == "dpt":
        logger.info("Running Diffusion Pseudotime (DPT)")
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata)
        adata.obs["pseudotime"] = adata.obs["dpt_pseudotime"]
    else:  # dtflow
        logger.info("Running DTFlow pseudotime")
        import dtflow

        adata.obs["pseudotime"] = dtflow.run(adata)

    if run_bltsa:
        logger.info("Running BLTSA trajectory refinement")
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            from scipy import sparse

            pandas2ri.activate()

            ro.globalenv["expression_matrix"] = (
                adata.X if not sparse.issparse(adata.X) else adata.X.toarray()
            )
            bltsa_path = os.environ.get("BLTSA_PATH", "/opt/BLTSA/BLTSA.R")
            logger.info("Loading BLTSA from %s", bltsa_path)

            ro.r(f'''
                source("{bltsa_path}")
                result <- BLTSA_trajectory(expression_matrix)
                pseudotime <- result$pseudotime
            ''')

            import numpy as np

            adata.obs["bltsa_pt"] = np.array(ro.r["pseudotime"])
        except Exception as exc:
            logger.warning("BLTSA failed (%s), skipping", exc)

    # Save result
    _save(adata, outfile)
    logger.info("Saved results to %s", outfile)


def _read_any(path):
    path = p.Path(path)
    if path.suffix in {".h5ad"}:
        return ad.read_h5ad(path)
    elif path.suffix in {".h5"}:
        return sc.read_10x_h5(path)
    else:  # assume 10x mtx directory
        return sc.read_10x_mtx(path, gex_only=True)


def _save(adata, out):
    out = p.Path(out)
    if out.suffix == ".h5ad":
        adata.write(out)
    else:
        adata.obs[["pseudotime"]].to_csv(out, sep="\t")


# scRN_AI v1.0.0
# Any usage is subject to this software's license.
