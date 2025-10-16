"""
Large-scale trajectory workflow with STAVIA (VIA 2.0).
"""
import anndata as ad
import pathlib as p
from VIA import VIA

def run(infile, root_cell, outfile):
    adata = ad.read_h5ad(infile)

    via = VIA(adata, jac_std_global=0.15, root_user=root_cell)
    via.run_VIA()

    # Store outputs
    adata.obs["via_pseudotime"] = via.pseudotime
    adata.obs["via_cluster"] = via.labels
    _save(adata, outfile)

def _save(adata, out):
    out = p.Path(out)
    if out.suffix == ".h5ad":
        adata.write(out)
    else:
        adata.obs.to_csv(out, sep="\t")