"""
Large-scale trajectory workflow with STAVIA (VIA 2.0).
"""
import anndata as ad
import pathlib as p
import scanpy as sc


def run(infile, root_cell, outfile):
    adata = ad.read_h5ad(infile)

    try:
        from pyVIA.core import VIA
    except ImportError:
        raise ImportError(
            "pyVIA is required for large-scale trajectory analysis. "
            "Install with: pip install pyVIA"
        )

    # Ensure PCA is computed
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata, n_comps=50)

    via = VIA(
        data=adata.obsm['X_pca'],
        true_label=adata.obs['leiden'].values if 'leiden' in adata.obs else None,
        root_user=root_cell,
        jac_std_global=0.15,
    )
    via.run_VIA()

    # Store outputs
    adata.obs["via_pseudotime"] = via.single_cell_pt_markov
    adata.obs["via_cluster"] = via.labels
    _save(adata, outfile)

def _save(adata, out):
    out = p.Path(out)
    if out.suffix == ".h5ad":
        adata.write(out)
    else:
        adata.obs.to_csv(out, sep="\t")