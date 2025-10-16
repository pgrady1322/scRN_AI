"""
Small-scale pseudotime workflow:
  • Loads raw counts → AnnData
  • Normalises + computes neighbours
  • Diffusion pseudotime (Scanpy) or DTFLOW
  • Optionally BLTSA branching inference
"""
import scanpy as sc
import anndata as ad
import pathlib as p

def run(infile, species, method, bltsa, outfile):
    adata = _read_any(infile)
    adata.uns["species"] = species

    # Basic QC & normalisation
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
    sc.pp.scale(adata)

    # Neighbours
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata)

    if method.lower() == "dpt":
        sc.tl.diffusion_map(adata)
        sc.tl.dpt(adata)
        adata.obs["pseudotime"] = adata.obs["dpt_pseudotime"]
    else:         # dtflow
        import dtflow
        adata.obs["pseudotime"] = dtflow.run(adata)         # quick wrapper

    if bltsa:
        from BLTSA import bltsa
        adata.obs["bltsa_pt"] = bltsa.run(adata.X)

    # Save result
    _save(adata, outfile)

def _read_any(path):
    path = p.Path(path)
    if path.suffix in {".h5ad"}:
        return ad.read_h5ad(path)
    elif path.suffix in {".h5"}:
        return ad.read_10x_h5(path)
    else:     # assume 10x mtx directory
        return sc.read_10x_mtx(path, gex_only=True)

def _save(adata, out):
    out = p.Path(out)
    if out.suffix == ".h5ad":
        adata.write(out)
    else:
        adata.obs[["pseudotime"]].to_csv(out, sep="\t")