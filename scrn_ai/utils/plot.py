import scanpy as sc
def umap_pseudotime(adata, key="pseudotime"):
    if "X_umap" not in adata.obsm.keys():
        sc.tl.umap(adata)
    sc.pl.umap(adata, color=key, cmap="viridis", show=False)