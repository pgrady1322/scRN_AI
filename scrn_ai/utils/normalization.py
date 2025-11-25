"""
Utility: normalise an AnnData file in-place and write new file.
Supports log1p, scran, sctransform, size-factor.
"""
import scanpy as sc, anndata as ad, pathlib as p, numpy as np

def run(infile, outfile, method):
    adata = ad.read_h5ad(infile)
    m = method.lower()

    if m == "log1p":
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif m == "size_factor":
        sc.pp.normalize_total(adata, target_sum=1e6)
    elif m in {"scran", "sctransform"}:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri; pandas2ri.activate()

        ro.globalenv["mat"] = sc.get.obs_df(adata, "counts") if "counts" in adata.layers else adata.to_df()
        if m == "scran":
            ro.r('''
                suppressMessages(library(scran))
                sf <- computeSumFactors(as.matrix(mat))
                mat_norm <- t(t(as.matrix(mat)) / sf)
            ''')
        else:  # sctransform
            ro.r('''
                suppressMessages(library(sctransform))
                mat_norm <- vst(as.matrix(mat))$y
            ''')
        mat_norm = np.array(ro.r["mat_norm"])
        adata.X = mat_norm
    else:
        raise ValueError(f"Unknown method: {method}")

    adata.write_h5ad(outfile)