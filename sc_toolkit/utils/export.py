import anndata as ad, pathlib as p

def run(infile, outdir, fmt):
    fmt = fmt.lower()
    outdir = p.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(infile)

    if fmt == "loom":
        from loompy import create
        loom_path = outdir / "adata.loom"
        adata.write_loom(loom_path)
    elif fmt == "mtx":
        from scipy import sparse, io
        io.mmwrite(outdir / "matrix.mtx", sparse.csr_matrix(adata.X))
        adata.obs_names.to_series().to_csv(outdir / "barcodes.tsv", index=False, header=False)
        adata.var_names.to_series().to_csv(outdir / "features.tsv", index=False, header=False)
    elif fmt == "csv":
        adata.to_df().to_csv(outdir / "matrix.csv")
        adata.obs.to_csv(outdir / "obs.csv")
        adata.var.to_csv(outdir / "var.csv")
    else:
        raise ValueError("format must be loom|mtx|csv")