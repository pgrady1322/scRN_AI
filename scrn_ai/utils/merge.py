import anndata as ad
import pathlib as p


def run(infiles, outfile):
    """Merge multiple AnnData files into one."""
    adatas = [ad.read_h5ad(f) for f in infiles]
    merged = ad.concat(adatas, join="outer")
    merged.obs_names_make_unique()
    p.Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(outfile)