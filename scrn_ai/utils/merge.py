import scanpy as sc, anndata as ad, pathlib as p
from functools import reduce

def run(infiles, outfile):
    adatas = [ad.read_h5ad(f) for f in infiles]
    merged = reduce(lambda x, y: x.concatenate(y, join="outer"), adatas)
    merged.write_h5ad(outfile)