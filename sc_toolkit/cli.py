import click
from importlib import import_module

@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def main():
    """Single-cell analysis toolbox."""
    pass

# === workflow commands =================================================
@main.command()
@click.option("--infile", "-i", type=click.Path(exists=True), required=True,
              help="Raw count matrix (.mtx, .h5, or .h5ad).")
@click.option("--species", "-s", required=True,
              help="Species name tag stored in AnnData.uns.")
@click.option("--method", "-m", default="dpt",
              type=click.Choice(["dpt", "dtflow"], case_sensitive=False))
@click.option("--bltsa", is_flag=True, help="Run BLTSA branching inference.")
@click.option("--outfile", "-o", required=True, type=click.Path())
def small(infile, species, method, bltsa, outfile):
    """Smaller-scale pseudotime analysis workflow."""
    mod = import_module(".small", __package__)
    mod.run(infile, species, method, bltsa, outfile)

@main.command()
@click.option("--infile", "-i", type=click.Path(exists=True), required=True,
              help="Merged AnnData (.h5ad) for large datasets.")
@click.option("--root_cell", type=str, default=None,
              help="Cell ID to use as VIA root.")
@click.option("--outfile", "-o", required=True, type=click.Path())
def large(infile, root_cell, outfile):
    """Large-scale trajectory mapping with STAVIA (VIA 2.0)."""
    mod = import_module(".large", __package__)
    mod.run(infile, root_cell, outfile)

# === utility commands ==================================================
@main.command()
@click.option("--infile", "-i", type=click.Path(exists=True), required=True)
@click.option("--outfile", "-o", type=click.Path(), required=True)
@click.option("--method", "-m", default="log1p",
              type=click.Choice(["log1p", "scran", "sctransform", "size_factor"]))
def ad_norm(infile, outfile, method):
    """Normalise counts within an AnnData object."""
    from .utils import normalization as norm
    norm.run(infile, outfile, method)

@main.command()
@click.option("--infiles", "-i", multiple=True, required=True,
              type=click.Path(exists=True))
@click.option("--outfile", "-o", type=click.Path(), required=True)
def ad_merge(infiles, outfile):
    """Merge multiple AnnData files (different batches/species/timepoints)."""
    from .utils import merge
    merge.run(infiles, outfile)

@main.command()
@click.option("--infile", "-i", type=click.Path(exists=True), required=True)
@click.option("--outdir", "-o", type=click.Path(), required=True)
@click.option("--format", "-f",
              type=click.Choice(["loom", "mtx", "csv"], case_sensitive=False),
              required=True)
def ad_export(infile, outdir, format):
    """Export AnnData to loom / mtx / csv."""
    from .utils import export
    export.run(infile, outdir, format)

if __name__ == "__main__":
    main()