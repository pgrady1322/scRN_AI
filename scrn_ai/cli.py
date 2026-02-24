import click
from importlib import import_module

@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def main():
    """Single-cell analysis toolbox."""
    pass

# === Main workflow commands (README-aligned) ===========================

@main.command()
@click.option("--input", "-i", type=click.Path(exists=True), required=True,
              help="Input file (.mtx, .h5ad, .loom, .csv).")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output .h5ad file path.")
@click.option("--min-genes", type=int, default=200,
              help="Minimum genes per cell (default: 200).")
@click.option("--min-cells", type=int, default=3,
              help="Minimum cells per gene (default: 3).")
@click.option("--max-genes", type=int, default=None,
              help="Maximum genes per cell (filter high outliers).")
@click.option("--max-mito-pct", type=float, default=None,
              help="Maximum mitochondrial percentage (e.g., 20.0 for 20%).")
def preprocess(input, output, min_genes, min_cells, max_genes, max_mito_pct):
    """
    Preprocess raw scRNA-seq data with QC filtering.
    
    Performs quality control filtering:
    - Filters cells by gene count thresholds
    - Filters genes by cell count thresholds
    - Optional: Filters high mitochondrial content cells
    
    Accepts multiple input formats: .mtx, .h5ad, .loom, .csv
    """
    from .workflows import preprocess as pp
    pp.run(input, output, min_genes, min_cells, max_genes, max_mito_pct)

@main.command()
@click.option("--input", "-i", type=click.Path(exists=True), required=True,
              help="Input .h5ad file.")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output .h5ad file.")
@click.option("--method", "-m", default="seurat",
              type=click.Choice(["seurat", "jmp", "log1p", "scran", "sctransform"], 
                               case_sensitive=False),
              help="Normalization method (default: seurat).")
@click.option("--algorithm", "-a", default="LogNormalize",
              type=click.Choice(["LogNormalize", "SCTransform", "TMM", "RLE", "UpperQuartile"],
                               case_sensitive=False),
              help="Specific algorithm within method (default: LogNormalize).")
@click.option("--scale-factor", type=float, default=10000,
              help="Scaling factor for normalization (default: 10000).")
def normalize(input, output, method, algorithm, scale_factor):
    """
    Normalize count data using various methods.
    
    Methods:
      - seurat: LogNormalize or SCTransform (Seurat R package)
      - jmp: TMM, RLE, or UpperQuartile (edgeR)
      - log1p: Simple log(x+1) transformation
      - scran: Deconvolution-based size factors
      - sctransform: Variance stabilizing transformation
    """
    from .workflows import normalization as norm
    norm.run(input, output, method, algorithm, scale_factor)

@main.command()
@click.option("--input", "-i", type=click.Path(exists=True), required=True,
              help="Input normalized .h5ad file.")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output image file (.png, .pdf, etc.).")
@click.option("--color-by", "-c", default="leiden",
              help="Observation key to color by (default: leiden).")
@click.option("--n-neighbors", type=int, default=15,
              help="Number of neighbors for UMAP (default: 15).")
@click.option("--min-dist", type=float, default=0.1,
              help="Minimum distance for UMAP (default: 0.1).")
@click.option("--cell-types", type=click.Path(exists=True), default=None,
              help="Optional: CSV with cell type annotations to overlay.")
def umap(input, output, color_by, n_neighbors, min_dist, cell_types):
    """
    Generate UMAP visualization for dimensional reduction.
    
    Computes UMAP embedding and generates visualization plot.
    Can overlay cell type annotations if provided.
    """
    from .workflows import visualization as viz
    viz.run_umap(input, output, color_by, n_neighbors, min_dist, cell_types)

@main.command()
@click.option("--input", "-i", type=click.Path(exists=True), required=True,
              help="Input normalized .h5ad file.")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output directory or .h5ad file.")
@click.option("--method", "-m", default="dpt",
              type=click.Choice(["dpt", "diffusion", "bltsa", "via"], case_sensitive=False),
              help="Pseudotime method (default: dpt).")
@click.option("--scale", type=click.Choice(["small", "large"], case_sensitive=False),
              default="small",
              help="Dataset scale: small (DPT/BLTSA) or large (VIA/STAVIA).")
@click.option("--root-cell", type=str, default=None,
              help="Root cell ID for pseudotime calculation.")
def pseudotime(input, output, method, scale, root_cell):
    """
    Perform pseudotime trajectory analysis.
    
    Methods:
      - dpt: Diffusion Pseudotime (Scanpy)
      - diffusion: Diffusion maps
      - bltsa: Branching trajectory inference
      - via: VIA/STAVIA for large-scale data
      
    Scale:
      - small: DPT, BLTSA for <50k cells
      - large: VIA/STAVIA for >50k cells
    """
    from .workflows import pseudotime as pt
    pt.run(input, output, method, scale, root_cell)

@main.command()
@click.option("--input", "-i", type=click.Path(exists=True), required=True,
              help="Input .h5ad file with normalized expression data.")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output directory for cell type annotations.")
@click.option("--timing", type=click.Choice(["pre_analysis", "post_analysis", "both"], case_sensitive=False),
              default="pre_analysis",
              help="When to run: pre_analysis (before UMAP), post_analysis (after trajectories), or both.")
@click.option("--model", "-m", default="gpt-4",
              type=click.Choice(["gpt-4", "gpt-4-turbo", "gpt-3.5-turbo"], case_sensitive=False),
              help="OpenAI model to use (default: gpt-4).")
@click.option("--confidence-threshold", type=float, default=0.7,
              help="Minimum confidence score to accept (0.0-1.0, default: 0.7).")
@click.option("--marker-genes", type=click.Path(exists=True), default=None,
              help="Optional: Path to custom marker gene CSV file.")
@click.option("--n-markers", type=int, default=10,
              help="Number of top marker genes per cluster (default: 10).")
@click.option("--max-clusters", type=int, default=50,
              help="Maximum number of clusters to annotate (default: 50).")
@click.option("--species", default="human",
              help="Species (human, mouse, etc., default: human).")
@click.option("--tissue", default=None,
              help="Tissue type if known (e.g., brain, blood).")
@click.option("--cluster-key", default="leiden",
              help="Key in adata.obs containing cluster assignments (default: leiden).")
def aitype(input, output, timing, model, confidence_threshold, marker_genes, 
           n_markers, max_clusters, species, tissue, cluster_key):
    """
    AI-powered cell type identification using ChatGPT.
    
    Uses OpenAI's GPT models to automatically identify cell types based on
    marker genes. Can run before analysis (pre_analysis) to guide UMAP/clustering,
    or after analysis (post_analysis) to annotate trajectories.
    
    Requires OPENAI_API_KEY environment variable to be set.
    
    Timing:
      - pre_analysis: Run after normalization, before UMAP/pseudotime
      - post_analysis: Run after trajectories to annotate final clusters
      - both: Run at both stages for validation
    
    Example:
      scrn_ai aitype -i normalized.h5ad -o cell_types/ --timing pre_analysis
    """
    from .workflows import aitype as at
    at.run(
        input_path=input,
        output_dir=output,
        timing=timing,
        model=model,
        confidence_threshold=confidence_threshold,
        marker_genes=marker_genes,
        n_markers=n_markers,
        max_clusters=max_clusters,
        species=species,
        tissue=tissue,
        cluster_key=cluster_key
    )

# === Additional workflow commands =====================================
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
    """Smaller-scale pseudotime analysis workflow (advanced parameters)."""
    mod = import_module(".small", __package__)
    mod.run(infile, species, method, bltsa, outfile)

@main.command()
@click.option("--infile", "-i", type=click.Path(exists=True), required=True,
              help="Merged AnnData (.h5ad) for large datasets.")
@click.option("--root_cell", type=str, default=None,
              help="Cell ID to use as VIA root.")
@click.option("--outfile", "-o", required=True, type=click.Path())
def large(infile, root_cell, outfile):
    """Large-scale trajectory mapping with STAVIA (VIA 2.0) - advanced parameters."""
    mod = import_module(".large", __package__)
    mod.run(infile, root_cell, outfile)

# === Utility commands ==================================================
@main.command()
@click.option("--infile", "-i", type=click.Path(exists=True), required=True)
@click.option("--outfile", "-o", type=click.Path(), required=True)
@click.option("--method", "-m", default="log1p",
              type=click.Choice(["log1p", "scran", "sctransform", "size_factor"]))
def ad_norm(infile, outfile, method):
    """Normalize counts within an AnnData object (basic methods)."""
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
@click.option("--format", "-f", "fmt",
              type=click.Choice(["loom", "mtx", "csv"], case_sensitive=False),
              required=True)
def ad_export(infile, outdir, fmt):
    """Export AnnData to loom / mtx / csv."""
    from .utils import export
    export.run(infile, outdir, fmt)

if __name__ == "__main__":
    main()