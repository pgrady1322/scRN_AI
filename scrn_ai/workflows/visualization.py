"""
UMAP visualization workflow.

Generates UMAP embeddings and visualizations with optional cell type overlays.
"""
import scanpy as sc
import anndata as ad
import pathlib as p
import matplotlib.pyplot as plt


def run_umap(input_file, output_file, color_by="leiden", n_neighbors=15, min_dist=0.1, cell_types=None):
    """
    Generate UMAP visualization for dimensional reduction.
    
    Parameters
    ----------
    input_file : str
        Path to input normalized .h5ad file
    output_file : str
        Path to output image file (.png, .pdf, etc.)
    color_by : str
        Observation key to color by (default: leiden)
    n_neighbors : int
        Number of neighbors for UMAP (default: 15)
    min_dist : float
        Minimum distance for UMAP (default: 0.1)
    cell_types : str, optional
        Path to CSV file with cell type annotations
    """
    print(f"[umap] Loading data from {input_file}")
    adata = ad.read_h5ad(input_file)
    
    # Load cell types if provided
    if cell_types is not None:
        print(f"[umap] Loading cell type annotations from {cell_types}")
        import pandas as pd
        ct_df = pd.read_csv(cell_types, index_col=0)
        # Merge cell types into adata.obs
        if 'cell_type' in ct_df.columns:
            adata.obs['cell_type'] = ct_df['cell_type']
            color_by = 'cell_type'
        elif 'annotation' in ct_df.columns:
            adata.obs['cell_type'] = ct_df['annotation']
            color_by = 'cell_type'
    
    # Compute neighbors if not already done
    if 'neighbors' not in adata.uns:
        print("[umap] Computing PCA...")
        if 'X_pca' not in adata.obsm:
            sc.pp.pca(adata, n_comps=50)
        
        print(f"[umap] Computing neighbors (n_neighbors={n_neighbors})...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    
    # Compute UMAP if not already done
    if 'X_umap' not in adata.obsm:
        print(f"[umap] Computing UMAP (min_dist={min_dist})...")
        sc.tl.umap(adata, min_dist=min_dist)
    
    # Compute leiden clustering if coloring by leiden and not present
    if color_by == 'leiden' and 'leiden' not in adata.obs:
        print("[umap] Computing Leiden clustering...")
        sc.tl.leiden(adata)
    
    # Generate plot
    print(f"[umap] Generating plot colored by '{color_by}'...")
    fig = sc.pl.umap(adata, color=color_by, show=False, return_fig=True)
    
    # Save figure
    output_path = p.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[umap] Saving to {output_file}")
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print("[umap] Complete!")
