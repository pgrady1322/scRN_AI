# Phase 1 Implementation: Command Alignment

This document describes the Phase 1 changes that align the CLI commands with the README documentation.

## Overview

Phase 1 focused on making the actual command-line interface match what is documented in the README. This included:
1. Renaming/aliasing existing commands
2. Creating new wrapper commands
3. Unifying the pseudotime interface
4. Adding preprocessing and UMAP commands

## Changes Made

### 1. New Commands Added

#### `preprocess` Command
**Purpose**: Preprocessing and QC filtering for raw scRNA-seq data  
**Location**: `scrn_ai/workflows/preprocess.py`

```bash
scrn_ai preprocess --input data.h5ad --output processed.h5ad \
    --min-genes 200 --min-cells 3 --max-mito-pct 20
```

**Features**:
- Supports multiple input formats (.mtx, .h5ad, .loom, .csv)
- QC filtering (min/max genes, mitochondrial content)
- Automatic format detection

#### `normalize` Command
**Purpose**: Unified normalization interface  
**Location**: `scrn_ai/workflows/normalization.py`

```bash
scrn_ai normalize --input processed.h5ad --output normalized.h5ad \
    --method seurat --algorithm LogNormalize
```

**Methods Supported**:
- **Seurat**: LogNormalize, SCTransform (via R/Seurat)
- **JMP**: TMM, RLE, UpperQuartile (via R/edgeR)
- **Basic**: log1p, scran, sctransform

**Note**: Legacy `ad_norm` command still available for backward compatibility.

#### `umap` Command
**Purpose**: UMAP visualization generation  
**Location**: `scrn_ai/workflows/visualization.py`

```bash
scrn_ai umap --input normalized.h5ad --output umap.png \
    --color-by leiden --cell-types annotations.csv
```

**Features**:
- Automatic PCA and neighbor computation
- Optional cell type overlay from CSV
- Configurable UMAP parameters

#### `pseudotime` Command
**Purpose**: Unified pseudotime analysis interface  
**Location**: `scrn_ai/workflows/pseudotime.py`

```bash
scrn_ai pseudotime --input normalized.h5ad --output results.h5ad \
    --method dpt --scale small --root-cell CELL_001
```

**Methods**:
- `dpt`: Diffusion Pseudotime
- `diffusion`: Diffusion maps
- `bltsa`: BLTSA trajectory inference
- `via`: VIA/STAVIA (for large-scale)

**Scales**:
- `small`: For datasets <50k cells (uses DPT, BLTSA)
- `large`: For datasets >50k cells (uses VIA/STAVIA)

### 2. Legacy Commands (Preserved)

The following commands remain for backward compatibility:

- `small`: Original small-scale pseudotime workflow
- `large`: Original large-scale VIA/STAVIA workflow
- `ad_norm`: Original normalization utility
- `ad_merge`: Merge multiple AnnData files
- `ad_export`: Export AnnData to various formats

All legacy commands are marked with `[LEGACY]` in their help text.

## File Structure

```
scrn_ai/
├── cli.py                      # Updated CLI with new commands
├── workflows/                  # NEW: Workflow modules
│   ├── __init__.py
│   ├── preprocess.py          # NEW: Preprocessing workflow
│   ├── normalization.py       # NEW: Normalization workflow
│   ├── visualization.py       # NEW: UMAP visualization
│   └── pseudotime.py          # NEW: Unified pseudotime
├── small.py                   # Legacy small-scale workflow
├── large.py                   # Legacy large-scale workflow
└── utils/                     # Existing utilities
    ├── normalization.py       # Original ad_norm implementation
    ├── merge.py
    ├── export.py
    └── plot.py
```

## Dependencies Added

Updated `env.yml` with:
- `r-seurat`: For Seurat normalization
- `bioconductor-edger`: For JMP/edgeR normalization
- `bioconductor-scran`: For scran normalization
- `r-sctransform`: For SCTransform
- `loompy`: For loom file support

## Testing

Run the test suite to verify all commands are accessible:

```bash
python test_phase1.py
```

This will test:
- ✅ All new commands show help properly
- ✅ Legacy commands still work
- ✅ Command-line options are recognized

## Migration Guide

### From Legacy Commands to New Commands

**Old way (ad_norm)**:
```bash
scrn_ai ad_norm -i input.h5ad -o output.h5ad -m log1p
```

**New way (normalize)**:
```bash
scrn_ai normalize --input input.h5ad --output output.h5ad --method log1p
```

**Old way (small)**:
```bash
scrn_ai small -i data.h5ad -s mouse -m dpt -o results.h5ad
```

**New way (pseudotime)**:
```bash
scrn_ai pseudotime --input data.h5ad --output results.h5ad \
    --method dpt --scale small
```

## What's Next: Phase 2

Phase 2 will implement the AItyping module:
- OpenAI API integration
- Marker gene detection
- Pre/post-analysis cell typing
- Confidence scoring
- Cell type annotation output

## Command Reference

### Quick Command Comparison

| Task | Legacy Command | New Command |
|------|----------------|-------------|
| Preprocess | ❌ N/A | `preprocess` |
| Normalize | `ad_norm` | `normalize` |
| UMAP | ❌ N/A | `umap` |
| Pseudotime (small) | `small` | `pseudotime --scale small` |
| Pseudotime (large) | `large` | `pseudotime --scale large` |
| Merge | `ad_merge` | `ad_merge` (unchanged) |
| Export | `ad_export` | `ad_export` (unchanged) |

### Full Command List

```bash
# Show all commands
scrn_ai --help

# New workflow commands
scrn_ai preprocess --help
scrn_ai normalize --help
scrn_ai umap --help
scrn_ai pseudotime --help

# Legacy workflow commands
scrn_ai small --help       # [LEGACY]
scrn_ai large --help       # [LEGACY]
scrn_ai ad_norm --help     # [LEGACY]

# Utility commands
scrn_ai ad_merge --help
scrn_ai ad_export --help
```

## Notes

- All new commands use `--input` and `--output` (long form) for consistency
- Legacy commands preserve their original short-form options (`-i`, `-o`)
- R integration is optional - commands fall back to Python implementations if R packages unavailable
- VIA/STAVIA requires `pyVIA` (installed via pip in env.yml)
