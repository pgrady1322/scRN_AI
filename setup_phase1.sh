#!/bin/bash
# Quick installation and test script for Phase 1 changes

echo "=================================="
echo "Phase 1: Command Alignment Setup"
echo "=================================="

# Check if in correct directory
if [ ! -f "scrn_ai/cli.py" ]; then
    echo "Error: Must run from scRN_AI repository root"
    exit 1
fi

# Create conda environment (optional - user can skip)
echo ""
echo "Step 1: Creating conda environment (optional - press Ctrl+C to skip)"
echo "---------------------------------------------------------------"
read -p "Create new conda environment? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    conda env create -f env.yml -n scrn_ai
    echo "Activate with: conda activate scrn_ai"
fi

# Test commands
echo ""
echo "Step 2: Testing Phase 1 commands"
echo "---------------------------------------------------------------"
python test_phase1.py

# Show example usage
echo ""
echo "=================================="
echo "Phase 1 Implementation Complete!"
echo "=================================="
echo ""
echo "New commands available:"
echo "  - scrn_ai preprocess    (QC filtering)"
echo "  - scrn_ai normalize     (Unified normalization)"
echo "  - scrn_ai umap          (Visualization)"
echo "  - scrn_ai pseudotime    (Trajectory analysis)"
echo ""
echo "Legacy commands preserved:"
echo "  - scrn_ai small         (Legacy pseudotime)"
echo "  - scrn_ai large         (Legacy VIA)"
echo "  - scrn_ai ad_norm       (Legacy normalization)"
echo "  - scrn_ai ad_merge      (Merge AnnData)"
echo "  - scrn_ai ad_export     (Export formats)"
echo ""
echo "Try: python -m scrn_ai.cli --help"
echo ""
echo "See PHASE1_IMPLEMENTATION.md for full documentation"
echo "=================================="
