#!/usr/bin/env python3
"""
Quick reference test for scrn_ai Phase 1 commands
Tests all commands are accessible and have proper help text
"""

import subprocess
import sys

def run_command(cmd):
    """Run command and return True if successful"""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=5
        )
        return result.returncode == 0
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return False

def main():
    print("=" * 60)
    print("scrn_ai Phase 1 - Quick Reference Test")
    print("=" * 60)
    
    tests = [
        ("Main help", "scrn_ai --help"),
        ("Preprocess", "scrn_ai preprocess --help"),
        ("Normalize", "scrn_ai normalize --help"),
        ("UMAP", "scrn_ai umap --help"),
        ("Pseudotime", "scrn_ai pseudotime --help"),
        ("AItyping", "scrn_ai aitype --help"),  # NEW: Phase 2
        ("Merge", "scrn_ai ad-merge --help"),
        ("Export", "scrn_ai ad-export --help"),
        ("Basic norm", "scrn_ai ad-norm --help"),
        ("Small workflow", "scrn_ai small --help"),
        ("Large workflow", "scrn_ai large --help"),
    ]
    
    passed = 0
    failed = 0
    
    for name, cmd in tests:
        print(f"\n Testing: {name}")
        print(f"  Command: {cmd}")
        if run_command(cmd):
            print(f"  ✅ PASS")
            passed += 1
        else:
            print(f"  ❌ FAIL")
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Results: {passed}/{len(tests)} passed, {failed}/{len(tests)} failed")
    print("=" * 60)
    
    if failed == 0:
        print("\n🎉 All commands working correctly (Phase 1 + Phase 2)!")
        print("\n📋 Command Reference:")
        print("  Phase 1 Commands:")
        print("    scrn_ai preprocess   - QC filtering and preprocessing")
        print("    scrn_ai normalize    - Unified normalization (Seurat/JMP/basic)")
        print("    scrn_ai umap         - UMAP visualization")
        print("    scrn_ai pseudotime   - Trajectory analysis (DPT/BLTSA/VIA)")
        print("\n  Phase 2 Commands:")
        print("    scrn_ai aitype       - AI-powered cell type identification ✨")
        print("\n  Utility Commands:")
        print("    scrn_ai ad-merge     - Merge AnnData files")
        print("    scrn_ai ad-export    - Export to loom/mtx/csv")
        print("    scrn_ai ad-norm      - Basic normalization")
        print("    scrn_ai small        - Advanced small-scale workflow")
        print("    scrn_ai large        - Advanced large-scale workflow")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed. Please check installation.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
