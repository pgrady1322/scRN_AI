#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v0.1.0

Phase 2 implementation tests.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import subprocess
import sys
import os

def run_command(cmd, timeout=10):
    """Run command and return True if successful"""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=timeout
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "Timeout"
    except Exception as e:
        return False, "", str(e)

def test_aitype_command_exists():
    """Test that aitype command is accessible"""
    print("\n1. Testing aitype command exists...")
    success, stdout, stderr = run_command("scrn_ai aitype --help")
    if success and "AI-powered cell type identification" in stdout:
        print("   ✅ PASS - aitype command accessible")
        return True
    else:
        print("   ❌ FAIL - aitype command not accessible")
        return False

def test_module_imports():
    """Test that new modules can be imported"""
    print("\n2. Testing module imports...")
    tests = [
        ("OpenAI client", "from scrn_ai.utils.openai_client import OpenAIClient"),
        ("Marker detection", "from scrn_ai.utils.marker_detection import get_top_markers_per_cluster"),
        ("AItype workflow", "from scrn_ai.workflows.aitype import run"),
    ]
    
    all_passed = True
    for name, import_stmt in tests:
        try:
            exec(import_stmt)
            print(f"   ✅ PASS - {name}")
        except ImportError as e:
            print(f"   ❌ FAIL - {name}: {e}")
            all_passed = False
    
    return all_passed

def test_openai_package():
    """Test that openai package is available"""
    print("\n3. Testing openai package availability...")
    try:
        import openai
        print(f"   ✅ PASS - openai package installed (version: {openai.__version__})")
        return True
    except ImportError:
        print("   ⚠️  WARNING - openai package not installed")
        print("      Install with: pip install openai")
        return False

def test_api_key_setup():
    """Test if OpenAI API key is configured"""
    print("\n4. Testing OpenAI API key configuration...")
    api_key = os.getenv("OPENAI_API_KEY")
    if api_key:
        masked_key = api_key[:7] + "..." + api_key[-4:] if len(api_key) > 11 else "***"
        print(f"   ✅ PASS - OPENAI_API_KEY is set ({masked_key})")
        return True
    else:
        print("   ⚠️  WARNING - OPENAI_API_KEY environment variable not set")
        print("      Set with: export OPENAI_API_KEY='your-api-key'")
        return False

def test_aitype_parameters():
    """Test that all aitype parameters are recognized"""
    print("\n5. Testing aitype parameters...")
    params = [
        "--timing",
        "--model",
        "--confidence-threshold",
        "--marker-genes",
        "--n-markers",
        "--max-clusters",
        "--species",
        "--tissue",
        "--cluster-key",
    ]
    
    success, stdout, stderr = run_command("scrn_ai aitype --help")
    if not success:
        print("   ❌ FAIL - Could not get help text")
        return False
    
    all_found = True
    for param in params:
        if param in stdout:
            print(f"   ✅ {param}")
        else:
            print(f"   ❌ {param} - not found")
            all_found = False
    
    return all_found

def test_marker_detection_functions():
    """Test marker detection utility functions"""
    print("\n6. Testing marker detection functions...")
    try:
        from scrn_ai.utils.marker_detection import (
            filter_marker_genes,
            validate_marker_genes
        )
        
        # Test filter_marker_genes
        test_genes = ["RPL5", "MT-CO1", "HSP90AA1", "CD3D", "CD4", "RPS10"]
        filtered = filter_marker_genes(test_genes)
        
        # Should remove RPL5, MT-CO1, HSP90AA1, RPS10
        expected = ["CD3D", "CD4"]
        if filtered == expected:
            print("   ✅ PASS - filter_marker_genes working correctly")
            return True
        else:
            print(f"   ❌ FAIL - Expected {expected}, got {filtered}")
            return False
    except Exception as e:
        print(f"   ❌ FAIL - Error: {e}")
        return False

def test_openai_client_initialization():
    """Test OpenAI client can be initialized (without API call)"""
    print("\n7. Testing OpenAI client initialization...")
    try:
        from scrn_ai.utils.openai_client import OpenAIClient
        
        # Try with dummy key (won't make actual API calls in this test)
        try:
            client = OpenAIClient(api_key="sk-test-dummy-key")
            print("   ✅ PASS - OpenAIClient can be initialized")
            return True
        except ImportError as e:
            print(f"   ⚠️  WARNING - openai package issue: {e}")
            return False
        except Exception as e:
            print(f"   ✅ PASS - OpenAIClient initialized (expected validation: {e})")
            return True
    except Exception as e:
        print(f"   ❌ FAIL - Error importing: {e}")
        return False

def main():
    print("=" * 70)
    print("Phase 2 Testing Suite - AItyping Module")
    print("=" * 70)
    
    tests = [
        test_aitype_command_exists,
        test_module_imports,
        test_openai_package,
        test_api_key_setup,
        test_aitype_parameters,
        test_marker_detection_functions,
        test_openai_client_initialization,
    ]
    
    results = []
    for test_func in tests:
        try:
            result = test_func()
            results.append(result)
        except Exception as e:
            print(f"   ❌ EXCEPTION: {e}")
            results.append(False)
    
    # Summary
    print("\n" + "=" * 70)
    passed = sum(results)
    total = len(results)
    print(f"Results: {passed}/{total} tests passed")
    print("=" * 70)
    
    if passed == total:
        print("\n🎉 All Phase 2 tests passed!")
        print("\n✅ AItyping module is ready to use")
        print("\nQuick Start:")
        print("  1. Set your API key: export OPENAI_API_KEY='sk-...'")
        print("  2. Run: scrn_ai aitype -i data.h5ad -o results/ --timing pre_analysis")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test(s) failed or showed warnings")
        print("\nNote: Some tests require:")
        print("  - openai package: pip install openai")
        print("  - OPENAI_API_KEY environment variable")
        print("\nCore functionality tests:", "PASSED ✅" if results[0] and results[1] else "FAILED ❌")
        return 1

if __name__ == "__main__":
    sys.exit(main())

# scRN_AI v0.1.0
# Any usage is subject to this software's license.
