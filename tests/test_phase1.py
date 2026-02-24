#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v0.1.0

Phase 1 implementation tests.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import subprocess
import sys

def test_command(cmd_args, description):
    """Test a command and report results."""
    print(f"\n{'='*70}")
    print(f"Testing: {description}")
    print(f"Command: {' '.join(cmd_args)}")
    print(f"{'='*70}")
    
    try:
        result = subprocess.run(
            cmd_args,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        if result.returncode == 0:
            print("✅ PASS - Command accessible")
            if result.stdout:
                print("\nOutput preview:")
                print(result.stdout[:500])
        else:
            print("❌ FAIL - Command returned non-zero exit code")
            if result.stderr:
                print("\nError:")
                print(result.stderr[:500])
        
        return result.returncode == 0
        
    except subprocess.TimeoutExpired:
        print("⚠️  TIMEOUT - Command took too long")
        return False
    except FileNotFoundError:
        print("❌ FAIL - Command not found")
        return False
    except Exception as e:
        print(f"❌ FAIL - Exception: {e}")
        return False

def main():
    """Run all command tests."""
    print("\n" + "="*70)
    print("Phase 1 Command Alignment Test Suite")
    print("="*70)
    
    tests = [
        # Main help
        (["python", "-m", "scrn_ai.cli", "--help"], "Main CLI help"),
        
        # New commands
        (["python", "-m", "scrn_ai.cli", "preprocess", "--help"], "preprocess command"),
        (["python", "-m", "scrn_ai.cli", "normalize", "--help"], "normalize command"),
        (["python", "-m", "scrn_ai.cli", "umap", "--help"], "umap command"),
        (["python", "-m", "scrn_ai.cli", "pseudotime", "--help"], "pseudotime command"),
        
        # Legacy commands (should still work)
        (["python", "-m", "scrn_ai.cli", "small", "--help"], "small command (legacy)"),
        (["python", "-m", "scrn_ai.cli", "large", "--help"], "large command (legacy)"),
        (["python", "-m", "scrn_ai.cli", "ad_norm", "--help"], "ad_norm command (legacy)"),
        (["python", "-m", "scrn_ai.cli", "ad_merge", "--help"], "ad_merge command"),
        (["python", "-m", "scrn_ai.cli", "ad_export", "--help"], "ad_export command"),
    ]
    
    results = []
    for cmd_args, description in tests:
        passed = test_command(cmd_args, description)
        results.append((description, passed))
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    passed_count = sum(1 for _, passed in results if passed)
    total_count = len(results)
    
    for description, passed in results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status} - {description}")
    
    print(f"\n{passed_count}/{total_count} tests passed")
    print("="*70)
    
    return 0 if passed_count == total_count else 1

if __name__ == "__main__":
    sys.exit(main())

# scRN_AI v0.1.0
# Any usage is subject to this software's license.
