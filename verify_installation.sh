#!/bin/bash
#
# Quick installation verification for scRN_AI
# Run this script to verify your installation is working correctly
#

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║         scRN_AI Installation Verification Test            ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo "❌ Python not found. Please install Python 3.8+ first."
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
echo "✓ Python version: $PYTHON_VERSION"

# Check if scrn_ai is installed
if ! command -v scrn_ai &> /dev/null; then
    echo "⚠️  scrn_ai command not found."
    echo "   Installing package in editable mode..."
    pip install -e . > /dev/null 2>&1
    
    if [ $? -eq 0 ]; then
        echo "✓ Package installed successfully"
    else
        echo "❌ Installation failed. Please run: pip install -e ."
        exit 1
    fi
else
    echo "✓ scrn_ai command found"
fi

echo ""
echo "Running comprehensive tests..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Run the quick test
python quick_test.py

# Check exit code
if [ $? -eq 0 ]; then
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║                  ✅ ALL TESTS PASSED                       ║"
    echo "║                                                            ║"
    echo "║  Your scRN_AI installation is working correctly!          ║"
    echo "║                                                            ║"
    echo "║  Next steps:                                               ║"
    echo "║  • See README.md for usage examples                        ║"
    echo "║  • See TESTING.md for detailed testing guide               ║"
    echo "║  • Try: scrn_ai --help                                  ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    exit 0
else
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║                  ❌ SOME TESTS FAILED                      ║"
    echo "║                                                            ║"
    echo "║  See TESTING.md for troubleshooting steps                  ║"
    echo "║                                                            ║"
    echo "║  Common fixes:                                             ║"
    echo "║  • pip install --force-reinstall -e .                      ║"
    echo "║  • Ensure virtual environment is activated                 ║"
    echo "║  • Check Python version (requires 3.8+)                    ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    exit 1
fi
