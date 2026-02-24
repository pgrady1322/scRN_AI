# Testing Guide for scRN_AI

This guide explains how to verify your scRN_AI installation and test all modules.

---

## Quick Installation Verification

After installing scRN_AI, run this simple test to verify everything is working:

```bash
# Navigate to the scRN_AI directory
cd /path/to/scRN_AI

# Run the quick test
python quick_test.py
```

**Expected Output**:
```
============================================================
scrn_ai Phase 1 - Quick Reference Test
============================================================

 Testing: Main help
  Command: scrn_ai --help
  ✅ PASS

 Testing: Preprocess
  Command: scrn_ai preprocess --help
  ✅ PASS

[... 9 more tests ...]

============================================================
Results: 11/11 passed, 0/11 failed
============================================================

🎉 All commands working correctly (Phase 1 + Phase 2)!
```

If you see **11/11 passed**, your installation is correct! ✅

---

## Detailed Testing Options

### 1. Quick Test (Recommended)

**Purpose**: Verify all CLI commands are accessible  
**Time**: ~5 seconds  
**File**: `quick_test.py`

```bash
python quick_test.py
```

**What it tests**:
- ✅ All 11 CLI commands respond to `--help`
- ✅ Entry points are correctly installed
- ✅ Basic command structure is intact

---

### 2. Phase 1 Test Suite

**Purpose**: Test Phase 1 workflow modules in detail  
**Time**: ~10 seconds  
**File**: `test_phase1.py`

```bash
python test_phase1.py
```

**What it tests**:
- ✅ Main CLI help accessible
- ✅ Preprocess command
- ✅ Normalize command
- ✅ UMAP command
- ✅ Pseudotime command
- ✅ Utility commands (merge, export, norm)
- ✅ Legacy commands (small, large)

---

### 3. Phase 2 Test Suite

**Purpose**: Test Phase 2 AI typing modules  
**Time**: ~10 seconds  
**File**: `test_phase2.py`

```bash
python test_phase2.py
```

**What it tests**:
- ✅ AItyping command exists
- ✅ OpenAI client module imports
- ✅ Marker detection module imports
- ✅ AItyping workflow imports
- ⚠️  OpenAI package availability
- ⚠️  API key configuration (optional)

**Note**: Some tests will show warnings if OpenAI API key is not set, but this is normal for basic verification.

---

### 4. Module Import Test

**Purpose**: Verify all Python modules can be imported  
**Time**: ~2 seconds

```bash
python -c "
from scrn_ai.utils.openai_client import OpenAIClient
from scrn_ai.utils.marker_detection import get_top_markers_per_cluster
from scrn_ai.workflows.aitype import run
from scrn_ai.workflows.preprocess import run as preprocess_run
from scrn_ai.workflows.normalization import run as normalize_run
from scrn_ai.workflows.visualization import run as viz_run
from scrn_ai.workflows.pseudotime import run as pseudo_run
print('✅ All modules import successfully!')
"
```

---

## Testing in Different Environments

### Testing in Virtual Environment

```bash
# Activate your virtual environment
source .venv/bin/activate  # macOS/Linux
# or
.venv\Scripts\activate     # Windows

# Run tests
python quick_test.py
```

### Testing in Conda Environment

```bash
# Activate conda environment
conda activate scrn_ai

# Run tests
python quick_test.py
```

### Testing in Docker Container

```bash
# Build the Docker image
docker build -t scrn_ai:0.1 .

# Run tests inside container
docker run --rm scrn_ai:0.1 python /app/quick_test.py
```

---

## Installation Verification Checklist

Use this checklist to verify your installation:

### Basic Installation ✓

- [ ] Python 3.8+ installed (`python --version`)
- [ ] Package installed (`pip list | grep sc-toolkit`)
- [ ] `scrn_ai` command available (`which scrn_ai`)
- [ ] Quick test passes (`python quick_test.py`)

### Phase 1 Modules ✓

- [ ] Preprocess command works (`scrn_ai preprocess --help`)
- [ ] Normalize command works (`scrn_ai normalize --help`)
- [ ] UMAP command works (`scrn_ai umap --help`)
- [ ] Pseudotime command works (`scrn_ai pseudotime --help`)

### Phase 2 Modules ✓

- [ ] AItyping command works (`scrn_ai aitype --help`)
- [ ] OpenAI package installed (`python -c "import openai"`)
- [ ] OpenAI client imports (`python -c "from scrn_ai.utils.openai_client import OpenAIClient"`)
- [ ] Marker detection imports (`python -c "from scrn_ai.utils.marker_detection import get_top_markers_per_cluster"`)

### Optional (For Production Use)

- [ ] OpenAI API key set (`echo $OPENAI_API_KEY`)
- [ ] R installed (for Seurat/JMP normalization) (`R --version`)
- [ ] Test data available (`.h5ad` files for real testing)

---

## Troubleshooting Test Failures

### Issue: "scrn_ai: command not found"

**Cause**: Package not installed or virtual environment not activated

**Solution**:
```bash
# Ensure you're in the project directory
cd /path/to/scRN_AI

# Install package in editable mode
pip install -e .

# OR activate your environment
source .venv/bin/activate
```

---

### Issue: "ModuleNotFoundError: No module named 'scrn_ai'"

**Cause**: Package not installed in current Python environment

**Solution**:
```bash
# Check which Python you're using
which python

# Install package
pip install -e /path/to/scRN_AI

# Verify installation
pip list | grep sc-toolkit
```

---

### Issue: "No module named 'openai'"

**Cause**: OpenAI package not installed (Phase 2 requirement)

**Solution**:
```bash
# Install openai package
pip install openai>=1.0.0

# Verify
python -c "import openai; print(openai.__version__)"
```

---

### Issue: Test shows "X/11 failed"

**Cause**: Some commands are not accessible

**Solution**:
```bash
# Reinstall package
pip install --force-reinstall -e .

# Clear Python cache
find . -type d -name __pycache__ -exec rm -rf {} +

# Try again
python quick_test.py
```

---

### Issue: "OPENAI_API_KEY environment variable not set"

**Cause**: OpenAI API key not configured (this is a WARNING, not an error)

**Solution** (optional, only needed for actual AI typing):
```bash
# Set API key for current session
export OPENAI_API_KEY="your-api-key-here"

# OR add to your shell profile for persistence
echo 'export OPENAI_API_KEY="your-api-key-here"' >> ~/.bashrc
# or ~/.zshrc for zsh
```

---

## Creating Test Data

If you want to test with actual data but don't have any, you can generate synthetic test data:

```python
import scanpy as sc
import numpy as np

# Generate synthetic data
adata = sc.AnnData(
    X=np.random.negative_binomial(5, 0.3, (1000, 2000))  # 1000 cells, 2000 genes
)
adata.obs_names = [f"Cell_{i}" for i in range(1000)]
adata.var_names = [f"Gene_{i}" for i in range(2000)]

# Add some metadata
adata.obs['n_counts'] = adata.X.sum(axis=1)
adata.var['n_cells'] = (adata.X > 0).sum(axis=0)

# Save
adata.write('test_data.h5ad')
print("✅ Created test_data.h5ad with 1000 cells and 2000 genes")
```

Then test the pipeline:
```bash
# Preprocess
scrn_ai preprocess \
    --input test_data.h5ad \
    --output processed.h5ad \
    --min-genes 100 \
    --min-cells 3

# Normalize
scrn_ai normalize \
    --input processed.h5ad \
    --output normalized.h5ad \
    --method log1p

# UMAP
scrn_ai umap \
    --input normalized.h5ad \
    --output umap_plot.png
```

---

## Automated Testing (For Developers)

If you're contributing to scRN_AI, run all tests before committing:

```bash
# Run all test suites
python quick_test.py && \
python test_phase1.py && \
python test_phase2.py

# Check for import errors
python -c "
import scrn_ai
from scrn_ai.workflows import preprocess, normalization, visualization, pseudotime, aitype
from scrn_ai.utils import openai_client, marker_detection
print('✅ All modules import successfully')
"

# Verify entry points
scrn_ai --help
scrn_ai preprocess --help
scrn_ai normalize --help
scrn_ai umap --help
scrn_ai pseudotime --help
scrn_ai aitype --help
```

---

## Continuous Integration (CI/CD)

For automated testing in CI/CD pipelines:

### GitHub Actions Example

```yaml
name: Test scRN_AI

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.11'
    
    - name: Install dependencies
      run: |
        pip install -e .
    
    - name: Run quick test
      run: |
        python quick_test.py
    
    - name: Run Phase 1 tests
      run: |
        python test_phase1.py
    
    - name: Run Phase 2 tests
      run: |
        python test_phase2.py
      continue-on-error: true  # OpenAI API key may not be set
```

---

## Expected Test Results

### All Tests Passing ✅

```
Quick Test:        11/11 passed
Phase 1 Tests:     7/10 passed (3 utility commands may use different invocation)
Phase 2 Tests:     2/7 passed (warnings for API key are normal)
Module Imports:    All successful
```

### Acceptable Warnings ⚠️

These warnings are **normal** and don't indicate problems:

- ⚠️  `OPENAI_API_KEY environment variable not set` - Only needed for actual AI typing
- ⚠️  `openai package issue` - If openai is installed but API key not set
- ⚠️  `RuntimeWarning: 'scrn_ai.cli' found in sys.modules` - Harmless import warning

---

## Getting Help

If tests fail and you can't resolve the issue:

1. **Check versions**:
   ```bash
   python --version  # Should be 3.8+
   pip --version
   pip list | grep -E "scanpy|anndata|openai|click"
   ```

2. **Check installation**:
   ```bash
   pip show sc-toolkit
   which scrn_ai
   ```

3. **Check environment**:
   ```bash
   echo $PATH
   echo $PYTHONPATH
   ```

4. **Review full output**:
   ```bash
   python quick_test.py > test_output.txt 2>&1
   cat test_output.txt
   ```

5. **Create an issue**: https://github.com/[your-org]/scRN_AI/issues
   - Include Python version
   - Include test output
   - Include `pip list` output

---

## Summary

**Quick Verification** (30 seconds):
```bash
cd scRN_AI
pip install -e .
python quick_test.py
```

**Expected**: `11/11 passed` ✅

**Full Verification** (2 minutes):
```bash
python quick_test.py
python test_phase1.py
python test_phase2.py
```

**Everything working?** You're ready to analyze single-cell data! 🎉

For usage examples, see the main [README.md](README.md).
