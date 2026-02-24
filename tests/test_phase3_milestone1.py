#!/usr/bin/env python3
"""
Comprehensive test suite for Phase 3 Milestone 1 - Configuration Parser

Tests all configuration parser functionality including:
- Loading defaults
- Loading YAML config files
- Environment variable substitution
- CLI overrides
- Validation
- Error handling
"""

import os
import sys
import tempfile
import yaml
from pathlib import Path

# Color codes for output
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
RESET = '\033[0m'

def print_test(name):
    """Print test name."""
    print(f"\n{BLUE}Testing:{RESET} {name}")

def print_pass(msg):
    """Print success message."""
    print(f"  {GREEN}✓{RESET} {msg}")

def print_fail(msg):
    """Print failure message."""
    print(f"  {RED}✗{RESET} {msg}")

def print_warning(msg):
    """Print warning message."""
    print(f"  {YELLOW}⚠{RESET} {msg}")

# Test counter
tests_passed = 0
tests_failed = 0

def run_test(test_func):
    """Run a test function and track results."""
    global tests_passed, tests_failed
    try:
        test_func()
        tests_passed += 1
        return True
    except Exception as e:
        tests_failed += 1
        print_fail(f"Error: {str(e)}")
        return False

# ============================================================================
# TEST SUITE
# ============================================================================

def test_import():
    """Test that ConfigParser can be imported."""
    print_test("Import ConfigParser")
    from scrn_ai.config import ConfigParser
    print_pass("ConfigParser imported successfully")

def test_load_defaults():
    """Test loading default configuration."""
    print_test("Load default configuration")
    from scrn_ai.config import ConfigParser
    
    parser = ConfigParser()
    
    # Check key default values
    assert parser.get('preprocessing.min_genes_per_cell') == 200
    print_pass(f"Default min_genes_per_cell: {parser.get('preprocessing.min_genes_per_cell')}")
    
    assert parser.get('preprocessing.min_cells_per_gene') == 3
    print_pass(f"Default min_cells_per_gene: {parser.get('preprocessing.min_cells_per_gene')}")
    
    assert parser.get('normalization.method') == 'seurat'
    print_pass(f"Default normalization method: {parser.get('normalization.method')}")
    
    assert parser.get('normalization.algorithm') == 'LogNormalize'
    print_pass(f"Default normalization algorithm: {parser.get('normalization.algorithm')}")
    
    assert parser.get('analysis.run_umap') is True
    print_pass(f"Default run_umap: {parser.get('analysis.run_umap')}")
    
    assert parser.get('analysis.umap_n_neighbors') == 15
    print_pass(f"Default UMAP n_neighbors: {parser.get('analysis.umap_n_neighbors')}")

def test_load_sample_config():
    """Test loading the sample configuration file."""
    print_test("Load sample configuration file")
    from scrn_ai.config import ConfigParser
    
    parser = ConfigParser('examples/sample_config.yaml')
    
    # Check loaded values
    assert parser.get('input.matrix_path') == './data/input/dataset.h5ad'
    print_pass(f"Input path: {parser.get('input.matrix_path')}")
    
    assert parser.get('normalization.method') == 'seurat'
    print_pass(f"Normalization method: {parser.get('normalization.method')}")
    
    assert parser.get('output.results_dir') == './output'
    print_pass(f"Results directory: {parser.get('output.results_dir')}")

def test_config_sections():
    """Test getting configuration sections."""
    print_test("Access configuration sections")
    from scrn_ai.config import ConfigParser
    
    parser = ConfigParser('examples/sample_config.yaml')
    
    # Test section getters
    input_config = parser.get_input_config()
    assert isinstance(input_config, dict)
    assert 'matrix_path' in input_config
    print_pass(f"Input config section: {len(input_config)} keys")
    
    preprocess_config = parser.get_preprocessing_config()
    assert isinstance(preprocess_config, dict)
    print_pass(f"Preprocessing config section: {len(preprocess_config)} keys")
    
    norm_config = parser.get_normalization_config()
    assert isinstance(norm_config, dict)
    assert norm_config['method'] == 'seurat'
    print_pass(f"Normalization config section: {len(norm_config)} keys")
    
    analysis_config = parser.get_analysis_config()
    assert isinstance(analysis_config, dict)
    print_pass(f"Analysis config section: {len(analysis_config)} keys")
    
    output_config = parser.get_output_config()
    assert isinstance(output_config, dict)
    print_pass(f"Output config section: {len(output_config)} keys")
    
    exec_config = parser.get_execution_config()
    assert isinstance(exec_config, dict)
    print_pass(f"Execution config section: {len(exec_config)} keys")

def test_env_var_substitution():
    """Test environment variable substitution."""
    print_test("Environment variable substitution")
    from scrn_ai.config import ConfigParser
    
    # Create temp config with env vars
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        config_data = {
            'input': {
                'matrix_path': '${TEST_DATA_DIR}/input.h5ad'
            },
            'output': {
                'results_dir': '${TEST_OUTPUT_DIR:-./default_output}'
            }
        }
        yaml.dump(config_data, f)
        temp_config = f.name
    
    try:
        # Set environment variable
        os.environ['TEST_DATA_DIR'] = '/test/data/path'
        
        parser = ConfigParser(temp_config)
        
        # Check substitution with existing var
        assert parser.get('input.matrix_path') == '/test/data/path/input.h5ad'
        print_pass(f"Substituted ${{TEST_DATA_DIR}}: {parser.get('input.matrix_path')}")
        
        # Check substitution with default (var doesn't exist)
        assert parser.get('output.results_dir') == './default_output'
        print_pass(f"Used default value: {parser.get('output.results_dir')}")
        
        # Cleanup
        del os.environ['TEST_DATA_DIR']
    finally:
        os.unlink(temp_config)

def test_cli_overrides():
    """Test command-line argument overrides."""
    print_test("Command-line overrides")
    from scrn_ai.config import ConfigParser
    
    parser = ConfigParser()
    
    original_min_genes = parser.get('preprocessing.min_genes_per_cell')
    print_pass(f"Original min_genes: {original_min_genes}")
    
    # Apply CLI overrides
    cli_args = {
        'preprocessing.min_genes_per_cell': 1000,
        'normalization.method': 'jmp',
        'analysis.run_umap': False
    }
    
    parser.merge_cli_overrides(cli_args)
    
    # Check overrides were applied
    assert parser.get('preprocessing.min_genes_per_cell') == 1000
    print_pass(f"Overridden min_genes: {parser.get('preprocessing.min_genes_per_cell')}")
    
    assert parser.get('normalization.method') == 'jmp'
    print_pass(f"Overridden normalization method: {parser.get('normalization.method')}")
    
    assert parser.get('analysis.run_umap') is False
    print_pass(f"Overridden run_umap: {parser.get('analysis.run_umap')}")

def test_merge_behavior():
    """Test that user config properly merges with defaults."""
    print_test("Config merging with defaults")
    from scrn_ai.config import ConfigParser
    
    # Create minimal config
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        config_data = {
            'input': {
                'matrix_path': '/custom/path.h5ad'
            },
            'preprocessing': {
                'min_genes_per_cell': 500  # Override default
            }
            # Don't specify normalization - should use defaults
        }
        yaml.dump(config_data, f)
        temp_config = f.name
    
    try:
        parser = ConfigParser(temp_config)
        
        # Check override was applied
        assert parser.get('preprocessing.min_genes_per_cell') == 500
        print_pass(f"User override applied: min_genes = {parser.get('preprocessing.min_genes_per_cell')}")
        
        # Check defaults are still present for unspecified values
        assert parser.get('preprocessing.min_cells_per_gene') == 3  # Default
        print_pass(f"Default preserved: min_cells = {parser.get('preprocessing.min_cells_per_gene')}")
        
        assert parser.get('normalization.method') == 'seurat'  # Default
        print_pass(f"Default preserved: normalization method = {parser.get('normalization.method')}")
    finally:
        os.unlink(temp_config)

def test_to_dict():
    """Test converting config to dictionary."""
    print_test("Export configuration to dict")
    from scrn_ai.config import ConfigParser
    
    parser = ConfigParser()
    config_dict = parser.to_dict()
    
    assert isinstance(config_dict, dict)
    print_pass(f"Config dict has {len(config_dict)} top-level keys")
    
    assert 'input' in config_dict
    assert 'preprocessing' in config_dict
    assert 'normalization' in config_dict
    assert 'analysis' in config_dict
    assert 'output' in config_dict
    assert 'execution' in config_dict
    print_pass("All expected sections present in dict")

def test_file_not_found():
    """Test error handling for missing config file."""
    print_test("Error handling - missing file")
    from scrn_ai.config import ConfigParser
    
    try:
        parser = ConfigParser('/nonexistent/config.yaml')
        print_fail("Should have raised FileNotFoundError")
    except FileNotFoundError as e:
        print_pass(f"Correctly raised FileNotFoundError: {str(e)[:50]}...")

def test_invalid_yaml():
    """Test error handling for invalid YAML."""
    print_test("Error handling - invalid YAML")
    from scrn_ai.config import ConfigParser
    from scrn_ai.config.parser import ConfigValidationError
    
    # Create invalid YAML file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write("invalid: yaml: content:\n  - bad: [unclosed")
        temp_config = f.name
    
    try:
        parser = ConfigParser(temp_config)
        print_fail("Should have raised ConfigValidationError")
    except ConfigValidationError as e:
        print_pass(f"Correctly raised ConfigValidationError: {str(e)[:50]}...")
    finally:
        os.unlink(temp_config)

def test_nested_access():
    """Test nested configuration access with dotted notation."""
    print_test("Nested configuration access")
    from scrn_ai.config import ConfigParser
    
    parser = ConfigParser()
    
    # Test various nesting levels
    value1 = parser.get('preprocessing.min_genes_per_cell')
    assert value1 == 200
    print_pass(f"Level 2 access: {value1}")
    
    value2 = parser.get('analysis.run_umap')
    assert value2 is True
    print_pass(f"Boolean access: {value2}")
    
    # Test with default for non-existent key
    value3 = parser.get('nonexistent.key', default='default_value')
    assert value3 == 'default_value'
    print_pass(f"Default value for missing key: {value3}")

def test_config_files_exist():
    """Test that all required config files exist."""
    print_test("Configuration files existence")
    
    files = [
        'scrn_ai/config/__init__.py',
        'scrn_ai/config/parser.py',
        'scrn_ai/config/defaults.yaml',
        'scrn_ai/config/schema.yaml',
        'examples/sample_config.yaml'
    ]
    
    for file in files:
        path = Path(file)
        assert path.exists(), f"Missing file: {file}"
        print_pass(f"Found: {file}")

# ============================================================================
# MAIN TEST RUNNER
# ============================================================================

def main():
    """Run all tests."""
    print("=" * 70)
    print(f"{BLUE}Phase 3 Milestone 1 - Configuration Parser Test Suite{RESET}")
    print("=" * 70)
    
    tests = [
        test_import,
        test_config_files_exist,
        test_load_defaults,
        test_load_sample_config,
        test_config_sections,
        test_env_var_substitution,
        test_cli_overrides,
        test_merge_behavior,
        test_to_dict,
        test_nested_access,
        test_file_not_found,
        test_invalid_yaml,
    ]
    
    for test in tests:
        run_test(test)
    
    # Print summary
    print("\n" + "=" * 70)
    total = tests_passed + tests_failed
    pass_rate = (tests_passed / total * 100) if total > 0 else 0
    
    print(f"\n{BLUE}Test Summary:{RESET}")
    print(f"  Total tests:  {total}")
    print(f"  {GREEN}Passed:{RESET}       {tests_passed}")
    print(f"  {RED}Failed:{RESET}       {tests_failed}")
    print(f"  Pass rate:    {pass_rate:.1f}%")
    
    if tests_failed == 0:
        print(f"\n{GREEN}🎉 All tests passed! Configuration parser is working correctly.{RESET}")
        return 0
    else:
        print(f"\n{RED}❌ Some tests failed. Please review the errors above.{RESET}")
        return 1

if __name__ == '__main__':
    sys.exit(main())
