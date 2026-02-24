#!/usr/bin/env python3
"""
scRN_AI v1.0.0

Tests for YAML configuration parser.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import os

import pytest
import yaml

from scrn_ai.config import ConfigParser
from scrn_ai.config.parser import ConfigValidationError


class TestConfigParser:
    """Test cases for ConfigParser class."""

    def test_load_defaults_only(self):
        """Test loading with defaults only (no config file)."""
        parser = ConfigParser()

        # Should have default values
        assert parser.get("preprocessing.min_genes_per_cell") == 200
        assert parser.get("normalization.method") == "seurat"
        assert parser.get("analysis.run_umap") is True

    def test_load_config_file(self, tmp_path):
        """Test loading from YAML config file."""
        # Create temporary config file
        config_file = tmp_path / "test_config.yaml"
        config_data = {
            "input": {"matrix_path": "/path/to/data.h5ad"},
            "preprocessing": {
                "min_genes_per_cell": 500  # Override default
            },
            "output": {"results_dir": "./test_output"},
        }

        with open(config_file, "w") as f:
            yaml.dump(config_data, f)

        # Load config
        parser = ConfigParser(str(config_file))

        # Check overrides
        assert parser.get("preprocessing.min_genes_per_cell") == 500

        # Check defaults are still present
        assert parser.get("preprocessing.min_cells_per_gene") == 3
        assert parser.get("normalization.method") == "seurat"

    def test_env_var_substitution(self, tmp_path):
        """Test environment variable substitution."""
        # Set test environment variable
        os.environ["TEST_DATA_DIR"] = "/test/data"
        os.environ["TEST_VAR"] = "test_value"

        config_file = tmp_path / "test_config.yaml"
        config_data = {
            "input": {"matrix_path": "${TEST_DATA_DIR}/input.h5ad"},
            "output": {"results_dir": "${TEST_OUTPUT_DIR:-./default_output}"},
        }

        with open(config_file, "w") as f:
            yaml.dump(config_data, f)

        parser = ConfigParser(str(config_file))

        # Check substitution with existing var
        assert parser.get("input.matrix_path") == "/test/data/input.h5ad"

        # Check substitution with default value
        assert parser.get("output.results_dir") == "./default_output"

        # Cleanup
        del os.environ["TEST_DATA_DIR"]
        del os.environ["TEST_VAR"]

    def test_merge_cli_overrides(self):
        """Test command-line argument overrides."""
        parser = ConfigParser()

        # Override with CLI args
        cli_args = {
            "preprocessing.min_genes_per_cell": 1000,
            "normalization.method": "jmp",
            "analysis.run_umap": False,
        }

        parser.merge_cli_overrides(cli_args)

        # Check overrides were applied
        assert parser.get("preprocessing.min_genes_per_cell") == 1000
        assert parser.get("normalization.method") == "jmp"
        assert parser.get("analysis.run_umap") is False

    def test_get_config_sections(self, tmp_path):
        """Test getting specific config sections."""
        config_file = tmp_path / "test_config.yaml"
        config_data = {
            "input": {"matrix_path": "/data.h5ad"},
            "preprocessing": {"min_genes_per_cell": 300},
            "normalization": {"method": "log1p"},
            "output": {"results_dir": "./output"},
        }

        with open(config_file, "w") as f:
            yaml.dump(config_data, f)

        parser = ConfigParser(str(config_file))

        # Test section getters
        input_config = parser.get_input_config()
        assert input_config["matrix_path"] == "/data.h5ad"

        preprocess_config = parser.get_preprocessing_config()
        assert preprocess_config["min_genes_per_cell"] == 300

        norm_config = parser.get_normalization_config()
        assert norm_config["method"] == "log1p"

    def test_validation_required_fields(self, tmp_path):
        """Test validation catches missing required fields."""
        config_file = tmp_path / "invalid_config.yaml"
        config_data = {
            # Missing required 'input.matrix_path'
            "preprocessing": {"min_genes_per_cell": 200}
        }

        with open(config_file, "w") as f:
            yaml.dump(config_data, f)

        parser = ConfigParser(str(config_file))

        # Should raise validation error
        with pytest.raises(ConfigValidationError):
            parser.validate()

    def test_to_dict(self):
        """Test converting config to dictionary."""
        parser = ConfigParser()
        config_dict = parser.to_dict()

        assert isinstance(config_dict, dict)
        assert "input" in config_dict
        assert "preprocessing" in config_dict
        assert "normalization" in config_dict

    def test_nonexistent_config_file(self):
        """Test error handling for missing config file."""
        with pytest.raises(FileNotFoundError):
            ConfigParser("/nonexistent/config.yaml")

    def test_invalid_yaml(self, tmp_path):
        """Test error handling for invalid YAML."""
        config_file = tmp_path / "invalid.yaml"

        # Write invalid YAML
        with open(config_file, "w") as f:
            f.write("invalid: yaml: content:\n  - bad")

        with pytest.raises(ConfigValidationError):
            ConfigParser(str(config_file))


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
