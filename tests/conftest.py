#!/usr/bin/env python3
"""
scRN_AI v1.0.0

Shared pytest fixtures.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import pytest
import yaml
from click.testing import CliRunner

from scrn_ai.cli import main as cli_main


@pytest.fixture
def runner():
    """Click CliRunner for invoking CLI commands without subprocess."""
    return CliRunner()


@pytest.fixture
def cli():
    """The top-level Click group."""
    return cli_main


@pytest.fixture
def sample_config(tmp_path):
    """Write a minimal valid config YAML and return its path."""
    cfg = {
        "input": {"matrix_path": "/tmp/test.h5ad"},
        "preprocessing": {"min_genes_per_cell": 500},
        "normalization": {"method": "log1p"},
        "output": {"results_dir": str(tmp_path / "output")},
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.dump(cfg))
    return str(path)


# scRN_AI v1.0.0
# Any usage is subject to this software's license.
