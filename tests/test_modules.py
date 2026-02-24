#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

Unit tests for importability and core module behaviour.

Tests cover:
  - Package-level imports
  - Workflow module imports
  - Utility module imports
  - Marker-detection filtering logic
  - OpenAI client initialisation (no API calls)
  - Version metadata

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import pytest


# ── Package metadata ─────────────────────────────────────────────────────

class TestPackageMetadata:
    """Package exposes expected version and metadata."""

    def test_version(self):
        import scrn_ai
        assert scrn_ai.__version__ == "1.0.0"

    def test_author(self):
        import scrn_ai
        assert scrn_ai.__author__ == "Patrick Grady"


# ── Core imports ─────────────────────────────────────────────────────────

class TestCoreImports:
    """All top-level sub-packages can be imported."""

    def test_import_package(self):
        import scrn_ai  # noqa: F401

    def test_import_cli(self):
        from scrn_ai import cli  # noqa: F401

    def test_import_utils(self):
        from scrn_ai import utils  # noqa: F401

    def test_import_workflows(self):
        from scrn_ai import workflows  # noqa: F401


# ── Workflow module imports ──────────────────────────────────────────────

class TestWorkflowImports:
    """Each workflow module can be imported and exposes its main callable."""

    @pytest.mark.parametrize(
        "module,func",
        [
            ("preprocess", "run"),
            ("normalization", "run"),
            ("visualization", "run_umap"),
            ("pseudotime", "run"),
            ("aitype", "run"),
        ],
    )
    def test_workflow_import_and_callable(self, module, func):
        mod = __import__(f"scrn_ai.workflows.{module}", fromlist=[func])
        assert callable(getattr(mod, func))


# ── Utility module imports ───────────────────────────────────────────────

class TestUtilityImports:
    """Each utility module can be imported."""

    def test_openai_client(self):
        from scrn_ai.utils.openai_client import OpenAIClient  # noqa: F401

    def test_marker_detection(self):
        from scrn_ai.utils.marker_detection import get_top_markers_per_cluster  # noqa: F401

    def test_export(self):
        from scrn_ai.utils.export import run  # noqa: F401

    def test_merge(self):
        from scrn_ai.utils.merge import run  # noqa: F401

    def test_plot(self):
        from scrn_ai.utils.plot import (
            umap_pseudotime,
            qc_violin,
            qc_scatter,
            dotplot,
            stacked_violin,
            pseudotime_heatmap,
        )  # noqa: F401


# ── Config imports ───────────────────────────────────────────────────────

class TestConfigImports:
    """Config sub-package exposes ConfigParser and validation error."""

    def test_config_parser(self):
        from scrn_ai.config import ConfigParser  # noqa: F401

    def test_validation_error(self):
        from scrn_ai.config.parser import ConfigValidationError  # noqa: F401


# ── Marker detection unit tests ──────────────────────────────────────────

class TestMarkerDetection:
    """Unit tests for marker gene filtering logic."""

    def test_filter_removes_ribosomal(self):
        from scrn_ai.utils.marker_detection import filter_marker_genes

        genes = ["RPL5", "RPS10", "CD3D"]
        result = filter_marker_genes(genes)
        assert "RPL5" not in result
        assert "RPS10" not in result
        assert "CD3D" in result

    def test_filter_removes_mitochondrial(self):
        from scrn_ai.utils.marker_detection import filter_marker_genes

        genes = ["MT-CO1", "CD4"]
        result = filter_marker_genes(genes)
        assert "MT-CO1" not in result
        assert "CD4" in result

    def test_filter_removes_hsp(self):
        from scrn_ai.utils.marker_detection import filter_marker_genes

        genes = ["HSP90AA1", "CD8A"]
        result = filter_marker_genes(genes)
        assert "HSP90AA1" not in result
        assert "CD8A" in result

    def test_filter_combined(self):
        from scrn_ai.utils.marker_detection import filter_marker_genes

        genes = ["RPL5", "MT-CO1", "HSP90AA1", "CD3D", "CD4", "RPS10"]
        result = filter_marker_genes(genes)
        assert result == ["CD3D", "CD4"]

    def test_filter_keeps_all_when_disabled(self):
        from scrn_ai.utils.marker_detection import filter_marker_genes

        genes = ["RPL5", "MT-CO1", "HSP90AA1"]
        result = filter_marker_genes(
            genes,
            remove_ribosomal=False,
            remove_mitochondrial=False,
            remove_hsp=False,
        )
        assert result == genes

    def test_validate_marker_genes(self):
        """validate_marker_genes is importable and callable."""
        from scrn_ai.utils.marker_detection import validate_marker_genes

        assert callable(validate_marker_genes)


# ── OpenAI client ────────────────────────────────────────────────────────

class TestOpenAIClient:
    """OpenAIClient can be initialised without making API calls."""

    def test_init_with_dummy_key(self):
        from scrn_ai.utils.openai_client import OpenAIClient

        # Should not raise — no actual API call happens until query time
        try:
            client = OpenAIClient(api_key="sk-test-dummy-key-00000000")
            assert client is not None
        except ImportError:
            pytest.skip("openai package not installed")

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
