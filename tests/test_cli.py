#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

CLI command tests — validates every Click command is accessible and
prints help without errors.  Uses Click's CliRunner (no subprocess).

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import pytest


# ── Phase 1 commands ─────────────────────────────────────────────────────

class TestPhase1CLI:
    """All Phase 1 workflow + utility commands respond to --help."""

    @pytest.mark.parametrize(
        "cmd",
        ["preprocess", "normalize", "umap", "pseudotime"],
        ids=lambda c: f"help-{c}",
    )
    def test_workflow_help(self, runner, cli, cmd):
        result = runner.invoke(cli, [cmd, "--help"])
        assert result.exit_code == 0, result.output
        assert cmd in result.output.lower() or "--help" in result.output

    @pytest.mark.parametrize(
        "cmd",
        ["small", "large", "ad-norm", "ad-merge", "ad-export"],
        ids=lambda c: f"help-{c}",
    )
    def test_legacy_and_utility_help(self, runner, cli, cmd):
        result = runner.invoke(cli, [cmd, "--help"])
        assert result.exit_code == 0, result.output

    def test_main_help(self, runner, cli):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "Single-cell analysis toolbox" in result.output


# ── Phase 2 commands ─────────────────────────────────────────────────────

class TestPhase2CLI:
    """AI-typing command is registered and advertises the right options."""

    def test_aitype_help(self, runner, cli):
        result = runner.invoke(cli, ["aitype", "--help"])
        assert result.exit_code == 0
        assert "AI-powered cell type identification" in result.output

    @pytest.mark.parametrize(
        "flag",
        [
            "--timing",
            "--model",
            "--confidence-threshold",
            "--marker-genes",
            "--n-markers",
            "--max-clusters",
            "--species",
            "--tissue",
            "--cluster-key",
        ],
    )
    def test_aitype_options(self, runner, cli, flag):
        result = runner.invoke(cli, ["aitype", "--help"])
        assert flag in result.output, f"Missing option {flag}"


# ── Preprocess options ───────────────────────────────────────────────────

class TestPreprocessOptions:
    """Preprocess command exposes expected QC flags."""

    @pytest.mark.parametrize(
        "flag",
        ["--min-genes", "--min-cells", "--max-genes", "--max-mito-pct"],
    )
    def test_preprocess_options(self, runner, cli, flag):
        result = runner.invoke(cli, ["preprocess", "--help"])
        assert flag in result.output


# ── Normalize options ────────────────────────────────────────────────────

class TestNormalizeOptions:
    """Normalize command lists all supported methods."""

    @pytest.mark.parametrize(
        "method",
        ["seurat", "jmp", "log1p", "scran", "sctransform"],
    )
    def test_normalize_methods(self, runner, cli, method):
        result = runner.invoke(cli, ["normalize", "--help"])
        assert method in result.output


# ── Pseudotime options ───────────────────────────────────────────────────

class TestPseudotimeOptions:
    """Pseudotime command lists all analysis methods and scales."""

    @pytest.mark.parametrize("method", ["dpt", "diffusion", "bltsa", "via"])
    def test_pseudotime_methods(self, runner, cli, method):
        result = runner.invoke(cli, ["pseudotime", "--help"])
        assert method in result.output

    @pytest.mark.parametrize("scale", ["small", "large"])
    def test_pseudotime_scales(self, runner, cli, scale):
        result = runner.invoke(cli, ["pseudotime", "--help"])
        assert scale in result.output

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
