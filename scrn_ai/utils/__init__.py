#!/usr/bin/env python3
"""
scRN_AI v1.0.0

Utility subpackage — re-exports all utility modules.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

_SUBMODULES = [
    "cytetype_client",
    "export",
    "marker_detection",
    "merge",
    "normalization",
    "plot",
]

__all__ = list(_SUBMODULES)


def __getattr__(name: str):
    """Lazy-load utility modules so that ``from scrn_ai import utils``
    does not eagerly import scanpy / matplotlib / etc."""
    if name in _SUBMODULES:
        import importlib

        return importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# scRN_AI v1.0.0
# Any usage is subject to this software's license.
