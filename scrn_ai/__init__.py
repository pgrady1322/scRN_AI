#!/usr/bin/env python3
"""
scRN_AI v1.0.0

Package initialization and version metadata.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

__version__ = "1.0.0"
__author__ = "Patrick Grady"
__email__ = ""

__all__ = ["cli", "utils", "workflows"]


def __getattr__(name: str):
    """Lazy-load subpackages so that ``import scrn_ai`` does not pull in
    heavy scientific dependencies (scanpy, matplotlib, …) until they are
    actually needed."""
    if name in __all__:
        import importlib

        return importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
