#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v1.0.0

Package installation and dependency configuration.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

from setuptools import setup, find_packages

setup(
    name="scrn_ai",
    version="1.0.0",
    description="Single-cell RNA-seq analysis toolkit with AI-powered cell typing",
    author="Patrick Grady",
    author_email="",
    url="https://github.com/patrickgrady/scRN_AI",
    packages=find_packages(),
    include_package_data=True,
    python_requires=">=3.8",
    install_requires=[
        "click>=8.0",
        "scanpy>=1.9",
        "numpy>=1.20",
        "pandas>=1.3",
        "scipy>=1.7",
        "anndata>=0.8",
        "scikit-learn>=1.0",
        "matplotlib>=3.4",
        "seaborn>=0.11",
        "leidenalg>=0.8",
        "pyVIA>=0.1",
        "rpy2>=3.5",
        "loompy>=3.0",
        "openai>=1.0.0",  # AI-powered cell typing
    ],
    entry_points={
        "console_scripts": [
            "scrn_ai=scrn_ai.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)

# scRN_AI v1.0.0
# Any usage is subject to this software's license.
