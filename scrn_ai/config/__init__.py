"""
Configuration management for scrn_ai pipelines.

This module provides YAML configuration parsing, validation, and 
command-line override capabilities for automated workflow execution.
"""

from .parser import ConfigParser

__all__ = ['ConfigParser']
