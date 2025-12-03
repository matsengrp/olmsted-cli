"""Olmsted CLI - Command-line interface for Olmsted data processing."""

__version__ = "0.1.0"

# Re-export main API class for convenience
from .api import OlmstedData

# Re-export commonly used types
from .types import (
    OlmstedClone,
    OlmstedDataset,
    OlmstedNode,
    OlmstedOutput,
    OlmstedTree,
)

__all__ = [
    "OlmstedData",
    "OlmstedNode",
    "OlmstedTree",
    "OlmstedClone",
    "OlmstedDataset",
    "OlmstedOutput",
]
