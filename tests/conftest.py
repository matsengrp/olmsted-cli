"""Shared pytest fixtures and configuration."""

import os
import sys
from pathlib import Path

# Add bin to Python path for tests
bin_dir = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(bin_dir))