"""Shared pytest fixtures and configuration."""

import os
import sys
import datetime
import shutil
from pathlib import Path
import pytest

# No need to add to path since modules are in the package

@pytest.fixture(scope="session")
def test_session_dir(request):
    """Create a single output directory for the entire test session."""
    # Get the root directory
    cli_root = Path(__file__).parent.parent
    test_output_root = cli_root / "_test_output"
    test_output_root.mkdir(exist_ok=True)
    
    # Create a unique directory for this test session
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    pid = os.getpid()
    session_dir = test_output_root / f"run_{pid}_{timestamp}"
    session_dir.mkdir(exist_ok=True)
    
    # Clean up old test runs at the start (keep last 10)
    test_dirs = sorted(test_output_root.glob("run_*"))
    if len(test_dirs) > 10:
        for old_dir in test_dirs[:-10]:
            shutil.rmtree(old_dir, ignore_errors=True)
    
    return session_dir