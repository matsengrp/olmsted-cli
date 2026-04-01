"""
Version information for olmsted-cli.

Single source of truth for the package version and build hash.
Used by --version flag, generated_by metadata, and pyproject.toml.
"""

import subprocess

#: Package version (keep in sync with pyproject.toml)
__version__ = "0.2.0"


def get_git_hash() -> str:
    """Get the short git commit hash, or 'unknown' if not in a git repo."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return "unknown"


def version_string() -> str:
    """Full version string: 'version (hash)'."""
    return f"{__version__} ({get_git_hash()})"
