"""
Version information for olmsted-cli.

The version comes from the installed package metadata, which setuptools-scm
derives from the latest git tag at build time (see pyproject.toml). The git
tag is the single source of truth — there is no hardcoded version here.
Used by the --version flag and the generated_by output metadata.
"""

import subprocess
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _package_version

try:
    __version__ = _package_version("olmsted-cli")
except PackageNotFoundError:
    # Running from a source tree that was never installed (no dist metadata).
    __version__ = "0+unknown"


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
