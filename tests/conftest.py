"""Shared pytest fixtures and configuration."""

import os
import sys
import datetime
import shutil
from pathlib import Path
import pytest
import json
import difflib
from typing import Any, Dict, List, Union

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


def normalize_json(obj: Any) -> Any:
    """Recursively sort dictionaries by key to enable comparison."""
    if isinstance(obj, dict):
        return {k: normalize_json(v) for k, v in sorted(obj.items())}
    elif isinstance(obj, list):
        return [normalize_json(item) for item in obj]
    else:
        return obj


def json_diff(obj1: Any, obj2: Any, path: str = "") -> List[str]:
    """Generate detailed differences between two JSON objects."""
    differences = []
    
    if type(obj1) != type(obj2):
        differences.append(f"{path}: Type mismatch - {type(obj1).__name__} vs {type(obj2).__name__}")
        return differences
    
    if isinstance(obj1, dict):
        keys1 = set(obj1.keys())
        keys2 = set(obj2.keys())
        
        # Check for missing/extra keys
        if keys1 - keys2:
            differences.append(f"{path}: Extra keys in first object: {keys1 - keys2}")
        if keys2 - keys1:
            differences.append(f"{path}: Extra keys in second object: {keys2 - keys1}")
        
        # Compare common keys
        for key in keys1 & keys2:
            new_path = f"{path}.{key}" if path else key
            differences.extend(json_diff(obj1[key], obj2[key], new_path))
    
    elif isinstance(obj1, list):
        if len(obj1) != len(obj2):
            differences.append(f"{path}: List length mismatch - {len(obj1)} vs {len(obj2)}")
        
        # Compare elements up to the shorter length
        for i in range(min(len(obj1), len(obj2))):
            new_path = f"{path}[{i}]"
            differences.extend(json_diff(obj1[i], obj2[i], new_path))
    
    elif obj1 != obj2:
        differences.append(f"{path}: Value mismatch - {repr(obj1)} vs {repr(obj2)}")
    
    return differences


def format_json_diff(file1: Union[str, Path], file2: Union[str, Path]) -> str:
    """Create a formatted diff output for two JSON files."""
    with open(file1) as f:
        data1 = json.load(f)
    with open(file2) as f:
        data2 = json.load(f)
    
    # Normalize for comparison
    norm1 = normalize_json(data1)
    norm2 = normalize_json(data2)
    
    # Get structural differences
    differences = json_diff(norm1, norm2)
    
    output = []
    output.append(f"\n{'='*60}")
    output.append(f"JSON Comparison: {Path(file1).name}")
    output.append(f"File 1: {file1}")
    output.append(f"File 2: {file2}")
    output.append(f"{'='*60}\n")
    
    if differences:
        output.append("Structural differences found:")
        for diff in differences[:10]:  # Limit to first 10 differences
            output.append(f"  - {diff}")
        if len(differences) > 10:
            output.append(f"  ... and {len(differences) - 10} more differences")
        output.append("")
    
    # Also show a unified diff of the pretty-printed JSON
    json1_str = json.dumps(norm1, indent=2, sort_keys=True)
    json2_str = json.dumps(norm2, indent=2, sort_keys=True)
    
    diff_lines = list(difflib.unified_diff(
        json1_str.splitlines(keepends=True),
        json2_str.splitlines(keepends=True),
        fromfile=str(file1),
        tofile=str(file2),
        n=3
    ))
    
    if diff_lines:
        output.append("Unified diff (first 100 lines):")
        output.extend(line.rstrip() for line in diff_lines[:100])
        if len(diff_lines) > 100:
            output.append(f"... and {len(diff_lines) - 100} more lines")
    
    return "\n".join(output)


class JSONAssertions:
    """Helper class for JSON assertions with detailed error output."""
    
    @staticmethod
    def assert_json_files_equal(file1: Union[str, Path], file2: Union[str, Path], message: str = ""):
        """Assert two JSON files are equal, showing differences if not."""
        with open(file1) as f:
            data1 = json.load(f)
        with open(file2) as f:
            data2 = json.load(f)
        
        norm1 = normalize_json(data1)
        norm2 = normalize_json(data2)
        
        if norm1 != norm2:
            diff_output = format_json_diff(file1, file2)
            if message:
                pytest.fail(f"{message}\n{diff_output}")
            else:
                pytest.fail(diff_output)
    
    @staticmethod
    def assert_json_equal(obj1: Any, obj2: Any, message: str = ""):
        """Assert two JSON objects are equal, showing differences if not."""
        norm1 = normalize_json(obj1)
        norm2 = normalize_json(obj2)
        
        if norm1 != norm2:
            differences = json_diff(norm1, norm2)
            
            output = ["\nJSON objects are not equal:"]
            if message:
                output.insert(0, message)
            
            output.append("\nStructural differences:")
            for diff in differences[:20]:
                output.append(f"  - {diff}")
            if len(differences) > 20:
                output.append(f"  ... and {len(differences) - 20} more differences")
            
            # Show a sample of the differences in unified diff format
            json1_str = json.dumps(norm1, indent=2, sort_keys=True)
            json2_str = json.dumps(norm2, indent=2, sort_keys=True)
            
            diff_lines = list(difflib.unified_diff(
                json1_str.splitlines(keepends=True),
                json2_str.splitlines(keepends=True),
                fromfile="object1",
                tofile="object2",
                n=3
            ))[:50]
            
            if diff_lines:
                output.append("\nUnified diff (first 50 lines):")
                output.extend(line.rstrip() for line in diff_lines)
            
            pytest.fail("\n".join(output))


@pytest.fixture
def json_assertions():
    """Provide JSON assertion helpers to tests."""
    return JSONAssertions()