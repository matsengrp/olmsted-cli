# Enhanced JSON Diff Output for pytest

## Overview

The olmsted-cli test suite now includes enhanced JSON comparison utilities that provide detailed diff output when JSON assertions fail. This makes it much easier to identify exactly what differs between expected and actual JSON output.

## Features

1. **Structural Differences**: Clear reporting of:
   - Missing or extra keys in objects
   - List length mismatches
   - Value differences with exact paths (e.g., `clones[0].size: 10 vs 11`)
   - Type mismatches

2. **Unified Diff Format**: Traditional unified diff view showing:
   - Side-by-side comparison of JSON content
   - Context lines around changes
   - Clear +/- indicators for additions/removals

3. **File and Object Comparison**: Support for both:
   - Direct JSON object comparison
   - JSON file comparison with file paths in output

## Usage

### In Test Classes

```python
class TestExample:
    @pytest.fixture(autouse=True)
    def setup_and_teardown(self, json_assertions):
        self.json_assertions = json_assertions
        
    def test_json_objects(self):
        expected = {"key": "value1"}
        actual = {"key": "value2"}
        
        # This will show detailed differences if they don't match
        self.json_assertions.assert_json_equal(expected, actual)
        
    def test_json_files(self):
        # This will show file paths and differences
        self.json_assertions.assert_json_files_equal(
            "expected.json", 
            "actual.json",
            "Custom error message"
        )
```

### Direct Usage

```python
from conftest import JSONAssertions

# Compare objects
JSONAssertions.assert_json_equal(obj1, obj2, "Objects differ")

# Compare files
JSONAssertions.assert_json_files_equal(file1, file2, "Files differ")
```

### In Directory Comparisons

The `compare_directories()` function now uses the enhanced diff output automatically when comparing JSON files in directories.

## Example Output

When a JSON comparison fails, you'll see output like:

```
JSON objects are not equal:

Structural differences:
  - : Extra keys in second object: {'extra_field'}
  - clones: List length mismatch - 2 vs 3
  - clones[0].size: Value mismatch - 10 vs 11
  - metadata.date: Value mismatch - '2024-01-01' vs '2024-01-02'

Unified diff (first 50 lines):
--- object1
+++ object2
@@ -2,22 +2,27 @@
   "clones": [
     {
       "id": "clone1",
-      "size": 10,
+      "size": 11,
       "v_gene": "IGHV1-2"
     },
     ...
```

## Implementation Details

The implementation is in `tests/conftest.py` and includes:

- `normalize_json()`: Recursively sorts dictionaries for consistent comparison
- `json_diff()`: Generates structural differences between JSON objects
- `format_json_diff()`: Creates formatted diff output for files
- `JSONAssertions`: Helper class providing assertion methods
- `json_assertions`: pytest fixture for easy test integration

## Benefits

1. **Faster Debugging**: No need to manually compare JSON files or print objects
2. **Clear Error Messages**: Exact paths to differences make fixes easier
3. **Complete Context**: Both structural and content differences are shown
4. **File Information**: When comparing files, paths are included in output
5. **Automatic Integration**: Works seamlessly with existing pytest infrastructure