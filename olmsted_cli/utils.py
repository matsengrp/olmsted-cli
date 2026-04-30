"""General-purpose utilities with minimal project dependencies.

This module contains pure utility functions and classes that depend only on
the standard library, tqdm, and ``constants.py``.  Any module in the project
can import from here without creating circular dependencies.

Higher-level utilities that depend on other project modules (schemas,
field_metadata, build_config) live in ``process_utils.py``.
"""

import argparse
import gzip
import uuid

from tqdm import tqdm

from .constants import VERBOSITY_HELP


def open_maybe_gzip(path, mode="rt"):
    """Open a file for reading, transparently handling ``.gz``.

    Returns a context-managerable file handle. Path inspection is by
    extension (``.gz``), not magic bytes — keeps the helper trivial and
    matches what ``write_olmsted_json`` does on the write side.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


# Module-level VerbosePrinter instance.  Call set_verbosity() early in each
# command's main() to configure the level.  All modules import ``vprint``
# from here instead of using naked print().


def set_verbosity(level=1):
    """Set the global verbosity level.

    Call this once at the start of each CLI command after parsing args::

        set_verbosity(args.verbose)

    All modules that ``from .utils import vprint`` will then respect the level
    because this mutates the existing object rather than replacing it.
    """
    vprint.level = level


# Constants for infinity handling
inf = float("inf")
neginf = float("-inf")


# Verbosity-aware printing
class VerbosePrinter:
    """
    Handle verbosity-aware printing for Olmsted CLI tools.

    This class provides a clean interface for printing messages at different
    verbosity levels without scattering if-statements throughout the code.

    Verbosity levels:
        0: Errors only (quiet mode)
        1: Normal status messages (default)
        2: Verbose output with detailed information
        3: Debug output with extensive diagnostic information

    Usage:
        vprint = VerbosePrinter(args.verbose)
        vprint.error("Something went wrong!")  # Always shown
        vprint.status("Processing file...")     # Level 1+
        vprint.verbose("Command arguments:")    # Level 2+
        vprint.debug(f"Mutation count: {n}")   # Level 3+

        # Or use the generic print with custom min_level
        vprint.print("Custom message", min_level=2)
    """

    def __init__(self, level=1):
        """
        Initialize the VerbosePrinter.

        Args:
            level: Verbosity level (0=quiet, 1=normal, 2=verbose, 3=debug)
        """
        self.level = level

    def print(self, *args, min_level=1, **kwargs):
        """
        Print if current verbosity level >= min_level.

        Args:
            *args: Arguments to pass to print()
            min_level: Minimum verbosity level required to print
            **kwargs: Keyword arguments to pass to print()
        """
        if self.level >= min_level:
            print(*args, **kwargs)

    def error(self, *args, **kwargs):
        """
        Always print errors (level 0+).

        Args:
            *args: Arguments to pass to print()
            **kwargs: Keyword arguments to pass to print()
        """
        print(*args, **kwargs)

    def status(self, *args, **kwargs):
        """
        Print status messages (level 1+).

        Args:
            *args: Arguments to pass to print()
            **kwargs: Keyword arguments to pass to print()
        """
        self.print(*args, min_level=1, **kwargs)

    def verbose(self, *args, **kwargs):
        """
        Print verbose messages (level 2+).

        Args:
            *args: Arguments to pass to print()
            **kwargs: Keyword arguments to pass to print()
        """
        self.print(*args, min_level=2, **kwargs)

    def debug(self, *args, **kwargs):
        """
        Print debug messages (level 3+).

        Args:
            *args: Arguments to pass to print()
            **kwargs: Keyword arguments to pass to print()
        """
        self.print(*args, min_level=3, **kwargs)

    def progress(self, iterable, min_level=1, **tqdm_kwargs):
        """
        Wrap an iterable with tqdm, respecting verbosity level.

        At levels below min_level, returns the bare iterable (no progress bar).
        At min_level and above, returns a tqdm-wrapped iterable.

        Args:
            iterable: The iterable to wrap.
            min_level: Minimum verbosity level to show progress bar (default: 1).
            **tqdm_kwargs: Keyword arguments passed to tqdm.

        Returns:
            The iterable, optionally wrapped with tqdm.
        """
        if self.level >= min_level:
            return tqdm(iterable, **tqdm_kwargs)
        return iterable


# Initialize default vprint (normal verbosity).
# Commands override this via set_verbosity() after parsing args.
vprint = VerbosePrinter(1)


def add_verbosity_args(parser):
    """Add standard -v/--verbose and -q/--quiet arguments to an argparser.

    Usage:
        parser = argparse.ArgumentParser(...)
        add_verbosity_args(parser)
        args = parser.parse_args()
        vprint = VerbosePrinter(args.verbose)
    """
    parser.add_argument(
        "-v", "--verbose",
        type=int,
        choices=[0, 1, 2, 3],
        default=1,
        help=VERBOSITY_HELP,
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Quiet mode — errors only (equivalent to -v 0)",
    )


def resolve_verbosity(args):
    """Resolve verbosity from args, handling -q flag."""
    if getattr(args, "quiet", False):
        args.verbose = 0
    return args.verbose


# Data extraction utilities
def get_optional_int(row, key, default=0):
    """
    Extract integer from row dictionary, returning default if missing or empty.

    Args:
        row (dict): Dictionary containing data (typically a CSV row)
        key (str): Key to extract from the dictionary
        default (int): Value to return if key is missing or empty (default: 0)

    Returns:
        int: The integer value from row[key], or default if missing/empty
    """
    value = row.get(key)
    return int(value) if value else default


# General utility functions
def comp(f, g):
    """
    Function composition: comp(f, g)(x) == f(g(x))
    """

    def h(*args, **kw_args):
        return f(g(*args, **kw_args))

    return h


def strip_ns(a):
    # Handle namespace stripping for both : and / separators
    return str(a).split(":")[-1].split("/")[-1]


def dict_subset(d, keys):
    return {k: d[k] for k in keys if k in d}


def merge(d, d2):
    """
    Merge d2 into d, returning a new dict (non-mutating).
    """
    d = d.copy()
    d.update(d2)
    return d


def get_in(d, path):
    """
    Retrieve value from nested dictionary using a path list.

    Args:
        d: Dictionary to traverse
        path: List of keys representing path to value

    Returns:
        Value at path or empty dict if path doesn't exist
    """
    return (
        d
        if len(path) == 0
        else get_in(d.get(path[0]) if isinstance(d, dict) else {}, path[1:])
    )


def clean_record(d):
    """
    Clean a record by removing namespaces and handling special values.

    Args:
        d: Data to clean (dict, list, or value)

    Returns:
        Cleaned data
    """
    if isinstance(d, list):
        return list(map(clean_record, d))
    elif isinstance(d, dict):
        return {strip_ns(k): clean_record(v) for k, v in d.items()}
    # can't have infinity in json
    elif d == inf or d == neginf:
        return None
    else:
        return d


def translate_dna_to_aa(dna_sequence):
    """
    Translate DNA sequence to amino acid sequence.
    Uses standard genetic code, handles ambiguous bases.
    """
    if not dna_sequence:
        return ""

    # Standard genetic code
    codon_table = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }

    aa_sequence = ""
    # Process in chunks of 3 nucleotides
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i : i + 3].upper()
        # Handle ambiguous bases by using 'X' for unknown amino acids
        if len(codon) == 3 and codon in codon_table:
            aa_sequence += codon_table[codon]
        else:
            aa_sequence += "X"  # Unknown amino acid for ambiguous codons

    return aa_sequence


# Key renaming utilities

def rename_keys(record, mapping, to_keep=None):
    """
    Rename keys in a record based on a mapping dictionary.

    Args:
        record: Dictionary to modify
        mapping: Dict mapping old keys to new keys
        to_keep: List of keys to keep with original name (copy, don't move)
    """
    if to_keep is None:
        to_keep = []

    for k in mapping.keys():
        if k in record:
            record[mapping[k]] = record.pop(k) if k not in to_keep else record[k]


def remap_list(lst, mapping):
    """Apply key renaming to all elements in a list."""
    for element in lst:
        rename_keys(element, mapping)


def remap_dict_values(d, mapping):
    """Apply key renaming to all values in a dictionary."""
    for v in d.values():
        rename_keys(v, mapping)


def try_del(d, attr):
    """Safely delete an attribute from a dictionary, ignoring errors."""
    try:
        del d[attr]
    except (KeyError, TypeError):
        pass


def listof(xs_str, f=None):
    """Split a colon-separated string and apply optional function to each element."""
    if f is None:
        f = lambda x: x
    return list(map(f, xs_str.split(":")))


def listofint(xs_str):
    """Split a colon-separated string and convert each element to int."""
    return listof(xs_str, int)


# JSON utility functions
def json_rep(x):
    """
    JSON serialization helper for non-standard types.

    Converts UUID objects to strings and other iterables to lists.
    Used as the 'default' parameter for json.dump().
    """
    if isinstance(x, uuid.UUID):
        return str(x)
    else:
        try:
            return list(x)
        except TypeError:
            raise


def natural_number(desc):
    """argparse type for positive integers."""
    def check(value):
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"{desc} must be a positive integer")
        return ivalue
    return check


def is_nullable_string(checker, instance):
    """JSON Schema type checker that treats None as a valid string."""
    return isinstance(instance, (str, type(None)))
