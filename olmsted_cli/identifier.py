"""
Identifier generation utility.

Single source of truth for all CLI-minted identifiers in the Olmsted
output JSON. Enforces the ``{datatype}-{uuid}`` shape at the signature
level so prefixless or format-origin-prefixed idents are impossible to
produce.

Convention (see ARCHITECTURE.md#identifier-conventions):

- ``*_id`` fields are reserved for input-derived identifiers. When
  synthesis is unavoidable, use the same ``{datatype}-{uuid}`` shape
  this minter produces.
- ``ident`` fields are always CLI-minted via ``IdentMinter.mint``.

Depends only on ``constants`` for the closed-set ``IdentDatatype`` type
alias — sits alongside ``constants.py``, ``types.py``, and ``utils.py``
near the bottom of the dependency hierarchy.
"""

from __future__ import annotations

import hashlib
import uuid as _uuid
from typing import Optional

from .constants import IDENT_DATATYPES, IdentDatatype


def deterministic_uuid(seed_base, counter: Optional[int] = None) -> str:
    """Generate a deterministic UUID-shaped string from a seed and counter.

    Uses MD5 of ``f"{seed_base}_{counter}"`` (or just ``seed_base`` when
    no counter is given), formatted as a canonical UUID string. Used to
    make golden data reproducible under ``--seed``.
    """
    seed_str = f"{seed_base}_{counter}" if counter is not None else str(seed_base)
    h = hashlib.md5(seed_str.encode()).hexdigest()
    return f"{h[:8]}-{h[8:12]}-{h[12:16]}-{h[16:20]}-{h[20:32]}"


class IdentMinter:
    """Generate ``{datatype}-{uuid}`` identifiers.

    Deterministic when ``seed`` is provided (for reproducible output
    under ``--seed``); random otherwise.

    The ``datatype`` argument is required on every call — there is no
    way to mint a bare or prefixless uuid through this class. Accepted
    datatypes are closed-set (see ``constants.IdentDatatype``); new
    object types whose ident is CLI-minted must register there.
    """

    def __init__(self, seed: Optional[int] = None):
        self._seed = seed
        self._counter = 0

    def mint(self, datatype: IdentDatatype) -> str:
        """Return a new ``{datatype}-{uuid}`` identifier.

        ``datatype`` is type-checked as ``Literal["dataset", "clone",
        "tree", "sample", "subject"]``; the runtime check catches
        anything slipping past static type-checking (or callers using
        ``Any``).
        """
        if datatype not in IDENT_DATATYPES:
            raise ValueError(
                f"datatype must be one of {IDENT_DATATYPES}; got {datatype!r}"
            )
        self._counter += 1
        if self._seed is not None:
            uuid_str = deterministic_uuid(self._seed, self._counter)
        else:
            uuid_str = str(_uuid.uuid4())
        return f"{datatype}-{uuid_str}"
