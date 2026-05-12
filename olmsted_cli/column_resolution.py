"""Resolve which CSV column supplies each canonical role.

Each role (sample, family, tree) accepts a small family of equivalent
column names — e.g. ``sample``, ``sample_id``, ``sample_name``. When more
than one is present, prefer ``_id`` > bare > ``_name``. Users can override
any role explicitly with ``--<role>-col`` flags.

``sample`` and ``family`` are required for PCP/tree CSV inputs; ``tree``
is optional and falls back to a single-tree-per-family interpretation
when absent.

Multi-variant value disagreement (e.g. both ``sample`` and ``sample_id``
present with different values on the same row) is a value-level concern
handled by :func:`check_row_role_conflicts`.
"""

from typing import Literal, Tuple

#: Closed set of role names recognized by the resolver. Used as a Literal
#: type alias and as a runtime tuple for assertions.
Role = Literal["sample", "family", "tree"]

#: Canonical column-name variants for each role, in preference order.
ROLE_VARIANTS: dict = {
    "sample": ("sample_id", "sample", "sample_name"),
    "family": ("family_id", "family", "family_name"),
    "tree": ("tree_id", "tree", "tree_name"),
}

REQUIRED_ROLES: Tuple[Role, ...] = ("sample", "family")
OPTIONAL_ROLES: Tuple[Role, ...] = ("tree",)


class RoleColumnNotFound(Exception):
    """A required role has no matching column, or an override names a missing one."""


class RoleColumnConflict(Exception):
    """Multiple variants of the same role disagree on a row's values."""


def resolve_role_columns(
    fieldnames,
    *,
    sample_override=None,
    family_override=None,
    tree_override=None,
    required_roles: Tuple[Role, ...] = REQUIRED_ROLES,
):
    """Map each role to its actual CSV column.

    Args:
        fieldnames: iterable of column names from the CSV header.
        sample_override: explicit column name for the sample role, or None.
        family_override: explicit column name for the family role, or None.
        tree_override: explicit column name for the tree role, or None.
        required_roles: roles that must resolve to a real column (raise
            if absent). Defaults to ``("sample", "family")``. Pass
            ``("family",)`` to make sample optional (e.g. trees-CSV
            where sample may be implicit from the companion PCP CSV).

    Returns:
        Dict ``{"sample": str | None, "family": str | None, "tree": str | None}``.
        Roles not in ``required_roles`` resolve to ``None`` when absent.

    Raises:
        RoleColumnNotFound: a required role has no matching column, or
            an override names a column not in ``fieldnames``.
    """
    fieldnames = list(fieldnames)
    fieldset = set(fieldnames)
    overrides = {
        "sample": sample_override,
        "family": family_override,
        "tree": tree_override,
    }
    resolved = {}
    for role, variants in ROLE_VARIANTS.items():
        override = overrides[role]
        if override is not None:
            if override not in fieldset:
                raise RoleColumnNotFound(
                    f"--{role}-col '{override}' not found in CSV columns. "
                    f"Available: {fieldnames}"
                )
            resolved[role] = override
            continue
        chosen = next((v for v in variants if v in fieldset), None)
        if chosen is None and role in required_roles:
            raise RoleColumnNotFound(
                f"No column found for required role '{role}'. "
                f"Expected one of {list(variants)} (or pass --{role}-col)."
            )
        resolved[role] = chosen
    return resolved


def find_present_variants(fieldnames):
    """List which variants per role are present in the CSV header.

    Args:
        fieldnames: iterable of column names.

    Returns:
        Dict ``{role: list[str]}`` in preference order. Roles with zero
        matches map to ``[]``.
    """
    fieldset = set(fieldnames)
    return {
        role: [v for v in variants if v in fieldset]
        for role, variants in ROLE_VARIANTS.items()
    }


def check_row_role_conflicts(row, present_variants):
    """Raise if any role's multiple present variants disagree on this row.

    Args:
        row: dict-like row from ``csv.DictReader``.
        present_variants: dict from :func:`find_present_variants`.

    Raises:
        RoleColumnConflict: two columns for the same role carry different
            non-empty values on this row.
    """
    for role, cols in present_variants.items():
        if len(cols) < 2:
            continue
        anchor = cols[0]
        anchor_val = row.get(anchor, "")
        for other in cols[1:]:
            other_val = row.get(other, "")
            if anchor_val != other_val:
                raise RoleColumnConflict(
                    f"Column conflict for role '{role}': "
                    f"{anchor!r}={anchor_val!r} vs {other!r}={other_val!r} "
                    f"in the same row. Pass --{role}-col to disambiguate."
                )
