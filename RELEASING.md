# Releasing olmsted-cli

olmsted-cli is published to [PyPI](https://pypi.org/project/olmsted-cli/) by a
**tag-triggered GitHub Actions workflow** using **PyPI Trusted Publishing
(OIDC)** — there are no stored API tokens. The package version is derived from
the git tag by `setuptools-scm`, so **the tag is the only place a version is
declared**.

- Workflow: [`.github/workflows/release.yml`](.github/workflows/release.yml)
- Version config: `[tool.setuptools_scm]` in [`pyproject.toml`](pyproject.toml)

## TL;DR — cut a release

```bash
# 0. Make sure main is green and you're on an up-to-date main.
git switch main && git pull

# 1. (Recommended) dry-run to TestPyPI first — see "Smoke-test" below.

# 2. Tag the release commit on main and push the tag.
git tag v0.4.0          # vX.Y.Z, semver, no leading zeros stripped
git push origin v0.4.0

# 3. Watch the run: Actions → "Release". The tag push builds the
#    wheel + sdist and publishes to PyPI automatically.

# 4. Verify.
pipx install olmsted-cli==0.4.0   # or: pip install olmsted-cli==0.4.0
olmsted --version                 # -> olmsted-cli 0.4.0 (<hash>)
```

That's it — no manual `build`/`twine`, no credentials on your machine.

## How versioning works

`setuptools-scm` computes the version from `git describe` at build time:

- **On a tag** `vX.Y.Z` with a clean tree → version is exactly `X.Y.Z`.
- **Between tags** (commits after the last tag) → a development version
  derived from the last tag plus the commit distance (e.g.
  `0.4.0.post3.dev…`). `local_scheme = "no-local-version"` strips the
  `+g<hash>` local segment so the version is always PyPI-acceptable.
- **No tags at all / un-installed source tree** → `version.py` falls back to
  `0+unknown`.

The runtime `__version__` (used by `olmsted --version` and the
`generated_by.version` field in output metadata) reads from the *installed
package metadata*, so it always matches the version the artifact was built at.
Do **not** add a hardcoded version string anywhere — the tag is the source of
truth.

### Choosing the version number

Plain [semver](https://semver.org/) on `vMAJOR.MINOR.PATCH`:

- **PATCH** — bug fixes, no output-schema or CLI-surface changes.
- **MINOR** — backward-compatible features / additive output fields.
- **MAJOR** — breaking changes to the CLI surface or the output schema.

While the package is `Development Status :: 4 - Beta`, output-schema changes may
still land in MINOR releases. Promote to `5 - Production/Stable` + `1.0.0` (edit
the classifier in `pyproject.toml`) once the schema has settled and CI enforces
the suite — at that point breaking output changes require a MAJOR bump.

## Smoke-test to TestPyPI first (recommended)

The workflow has a manual path that publishes to **TestPyPI** instead of PyPI,
so you can validate the artifact before cutting the real tag:

1. Actions → **Release** → **Run workflow** → `target: testpypi`.
2. After it succeeds, install from TestPyPI in a clean environment:
   ```bash
   pipx run --spec \
     --index-url https://test.pypi.org/simple/ \
     --pip-args "--extra-index-url https://pypi.org/simple/" \
     olmsted-cli olmsted --version
   ```
   (the extra index lets runtime deps like `ete3`/`six` resolve from real PyPI).

TestPyPI rejects re-uploads of an existing version, so a dirty-tree dev version
(`…dev…`) is fine for repeated dry-runs; only clean tags produce reusable
release numbers.

## One-time setup (PyPI project owner)

Trusted publishing must be configured once per index before the first publish.
This is done in the **web UI** by someone with owner rights — it cannot be
scripted from here.

1. **Create the GitHub environments.** Repo → Settings → Environments → add
   `pypi` and `testpypi` (names must match `environment:` in `release.yml`).
   Optionally add required reviewers to `pypi` to gate production publishes.

2. **Register the trusted publisher** on each index:
   - PyPI: <https://pypi.org/manage/account/publishing/>
   - TestPyPI: <https://test.pypi.org/manage/account/publishing/>

   Use these values (identical except the environment):

   | Field | Value |
   |-------|-------|
   | PyPI Project Name | `olmsted-cli` |
   | Owner | `matsengrp` |
   | Repository name | `olmsted-cli` |
   | Workflow filename | `release.yml` |
   | Environment name | `pypi` (or `testpypi`) |

   For a brand-new project name, use the **"pending publisher"** form — it
   reserves the name and lets the first workflow run create the project, so no
   manual `twine upload` is ever needed.

## Troubleshooting

- **Version came out as `0.4.0.post…dev…` instead of `0.4.0`.** The build ran
  on an untagged or dirty commit. Tag the exact release commit on a clean tree
  and push the tag.
- **`fetch-depth` / "no tags found" in CI.** The workflow checks out with
  `fetch-depth: 0` so `setuptools-scm` can see tags; don't remove it.
- **Publish step skipped.** `publish-pypi` only runs for `push` events on
  `refs/tags/v*`; a `workflow_dispatch` run goes to `publish-testpypi` instead.
- **`ModuleNotFoundError` on a fresh install.** A runtime dependency isn't
  declared. Reproduce in a *clean* venv (`python -m venv` → `pip install` the
  wheel), not an editable dev env, which masks missing deps. (This is how the
  undeclared `ete3 → six` dependency was found.)
