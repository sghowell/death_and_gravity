# P8: complete snapshot collection and regression audit

This is verification infrastructure, not a change to the mathematical
action, a frozen checkpoint, a dependency or a certificate.

## Diagnosed delay

The initial S6.78--S6.81 complete snapshots became CPU-bound during
pytest importlib collection, before running tests. A separate Python
stack diagnostic showed recursive _NamespacePath._recalculate calls
through pytest's parent-package discovery. This archive's deeply nested
empty namespace directories caused repeated dynamic parent-path
recalculation. A separate deep focused science invocation exhibited
the same stack.

The four incomplete smaller full runs and two collection diagnostics
were interrupted and superseded. Their interruption results are not
test passes or mathematical counterexamples. An additional S6.82 native
ordinary invocation was restarted with the explicit namespace setup;
its interruption was in SymPy, not evidence that its entire elapsed
time was collection. That run is likewise not counted as passed.

## Exact snapshot and narrow fix

[scripts/p8_snapshot_regression.py](../scripts/p8_snapshot_regression.py)
captures every existing test_*.py and *_test.py file below the explicit
P8 root, using rg. It records the exact resolved path set and every
test file's SHA-256. Collection must contain every captured file and
no extra file. Missing or modified files fail at collection or at
session finish. Only subsequently added test files/directories outside
the captured set are deferred to a later snapshot.

Before importlib collection, the runner constructs static namespace
packages for the exact captured-test ancestor directories relative to
the repository's explicit pytest root. It accepts only valid namespace
identifiers and directories without __init__.py, and refuses any
already loaded conflicting module. It executes no package code and
does not preload or replace any scientific/test module.

The namespace __path__ and spec search path are ordinary one-element
lists naming the same physical directory. Pytest still imports the
same test/conftest modules, applies assertion rewriting, and assigns
their ordinary importlib-qualified names and file origins. Only the
recursive dynamic parent-path bookkeeping is removed. The context
restores its created namespace modules and parent attributes on exit.
The complete collection has 230 such ancestor directories.

## Verification of the collector

Seventeen helper tests pass in 0.98 seconds. They cover both test-name
patterns, captured-file mutation, missing/extra files, ignoring only
later additions, optional tracing without selection changes, namespace
metadata and cleanup, regular/existing/ambiguous package refusal, and
independent conftests plus duplicate test basenames forty directories
deep. Ruff passes for the runner and tests.

The complete 409-file snapshot collects all 6334 tests successfully.
Its path-list SHA-256 is

    70232b58ca056dedb388415178d3961bb0155c7e6d1667491ae12d382e3ad648

This snapshot was captured before S6.83 tests were added. It includes
every P8 checkpoint through S6.82, without dropping a frozen test.

The complete regression passes 6334 tests in 2282.15 seconds, with
all captured-file content hashes unchanged through session finish.
The adapter counters are 24110 domain fallbacks, 7340 exact descents
and 4 mixed fallbacks. This superseding full snapshot covers all
S6.78--S6.82 checkpoints, not just the newly added tests.

## Scientific arithmetic boundary

The pre-existing opt-in exact Gaussian GCD adapter is unchanged. Its
128 original normalized-tuple comparisons pass; every actual descent
checks both exact cofactor identities, with other cases falling back
to the original implementation. The complete suite retains its
per-test SymPy seed.

Native focused replays may use the same static namespace setup without
entering the GCD-adapter context. They explicitly assert the original
SymPy GCD method before and after pytest. The S6.83 science replay in
this mode passes all 31 tests in 3.30 seconds. This is a collection
optimization, not replacement of the native mathematics.

The full CLI is:

    .venv/bin/python -u scripts/p8_snapshot_regression.py

Use --collect-only for the complete collection diagnostic and
--trace-collection to log captured test files as collection starts.
