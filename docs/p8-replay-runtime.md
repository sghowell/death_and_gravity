# Replaying long P8 certificate chains

Use the read-only entry point from the repository's Python environment:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/path/to/checkpoint
    .venv/bin/python scripts/p8_replay.py cli p8_package_name.verify
    .venv/bin/python scripts/p8_replay.py full

The first two commands retain unmodified scientific SymPy and check its
GCD implementation identity before and after execution. Only full mode
delegates to the existing independently audited exact-GCD snapshot runner.
Its collection, hash and adapter-restoration checks are unchanged.
Full mode also accepts --collect-only and --trace-collection.

The wrapper sets the recursion allowance to 4000 for pytest imports of the
long certificate chain and disables the integer-to-string digit cap for
trusted, bounded local exact-arithmetic results. This is not an algorithm
change, approximation, test exclusion or relaxation of a proof gate.
Do not reuse this unlimited-formatting setting for an untrusted public
integer-input service.

These allowances address two observed interpreter failures: S6.131's
ordinary pytest import reached the default recursion limit, and S6.132's
native report contained exact rational integers beyond Python's default
4300-digit formatting cap. The failed runs are not passes. Their frozen
scientific sources and reports are not edited to accommodate the runner;
fresh verification is required.

Ordinary mode snapshots the checkpoint's immediate default-pattern test
files. CLI mode performs the certificate's --check operation, never report
generation or overwrite. Full mode captures the complete P8 test snapshot.
Checkpoint paths outside problems/P8, missing/empty targets, external test
symlinks and arbitrary non-certificate module names are rejected.

Bootstrap unit tests are separate from scientific P8 certification:

    .venv/bin/python -m pytest scripts/tests/test_p8_replay.py -q

Their mocks check delegation and argument/adapter restoration; they do not
substitute for actual scientific or native-certificate replays.
