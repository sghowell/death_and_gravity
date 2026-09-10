"""Read-only P8 replay entry point with explicit interpreter allowances.

Only full mode enables the separately audited exact GCD regression adapter.
Ordinary and CLI modes assert that scientific SymPy remains unmodified.
"""

import argparse
import os
import re
import runpy
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RECURSION_ALLOWANCE = 4000


def configure_runtime():
    # Inputs are bounded, trusted local exact-arithmetic certificates, not
    # untrusted integer strings accepted from a public service.
    sys.setrecursionlimit(RECURSION_ALLOWANCE)
    sys.set_int_max_str_digits(0)
    print(
        "P8 replay runtime: recursion",
        sys.getrecursionlimit(),
        "integer-string digits",
        sys.get_int_max_str_digits(),
        "(0 means unlimited exact formatting)",
        flush=True,
    )


def prepare_imports(repo):
    sys.path.insert(0, str(repo / "scripts"))
    for source in sorted((repo / "problems/P8").rglob("src")):
        sys.path.insert(0, str(source))
    from sympy.core.random import seed

    seed(0)


def checkpoint_tests(repo, target):
    path = Path(target)
    if not path.is_absolute():
        path = repo / path
    path = path.resolve()
    tests = path if path.name == "tests" else path / "tests"
    tests = tests.resolve()
    p8 = (repo / "problems/P8").resolve()
    if not tests.is_relative_to(p8) or tests == p8 or not tests.is_dir():
        raise ValueError(
            "Require an existing checkpoint tests directory inside problems/P8"
        )
    files = sorted(
        {
            p.resolve()
            for pattern in ("test_*.py", "*_test.py")
            for p in tests.glob(pattern)
            if p.is_file()
        }
    )
    if not files or any(p.parent != tests for p in files):
        raise ValueError(
            "Require nonempty, direct, in-directory default-pattern test files"
        )
    return tests, files


def certificate_module(value):
    if not re.fullmatch(r"p8_[a-z0-9_]+\.verify", value):
        raise ValueError("Require a P8 certificate module such as p8_affine.verify")
    return value


def native_adapter():
    import p8_exact_regression

    return p8_exact_regression


def ordinary(repo, target):
    import pytest
    from p8_snapshot_regression import CompleteSnapshot, static_snapshot_namespaces

    tests, files = checkpoint_tests(repo, target)
    snapshot = CompleteSnapshot(tests, files)
    adapter = native_adapter()
    assert adapter.PolyElement._gcd is adapter.ORIGINAL_GCD
    try:
        with static_snapshot_namespaces(snapshot, repo) as count:
            print("Native ordinary static namespace ancestors", count, flush=True)
            return pytest.main(
                [
                    str(tests),
                    "-q",
                    "--rootdir",
                    str(repo),
                    "--import-mode=importlib",
                    "-p",
                    "no:faulthandler",
                ],
                plugins=[snapshot],
            )
    finally:
        assert adapter.PolyElement._gcd is adapter.ORIGINAL_GCD


def cli(module):
    module = certificate_module(module)
    adapter = native_adapter()
    assert adapter.PolyElement._gcd is adapter.ORIGINAL_GCD
    old_argv = sys.argv
    try:
        sys.argv = [module, "--check"]
        runpy.run_module(module, run_name="__main__")
        return 0
    finally:
        sys.argv = old_argv
        assert adapter.PolyElement._gcd is adapter.ORIGINAL_GCD


def full(repo, *, collect_only=False, trace_collection=False):
    script = repo / "scripts/p8_snapshot_regression.py"
    old_argv = sys.argv
    try:
        sys.argv = [str(script)]
        if collect_only:
            sys.argv.append("--collect-only")
        if trace_collection:
            sys.argv.append("--trace-collection")
        # The existing full runner owns adapter checks, snapshot capture,
        # collection completeness, source hashes and adapter restoration.
        runpy.run_path(str(script), run_name="__main__")
        return 0
    finally:
        sys.argv = old_argv


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    modes = result.add_subparsers(dest="mode", required=True)
    normal = modes.add_parser("ordinary", help="Unmodified-SymPy checkpoint tests")
    normal.add_argument("checkpoint")
    native = modes.add_parser("cli", help="Unmodified-SymPy native certificate --check")
    native.add_argument("module")
    complete = modes.add_parser("full", help="Existing complete P8 snapshot regression")
    complete.add_argument("--collect-only", action="store_true")
    complete.add_argument("--trace-collection", action="store_true")
    return result


def main(argv=None):
    arguments = parser()
    args = arguments.parse_args(argv)
    repo = REPO.resolve()
    # Validate targets before configuring/importing any scientific module.
    try:
        if args.mode == "ordinary":
            checkpoint_tests(repo, args.checkpoint)
        elif args.mode == "cli":
            certificate_module(args.module)
    except ValueError as exc:
        arguments.error(str(exc))
    configure_runtime()
    os.chdir(repo)
    prepare_imports(repo)
    if args.mode == "ordinary":
        return ordinary(repo, args.checkpoint)
    if args.mode == "cli":
        return cli(args.module)
    return full(
        repo, collect_only=args.collect_only, trace_collection=args.trace_collection
    )


if __name__ == "__main__":
    raise SystemExit(main())
