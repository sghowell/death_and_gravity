"""Replay bootstrap tests; no scientific certificate is stubbed as verified."""

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location(
    "p8_replay_under_test", SCRIPTS / "p8_replay.py"
)
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


@pytest.fixture
def checkpoint(tmp_path):
    tests = tmp_path / "problems/P8/checkpoint/tests"
    tests.mkdir(parents=True)
    (tests / "test_example.py").write_text("def test_example():\n    assert True\n")
    (tests / "example_test.py").write_text("def test_second():\n    assert True\n")
    (tests / "helper.py").write_text("VALUE = 1\n")
    return tmp_path, tests


@pytest.fixture
def adapter(monkeypatch):
    original = object()
    fake = SimpleNamespace(
        ORIGINAL_GCD=original, PolyElement=SimpleNamespace(_gcd=original)
    )
    monkeypatch.setattr(replay, "native_adapter", lambda: fake)
    return fake


def test_runtime_allows_large_exact_integers_without_rounding():
    old_recursion = sys.getrecursionlimit()
    old_digits = sys.get_int_max_str_digits()
    try:
        replay.configure_runtime()
        assert sys.getrecursionlimit() == 4000
        assert sys.get_int_max_str_digits() == 0
        exact = 10**5000 + 7
        printed = str(exact)
        assert len(printed) == 5001
        assert printed[0] == "1" and printed[-1] == "7"
        assert int(printed) == exact
    finally:
        sys.set_int_max_str_digits(old_digits)
        sys.setrecursionlimit(old_recursion)


@pytest.mark.parametrize(
    "target", ("problems/P8/checkpoint", "problems/P8/checkpoint/tests")
)
def test_checkpoint_resolves_both_default_patterns_only(checkpoint, target):
    repo, tests = checkpoint
    actual, files = replay.checkpoint_tests(repo, target)
    assert actual == tests
    assert {p.name for p in files} == {"test_example.py", "example_test.py"}


@pytest.mark.parametrize("target", (".", "problems/P4", "../", "problems/P8/missing"))
def test_foreign_or_missing_checkpoint_rejected(checkpoint, target):
    repo, _ = checkpoint
    with pytest.raises(ValueError):
        replay.checkpoint_tests(repo, target)


def test_empty_checkpoint_rejected(checkpoint):
    repo, _ = checkpoint
    (repo / "problems/P8/empty/tests").mkdir(parents=True)
    with pytest.raises(ValueError):
        replay.checkpoint_tests(repo, "problems/P8/empty")


def test_symlinked_external_test_rejected(checkpoint):
    repo, tests = checkpoint
    outside = repo / "outside.py"
    outside.write_text("def test_external():\n    assert True\n")
    (tests / "test_external.py").symlink_to(outside)
    with pytest.raises(ValueError):
        replay.checkpoint_tests(repo, tests)


@pytest.mark.parametrize(
    "module", ("p8_affine.verify", "p8_vacuum_fermion_four_scalar.verify")
)
def test_valid_module_names(module):
    assert replay.certificate_module(module) == module


@pytest.mark.parametrize(
    "module",
    (
        "os",
        "p8_affine",
        "p8_.verify",
        "p8_X.verify",
        "../p8_a.verify",
        "p8_a.verify.extra",
        "p8_a.verify;pwd",
        "",
    ),
)
def test_arbitrary_module_or_shell_text_rejected(module):
    with pytest.raises(ValueError):
        replay.certificate_module(module)


def test_cli_is_read_only_check_and_restores_arguments(monkeypatch, adapter):
    original_argv = sys.argv
    calls = []

    def run(module, *, run_name):
        calls.append((module, run_name, list(sys.argv)))
        assert adapter.PolyElement._gcd is adapter.ORIGINAL_GCD

    monkeypatch.setattr(replay.runpy, "run_module", run)
    assert replay.cli("p8_affine.verify") == 0
    assert calls == [("p8_affine.verify", "__main__", ["p8_affine.verify", "--check"])]
    assert sys.argv is original_argv


def test_cli_rejects_an_algorithm_mutation(monkeypatch, adapter):
    def run(*args, **kwargs):
        adapter.PolyElement._gcd = object()

    monkeypatch.setattr(replay.runpy, "run_module", run)
    with pytest.raises(AssertionError):
        replay.cli("p8_affine.verify")


def test_cli_restores_arguments_on_failure(monkeypatch, adapter):
    original_argv = sys.argv

    def run(*args, **kwargs):
        raise RuntimeError("diagnostic failure")

    monkeypatch.setattr(replay.runpy, "run_module", run)
    with pytest.raises(RuntimeError, match="diagnostic"):
        replay.cli("p8_affine.verify")
    assert sys.argv is original_argv


@pytest.mark.parametrize("exit_code", (0, 1, 5))
def test_ordinary_preserves_snapshot_plugin_and_exit_code(
    checkpoint, monkeypatch, adapter, exit_code
):
    repo, tests = checkpoint
    calls = []

    def run(arguments, *, plugins):
        calls.append(arguments)
        assert len(plugins) == 1
        assert plugins[0].root == tests
        assert {p.name for p in plugins[0].paths} == {
            "test_example.py",
            "example_test.py",
        }
        plugins[0].check_hashes()
        assert adapter.PolyElement._gcd is adapter.ORIGINAL_GCD
        return exit_code

    monkeypatch.setattr(pytest, "main", run)
    assert replay.ordinary(repo, "problems/P8/checkpoint") == exit_code
    assert "--import-mode=importlib" in calls[0]
    assert calls[0][0] == str(tests)
    assert "--exact" not in calls[0]


@pytest.mark.parametrize(
    "collect,trace", ((False, False), (True, False), (False, True), (True, True))
)
def test_full_delegates_to_existing_runner_and_restores_arguments(
    tmp_path, monkeypatch, collect, trace
):
    original_argv = sys.argv
    calls = []

    def run(path, *, run_name):
        calls.append((path, run_name, list(sys.argv)))

    monkeypatch.setattr(replay.runpy, "run_path", run)
    assert replay.full(tmp_path, collect_only=collect, trace_collection=trace) == 0
    path = str(tmp_path / "scripts/p8_snapshot_regression.py")
    expected = (
        [path]
        + (["--collect-only"] if collect else [])
        + (["--trace-collection"] if trace else [])
    )
    assert calls == [(path, "__main__", expected)]
    assert sys.argv is original_argv


@pytest.mark.parametrize(
    "args", ([], ["unknown"], ["ordinary"], ["cli"], ["full", "unexpected"])
)
def test_incomplete_or_extra_cli_arguments_rejected(args):
    with pytest.raises(SystemExit) as error:
        replay.parser().parse_args(args)
    assert error.value.code == 2
