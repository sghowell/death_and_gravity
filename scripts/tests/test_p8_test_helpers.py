"""Import actual sibling helpers; never stub or skip a scientific test."""

import importlib
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from p8_snapshot_regression import CompleteSnapshot
from p8_test_helpers import sibling_test_imports


@pytest.fixture
def snapshot(tmp_path):
    root = tmp_path / "fixture"
    root.mkdir()
    (root / "_p8_fixture_helper.py").write_text("VALUE = 37\n")
    (root / "test_actual.py").write_text(
        "from _p8_fixture_helper import VALUE\n"
        "def test_actual():\n    assert VALUE == 37\n"
    )
    return CompleteSnapshot.capture(root)


def test_importlib_replays_actual_helper_without_test_selection_change(snapshot):
    code = (
        "import sys, pytest\n"
        f"sys.path.insert(0, {str(SCRIPTS)!r})\n"
        "from p8_snapshot_regression import CompleteSnapshot\n"
        "from p8_test_helpers import sibling_test_imports\n"
        f"snapshot = CompleteSnapshot.capture({str(snapshot.root)!r})\n"
        "paths, hashes = snapshot.paths, dict(snapshot.hashes)\n"
        "with sibling_test_imports(snapshot) as count:\n"
        "    assert count == 1\n"
        "    result = pytest.main([str(snapshot.root), '-q', '--import-mode=importlib', '-p', 'no:faulthandler'], plugins=[snapshot])\n"
        "assert snapshot.paths == paths and snapshot.hashes == hashes\n"
        "raise SystemExit(result)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "1 passed" in result.stdout
    assert "all 1 captured files present and unchanged" in result.stdout


def test_paths_restored_and_helper_body_unchanged(snapshot):
    before = list(sys.path)
    path = snapshot.root / "_p8_fixture_helper.py"
    body = path.read_bytes()
    with sibling_test_imports(snapshot) as count:
        assert count == 1
        assert importlib.util.find_spec("_p8_fixture_helper").origin == str(path)
    assert sys.path == before
    assert path.read_bytes() == body


def test_helper_mutation_fails_and_paths_restored(snapshot):
    before = list(sys.path)
    with (
        pytest.raises(RuntimeError, match="helper changed"),
        sibling_test_imports(snapshot),
    ):
        (snapshot.root / "_p8_fixture_helper.py").write_text("VALUE = 99\n")
    assert sys.path == before


@pytest.mark.parametrize("case", ("existing", "duplicate", "symlink", "shadow"))
def test_ambiguous_helper_rejected_before_path_mutation(snapshot, monkeypatch, case):
    root = snapshot.root
    if case == "existing":
        monkeypatch.setitem(sys.modules, "_p8_fixture_helper", SimpleNamespace())
    elif case == "duplicate":
        other = root / "other"
        other.mkdir()
        (other / "_p8_fixture_helper.py").write_text("VALUE = 99\n")
        (other / "test_other.py").write_text("def test_other():\n    pass\n")
        snapshot = CompleteSnapshot.capture(root)
    elif case == "symlink":
        outside = root.parent / "outside.py"
        outside.write_text("VALUE = 99\n")
        (root / "_p8_other_helper.py").symlink_to(outside)
    else:
        other = root.parent / "other"
        other.mkdir()
        (other / "_p8_fixture_helper.py").write_text("VALUE = 99\n")
        monkeypatch.syspath_prepend(str(other))
    before = list(sys.path)
    with pytest.raises(ValueError), sibling_test_imports(snapshot):
        pytest.fail("Ambiguous helper accepted")
    assert sys.path == before


def test_wrong_import_origin_rejected(snapshot, monkeypatch):
    with (
        pytest.raises(RuntimeError, match="Wrong sibling test helper"),
        sibling_test_imports(snapshot),
    ):
        monkeypatch.setitem(
            sys.modules, "_p8_fixture_helper", SimpleNamespace(__file__=__file__)
        )
