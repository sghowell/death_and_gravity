"""Collection isolation must never silently exclude a captured P8 test."""

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPTS=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SCRIPTS))
from p8_snapshot_regression import CompleteSnapshot, static_snapshot_namespaces


@pytest.fixture
def two_files(tmp_path):
    first=tmp_path/"test_existing.py"
    second=tmp_path/"existing_test.py"
    first.write_text("def test_one():\n    assert 2+2==4\n")
    second.write_text("def test_two():\n    assert 3*3==9\n")
    return tmp_path,(first,second)


def test_capture_both_default_patterns_and_exact_hashes(two_files):
    root,files=two_files
    snapshot=CompleteSnapshot.capture(root)
    assert snapshot.paths==frozenset(files)
    assert set(snapshot.hashes)==set(files)
    snapshot.check_hashes()


def test_later_files_and_directories_are_not_imported(two_files):
    root,files=two_files
    snapshot=CompleteSnapshot.capture(root)
    later=root/"test_later.py"
    later.write_text("raise RuntimeError('Must not import later file')\n")
    newdir=root/"later"
    newdir.mkdir()
    (newdir/"conftest.py").write_text("raise RuntimeError('Must not import later conftest')\n")
    (newdir/"test_later.py").write_text("raise RuntimeError('Must not collect later directory')\n")
    assert snapshot.pytest_ignore_collect(later,None) is True
    assert snapshot.pytest_ignore_collect(newdir,None) is True
    assert snapshot.pytest_ignore_collect(files[0],None) is None
    code=("import sys,pytest\nfrom pathlib import Path\n"
          f"sys.path.insert(0,{str(SCRIPTS)!r})\n"
          "from p8_snapshot_regression import CompleteSnapshot\n"
          f"snapshot=CompleteSnapshot({str(root)!r},{[str(p) for p in files]!r})\n"
          f"raise SystemExit(pytest.main([{str(root)!r},'-q','--import-mode=importlib','-p','no:faulthandler'],plugins=[snapshot]))\n")
    result=subprocess.run([sys.executable,"-u","-c",code],check=False,capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    assert "2 passed" in result.stdout
    assert "all 2 captured files present and unchanged" in result.stdout


@pytest.mark.parametrize("case",("missing","extra"))
def test_missing_or_extra_collected_file_is_an_error(two_files,case):
    root,files=two_files
    snapshot=CompleteSnapshot.capture(root)
    actual=files[:1] if case=="missing" else files+(root/"extra_test.py",)
    session=SimpleNamespace(items=[SimpleNamespace(path=path) for path in actual])
    with pytest.raises(pytest.UsageError,match="Incomplete test snapshot"):
        snapshot.pytest_collection_finish(session)


@pytest.mark.parametrize("when",("collection","finish"))
def test_mutated_snapshot_file_fails_at_both_boundaries(two_files,when):
    root,files=two_files
    snapshot=CompleteSnapshot.capture(root)
    files[0].write_text("def test_changed():\n    assert False\n")
    session=SimpleNamespace(items=[SimpleNamespace(path=path) for path in files])
    with pytest.raises(pytest.UsageError,match="Captured test file changed"):
        if when=="collection":
            snapshot.pytest_collection_finish(session)
        else:
            snapshot.pytest_sessionfinish(session)


def test_empty_and_outside_root_snapshots_are_rejected(two_files):
    root,_files=two_files
    for files in ([],[root.parent/"outside.py"]):
        with pytest.raises(ValueError):
            CompleteSnapshot(root,files)


def test_unrelated_paths_are_deferred_and_snapshot_ancestors_retained(two_files):
    root,files=two_files
    snapshot=CompleteSnapshot.capture(root)
    assert snapshot.pytest_ignore_collect(root,None) is None
    assert snapshot.pytest_ignore_collect(root.parent,None) is None
    assert snapshot.pytest_ignore_collect(root/"ordinary_module.py",None) is None
    assert snapshot.pytest_ignore_collect(files[1],None) is None


def test_optional_collection_trace_changes_no_snapshot_selection(two_files,capsys):
    root,files=two_files
    plain=CompleteSnapshot.capture(root)
    traced=CompleteSnapshot.capture(root,trace_collection=True)
    assert plain.paths==traced.paths
    assert plain.hashes==traced.hashes
    plain.pytest_collectstart(SimpleNamespace(path=files[0]))
    assert capsys.readouterr().out==""
    traced.pytest_collectstart(SimpleNamespace(path=files[0]))
    assert capsys.readouterr().out=="P8 collecting: test_existing.py\n"
    traced.pytest_collectstart(SimpleNamespace(path=root/"later_test.py"))
    assert capsys.readouterr().out==""
    with pytest.raises(TypeError):
        CompleteSnapshot(root,files,trace_collection=1)


@pytest.fixture
def nested_snapshot(tmp_path):
    root=tmp_path/"snapshot_fixture"
    shallow=root/"tests"
    deep=root.joinpath(*[f"level_{index}" for index in range(40)],"tests")
    shallow.mkdir(parents=True)
    deep.mkdir(parents=True)
    for index,directory in enumerate((shallow,deep)):
        (directory/"conftest.py").write_text(
            f"import pytest\n@pytest.fixture\ndef local_value():\n    return {index}\n")
        (directory/"test_same.py").write_text(
            f"def test_same(local_value):\n    assert local_value=={index}\n")
    return tmp_path,CompleteSnapshot.capture(root)


def test_static_namespace_metadata_paths_and_cleanup(nested_snapshot):
    repo,snapshot=nested_snapshot
    expected={".".join(parent.relative_to(repo).parts):parent
              for path in snapshot.paths for parent in path.parents
              if parent!=repo and parent.is_relative_to(repo)}
    hashes=dict(snapshot.hashes)
    with static_snapshot_namespaces(snapshot,repo) as count:
        assert count==len(expected)
        for name,path in expected.items():
            module=sys.modules[name]
            assert type(module.__path__) is list
            assert module.__path__==[str(path)]
            assert module.__spec__.submodule_search_locations==[str(path)]
            assert module.__spec__.origin is None
            assert module.__package__==name
            parent,_,child=name.rpartition(".")
            if parent:
                assert getattr(sys.modules[parent],child) is module
        assert snapshot.hashes==hashes
        snapshot.check_hashes()
    assert not any(name in sys.modules for name in expected)


def test_deep_importlib_collection_preserves_tests_conftests_and_origins(nested_snapshot):
    repo,snapshot=nested_snapshot
    code=("import sys,pytest\nfrom pathlib import Path\n"
          f"sys.path.insert(0,{str(SCRIPTS)!r})\n"
          "from p8_snapshot_regression import CompleteSnapshot,static_snapshot_namespaces\n"
          f"snapshot=CompleteSnapshot.capture({str(snapshot.root)!r})\n"
          f"with static_snapshot_namespaces(snapshot,{str(repo)!r}):\n"
          f"    result=pytest.main([{str(snapshot.root)!r},'--rootdir',{str(repo)!r},'-q','--import-mode=importlib','-p','no:faulthandler'],plugins=[snapshot])\n"
          "    for path in snapshot.paths:\n"
          f"        name='.'.join(path.relative_to(Path({str(repo)!r})).with_suffix('').parts)\n"
          "        assert Path(sys.modules[name].__file__).resolve()==path\n"
          "raise SystemExit(result)\n")
    result=subprocess.run([sys.executable,"-u","-c",code],check=False,
                          capture_output=True,text=True,timeout=30)
    assert result.returncode==0,result.stdout+result.stderr
    assert "2 passed" in result.stdout
    assert "all 2 captured files present and unchanged" in result.stdout


@pytest.mark.parametrize("case",("regular","existing","invalid_name","outside","same_root"))
def test_namespace_guard_rejects_without_partial_mutation(nested_snapshot,case,monkeypatch):
    repo,snapshot=nested_snapshot
    before=set(sys.modules)
    if case=="regular":
        (snapshot.root/"__init__.py").write_text("raise RuntimeError('Must not execute')\n")
    elif case=="existing":
        sentinel=SimpleNamespace(__path__=[str(snapshot.root)])
        monkeypatch.setitem(sys.modules,"snapshot_fixture",sentinel)
        before.add("snapshot_fixture")
    elif case=="invalid_name":
        invalid=snapshot.root/"not-an-identifier"
        invalid.mkdir()
        path=invalid/"test_bad.py"
        path.write_text("def test_bad():\n    pass\n")
        snapshot=CompleteSnapshot.capture(snapshot.root)
    elif case=="outside":
        repo=repo/"unrelated"
    elif case=="same_root":
        repo=snapshot.root
    with pytest.raises(ValueError), static_snapshot_namespaces(snapshot,repo):
        pytest.fail("Guard allowed an ambiguous namespace")
    assert set(sys.modules)==before
    if case=="existing":
        assert sys.modules["snapshot_fixture"] is sentinel


def test_namespace_cleanup_after_exception(nested_snapshot):
    repo,snapshot=nested_snapshot
    with pytest.raises(RuntimeError,match="deliberate"), static_snapshot_namespaces(snapshot,repo):
        raise RuntimeError("deliberate")
    assert not any(name=="snapshot_fixture" or name.startswith("snapshot_fixture.")
                   for name in sys.modules)
