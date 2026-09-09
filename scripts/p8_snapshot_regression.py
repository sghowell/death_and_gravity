"""Complete, immutable P8 test snapshot with normal lazy directory collection.

The snapshot includes every existing default-pattern Python test file. Only
new files/directories absent from that snapshot are ignored. Missing,
changed or unexpectedly collected test files are errors, not exclusions.
The existing exact GCD adapter and per-test SymPy seed are unchanged.
"""

import argparse
import hashlib
import importlib.machinery
import importlib.util
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

import p8_exact_regression as runner
import pytest
from sympy.core.random import seed

REPO=Path(__file__).resolve().parents[1]


@contextmanager
def static_snapshot_namespaces(snapshot,repo):
    """Pin only empty namespace ancestors to their captured physical paths.

    Pytest's importlib collector otherwise creates dynamic NamespacePath
    objects whose repeated parent recalculation is exponential in this
    archive's depth. No test/source module is imported here; pytest still
    loads and rewrites those modules normally. Refuse regular packages,
    existing modules and ambiguous paths rather than replace their behavior.
    """
    repo=Path(repo).resolve()
    if snapshot.root==repo or not snapshot.root.is_relative_to(repo):
        raise ValueError("Require the snapshot strictly inside its explicit import root")
    ancestors={parent for path in snapshot.paths for parent in path.parents
               if parent!=repo and parent.is_relative_to(repo)}
    planned=[]
    for path in sorted(ancestors,key=lambda value:(len(value.parts),str(value))):
        parts=path.relative_to(repo).parts
        if not all(part.isidentifier() for part in parts):
            raise ValueError(f"Noncanonical namespace path: {path}")
        if (path/"__init__.py").exists():
            raise ValueError(f"Cannot replace a regular package: {path}")
        name=".".join(parts)
        if name in sys.modules:
            raise ValueError(f"Cannot replace an existing module: {name}")
        planned.append((name,path))
    created=[]
    try:
        for name,path in planned:
            spec=importlib.machinery.ModuleSpec(name,loader=None,is_package=True)
            spec.submodule_search_locations=[str(path)]
            module=importlib.util.module_from_spec(spec)
            # Keep ordinary lists, not recursively dynamic _NamespacePath.
            assert type(module.__path__) is list
            sys.modules[name]=module
            parent_name,_,child=name.rpartition(".")
            parent=sys.modules[parent_name] if parent_name else None
            if parent is not None:
                setattr(parent,child,module)
            created.append((name,module,parent,child))
        yield len(created)
    finally:
        for name,module,parent,child in reversed(created):
            if parent is not None and getattr(parent,child,None) is module:
                delattr(parent,child)
            if sys.modules.get(name) is module:
                del sys.modules[name]


class CompleteSnapshot:
    def __init__(self,root,paths,*,trace_collection=False):
        if type(trace_collection) is not bool:
            raise TypeError("Require a native collection-trace flag")
        self.trace_collection=trace_collection
        self.root=Path(root).resolve()
        self.paths=frozenset(Path(path).resolve() for path in paths)
        if not self.paths or any(not path.is_relative_to(self.root) for path in self.paths):
            raise ValueError("Require a nonempty test snapshot inside its explicit root")
        self.hashes={path:hashlib.sha256(path.read_bytes()).hexdigest() for path in self.paths}
        self.ancestors={parent for path in self.paths for parent in path.parents
                        if parent==self.root or parent.is_relative_to(self.root)}

    @classmethod
    def capture(cls,root,*,trace_collection=False):
        root=Path(root).resolve()
        result=subprocess.run(["rg","--files",str(root),"-g","test_*.py","-g","*_test.py"],
                              check=True,capture_output=True,text=True)
        return cls(root,result.stdout.splitlines(),trace_collection=trace_collection)

    def pytest_collectstart(self,collector):
        if self.trace_collection:
            path=collector.path.resolve()
            if path in self.paths:
                print("P8 collecting:",path.relative_to(self.root),flush=True)

    def pytest_ignore_collect(self,collection_path,config):
        path=collection_path.resolve()
        if not path.is_relative_to(self.root):
            return None
        if path.is_dir():
            return True if path not in self.ancestors else None
        if path.match("test_*.py") or path.match("*_test.py"):
            return True if path not in self.paths else None
        return None

    def pytest_collection_finish(self,session):
        actual={item.path.resolve() for item in session.items}
        if actual!=self.paths:
            missing=sorted(str(path) for path in self.paths-actual)
            extra=sorted(str(path) for path in actual-self.paths)
            raise pytest.UsageError(f"Incomplete test snapshot: missing={missing}, extra={extra}")
        self.check_hashes()
        print("P8 collection complete:",len(session.items),"tests; all",len(self.paths),
              "captured files present and unchanged",flush=True)

    def check_hashes(self):
        for path,expected in self.hashes.items():
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
                raise pytest.UsageError(f"Captured test file changed: {path}")

    def pytest_sessionfinish(self,session):
        self.check_hashes()

    def pytest_runtest_setup(self,item):
        seed(0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-collection",action="store_true")
    parser.add_argument("--collect-only",action="store_true")
    args=parser.parse_args()
    seed(0)
    snapshot=CompleteSnapshot.capture(REPO/"problems/P8",trace_collection=args.trace_collection)
    relative=sorted(str(path.relative_to(REPO)) for path in snapshot.paths)
    payload="\n".join(relative)+"\n"
    print("P8 complete test-file snapshot:",len(relative),"files; SHA256",
          hashlib.sha256(payload.encode()).hexdigest(),flush=True)
    print("Exact GCD adapter self-check:",runner.self_check(),flush=True)
    arguments=[str(snapshot.root),"-q","--full-trace","-p","no:faulthandler",
               "--rootdir",str(REPO)]
    if args.collect_only:
        arguments.extend(["--collect-only","-qq"])
    with static_snapshot_namespaces(snapshot,REPO) as namespaces, runner.exact_runner() as counts:
        print("Static snapshot namespace ancestors:",namespaces,flush=True)
        try:
            result=pytest.main(runner.pytest_arguments(arguments),plugins=[snapshot])
        finally:
            print("Exact GCD runner counters:",dict(counts),flush=True)
    return result


if __name__=="__main__":
    raise SystemExit(main())
