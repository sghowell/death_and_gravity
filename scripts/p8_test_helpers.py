"""Expose frozen sibling test helpers without changing importlib collection."""

import hashlib
import importlib.util
import sys
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def sibling_test_imports(snapshot):
    """Temporarily expose private helpers beside captured tests, with guards.

    Pytest importlib mode deliberately does not put test directories on
    sys.path. A few archived tests use absolute sibling-helper imports.
    Preserve those files rather than rewriting their frozen import lines.
    Test selection remains entirely owned by CompleteSnapshot.
    """
    helpers = {}
    for directory in sorted({path.parent for path in snapshot.paths}):
        for candidate in sorted(directory.glob("_*.py")):
            if candidate.name.startswith("__"):
                continue
            path = candidate.resolve()
            name = candidate.stem
            if not name.isidentifier() or path.parent != directory:
                raise ValueError(f"Noncanonical sibling test helper: {candidate}")
            if name in helpers:
                raise ValueError(f"Ambiguous sibling test helper: {name}")
            if name in sys.modules:
                raise ValueError(f"Cannot replace an existing helper module: {name}")
            spec = importlib.util.find_spec(name)
            if spec is not None and (
                spec.origin is None or Path(spec.origin).resolve() != path
            ):
                raise ValueError(f"Sibling test helper would shadow a module: {name}")
            helpers[name] = path

    hashes = {
        path: hashlib.sha256(path.read_bytes()).hexdigest() for path in helpers.values()
    }
    previous_path = list(sys.path)
    try:
        for directory in sorted({str(path.parent) for path in helpers.values()}):
            if directory not in sys.path:
                sys.path.append(directory)
        yield len(helpers)
    finally:
        sys.path[:] = previous_path
        for name, path in helpers.items():
            module = sys.modules.get(name)
            if module is not None:
                origin = getattr(module, "__file__", None)
                if origin is None or Path(origin).resolve() != path:
                    raise RuntimeError(f"Wrong sibling test helper imported: {name}")
        for path, expected in hashes.items():
            if (
                not path.is_file()
                or hashlib.sha256(path.read_bytes()).hexdigest() != expected
            ):
                raise RuntimeError(f"Captured sibling test helper changed: {path}")
