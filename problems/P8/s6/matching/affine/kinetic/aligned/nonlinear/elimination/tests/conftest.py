"""Discover the existing source roots without changing any frozen package."""
import sys
from pathlib import Path

P8 = next(path for path in Path(__file__).resolve().parents if path.name == "P8")
for source in sorted(P8.rglob("src")):
    sys.path.insert(0, str(source))
