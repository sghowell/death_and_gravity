"""Discover local packages without changing a frozen ancestor."""
import sys
from pathlib import Path

for source in sorted(Path(__file__).resolve().parents[6].rglob("src")):
    sys.path.insert(0, str(source))
