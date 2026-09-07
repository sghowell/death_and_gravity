"""Expose local P8 packages without changing any frozen ancestor."""
import sys
from pathlib import Path

P8 = Path(__file__).resolve().parents[4]
for source in sorted(P8.rglob("src")):
    sys.path.insert(0, str(source))
