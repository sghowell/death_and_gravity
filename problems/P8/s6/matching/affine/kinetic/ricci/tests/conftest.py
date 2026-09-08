"""Discover local scientific packages without editing frozen ancestors."""
import sys
from pathlib import Path

P8 = Path(__file__).resolve().parents[6]
for source in sorted(P8.rglob("src")):
    sys.path.insert(0, str(source))
