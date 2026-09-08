"""Discover the P8 packages without editing a frozen ancestor."""
import sys
from pathlib import Path

for source in sorted(Path(__file__).resolve().parents[7].rglob("src")):
    sys.path.insert(0, str(source))
