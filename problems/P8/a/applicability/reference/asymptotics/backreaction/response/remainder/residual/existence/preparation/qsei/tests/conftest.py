"""Import only this checkpoint, the immutable ancestors and A.9 sibling."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(path/"src") for path in (ROOT, *ROOT.parents)]
sys.path.insert(0, str(ROOT.parents[2]/"qsei"/"src"))
