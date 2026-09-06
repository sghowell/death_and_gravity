import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(p/"src") for p in (ROOT, *ROOT.parents)]
sys.path.insert(0, str(ROOT.parents[2]/"qsei"/"src"))
