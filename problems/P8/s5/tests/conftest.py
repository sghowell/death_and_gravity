import sys
from pathlib import Path

S5 = Path(__file__).resolve().parents[1]
for source in (S5/"src", S5.parent/"src"):
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
