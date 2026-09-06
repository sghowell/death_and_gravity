import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P8 = ROOT.parents[7]
for source in (ROOT/"src", *(parent/"src" for parent in ROOT.parents[:7]),
               P8/"src", P8/"s5"/"src", P8/"s5"/"physical"/"src"):
    sys.path.insert(0, str(source))
