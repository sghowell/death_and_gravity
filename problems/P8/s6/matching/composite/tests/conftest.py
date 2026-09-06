import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P8 = ROOT.parents[2]
MATCHING = ROOT.parent
MATTER = P8/"s5"/"matter"
for source in (ROOT/"src", MATCHING/"bimetric"/"general"/"monotonic"/"src",
               MATCHING/"bimetric"/"general"/"src", MATCHING/"bimetric"/"src",
               MATCHING/"src", P8/"s6"/"src", P8/"src", P8/"s5"/"src",
               P8/"s5"/"physical"/"src", P8/"s5"/"control"/"src", P8/"s5"/"scattering"/"src",
               MATTER/"src", MATTER/"physical"/"src", MATTER/"physical"/"control"/"src",
               MATTER/"physical"/"control"/"tree"/"src", MATTER/"physical"/"control"/"tree"/"loops"/"src",
               MATTER/"physical"/"control"/"tree"/"loops"/"tensor"/"src",
               MATTER/"physical"/"control"/"tree"/"loops"/"tensor"/"scalars"/"src"):
    sys.path.insert(0, str(source))
