"""Standalone access to unchanged A17 and its photon/general-clock lineage."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FOCUSING_ROOT = ROOT.parent
MAXWELL_ROOT = FOCUSING_ROOT.parent
A_ROOT = MAXWELL_ROOT.parents[1]
CLOCK_ROOT = A_ROOT/"applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei"
sys.path[:0] = [str(ROOT/"src"), str(FOCUSING_ROOT/"src"), str(MAXWELL_ROOT/"src"),
                str(CLOCK_ROOT.parents[2]/"qsei"/"src"),
                *(str(path/"src") for path in (CLOCK_ROOT, *CLOCK_ROOT.parents))]
