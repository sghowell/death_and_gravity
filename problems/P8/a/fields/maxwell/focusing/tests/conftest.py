"""Standalone access to the pinned Maxwell and general-clock lineage."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAXWELL_ROOT = ROOT.parent
A_ROOT = MAXWELL_ROOT.parents[1]
CLOCK_ROOT = A_ROOT/"applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei"
sys.path[:0] = [str(ROOT/"src"), str(MAXWELL_ROOT/"src"),
                str(CLOCK_ROOT.parents[2]/"qsei"/"src"),
                *(str(path/"src") for path in (CLOCK_ROOT, *CLOCK_ROOT.parents))]
