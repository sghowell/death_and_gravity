"""Read-only first parent-matching obstruction and escape-control replay."""

import argparse
import hashlib
import json
from pathlib import Path

from p8_uv.verify import check_prior as uv_check_prior

from . import background, conformal, convexity, escape, frames, horndeski

ROOT = Path(__file__).resolve().parents[2]
S6 = ROOT.parent
REPORT = ROOT/"certificates"/"parent-tests.json"
PRIOR_SHA = "6a7ab2a7bf6719d3829795779a876dab0694a5b99383726c14f0880515d2d7d4"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    path = S6/"certificates"/"vacuum-extension.json"
    if sha(path) != PRIOR_SHA:
        raise ValueError("Published S6.1 checkpoint changed")
    prior = json.loads(path.read_text())
    for relative, expected in prior["source_sha256"].items():
        if sha(S6/relative) != expected:
            raise ValueError(f"Published S6.1 input changed: {relative}")
    uv_check_prior()
    return PRIOR_SHA


def build_report():
    previous = prior_checks()
    sources = sorted(ROOT.glob("src/p8_matching/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S6.2.parents", "related_claims": ["P8-S6.2.frames"], "date": "2026-09-04",
            "status": "FIRST_MATCHING_ANSATZES_OBSTRUCTED; GENERAL_UV_MATCHING_OPEN",
            "prior_vacuum_extension_sha256": previous,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
            "stable_algebraic_parent": convexity.derive(),
            "minimal_parent": background.derive(),
            "nonminimal_parent_error_bounds": conformal.derive(),
            "nonminimal_escape_control_not_P8_witness": escape.derive(),
            "Einstein_parent_global_tensor_matching": frames.derive(),
            "exact_quartic_Horndeski_map": horndeski.derive(),
            "written_lemmas": ["notes/parent-tests.md", "notes/frame-tests.md"],
            "tested_parent_hypotheses": "FORMULATION.md; selected ansatzes, not universal UV assumptions",
            "not_established": ["exclusion of all parents for the smooth splice or D witness",
                                "a general D-row or DHOST UV no-go", "finite-gravity positivity verdict",
                                "quantum-stress or higher-derivative parent classification",
                                "a controlled alternative matching construction", "M1 interaction control",
                                "all-orders or nonlinear stability", "P8(a) or full P8 completion"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("S6.2 report differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.2: first parent-matching obstructions and nonminimal escape control replay passed; general UV matching OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
