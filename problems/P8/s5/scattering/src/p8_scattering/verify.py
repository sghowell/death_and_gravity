"""Read-only S5.3b hard finite-time tree bound and S5.4 applicability replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_control.verify import check_inputs
from p8_physical.vertices import tensor_basis

from . import amplitude, majorant, uv_audit
from .frequency import Leg

ROOT = Path(__file__).resolve().parents[2]
CONTROL = ROOT.parent/"control"
REPORT = ROOT/"certificates"/"tree-control.json"
PRIOR_SHA = "e182824e06c326707239eedd5ee562dad39a1b04f4192fc79e0f59977cf276be"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    report = CONTROL/"certificates"/"adiabatic.json"
    if sha(report) != PRIOR_SHA:
        raise ValueError("Oscillator checkpoint changed")
    for relative, expected in json.loads(report.read_text())["source_sha256"].items():
        if sha(CONTROL/relative) != expected:
            raise ValueError(f"Oscillator input changed: {relative}")
    check_inputs()
    return PRIOR_SHA


def fixture(kinds, scale):
    waves = tuple(tuple(scale*k for k in wave) for wave in ((0, 3, 4), (0, 3, -4), (-4, -3, 0), (4, -3, 0)))
    Es, Et = sp.symbols("Es Et", real=True)
    return tuple(Leg(k, kind, sign*(Es if kind == "s" else Et), tensor_basis(k)[0] if kind == "t" else None)
                 for k, kind, sign in zip(waves, kinds, (1, 1, -1, -1)))


@cache
def frozen_example(kinds, scale, point, chart):
    legs = fixture(kinds, scale)
    result = amplitude.tree(legs, point, chart)
    Es, Et = sp.symbols("Es Et", real=True)
    squares = {Es: amplitude.omega_squared("s", legs[0].wave, point, chart),
               Et: amplitude.omega_squared("t", legs[0].wave, point, chart)}
    if any(v <= 0 for v in squares.values()):
        raise ValueError("Frozen external oscillator is not positive")
    def on_shell(value):
        # Reducing even powers first keeps the all-scalar case rational.
        value = sp.cancel(value).subs({Es**2: squares[Es], Et**2: squares[Et]})
        return sp.simplify(value.subs({key: sp.sqrt(v) for key, v in squares.items()}))
    shells = [on_shell(r) for r in result["external_shell_residuals"]]
    if any(r != 0 for r in shells):
        raise ValueError("Frozen external shell check failed")
    channels = []
    for channel in result["channels"]:
        entries = []
        for entry in channel["entries"]:
            inverse = on_shell(entry["inverse_propagator"])
            if inverse == 0:
                raise ValueError("Frozen internal pole")
            entries.append({"kind": entry["kind"], "inverse_propagator": str(inverse),
                            "contribution": str(on_shell(entry["contribution"]))})
        channels.append({"partition": [list(pair) for pair in channel["partition"]], "entries": entries})
    return {"kinds": list(kinds), "scale": scale, "wavevectors": [list(l.wave) for l in legs],
            "chart": chart, "x": str(point), "external_energy_squares": {str(k): str(v) for k, v in squares.items()},
            "external_shell_residuals": list(map(str, shells)),
            "contact": str(on_shell(result["contact"])), "channels": channels,
            "total": str(on_shell(result["total"])),
            "scope": result["scope"]}


@cache
def build_report():
    prior = prior_checks()
    free = amplitude.free_redefinition_check()
    if free["residual"] != 0:
        raise ValueError("Free-field redefinition cancellation failed")
    examples = [frozen_example(kinds, scale, point, chart) for kinds, scale, point, chart in
                (("ssss", 10, 0, "gamma"), ("ssss", 10, sp.Rational(3, 5), "unitary"),
                 ("stst", 10, 0, "gamma"), ("tttt", 10, 0, "gamma"),
                 ("ssss", 10**8, 0, "gamma"))]
    sources = sorted(ROOT.glob("src/p8_scattering/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claims": ["P8-S5.3b.D", "P8-S5.4.audit"], "date": "2026-09-04",
            "status": "FINITE_TIME_HARD_TREE_CRITERION_VERIFIED; UV_ACCEPTANCE_CONTRACT_REQUIRED",
            "prior_oscillator_sha256": prior,
            "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in sources},
            "free_redefinition_cancellation": {"contact": str(free["contact"]),
                                               "exchange": list(map(str, free["exchange"])), "residual": "0"},
            "frozen_tree_regressions_not_used_as_control_proof": examples,
            "uniform_finite_time_tree_majorant": majorant.build(),
            "uv_applicability": uv_audit.derive(),
            "written_lemmas": ["notes/finite-time.md", "notes/frequency.md", "notes/uv-applicability.md"],
            "not_established": ["inclusive/forward or zero-mode amplitudes", "infinite-time S matrix",
                                "vacuum production or full Fock-space norm", "all-orders, loops or UV cutoff",
                                "chart-independent frozen scattering approximation", "global vacuum or accumulated long-time evolution",
                                "M1 interaction control", "nonlinear PDE stability", "UV positivity or completion"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("S5.3b/4 certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.3b: hard finite-time tree bound replay passed; S5.4 UV applicability audited, no UV verdict")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
