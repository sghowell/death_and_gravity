"""The certified EC centre data shared by the Theorem C centre computations.

The A3 shooting parameters (a*, mu*) of the certified EC solution lie in
m + (K(X) - m): the midpoint m stored in results/a3_midpoint.json plus the
stored Krawczyk offset K_minus_m (a3driver: K(X) inside int X proves a unique
zero in K(X)).  The midpoint alone (radius <= 1e-41) does NOT contain them:
    a*  - a_m  = 7.031e-12  +/- 8.87e-16,
    mu* - mu_m = 2.3757e-9  +/- 4.61e-14.
These are the balls of modecount.certified_centre and of every certified tube
build (Tube.extend_centre) and of scripts/s43_parts.certified_centre_balls.
Until 2026-10-16, tc_commutator._centre_expansion (Lemma D, the S4-3 centre
region) and tc_monotone.run_centre (the P4-3 centre tail) used the bare
midpoint; both now build their series here.
"""
from __future__ import annotations

import json
import os

from flint import arb

A3_RECORD = os.path.join(os.path.dirname(__file__), "..", "..", "..", "results",
                         "a3_midpoint.json")


def _ball(s):
    """Ball containing the decimal ball "[mid +/- rad]" printed by arb.  The
    float radius is rounded up by arb's mag conversion; the endpoint check
    proves the containment at the working precision (and widens otherwise)."""
    mid, rad = s.strip()[1:-1].split(" +/- ")
    b = arb(mid) + arb(0, float(rad))
    lo, hi = arb(mid) - arb(rad), arb(mid) + arb(rad)
    if not (b.contains(lo) and b.contains(hi)):
        b = b.union(lo).union(hi)
    return b


def certified_centre_balls(path=A3_RECORD):
    """(a, mu) balls containing the certified (a*, mu*), built at the ambient
    working precision: a3_midpoint.json midpoint + its Krawczyk offset."""
    aj = json.load(open(path))
    a = _ball(aj["a"]) + _ball(aj["K_minus_m"][1])
    mu = _ball(aj["mu"]) + _ball(aj["K_minus_m"][2])
    return a, mu


def certified_centre_expansion(prec=256, K=30):
    """The certified A2 centre expansion at the certified (a*, mu*) balls."""
    from . import centre
    from .arbseries import precision
    with precision(prec):
        a, mu = certified_centre_balls()
        ce = centre.centre_expansion(mu * (2 * a).exp(), nhat=(-a).exp(), K=K)
        if not ce.certify().ok:
            raise RuntimeError("centre certificate failed")
    return ce
