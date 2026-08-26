"""S2, item A3: the Krawczyk certificate for the Evans-Coleman matching (slow: ~3-5 min).

Runs the whole certified pipeline once at the Newton-refined midpoint (constants below,
produced by ``python -m p4.validated.a3driver``): A1 sonic Taylor model at x0 = -0.05 with
the y-tail certificate -> validated integration of the augmented 7D system to x_c ->
A2 centre Taylor model in mu with the y-tail certificate -> F(m), F'(X) -> Krawczyk.
"""
import json
import os

import numpy as np
import pytest
from flint import arb
from p4 import shoot
from p4.validated import a3driver, matching, tmint
from p4.validated.arbseries import precision

HERE = os.path.dirname(os.path.abspath(__file__))
# Newton-refined midpoint (from a3driver; see notes/s2-validated-shooting.md section 4)
X_C = -3.0
M_V0 = "0.112439401388092"
M_A = None            # filled from the driver output below
M_MU = None
HALF = (1e-9, 1e-8, 1e-8)


def _midpoint():
    if M_A is not None:
        return arb(M_V0), arb(M_A), arb(M_MU)
    p = os.path.join(HERE, "..", "results", "a3_midpoint.json")
    if os.path.exists(p):
        d = json.load(open(p))
        return arb(d["V0"]), arb(d["a"]), arb(d["mu"])
    return arb(a3driver.S1_V0), -arb(a3driver.S1_NHAT).log(), arb(a3driver.S1_NHAT) ** 2 * arb(a3driver.S1_WHAT)


def test_krawczyk_certificate():           # slow (~3-5 min): the full certified pipeline
    with precision(384):
        c, a_c, mu_c = _midpoint()
        res = a3driver.run(c, a_c, mu_c, X_C, *HALF)
        st, det = res["state"], res["det"]
        # (i) the point enclosure at x_c contains S1's double-precision trajectory
        sh = shoot.shoot(float(c), 1, x_end=X_C, delta=0.05)
        up = st.u_point()
        for i in range(3):
            assert (up[i] + arb(0, 1e-8)).contains(arb(float(sh.y_end[1 + i]))), i
        # (ii) enclosure widths
        assert max(st.pt.widths()) < 1e-20
        assert max(float(y.rad()) / abs(float(y)) for y in st.y_interval()) < 1e-2
        # (iii) the Krawczyk test passes: unique zero of F in the box
        assert res["ok"], (res["Kbox"], det["Fm"])
        for i, hw in enumerate(HALF):
            assert float(abs(res["Kbox"][i]).abs_upper()) < hw
        # (iv) A4: V = v e^x has exactly one zero on (-inf, 0] along the certified solution
        ok, info = matching.sign_certificate_v(st, res["centre"], a_c)
        assert ok, info
        assert -0.3 < info["zero_step"][0] < -0.2                     # S1: the zero is at x = -0.2509
