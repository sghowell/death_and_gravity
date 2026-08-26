"""Driver for A3: Newton refinement of (V0, a, mu) and the Krawczyk certificate.

    python -m p4.validated.a3driver [x_c] [w_V0] [w_a] [w_mu] [newton_iters]

Each iteration: sonic-side state at x0 = -0.05 for V0 in [c - w, c + w] (A1 Taylor
model + y-tail certificate), validated integration to x_c (``tmint``), centre side
(A2 Taylor model in mu + y-tail certificate), F(m), F'(m); a Newton step updates the
midpoint m = (c, a_c, mu_c).  The last iteration runs the Krawczyk test on the box
X = [c +/- w_V0] x [a_c +/- w_a] x [mu_c +/- w_mu].  Prints a JSON summary.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from flint import arb

from . import matching, tmint
from .arbseries import precision

S1_V0 = "0.112439401388092"
S1_NHAT, S1_WHAT = "1.2365999612", "5.82098013"


def run(c, a_c, mu_c, x_c, w, w_a, w_mu, verbose=False, hmax=0.02):
    t0 = time.time()
    st = matching.sonic_initial_state(c, w, x0=-0.05)
    it = tmint.Integrator(K=28, hmax=hmax, verbose=verbose)
    it.integrate(st, x_c)
    cs = matching.CentreSide.build(mu_c, w_mu)
    ok, Kbox, det = matching.krawczyk(st, cs, w_a, (c, a_c, mu_c))
    return dict(state=st, centre=cs, ok=ok, Kbox=Kbox, det=det, time=time.time() - t0)


def main(argv):
    x_c = float(argv[0]) if argv else -3.0
    w = float(argv[1]) if len(argv) > 1 else 1e-9
    w_a = float(argv[2]) if len(argv) > 2 else 1e-8
    w_mu = float(argv[3]) if len(argv) > 3 else 1e-8
    iters = int(argv[4]) if len(argv) > 4 else 2
    with precision(384):
        c = arb(S1_V0)
        a_c = -arb(S1_NHAT).log()
        mu_c = arb(S1_NHAT) ** 2 * arb(S1_WHAT)
        summary = []
        for k in range(iters):
            res = run(c, a_c, mu_c, x_c, w, w_a, w_mu)
            st, det = res["state"], res["det"]
            Fm, newton = det["Fm"], det["newton"]
            entry = dict(iteration=k, V0=c.str(30), a=a_c.str(30), mu=mu_c.str(30),
                         Fm=[f.str(5, radius=True) for f in Fm],
                         newton_step=[float(n) for n in newton], krawczyk_ok=res["ok"],
                         K_minus_m=[k_.str(5, radius=True) for k_ in res["Kbox"]],
                         halfwidths=[w, w_a, w_mu], steps=len(st.log), time=res["time"],
                         pt_width=max(st.pt.widths()), iv_u=max(st.iv.widths()[:4]),
                         iv_y=max(st.iv.widths()[4:]),
                         y_interval=[y.str(8, radius=True) for y in st.y_interval()],
                         u_point=[u.str(20, radius=True) for u in st.u_point()])
            print(json.dumps(entry, indent=1), flush=True)
            summary.append(entry)
            if res["ok"]:
                out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "results",
                                   "a3_midpoint.json")
                with open(os.path.normpath(out), "w") as fh:
                    json.dump(dict(V0=c.str(40), a=a_c.str(40), mu=mu_c.str(40), x_c=x_c,
                                   halfwidths=[w, w_a, w_mu], K_minus_m=entry["K_minus_m"]), fh, indent=1)
                break
            c = c + arb(float(newton[0]))
            a_c = a_c + arb(float(newton[1]))
            mu_c = mu_c + arb(float(newton[2]))
        return summary


if __name__ == "__main__":
    main(sys.argv[1:])
