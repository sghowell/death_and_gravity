"""Applicability checks, not a UV-completion exclusion or a positivity prior."""

from functools import cache

import sympy as sp
from p8.signs import even_rational_positive
from p8_s5 import nonlinear_d as m


@cache
def derive():
    F_vac = sp.factor(m.functions()["F"].subs(m.X, 0))
    positive = even_rational_positive(-F_vac, m.u)
    conformal_length = sp.integrate(1/m.d, (m.u, -sp.oo, sp.oo))
    k = sp.Symbol("k", positive=True)
    ratio = k/sp.sqrt(m.d)
    return {"constant_clock_F": str(F_vac), "minus_F_strictly_positive": positive,
            "flat_constant_clock_metric_equation_requires": "F(phi0,0)=0",
            "finite_constant_clock_Minkowski_vacuum": "absent for this explicit polynomial continuation to X=0",
            "conformal_time_length_over_tau": str(conformal_length),
            "fixed_comoving_k_over_Eref": str(ratio),
            "past_ratio_limit": str(sp.limit(ratio, m.u, -sp.oo)),
            "future_ratio_limit": str(sp.limit(ratio, m.u, sp.oo)),
            "tensor_Planck_mass": "M^2, nonzero; graviton is not removed",
            "positivity_verdict": "NOT_APPLIED: stationary scattering/vacuum matching and gravitational dispersion assumptions not frozen",
            "warning": "Neither absence of this vacuum nor finite conformal time excludes all possible UV completions"}
