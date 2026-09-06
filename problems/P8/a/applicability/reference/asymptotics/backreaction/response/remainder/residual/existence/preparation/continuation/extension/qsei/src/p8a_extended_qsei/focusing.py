"""Fresh actual-curvature index test; no cosmological incompleteness claim."""

from functools import cache

import sympy as sp
from p8a_extension import bounds as actual
from p8a_see_qsei.scattering import nonnegative

from . import sampling


def index_lower(duration):
    """Conditional lower J*T0 for an actual free-domain segment of this length.

    duration=tau/T0. Passing this necessary envelope check does not prove
    that a segment of a prescribed length fits at a prescribed point.
    """
    duration = nonnegative(duration, positive=True)
    if duration > sampling.calibration()["proper_span_upper"]:
        raise ValueError("duration exceeds the certified source-free envelope")
    return 3/duration-duration/8


@cache
def calibration():
    data, sampler = actual.calibration(), sampling.calibration()
    span = sampler["proper_span_upper"]
    geo = data["geometry"]
    ricci = 3*(geo["u_cap"]+geo["h_max"]**2)/geo["a_min"]**2
    q2 = 180*data["delta"]*sampler["rounded_absolute_coefficient"]
    lower = index_lower(span)
    expansion = 3*geo["h_max"]/geo["a_min"]
    margins = {"Ricci": sp.Rational(1, 4)-ricci,
               "positive_index_energy_factor": 3-span**2/8,
               "no_comoving_index_trigger": lower-expansion}
    if any(value <= 0 for value in margins.values()):
        raise ValueError("the extended comoving index screen failed")
    return {"proper_span_upper": span, "Ricci_times_T0_squared_cap": sp.Rational(1, 4),
            "expansion_times_T0_cap": expansion, "index_times_T0_lower": lower,
            "normal_jacobian_lower": (geo["a_min"]/geo["a_max"])**3,
            "Q2_over_T0_squared": q2,
            "Q2_over_available_duration_squared_lower": q2/span**2,
            "strict_margins": margins}


def identities():
    x = sp.Symbol("x", real=True)
    a = sp.Function("a", positive=True)(x)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    tau, s = sp.symbols("tau s", positive=True)
    ricci = 3*(sp.diff(h/a, x)/a+(h/a)**2)
    return {"actual_FK_Ricci": sp.simplify(ricci+3*(u+h**2)/a**2),
            "one_sided_index_Sobolev_weight": sp.integrate(tau-s, (s, 0, tau))-tau**2/2,
            "index_lower_and_duration": sp.cancel((3-tau**2/8)/tau-(3/tau-tau/8))}
