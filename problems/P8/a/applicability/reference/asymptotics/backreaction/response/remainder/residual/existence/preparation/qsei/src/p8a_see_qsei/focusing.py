"""Why this short exact-solution/QSEI slab does not close the focusing gate.

The index-form exclusion is only for the comoving normals to constant
cosmic-time surfaces and segments wholly inside the certified free half.
It is not a no-go for other hypersurfaces or for an unknown continuation.
"""

from functools import cache

import sympy as sp
from p8a_preparation import bounds as prior

from . import sampling
from .scattering import nonnegative


def index_lower(duration):
    """Lower J*T0 for g(0)=1,g(tau)=0; duration=tau/T0."""
    duration = nonnegative(duration, positive=True)
    if duration > sampling.calibration()["proper_span_upper"]:
        raise ValueError("the requested duration leaves the certified free-half domain")
    return 3/duration-duration/8


@cache
def calibration():
    data = prior.calibration()
    span = sampling.calibration()["proper_span_upper"]
    ricci_actual_cap = 3*(data["geometry"]["u_cap"]+sp.Rational(1, 4))/4
    lower = index_lower(span)
    expansion = sp.Rational(3, 4)
    q2 = 180*data["delta"]*sampling.calibration()["rounded_coefficient"]
    margins = {"Ricci": sp.Rational(1, 4)-ricci_actual_cap,
               "index_energy_factor": 3-span**2/8,
               "no_index_trigger": lower-expansion}
    if any(value <= 0 for value in margins.values()):
        raise ValueError("the short-slab focusing exclusion failed")
    return {"proper_span_upper": span, "Ricci_times_T0_squared_cap": sp.Rational(1, 4),
            "expansion_times_T0_cap": expansion, "index_times_T0_lower": lower,
            "normal_jacobian_lower": sp.Rational(8, 27),
            "Q2_over_T0_squared": q2, "Q2_over_available_duration_squared_lower": q2/span**2,
            "strict_margins": margins}


def identities():
    x = sp.Symbol("x", real=True)
    a = sp.Function("a", positive=True)(x)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    proper_h = h/a
    tau = sp.Symbol("tau", positive=True)
    s = sp.Symbol("s", nonnegative=True)
    return {
        "FK_timelike_Ricci": sp.simplify(3*(sp.diff(proper_h, x)/a+proper_h**2)+3*(u+h**2)/a**2),
        "one_sided_Sobolev_weight": sp.integrate(tau-s, (s, 0, tau))-tau**2/2,
        "index_lower_algebra": sp.simplify((3-tau**2/8)/tau-(3/tau-tau/8)),
    }
