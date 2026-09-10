"""Fixed-scale vacuum derivatives independently recover the earlier quadratic anchors."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_ms_mass import anchors

from . import massless


@cache
def data():
    d = massless.data()
    old = anchors.data()
    e = d["symbols"]["epsilon"]
    A = s.exp(s.EulerGamma * e) * s.gamma(1 + e)
    L = s.Symbol("log_mass_squared_over_mu_squared", real=True)
    m, mu = s.symbols("m mu", positive=True)
    checks = {}
    polys = {}
    for sector, expected in (
        ("scalar", -3 * L * L + 14 * L - 19),
        ("gauge", 6 * L * L - 16 * L + 18),
    ):
        v = d[sector + "_raw_vacuum_Gamma_rational_factor"]
        c = d[sector + "_proper_fermion_CT_rational_factor"]
        normal = A * A * v * s.exp(-2 * e * L) + A * c * s.exp(-e * L)
        finite = s.simplify(s.diff(normal, e, 2).subs(e, 0) / 2)
        polys[sector] = finite
        mass = (
            A * A * (4 - 4 * e) * (3 - 4 * e) * v + A * (4 - 2 * e) * (3 - 2 * e) * c
        ) / e**2
        key = (
            sector
            + "_paired_mass_reference_in_"
            + ("NY_squared" if sector == "scalar" else "NYaCf")
            + "_m_squared_over_Q_squared_units"
        )
        checks[sector + "_full_regulated_mass_derivative"] = s.simplify(mass - old[key])
        checks[sector + "_general_scale_finite_polynomial"] = s.simplify(
            finite - expected
        )
        literal = m**4 * finite.subs(L, s.log(m * m / (mu * mu)))
        coefficient = s.simplify(s.diff(literal, m, 2).subs(mu, m) / m**2)
        checks[sector + "_finite_fixed_mu_mass_derivative"] = coefficient - (
            -56 if sector == "scalar" else 40
        )
    bare = m**4 * (mu / m) ** (4 * e)
    counter = m**4 * (mu / m) ** (2 * e)
    checks["raw_vacuum_fixed_mu_homogeneity"] = s.simplify(
        m * m * s.diff(bare, m, 2) - (4 - 4 * e) * (3 - 4 * e) * bare
    )
    checks["counterterm_fixed_mu_homogeneity"] = s.simplify(
        m * m * s.diff(counter, m, 2) - (4 - 2 * e) * (3 - 2 * e) * counter
    )
    return {
        "general_scale_finite_massless_vacuum_polynomials": polys,
        "dictionary": "L=log(m^2/mu^2). Differentiate at fixed mu before setting mu=m. Scalar units NYm^4/Q^2 and gauge units NaC_Fm^4/Q^2.",
        "scope": "Checks the massless-exchange reference and its proper forest, not a derivative of the already mu=m finite constant or the full scalar OS-paired vacuum.",
        "checks": checks,
    }
