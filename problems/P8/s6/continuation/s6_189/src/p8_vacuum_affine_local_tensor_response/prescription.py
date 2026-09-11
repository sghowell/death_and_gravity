"""Literal fixed local Proca prescription and independent physical stress signs."""

from functools import cache

import sympy as s
from p8_vacuum_affine_proca_gaussian import stress
from p8_vacuum_affine_proca_gaussian.bridge import KAPPA, MASS
from p8_vector_state.wkb import u


@cache
def data():
    old = stress.prescription()
    full = old["mu_mass_specialization"]
    names = {str(x): x for x in full.free_symbols}
    R, C, E, m = s.symbols("R_old Weyl_squared Euler_density mass", real=True)
    changed = full.subs(
        {
            names["R_old"]: R,
            names["Ricci_squared"]: (C - E) / 2 + R * R / 3,
            names["Riemann_squared"]: 2 * C - E + R * R / 3,
            names["positive_mass"]: m,
        }
    )
    expected = (
        s.Rational(5, 2) * m**4
        + s.Rational(5, 3) * m * m * R
        - C / 30
        + E / 90
        - R * R / 18
    )
    t = s.Symbol("t", real=True)
    scale = s.Function("scale", positive=True)(t)
    lapse = s.Function("lapse", positive=True)(t)
    curvature = 6 * (
        s.diff(scale, t, 2) / (scale * lapse**2)
        + s.diff(scale, t) ** 2 / (scale**2 * lapse**2)
        - s.diff(scale, t) * s.diff(lapse, t) / (scale * lapse**3)
    )
    coefficients = (
        s.Rational(5, 2),
        s.Rational(5, 3) * curvature,
        -(curvature**2) / 18,
    )
    local = stress.local()["actual_local_coefficients"]
    checks = {
        "literal_covariant_finite_basis": s.expand(changed - expected),
        "same_actual_mass": MASS - 1000,
        "same_actual_kappa": KAPPA - s.Integer(10) ** 800,
    }
    signs = {}
    for n, term in enumerate(coefficients):
        L = lapse * scale**3 * term
        EL_N = s.diff(L, lapse) - s.diff(s.diff(L, s.diff(lapse, t)), t)
        EL_a = (
            s.diff(L, scale)
            - s.diff(s.diff(L, s.diff(scale, t)), t)
            + s.diff(s.diff(L, s.diff(scale, t, 2)), t, 2)
        )
        rho = s.simplify((-EL_N / scale**3).subs(lapse, 1).doit())
        pressure = s.simplify((EL_a / (3 * lapse * scale**2)).subs(lapse, 1).doit())
        values = {
            "energy": s.simplify(rho.subs(scale, (1 + t * t) ** 2).doit()),
            "pressure": s.simplify(pressure.subs(scale, (1 + t * t) ** 2).doit()),
        }
        signs[n] = values
        for label, value in values.items():
            checks[f"independent_covariant_sign_{n}_{label}"] = s.factor(
                value - local[n][label].subs(u, t)
            )
    return {
        "literal_fixed_finite_density_before_64_pi_squared": changed,
        "compact_variational_representative": expected - E / 90,
        "same_prescription": "The source-pinned mu=m finite local action is matched to fourth-order mode subtraction. It is not the full state-dependent determinant. Euler and box-R compact variations vanish; no infinite-time flux is discarded.",
        "independent_physical_local_stress_signs": signs,
        "physical_sign": "In +---, delta S=(1/2) integral sqrt(-g) T_ab delta g^ab, so rho=-EL_N/a^3 and pressure=EL_a/(3Na^2); R_old=-R_P8.",
        "fixed_scalar_profile": "The fixed QG1 scalar coefficient and original homogeneous M1 are independent of a compact unimodular TT gamma history with N1,u=t,X1. They are not changed to cancel this metric Hessian.",
        "checks": checks,
    }
