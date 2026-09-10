"""Explicit new field content and independently normalized SU(3) algebra."""

from functools import cache

import sympy as sp


@cache
def data():
    I = sp.I
    generators = [
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / 2,
        sp.Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]]) / 2,
        sp.diag(1, -1, 0) / 2,
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / 2,
        sp.Matrix([[0, 0, -I], [0, 0, 0], [I, 0, 0]]) / 2,
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / 2,
        sp.Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]]) / 2,
        sp.diag(1, 1, -2) / (2 * sp.sqrt(3)),
    ]
    y, phi, m = sp.symbols(
        "positive_Yukawa scalar_Phi positive_fermion_mass", positive=True
    )
    flavor_y = sp.diag(y, -y)
    exchange = sp.Matrix([[0, 1], [1, 0]])
    mass = sp.eye(2) * m + flavor_y * phi
    checks = {}
    for i, t in enumerate(generators):
        checks[f"generator_Hermitian_{i}"] = sp.expand(t - t.conjugate().T)
        checks[f"generator_traceless_{i}"] = sp.trace(t)
        for j, u in enumerate(generators):
            checks[f"fundamental_trace_{i}_{j}"] = sp.trace(t * u) - sp.Rational(
                int(i == j), 2
            )
    checks["fundamental_Casimir"] = sum(
        (t * t for t in generators), sp.zeros(3)
    ) - sp.Rational(4, 3) * sp.eye(3)
    checks["mass_flavor_exchange_parity"] = (
        exchange * mass.subs(phi, -phi) * exchange - mass
    )
    checks["opposite_Yukawa_trace"] = sp.trace(flavor_y)
    checks["active_Dirac_Yukawa_square_trace"] = (
        sp.trace(sp.kronecker_product(flavor_y, sp.eye(3)) ** 2) - 6 * y**2
    )
    return {
        "new_candidate": "GY14-unbroken: SU(3), fourteen vectorlike fundamental Dirac fermions, two with opposite Phi Yukawas and twelve inert spectators",
        "unchanged_scalar_tree_part": "canonical Phi,H potential m_Phi^2 Phi^2/2 + M H^2/2 + G H Phi^2/2 + L Phi^4/24; no gauge Higgs or direct H Yukawa",
        "light_scalar_mass_squared": 1,
        "gauge_generators": generators,
        "fundamental_Casimir": sp.Rational(4, 3),
        "active_Dirac_flavors": 2,
        "spectator_Dirac_flavors": 12,
        "all_Dirac_flavors": 14,
        "gauge_dimension": 8,
        "active_color_flavor_multiplicity": 6,
        "one_loop_gauge_b0": sp.Rational(5, 3),
        "below_fermion_threshold_one_loop_b0": 11,
        "flavor_Yukawa_matrix": flavor_y,
        "flavor_mass_matrix": mass,
        "parity_exchange": exchange,
        "perturbative_gauge_boson_mass_squared": [0] * 8,
        "scope": "New prospective boundary data, not a replacement of the frozen scalar action. Vectorlike fermions cancel the perturbative gauge anomaly. Scalar parity is a flavor exchange, not a massive chiral transformation. No full threshold matching, confinement gap, rolling action or common bounce parent is asserted.",
        "checks": checks,
    }
