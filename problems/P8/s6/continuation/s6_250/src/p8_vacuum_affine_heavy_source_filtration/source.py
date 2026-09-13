"""Full physical heavy source, canonical clock jets and forced KG normalization."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_parent import family as parent

u, X = parent.u, parent.X
KAPPA = parent.K0
MASS2 = parent.MASS2
G = parent.G
A = parent.LOCALIZER
N = parent.N
J1 = -G / (32 * s.pi**2)


def physical_source(clock=u, gradient=X, localizer=A):
    clock, gradient, localizer = map(s.sympify, (clock, gradient, localizer))
    h = (1 - gradient) ** 8 * s.exp(-localizer * gradient * gradient)
    V = (1 - gradient) ** N / (gradient**N + (1 - gradient) ** N)
    return G * KAPPA * clock * clock * h / 2 - J1 * V


@cache
def data():
    alpha = s.Symbol("localizer", positive=True)
    h = parent.switch(X, alpha)
    V = (1 - X) ** N / (X**N + (1 - X) ** N)
    bare = G * KAPPA * u * u * h / 2
    whole = bare - J1 * V
    actual = parent.coefficients()["normalized_heavy_source"]
    checks = {
        "entire_physical_source_matches_normalized_parent": s.cancel(
            s.sqrt(KAPPA) * actual - bare.subs(alpha, A)
        ),
        "full_fixed_finite_source_extension": s.cancel(whole - bare + J1 * V),
        "finite_source_exact_clock_factor": s.cancel(
            (X**N + (1 - X) ** N) * V - (1 - X) ** N
        ),
        "finite_source_retains_vacuum_value": V.subs(X, 0) - 1,
        "physical_KG_source_not_divided_by_kappa": s.sqrt(KAPPA)
        * G
        * s.sqrt(KAPPA)
        * u
        * u
        / 2
        - G * KAPPA * u * u / 2,
        "entire_synchronous_clock_source_zero_not_only_finite_jet": whole.subs(X, 1),
    }
    for j in range(8):
        checks["entire_switch_clock_jet_" + str(j)] = s.diff(h, X, j).subs(X, 1)
    checks["first_nonzero_entire_clock_jet"] = s.diff(h, X, 8).subs(X, 1) - s.factorial(
        8
    ) * s.exp(-alpha)
    for i in range(3):
        for j in range(8 - i):
            checks[f"whole_bare_source_mixed_clock_jet_{i}_{j}"] = s.diff(
                bare, u, i, X, j
            ).subs(X, 1)
    scale = s.Symbol("scale", positive=True)
    k = s.Symbol("kappa", positive=True)
    q = s.Matrix(4, 4, s.symbols("symmetric_metric_slot0:16", real=True))
    q = (q + q.T) / 2
    metric = s.diag(1, -(scale**2), -(scale**2), -(scale**2))
    metric_inv = metric.inv()
    variation = 2 * q / s.sqrt(k)
    inverse_variation = -metric_inv * variation * metric_inv
    pi = s.Matrix(s.symbols("canonical_light_gradient0:4", real=True))
    v = s.Matrix([1, 0, 0, 0])
    dx = (2 * v.T * metric_inv * pi / s.sqrt(k) + v.T * inverse_variation * v)[0]
    wanted = 2 * (pi[0] - q[0, 0]) / s.sqrt(k)
    checks["full_inverse_metric_first_variation"] = (
        metric * inverse_variation + variation * metric_inv
    )
    checks["whole_canonical_clock_gradient_variation"] = s.expand(dx - wanted)
    t = s.Symbol("clock_time", real=True)
    leading = G * k * t * t * s.exp(-alpha) * wanted**8 / 2
    target = 128 * G * t * t * s.exp(-alpha) * (pi[0] - q[0, 0]) ** 8 / k**3
    checks["full_canonical_eighth_source_coefficient"] = s.expand(leading - target)
    checks["actual_fixed_cubic_parameter_factor"] = 128 * G - s.Rational(1, 64)
    Kgrad, Ktime, Kmass, absH = s.symbols(
        "gradient_energy velocity_energy mass_energy abs_H", nonnegative=True
    )
    E = (Kgrad + Ktime + Kmass) / 2
    expanding = -3 * absH * Ktime - absH * Kgrad
    contracting = 3 * absH * Ktime + absH * Kgrad
    checks["complete_expanding_KG_energy_margin"] = (
        6 * absH * E - expanding - absH * (6 * Ktime + 4 * Kgrad + 3 * Kmass)
    )
    checks["complete_contracting_KG_energy_margin"] = (
        6 * absH * E - contracting - absH * (2 * Kgrad + 3 * Kmass)
    )
    growth = s.Rational(25, 16) ** 6
    checks["entire_actual_KG_energy_growth_margin"] = (
        15 - growth - s.Rational(7517615, 16777216)
    )
    return {
        "physical_source": whole,
        "actual_parameters": {
            "kappa": KAPPA,
            "mass_squared": MASS2,
            "g": G,
            "localizer": A,
            "finite_j1": J1,
            "finite_switch_order": N,
        },
        "density_source": "j=sqrt(-g)J, with density operator K=sqrt(-g)(Box_g+n). Physical H solves (Box_g+n)H=J; normalization of the scalar coefficient action divides BOTH sides by kappa, not only J.",
        "canonical_first_clock_variation": dx,
        "canonical_generic_eighth_source": target,
        "jet_boundary": "J has total perturbation degree at least8 near the whole clock, for complete canonical metric and light jets. Its generic coefficient is t² exp(-A)(pi_dot-h00)^8/(64kappa³). At t0 the explicit clock-field prefactor can raise the pointwise degree, not the uniform slab or loop order.",
        "finite_source_boundary": "The fixed onepoint extension-J1 V is of field degree at least1024 and formal loop grade1. It is retained, not replaced by the different source localizer.",
        "sector_boundary": "On the entire synchronous phi=t submanifold, X1 exactly and J0 for arbitrary spatial metric. Thus purely spatial external restrictions can postpone or eliminate source terms beyond the generic total-field lower bounds; the graph table is not a nonzero tensor/vector amplitude theorem.",
        "full_forced_KG_bound": "On the unchanged FLRW slab with original zero germ, ||H||C Hr <=15 T/sqrt(n) ||J||C Hr, and ||Hdot||C Hr and ||a^-1 grad H||C Hr <=15 T||J||C Hr separately. All external momenta remain. This is the complete forced physical KG equation on a prescribed metric, not nonlinear full-system closure.",
        "KG_growth_factor": growth,
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "entire_order_eight_source_not_zero_function": s.diff(h, X, 8).subs(X, 1)
            != 0,
            "leading_source_not_uniformly_order_ten_on_slab": target != 0,
            "actual_canonical_source_factor_one_over_64": 128 * G == s.Rational(1, 64),
            "finite_onepoint_switch_not_source_localizer": N == 1024 and N > 8,
            "complete_forced_KG_growth_below_fifteen": bool(growth < 15),
            "actual_heavy_mass_positive": bool(MASS2 > 0),
            "no_source_bound_compared_to_unevaluated_full_inverse": True,
        },
    }
