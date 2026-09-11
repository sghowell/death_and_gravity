"""Literal full-target clock-jet tree and correct free-matter charge reduction."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family
from p8_vacuum_analytic_affine_parent import dynamics

u = family.u
H = 4 * u / (1 + u * u)
h = (1 + u * u) ** 3
delta = 1 / (2 * h)
ell = 1 / (10 * (1 + u * u) ** 6)


def time(value):
    return s.diff(value, u)


def weighted(value):
    return time(value) + 3 * H * value


def euler(density, field, order=2):
    result = 0
    for j in range(order + 1):
        value = s.diff(density, s.diff(field, u, j))
        for _ in range(j):
            value = weighted(value)
        result += (-1) ** j * value
    return s.factor(result)


@cache
def clock_coefficients():
    X = family.X
    R = 1 + (X - 1) / h
    F = family.original.data()["original_retuned_tree_scalar"]
    actual = family.data()
    checks = {}
    for j in range(4):
        checks["literal_full_R_clock_jet_" + str(j)] = s.simplify(
            s.diff(actual["R"] - R, X, j).subs(X, 1)
        )
    for j in range(3):
        checks["literal_full_F_clock_jet_" + str(j)] = s.simplify(
            s.diff(actual["F"] - F, X, j).subs(X, 1)
        )
    B = -R / 2
    fphi = -s.diff(R, u) / 2
    omega_u = -s.diff(R, u) / (4 * R)
    omega_X = -s.diff(R, X) / (4 * R)
    U = R ** -s.Rational(3, 4)
    a = 2 * U * B / 3
    b = U * (4 * B * s.sqrt(X) * omega_u + 2 * s.sqrt(X) * fphi)
    f = U * (F + 6 * B * X * omega_u**2 + 6 * X * fphi * omega_u)
    speed = s.symbols("speed", positive=True)
    Is = (U * 12 * X * fphi * omega_X).subs(X, speed**2)
    IN = -Is.subs(speed, 1)
    INN = (s.diff(Is, speed) + 2 * Is).subs(speed, 1)

    def jets(expr):
        return [
            s.factor(expr.subs(X, 1)),
            s.factor(-2 * s.diff(expr, X).subs(X, 1)),
            s.factor((6 * s.diff(expr, X) + 4 * s.diff(expr, X, 2)).subs(X, 1)),
        ]

    rows = {name: jets(expr) for name, expr in (("a", a), ("b", b), ("f", f), ("U", U))}
    # Exact I(u,1)=0 boundary. Its lapse jets must be differentiated
    # before discarding the background value.
    rows["b"][1] -= IN
    rows["b"][2] -= INN
    rows["f"][1] -= s.diff(IN, u)
    rows["f"][2] -= s.diff(INN, u) - 2 * s.diff(IN, u)
    checks["primitive_first_lapse_zero"] = s.factor(IN)
    checks["primitive_second_lapse_not_dropped"] = s.factor(
        INN + 3 * s.diff(h, u) / h**3
    )
    return {
        "coefficient_clock_N_jets": rows,
        "primitive_N": IN,
        "primitive_NN": INN,
        "checks": checks,
    }


@cache
def quadratic():
    p = clock_coefficients()
    rows = p["coefficient_clock_N_jets"]
    eps, n, v, vd, sd = s.symbols("epsilon n v vdot matterdot", real=True)
    N = 1 + eps * n
    coeff = {
        name: sum(vals[j] * (eps * n) ** j / s.factorial(j) for j in range(3))
        for name, vals in rows.items()
    }
    vol = 1 + 3 * eps * v + s.Rational(9, 2) * eps**2 * v * v
    # Exactly the full homogeneous transformed light action plus original M1.
    literal = vol * (
        9 * coeff["a"] * (H + eps * vd) ** 2 / N
        + 3 * coeff["b"] * (H + eps * vd)
        + N * coeff["f"]
        + coeff["U"] * (ell + eps * sd) ** 2 / (2 * N)
    )
    raw = s.factor(s.diff(literal, eps, 2).subs(eps, 0) / 2)
    vv = s.factor(s.expand(raw).coeff(v, 1).coeff(vd, 1))
    # Compact weighted IBP of vv*v*v' and 3ell*v*sigma'.
    normalized = s.expand(raw - vv * v * vd - weighted(vv) * v * v / 2)
    B = s.factor(normalized.coeff(n, 1).coeff(vd, 1))
    C = s.factor(normalized.coeff(n, 2))
    w = s.factor(normalized.coeff(n, 1).coeff(sd, 1))
    theta = B / 6
    J = s.factor(C + 3 * theta * theta - w * w / 2)
    expected = (
        -3 * vd * vd
        + (J + w * w / 2 - 3 * theta * theta) * n * n
        + 6 * theta * n * vd
        + sd * sd / 2
        + w * n * sd
        + 3 * ell * v * sd
    )
    checks = {
        "literal_quadratic_and_full_weighted_boundary": s.factor(normalized - expected),
        "actual_lapse_pivot_matches_independent_full_Hamiltonian": s.factor(
            J - dynamics.clock_lapse()["J"]
        ),
        "actual_matter_mixed_normalization": s.factor(w - ell * (3 * delta - 1)),
        "actual_tree_trace_kinetic_normalization": s.factor(
            s.expand(normalized).coeff(vd, 2) + 3
        ),
        "actual_matter_kinetic_normalization": s.factor(
            s.expand(normalized).coeff(sd, 2) - s.Rational(1, 2)
        ),
        "fixed_matter_charge_background_Ward": s.factor(time(ell) + 3 * H * ell),
        "physical_scale_source_mass_cancels_with_boundary": s.factor(
            s.expand(normalized).coeff(v, 2)
        ),
        "background_constraint_no_n_v_contact": s.factor(
            s.expand(normalized).coeff(n, 1).coeff(v, 1)
        ),
    }
    return {
        "actual_raw_quadratic_density": raw,
        "boundary_reduced_quadratic_density": normalized,
        "theta": s.factor(theta),
        "J": J,
        "w": w,
        "ell": ell,
        "delta": delta,
        "checks": checks,
    }


@cache
def adapted():
    q = quadratic()
    eta, scale = (s.Function(name)(u) for name in ("eta", "invariant_scale"))
    n = s.diff(eta, u)
    vhat = scale + H * eta - delta * n
    vdot = s.diff(vhat, u)
    theta, J, w = q["theta"], q["J"], q["w"]
    # Fixed-charge Routh reduction, not on-shell velocity substitution in L.
    effective = (
        -3 * vdot**2
        + 6 * theta * n * vdot
        + (J - 3 * theta**2) * n * n
        - 3 * ell * w * n * vhat
        - s.Rational(9, 2) * ell**2 * vhat**2
    )
    currents = s.ImmutableMatrix([euler(effective, field) for field in (eta, scale)])
    coefficients = {
        str(i) + str(j): {
            order: s.factor(s.diff(currents[i], s.diff(field, u, order)))
            for order in range(5)
        }
        for i in range(2)
        for j, field in enumerate((eta, scale))
    }
    A = -6 * delta**2
    checks = {
        "actual_fourth_time_coefficient": s.factor(coefficients["00"][4] - A),
        "actual_third_time_scale_cross": s.factor(coefficients["01"][3] - 6 * delta),
        "actual_third_scale_time_cross": s.factor(coefficients["10"][3] + 6 * delta),
        "actual_second_scale_coefficient": s.factor(coefficients["11"][2] - 6),
        "strict_time_pivot_continuous_lower": 6 * s.Rational(32, 125) ** 2
        - s.Rational(6144, 15625),
        "strict_time_inverse_reciprocal": s.Rational(6144, 15625)
        * s.Rational(15625, 6144)
        - 1,
    }
    for ij, top in (("01", 3), ("10", 3), ("11", 2)):
        for order in range(top + 1, 5):
            checks[f"derivative_inventory_zero_{ij}_{order}"] = coefficients[ij][order]
    # Re-derive both charge-reduced currents independently in hat variables.
    n0, v0, s0 = (s.Function(name)(u) for name in ("hat_n", "hat_v", "sigma"))
    sd = s.diff(s0, u)
    vd = s.diff(v0, u)
    original = (
        -3 * vd**2
        + (J + w * w / 2 - 3 * theta * theta) * n0 * n0
        + 6 * theta * n0 * vd
        + sd * sd / 2
        + w * n0 * sd
        - 3 * ell * vd * s0
    )
    charge = sd + w * n0 + 3 * ell * v0
    matter_rate = -w * n0 - 3 * ell * v0
    metric = s.Matrix(
        [euler(original, field, 1).subs(sd, matter_rate) for field in (n0, v0)]
    )
    eff0 = (
        -3 * vd**2
        + 6 * theta * n0 * vd
        + (J - 3 * theta**2) * n0 * n0
        - 3 * ell * w * n0 * v0
        - s.Rational(9, 2) * ell**2 * v0 * v0
    )
    checks["literal_free_matter_charge_equation"] = s.factor(
        euler(original, s0, 1) + weighted(charge)
    )
    checks["fixed_charge_metric_Euler_not_wrong_substitution"] = s.ImmutableMatrix(
        (metric - s.Matrix([euler(eff0, field, 1) for field in (n0, v0)])).applyfunc(
            s.factor
        )
    )
    checks["actual_zero_charge_reconstruction"] = charge.subs(sd, matter_rate)
    mapped = metric.subs({n0: n, v0: vhat}, simultaneous=True).doit()
    physical_N = mapped[0] - delta * mapped[1]
    transformed = s.Matrix([-weighted(physical_N) + H * mapped[1], mapped[1]])
    checks["independent_physical_current_to_adapted_tree"] = s.ImmutableMatrix(
        (currents - transformed).applyfunc(s.factor)
    )
    return {
        "adapted_density": effective,
        "adapted_tree_currents": currents,
        "complete_derivative_coefficients": coefficients,
        "A": A,
        "A_absolute_lower": s.Rational(6144, 15625),
        "A_inverse_upper": s.Rational(15625, 6144),
        "source_reconstruction": {"n": n, "physical_v": scale + H * eta, "hat_v": vhat},
        "matter_reconstruction_rate": -3 * ell * vhat - w * n,
        "checks": checks,
    }
