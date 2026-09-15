"""Independent stress, covariance, fixed matching and unresolved-light tests."""

import pytest
import sympy as s
from p8_affine import verify as base
from p8_vacuum_affine_gaussian_metric_pole_matching import (
    audit,
    matching,
    poles,
    source,
    spectral,
)

ROWS = audit.residuals()
GATES = audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_each_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_each_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_unsupported_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("mass2", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
@pytest.mark.parametrize("energy_ratio", (5, 9, 17))
def test_literal_scalar_pair_stress_and_full_spin_norms(mass2, energy_ratio):
    energy = energy_ratio * mass2
    e = s.sqrt(energy) / 2
    k = s.sqrt(energy / 4 - mass2)
    eta = s.diag(1, -1, -1, -1)
    a = s.Matrix([e, 0, 0, k])
    b = s.Matrix([e, 0, 0, -k])
    tensor = -(a * b.T + b * a.T - eta * ((a.T * eta * b)[0] + mass2)) / 2
    assert tensor[0, :] == s.zeros(1, 4)
    B = tensor[1:4, 1:4]
    trace = s.simplify(s.trace(B) ** 2 / 3)
    tf = s.simplify((s.trace(B * B) - s.trace(B) ** 2 / 3) / 5)
    phase = s.sqrt(1 - 4 * mass2 / energy) / (32 * s.pi**2)
    assert s.simplify(phase * trace - spectral.density("scalar", 0, energy, mass2)) == 0
    assert s.simplify(phase * tf - spectral.density("scalar", 2, energy, mass2)) == 0


@pytest.mark.parametrize(
    "external,internal", ((1, 1), (1, 3), (2, 5), (s.Rational(3, 2), s.Rational(7, 2)))
)
@pytest.mark.parametrize(
    "cosine", (s.Rational(-2, 3), s.Rational(1, 5), s.Rational(4, 5))
)
def test_independent_literal_gravitational_scalar_pair_tree(external, internal, cosine):
    mu, nu, x = map(s.sympify, (external, internal, cosine))
    energy = 8 * (mu + nu)
    e = s.sqrt(energy) / 2
    pe = s.sqrt(energy / 4 - mu)
    pf = s.sqrt(energy / 4 - nu)
    sine = s.sqrt(1 - x * x)
    a = s.Matrix([e, 0, 0, pe])
    b = s.Matrix([e, 0, 0, -pe])
    c = s.Matrix([-e, -pf * sine, 0, -pf * x])
    d = s.Matrix([-e, pf * sine, 0, pf * x])
    eta = s.diag(1, -1, -1, -1)

    def stress(v, w, mass2):
        return v * w.T + w * v.T - eta * ((v.T * eta * w)[0] + mass2)

    A = stress(a, b, mu)
    B = stress(c, d, nu)
    contraction = s.trace(eta * A * eta * B) - s.trace(eta * A) * s.trace(eta * B) / 2
    assert (
        s.simplify(
            -contraction / (3 * energy)
            - spectral.scalar_tree(energy, x, mu, nu, s.Integer(3))
        )
        == 0
    )


@pytest.mark.parametrize("mu,nu", ((1, 1), (1, 3), (2, 5)))
@pytest.mark.parametrize("z", (s.Rational(-1, 2), s.Rational(1, 3), s.Integer(1)))
def test_independent_full_angular_scalar_sew(mu, nu, z):
    mu, nu = map(s.sympify, (mu, nu))
    ss = 8 * (mu + nu)
    x = s.Symbol("x", real=True)
    # The azimuthal mean of y^2 for axes of relative cosine z.
    mean_y2 = z * z * x * x + (1 - z * z) * (1 - x * x) / 2
    a = (ss * ss + 8 * mu * nu) / (4 * ss)
    b = -(ss - 4 * mu) * (ss - 4 * nu) / (4 * ss)
    average = s.integrate((a + b * x * x) * (a + b * mean_y2), (x, -1, 1)) / 2
    a0 = (ss + 2 * mu) * (ss + 2 * nu) / (6 * ss)
    a2 = -(ss - 4 * mu) * (ss - 4 * nu) / (6 * ss)
    assert s.factor(average - a0 * a0 - a2 * a2 * s.legendre(2, z) / 5) == 0


@pytest.mark.parametrize("species", ("scalar", "vector"))
@pytest.mark.parametrize("spin", (0, 2))
@pytest.mark.parametrize("order", (3, 4, 5, 6))
def test_exact_moments_at_additional_orders(species, spin, order):
    # Different variable r=4nu/sigma turns the moment into a beta integral.
    r = s.Symbol("r", positive=True)
    nu = s.Rational(3, 2)
    sigma = 4 * nu / r
    if species == "scalar":
        polynomial = (
            (sigma - 4 * nu) ** 2 / 3840 if spin == 2 else (sigma + 2 * nu) ** 2 / 384
        )
    else:
        polynomial = (
            (13 * sigma**2 + 56 * nu * sigma + 48 * nu**2) / 3840
            if spin == 2
            else (sigma**2 - 4 * nu * sigma + 12 * nu**2) / 384
        )
    polynomial = s.Poly(
        s.cancel(polynomial * (4 * nu / r**2) / sigma ** (order + 1)), r
    )
    actual = (
        sum(
            c
            * s.gamma(j + 1)
            * s.gamma(s.Rational(3, 2))
            / s.gamma(j + s.Rational(5, 2))
            for (j,), c in polynomial.terms()
        )
        / s.pi**2
    )
    assert (
        s.simplify(actual - spectral.moment(species, spin, order).subs(spectral.NU, nu))
        == 0
    )


@pytest.mark.parametrize(
    "omega,lam,X,Y",
    ((2, 1, 3, 5), (3, 2, 7, -1), (5, s.Rational(3, 2), -2, 4), (7, 3, 2, 9)),
)
def test_independent_covariance_linear_solve(omega, lam, X, Y):
    omega, lam, X, Y = map(s.sympify, (omega, lam, X, Y))
    u, v, c = s.symbols("u v c")
    matrix = s.Matrix([[u, c], [c, v]])
    J = s.Matrix([[0, 1], [-1, 0]])
    G = s.diag(X, Y)
    residual = (
        lam * matrix - omega * (J * matrix - matrix * J) - (J * G + (J * G).T) / 2
    )
    solved = s.solve(list(residual), (u, v, c))
    assert solved[u] == omega * (Y - X) / (lam * lam + 4 * omega * omega)
    assert solved[v] == -solved[u]
    assert solved[c] == lam * (Y - X) / (2 * (lam * lam + 4 * omega * omega))
    detector = s.diag(s.Integer(11), s.Integer(-3))
    actual = -s.trace(detector * matrix.subs(solved)) / 2
    assert actual == omega * (-3 - 11) * (Y - X) / (2 * (lam * lam + 4 * omega * omega))


@pytest.mark.parametrize(
    "matrix",
    (
        s.Matrix([[2, 1, 3, -1], [1, 4, 2, 5], [3, 2, -2, 1], [-1, 5, 1, 3]]),
        s.Matrix([[1, 2, 0, 4], [2, -3, 5, 0], [0, 5, 7, -2], [4, 0, -2, 6]]),
        s.diag(2, -1, 5, 7),
    ),
)
def test_full_nonunimodular_metric_volume_second_variation(matrix):
    eta = s.diag(1, -1, -1, -1)
    e = s.Symbol("epsilon")
    determinant = -(eta + 2 * e * matrix).det()
    actual = s.series(s.sqrt(determinant), e, 0, 3).removeO()
    A = eta * matrix
    expected = 1 + e * s.trace(A) + e * e * (s.trace(A) ** 2 / 2 - s.trace(A * A))
    assert s.expand(actual - expected) == 0
    full = source.parent.fixed_functions()
    volume = (
        source.KAPPA * full["F_full"].subs({source.parent.u: 0, source.parent.X: 0})
        + source.fixed_coefficients()["volume"]
    )
    assert s.cancel(volume) == 0
    assert s.expand(volume * actual) == 0


@pytest.mark.parametrize("sector", ("TT", "trace"))
def test_independent_fixed_local_metric_Euler_variation(sector):
    c = source.fixed_coefficients()
    t = s.Symbol("time")
    h = s.Function("h")(t)
    density = (
        c["R_old"] * s.diff(h, t) ** 2 / 4
        + c["Weyl_squared"] * s.diff(h, t, 2) ** 2 / 2
        if sector == "TT"
        else -c["R_old"] * s.diff(h, t) ** 2 / 2
        + 3 * c["R_old_squared"] * s.diff(h, t, 2) ** 2
    )
    euler = -s.diff(s.diff(density, s.diff(h, t)), t) + s.diff(
        s.diff(density, s.diff(h, t, 2)), t, 2
    )
    wanted = (
        -c["R_old"] * s.diff(h, t, 2) / 2 + c["Weyl_squared"] * s.diff(h, t, 4)
        if sector == "TT"
        else c["R_old"] * s.diff(h, t, 2) + 6 * c["R_old_squared"] * s.diff(h, t, 4)
    )
    assert s.expand(euler - wanted) == 0


def test_light_scalar_curvature_remains_unmatched_in_the_whole_kernel():
    packet = poles.data()
    cR, cW, cR2 = packet["unmatched_light_curvature_coefficients"]
    full = packet["all_three_cuts_with_unmatched_light_polynomial"]
    k = source.KAPPA
    p = poles.P
    assert all(len(value.atoms(s.Integral)) == 3 for value in full.values())
    assert s.diff(full[2], cR) == 2 * p / k
    assert s.diff(full[2], cW) == -4 * p * p / k
    assert s.diff(full[0], cR) == -4 * p / k
    assert s.diff(full[0], cR2) == -24 * p * p / k
    residue = packet["full_parameterized_Gaussian_massless_residue"]
    assert residue.has(cR)
    assert s.cancel(residue * (k + poles.delta_kappa() + 2 * cR) - k) == 0
    assert "NOT established" in packet["exact_covariant_low_momentum_structure"]
    assert "do NOT establish" in packet["rigorous_bounds"]


def test_fixed_H_Proca_Newton_is_not_the_unlicensed_Phi_extension():
    n = source.HEAVY_MASS2
    m = source.VECTOR_MASS2
    fixed = source.fixed_coefficients()
    wanted = (n * (s.log(n) - 1) + 5 * m) / (192 * s.pi**2)
    assert s.factor(fixed["R_old"] - wanted) == 0
    generic_Phi = source.scalar_density(s.Integer(1), s.Integer(0)).coeff(source.R) / (
        64 * s.pi**2
    )
    assert generic_Phi == -1 / (192 * s.pi**2)
    assert s.factor(poles.delta_kappa() - 2 * (wanted + generic_Phi)) == 1 / (
        96 * s.pi**2
    )
    assert (
        "not set by applying the H prescription" in source.data()["whole_source_scope"]
    )


def test_low_disk_margins_and_prescribed_residue_are_separate():
    packet = poles.data()
    assert all(
        value > 0 for value in packet["all_coarse_positive_arithmetic_margins"].values()
    )
    assert all(
        len(value.atoms(s.Integral)) == 2
        for value in packet["complete_radial_A2_and_H0"]
    )
    assert packet["positive_canonical_H_Proca_pole_residue"] == source.KAPPA / (
        source.KAPPA + poles.delta_kappa()
    )
    assert "H/Proca" in audit.observable()["established"]
    assert (
        "Light-scalar finite gravitational curvature matching"
        in audit.observable()["not_established"]
    )


def test_parent_warmup_does_not_change_the_packet_contract():
    before = base.serialize(audit.packets())
    audit.previous.packets()
    assert base.serialize(audit.packets()) == before
    assert source.data()["checks"] is not source.previous.data()["checks"]
    assert (len(ROWS), audit.scalar_entry_count(), len(GATES)) == (121, 368, 36)


def test_original_frontiers_and_rejected_packet_remain_unchanged():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 141 and len(audit.frontier()) == 9
    assert audit.qualifications() == audit.previous.qualifications()
    assert len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert audit.controls()["rejected_inputs"] == 73
    assert audit.require_parameters(audit.parameters()) == audit.parameters()


def test_full_four_point_and_original_P8_obligations_not_closed():
    scope = audit.observable()
    for phrase in (
        "tadpole",
        "LSZ",
        "Regge",
        "V/G/B/P8 closure",
        "full Phi-inclusive pole",
    ):
        assert phrase in scope["not_established"]
    assert "not inferred from the H prescription" in scope["domain"]
    assert "nor the complete scalar four-point" in matching.data()["matching_boundary"]
