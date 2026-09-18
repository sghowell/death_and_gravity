"""Independent original-source, auxiliary-limit, Ward, splitting and scope tests."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import (
    audit,
    bounds,
    conventions,
    kernel,
    measure,
    radiative,
    source,
)
from p8_vacuum_affine_radiative_state_soft_index import recoil

ROWS, GATES = audit.residuals(), audit.gates()
Q = s.Matrix([1, 0, 1, 0])
U = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
N = Q[1:4, 0]
PLUS = s.diag(0, 1, 0, -1) / s.sqrt(2)
CROSS = s.zeros(4)
CROSS[1, 3] = CROSS[3, 1] = 1 / s.sqrt(2)
POLS = (PLUS, (PLUS + s.I * CROSS) / s.sqrt(2))


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_analytic_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_rejected_inputs_and_scope(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "j,expected", tuple(enumerate((1, -s.Rational(1, 3), s.Rational(4, 15))))
)
def test_fixed_interval_endpoint_derivatives(j, expected):
    v, t = s.symbols("v t", real=True)
    integrand = 2 / (t + 1 - (t - 1) * v * v)
    assert s.integrate(s.diff(integrand, t, j).subs(t, 1), (v, 0, 1)) == expected


@pytest.mark.parametrize(
    "t",
    (
        s.S.One,
        1 + s.Rational(1, 10**30),
        s.Rational(29, 16),
        s.Rational(17, 8),
        s.Integer(7),
    ),
)
def test_massive_real_kernel_independent_integral_calibration(t):
    with mp.workdps(100):
        mt = mp.mpf(str(s.N(t, 110)))

        def derivative(v):
            D = mt + 1 - (mt - 1) * v * v
            return 16 * mt / D - 4 * (2 * mt * mt - 1) * (1 - v * v) / D**2

        integrated = mp.quad(derivative, [0, 1])
        actual = mp.mpf(str(s.N(kernel.massive_fprime(t), 100)))
        assert abs(integrated - actual) < mp.mpf("1e-70")
        assert 0 < actual < 56


@pytest.mark.parametrize(
    "bad",
    (
        None,
        True,
        False,
        "1",
        1.0,
        s.Float(1),
        [],
        {},
        s.Matrix([1]),
        s.I,
        s.oo,
        s.nan,
        s.Symbol("x"),
    ),
)
def test_exact_scalar_apis_reject_inexact_or_unresolved_values(bad):
    for call in (
        kernel.massive_fprime,
        kernel.massive_c,
        kernel.energy,
        bounds.coefficient_change_upper,
        measure.born_seed_upper,
        measure.full_tree_upper,
    ):
        with pytest.raises((TypeError, ValueError)):
            call(bad)


@pytest.mark.parametrize("bad", (0, -1, s.Rational(1, 2), 8))
def test_massive_real_invariant_domain(bad):
    with pytest.raises(ValueError):
        kernel.massive_fprime(bad)


@pytest.mark.parametrize("bad", (1, s.Rational(3, 2), 8))
def test_massive_same_orientation_domain(bad):
    with pytest.raises(ValueError):
        kernel.massive_c(bad)


@pytest.mark.parametrize(
    "R", (0, s.Rational(1, 10**100), s.Rational(1, 128), s.Rational(1, 8))
)
def test_original_parameter_continuity_bound(R):
    assert bounds.coefficient_upper() == 320 / s.Integer(10) ** 1200
    assert bounds.coefficient_change_upper(R) == 25000 * R / s.Integer(10) ** 1200


@pytest.mark.parametrize(
    "x", (s.Rational(1, 10**30), s.Rational(1, 128), s.Rational(1, 8))
)
def test_named_measure_integrals_and_original_parameters(x):
    L = -s.log(x)
    assert (
        measure.born_seed_upper(x)
        == 2000 * x * (1 + L) / s.Integer(10) ** 1600
        + 3000 * x * x * (L * L + L + s.Rational(1, 2)) / s.Integer(10) ** 2400
    )
    assert (
        measure.full_tree_upper(x)
        == 4 * s.Rational(1, 10) ** 1187 * x * (1 + L)
        + 3000 * x * x * (L * L + L + s.Rational(1, 2)) / s.Integer(10) ** 2400
    )
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2


def test_log_integrability_does_not_remove_a_leading_hard_loop_pole():
    x, e = s.symbols("x e", positive=True)
    assert s.expand(s.integrate(-s.log(e), (e, 0, x)) - x * (1 - s.log(x))) == 0
    assert (
        s.limit(
            s.integrate(1 / e, (e, s.Symbol("cut", positive=True), x)),
            s.Symbol("cut", positive=True),
            0,
        )
        == s.oo
    )
    assert "hard loop" in audit.observable()["not_established"]


def test_forward_limit_nonzero_phase_not_erased():
    value = conventions.data()["whole_forward_limit_phase"]
    assert value == -3281 * s.sqrt(2) * s.I / (6000 * s.pi)
    assert s.im(value) != 0


@pytest.mark.parametrize("pol", POLS)
def test_arbitrary_unequal_collinear_splitting(pol):
    E, states = radiative.calibration_states()
    base = states["two"]
    alpha = s.Rational(1, 10**40)
    split = [alpha * base[0], (1 - alpha) * base[0], base[1]]
    actual = radiative.coefficient(E, split, U, N, pol)
    expected = radiative.coefficient(E, base, U, N, pol)
    assert s.factor(actual - expected) == 0


@pytest.mark.parametrize(
    "h,j,pol",
    tuple(
        product(
            (s.Rational(1, 100), s.Rational(1, 10**10)),
            (s.Rational(1, 3), s.Rational(3, 4)),
            POLS,
        )
    ),
)
def test_null_Gram_bound_in_joint_soft_direction_approaches(h, j, pol):
    ni = s.Matrix([2 * h / (1 + h * h), (1 - h * h) / (1 + h * h), 0])
    hj = h * j
    nj = s.Matrix([0, (1 - hj * hj) / (1 + hj * hj), 2 * hj / (1 + hj * hj)])
    qi, qj = s.Matrix([1, *ni]), s.Matrix([1, *nj])
    di, dj = recoil.dot(qi, Q), recoil.dot(qj, Q)
    P = lambda a, b: (a.T * pol * b)[0]
    S = P(qi, qi) * dj / di + P(qj, qj) * di / dj - 2 * P(qi, qj)
    gap = 1 - ni.dot(nj)
    assert s.simplify(S * s.conjugate(S)) <= 4 * gap * gap


@pytest.mark.parametrize("pol", POLS)
def test_exact_soft_collinear_input_is_not_arbitrarily_assigned(pol):
    with pytest.raises(ValueError, match="soft-collinear"):
        radiative.coefficient(s.Rational(5, 4), [Q / 32], U, N, pol)


@pytest.mark.parametrize(
    "mutation",
    (
        "energy_none",
        "wrong_u",
        "wrong_n",
        "missing_quanta",
        "zero_ray",
        "bad_ray",
        "bad_tensor",
        "traceful",
        "non_TT",
        "nonunit",
    ),
)
def test_physical_coefficient_api_domain(mutation):
    E, states = radiative.calibration_states()
    args = [E, states["one"], U, N, PLUS]
    if mutation == "energy_none":
        args[0] = None
    elif mutation == "wrong_u":
        args[2] = [1, 0]
    elif mutation == "wrong_n":
        args[3] = [0, 0, 2]
    elif mutation == "missing_quanta":
        args[1] = None
    elif mutation == "zero_ray":
        args[1] = [s.zeros(4, 1)]
    elif mutation == "bad_ray":
        args[1] = [s.Matrix([s.Rational(1, 32), 0, 0, 0])]
    elif mutation == "bad_tensor":
        args[4] = None
    elif mutation == "traceful":
        args[4] = s.diag(0, 1, 0, 1) / s.sqrt(2)
    elif mutation == "non_TT":
        args[4] = s.diag(0, 1, -1, 0) / s.sqrt(2)
    else:
        args[4] = 2 * PLUS
    with pytest.raises((TypeError, ValueError)):
        radiative.coefficient(*args)


def mpv(v):
    val = s.N(v, 115)
    return mp.mpc(str(s.re(val)), str(s.im(val)))


def independent_auxiliary_coefficient(points, rays, A, auxiliary=None):
    allp = tuple(points) + tuple(rays)
    den = [mpv(recoil.dot(p, Q)).real for p in allp]
    pij = [[mpv((p.T * A * r)[0]) for r in allp] for p in allp]
    sb = sum(pij[i][i] / den[i] for i in range(len(allp)))
    real, imag = mp.mpc(0), mp.mpc(0)
    for i in range(len(allp)):
        for j in range(i + 1, len(allp)):
            zij = mpv(recoil.dot(allp[i], allp[j])).real
            Sij = (
                pij[i][i] * den[j] / den[i]
                + pij[j][j] * den[i] / den[j]
                - 2 * pij[i][j]
            )
            if zij == 0:
                assert abs(Sij) < mp.mpf("1e-95")
                continue  # Relative-collinear removable product, not d_j=0.
            magnitude = abs(zij)
            massive_pair = j < 4
            if massive_pair or auxiliary is not None:
                mi = mp.mpf(1) if i < 4 else auxiliary * mpv(allp[i][0]).real
                mj = mp.mpf(1) if j < 4 else auxiliary * mpv(allp[j][0]).real
                mass_product = mi * mj
                if massive_pair and magnitude == 1:
                    fp = mp.mpf(22) / 3
                else:
                    root = mp.sqrt(magnitude * magnitude - mass_product * mass_product)
                    angle = mp.acosh(magnitude / mass_product)
                    ratio = angle / root
                    derivative = (root - magnitude * angle) / root**3
                    fp = (
                        8 * magnitude * ratio
                        + 2
                        * (2 * magnitude * magnitude - mass_product * mass_product)
                        * derivative
                    )
                real -= fp * Sij / (16 * mp.pi**2)
            else:
                scale = mpv(allp[j][0]).real
                if i >= 4:
                    scale *= mpv(allp[i][0]).real
                real -= (mp.log(2 * magnitude / scale) + 1) * Sij / (4 * mp.pi**2)
            if allp[i][0] * allp[j][0] > 0:
                if massive_pair or auxiliary is not None:
                    imag += (
                        zij
                        * (2 * zij * zij - 3 * mass_product * mass_product)
                        / (zij * zij - mass_product * mass_product) ** mp.mpf("1.5")
                        * Sij
                        / (8 * mp.pi)
                    )
                else:
                    imag += 2 * Sij / (8 * mp.pi)
    phase = mp.mpf(0)
    for i, Dval in enumerate(den):
        if i < 4:
            phase += Dval * mp.log(abs(Dval))
        else:
            energy = mpv(allp[i][0]).real
            scale = energy if auxiliary is None else auxiliary * energy
            phase += Dval * mp.log(abs(Dval) / scale)
    real += sb * phase / (4 * mp.pi**2)
    imag -= sb * sum(Dval for p, Dval in zip(allp, den) if p[0] > 0) / (4 * mp.pi)
    return real + 1j * imag


@pytest.mark.parametrize(
    "label,pol,power", tuple(product(("one", "two", "four"), POLS, (8, 16, 30)))
)
def test_auxiliary_mass_limit_including_the_phase(label, pol, power):
    E, states = radiative.calibration_states()
    points, rays, _ = recoil.momenta(E, states[label], U)
    exact = radiative.coefficient(E, rays, U, N, pol)
    with mp.workdps(110):
        eps = mp.mpf(10) ** (-power)
        raw = independent_auxiliary_coefficient(points, rays, pol, auxiliary=eps)
        assert abs(raw - mpv(exact)) < 100 * eps * eps * (1 + abs(mp.log(eps)))


@pytest.mark.parametrize(
    "label,pol", tuple(product(("born", "one", "two", "split", "four"), POLS))
)
def test_original_coefficient_and_total_energy_difference_calibrations(label, pol):
    E, states = radiative.calibration_states()
    rays = states[label]
    actual = radiative.coefficient(E, rays, U, N, pol)
    born = radiative.coefficient(E, [], U, N, pol)
    with mp.workdps(100):
        value = mpv(actual)
        assert abs(value) < 320
        if rays:
            R = mpv(sum(ray[0] for ray in rays)).real
            assert abs(value - mpv(born)) < 25000 * R
        else:
            assert actual == born


def test_omitting_the_phase_leaves_a_nonzero_mass_log_coefficient():
    E, states = radiative.calibration_states()
    row = radiative.components(E, states["one"], U, N, PLUS)
    lost = s.factor(row["Sbar"] * sum(row["D"][4:]) / (4 * s.pi**2))
    assert lost != 0
    assert "delta Z_j" in radiative.data()["whole_normal_derivative_identity"]
