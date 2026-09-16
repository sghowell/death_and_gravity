"""Independent components, spectral integrals, Feynman sheets and scope checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_spectator_gravity_insertion import audit, bounds, matching, source
from p8_vacuum_affine_spectator_gravity_insertion import insertion as ins

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "energy,momentum,scale",
    (
        (s.Rational(5, 4), s.Rational(3, 4), 1),
        (s.Rational(5, 3), s.Rational(4, 3), 2),
        (s.Rational(13, 12), s.Rational(5, 12), 3),
    ),
)
@pytest.mark.parametrize("partition", range(3))
def test_independent_conserved_component_stress_projectors(
    energy, momentum, scale, partition
):
    eta = s.diag(1, -1, -1, -1)

    def dot(u, v):
        return (u.T * eta * v)[0]

    E, r = energy, momentum
    ps = tuple(
        scale * s.Matrix(v)
        for v in (
            (E, 0, 0, r),
            (E, 0, 0, -r),
            (-E, -r * s.Rational(3, 5), 0, -r * s.Rational(4, 5)),
            (-E, r * s.Rational(3, 5), 0, r * s.Rational(4, 5)),
        )
    )
    mu = dot(ps[0], ps[0])
    a, b, _ = [dot(ps[x] + ps[y], ps[x] + ps[y]) for x, y in ((0, 1), (0, 2), (0, 3))]
    left, right = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))[partition]
    aa, bb, cc = ins.channels(a, b, mu)[partition]
    x, y = (ps[i] for i in left)
    u, v = (ps[i] for i in right)
    tx = x * y.T + y * x.T - eta * (dot(x, y) + mu)
    tu = u * v.T + v * u.T - eta * (dot(u, v) + mu)
    assert tx * eta * (x + y) == s.zeros(4, 1)
    assert tu * eta * (x + y) == s.zeros(4, 1)
    tr0 = s.trace(eta * tx) * s.trace(eta * tu) / 3
    tr2 = s.trace(tx * eta * tu * eta) - tr0
    assert s.factor(tr0 - ins.numerator0(aa, mu)) == 0
    assert s.factor(tr2 - ins.numerator2(aa, bb, cc, mu)) == 0


@pytest.mark.parametrize("species", ("scalar", "vector"))
@pytest.mark.parametrize("spin", (0, 2))
@pytest.mark.parametrize("p_value", (2, -1, 3 + 2j, -1 - 3j))
def test_independent_radial_spectral_and_closed_integrals(species, spin, p_value):
    with mp.workdps(100):
        p = mp.mpc(p_value)
        fn = s.lambdify(
            (ins.P, ins.NU), ins.radial_closed(species, spin, ins.P, ins.NU), "mpmath"
        )
        weight = s.lambdify(
            source.spectral.V, source.spectral.weight(species, spin), "mpmath"
        )
        rho = s.lambdify(
            (ins.A, ins.NU),
            source.spectral.density(species, spin, ins.A, ins.NU),
            "mpmath",
        )
        radial = mp.quad(lambda v: p * weight(v) / (12 + p * (1 - v * v)), [0, 1])
        spectral = mp.quad(
            lambda sigma: (
                (64 if spin == 2 else 768)
                * mp.pi**2
                * p
                * rho(sigma, 3)
                / (sigma**3 * (sigma + p))
            ),
            [12, 24, mp.inf],
        )
        assert abs(fn(p, 3) - radial) < mp.mpf("1e-80")
        assert abs(radial - spectral) < mp.mpf("1e-80")


@pytest.mark.parametrize("species", ("scalar", "vector"))
@pytest.mark.parametrize("spin", (0, 2))
def test_Feynman_boundary_and_wrong_principal_product_negative_control(species, spin):
    with mp.workdps(100):
        fn = s.lambdify(
            (ins.P, ins.NU), ins.radial_closed(species, spin, ins.P, ins.NU), "mpmath"
        )
        rho = s.lambdify(
            (ins.A, ins.NU),
            source.spectral.density(species, spin, ins.A, ins.NU),
            "mpmath",
        )
        a = mp.mpf(18)
        p = -a - mp.j * mp.mpf("1e-50")
        expected = -(64 if spin == 2 else 768) * mp.pi**3 * rho(a, 3) / a**2
        assert abs(mp.im(fn(p, 3)) - expected) < mp.mpf("1e-45")
        ratio = mp.sqrt(p / (12 + p))
        correct = mp.atanh(ratio) / ((12 + p) * ratio)
        wrong = mp.atanh(ratio) / mp.sqrt(p * (12 + p))
        assert abs(correct + wrong) < mp.mpf("1e-80")
        assert mp.im(correct) > 0 and mp.im(wrong) < 0


@pytest.mark.parametrize("species", ("scalar", "vector"))
@pytest.mark.parametrize("spin", (0, 2))
def test_zero_momentum_subtracted_radial_value(species, spin):
    assert ins.radial_closed(species, spin, 0, 3) == 0


@pytest.mark.parametrize(
    "energy,transfer",
    (
        (s.Rational(25, 4), -s.Rational(9, 8)),
        (9, -1),
        (16, -6),
        (9, -s.Rational(1, 10**204)),
        (9, -s.Rational(1, 10**400)),
        (9, -5 + s.Rational(1, 10**400)),
    ),
)
def test_original_all_nonforward_signed_interval_domain(energy, transfer):
    assert bounds.original_known_rate_interval(energy, transfer) == (
        -s.Rational(3, 10**602),
        0,
    )
    c = 4 - energy - transfer
    ag = ins.gravity_born(energy, transfer, 1, source.KAPPA)
    z = 1 + 2 * transfer / (energy - 4)
    am = source.forward.matter_born(energy, z, source.HEAVY_MASS2, source.CUBIC)
    assert ag > 0 and am > 0
    assert 0 < ag / (ag + am) < 1 and c < 0


@pytest.mark.parametrize(
    "value", (True, False, 1.0, s.Float(1), "1", None, s.I, s.oo, s.nan, s.Symbol("x"))
)
def test_inexact_and_unsupported_domain_values_rejected(value):
    with pytest.raises((ValueError, TypeError)):
        bounds.require_physical(9, value)


@pytest.mark.parametrize(
    "energy,transfer", ((4, -1), (17, -1), (9, 0), (9, -5), (9, 1))
)
def test_physical_endpoints_and_wrong_energy_rejected(energy, transfer):
    with pytest.raises(ValueError):
        bounds.require_physical(energy, transfer)


def test_remaining_physical_coordinates_are_required_not_default_zero():
    with pytest.raises(TypeError):
        matching.remaining_local(9, -1, 1, source.KAPPA)
    with pytest.raises(TypeError):
        matching.remaining_local(9, -1, 1, source.KAPPA, 0, 0)


def test_known_reference_can_be_changed_by_unassigned_matching():
    a, b, mu, k = ins.A, ins.B, ins.MU, ins.K
    f = matching.remaining_local(
        a, b, mu, k, matching.ALPHA, matching.BETA, matching.DK
    )
    assert all(s.diff(f, x) != 0 for x in (matching.ALPHA, matching.BETA, matching.DK))
    assert s.factor(s.diff(f, matching.DK) + ins.gravity_born(a, b, mu, k) / k) == 0


def test_no_original_status_or_parameter_change():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.frontier() == audit.previous.frontier()
    assert audit.qualifications() == audit.previous.qualifications()
    assert audit.require_parameters(audit.parameters()) == audit.parameters()


def test_fixed_local_coefficients_are_prior_source_values():
    alpha, beta, dk = matching.fixed_shifts()
    ell = s.log(source.HEAVY_MASS2) + 2
    assert alpha == -ell / 40 and beta == -3 * ell / 10
    assert dk == source.poles.delta_kappa()
    assert source.gaussian.fixed_coefficients()["R_old"] * 2 == dk


def test_M1_scale_derivative_does_not_fix_a_local_value():
    scale = s.Symbol("scale_squared", positive=True)
    a, b, mu, k = ins.A, ins.B, ins.MU, ins.K
    q2 = sum(x * x for x, _, _ in ins.channels(a, b, mu))
    derivative = scale * s.diff(ins.massless_known(a, b, mu, k, scale), scale)
    assert s.factor(derivative - (q2 + 12 * mu**2) / (640 * s.pi**2 * k**2)) == 0


def test_exact_Gaussian_bound_not_unknown_higher_loop_bound():
    record = bounds.data()
    assert record["whole_complex_disk_radius"] == 16
    assert record["whole_Gaussian_Dyson_relative_remainder_bound"] < s.Rational(
        5, 10**1204
    )
    assert all(v > 0 for v in record["whole_positive_arithmetic_margins"].values())
    assert "independent higher-loop" in record["Gaussian_disk_boundary"]
