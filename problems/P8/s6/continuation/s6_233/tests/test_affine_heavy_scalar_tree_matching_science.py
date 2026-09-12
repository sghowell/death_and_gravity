"""Independent canonical model, literal exchanges and massive all-angle tree matching."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import audit, cut, matching, model


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_all_exact_residuals(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_all_proof_gates(name, value):
    assert value is True, name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_original_frontier_and_separately_named_model():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 89
    assert audit.frontier() == audit.previous.frontier()
    assert all(row["status"] != "COMPLETE" for row in audit.frontier())
    assert "SEPARATE_CLASSICAL_TWO_SCALAR_MODEL" in audit.ITEM["status"]
    assert audit.require_model("V2S-T1", 2) == ("V2S-T1", s.Integer(2))
    assert audit.require_model("V2S-T1", matching.ENERGY) == ("V2S-T1", matching.ENERGY)


def test_independent_full_kinetic_mass_and_potential_square():
    D, g = s.symbols("D g", positive=True)
    phi, H, pdot, Hdot = s.symbols("phi H pdot Hdot", real=True)
    M2 = D + 2
    C = -g * g * (3 / D - 2 / D**2)
    L = (
        (pdot * pdot + Hdot * Hdot) / 2
        - phi * phi / 2
        - M2 * H * H / 2
        + g * H * phi * phi / 2
        + C * phi**4 / 24
    )
    kinetic = s.hessian(L, (pdot, Hdot))
    V = -L.subs({pdot: 0, Hdot: 0})
    assert kinetic == s.eye(2)
    assert s.hessian(V, (phi, H)).subs({phi: 0, H: 0}) == s.diag(1, M2)
    valley = g * phi * phi / (2 * M2)
    rest = g * g * (D - 1) * phi**4 / (6 * D * D * (D + 2))
    assert s.cancel(V - phi * phi / 2 - M2 * (H - valley) ** 2 / 2 - rest) == 0
    coefficient = g * g * (D - 1) / (6 * D * D * (D + 2))
    assert (
        s.cancel(s.diff(V, phi).subs(H, valley) - phi - 4 * coefficient * phi**3) == 0
    )


def test_independent_global_quadratic_coercivity_certificate():
    D, g, phi, H = model.D, model.g, model.phi, model.H
    remainder = (
        phi * phi / 3
        + (3 * D + 5) * (H - 3 * g * phi * phi / (2 * (3 * D + 5))) ** 2 / 6
        + g * g * (9 * D - 10) * phi**4 / (24 * D * D * (3 * D + 5))
    )
    assert s.cancel(model.V - (phi * phi + H * H) / 6 - remainder) == 0
    assert model.GAP > 2
    assert (9 * model.GAP - 10) > 0
    assert model.data()["canonical_kinetic_matrix"] == s.eye(2)


def test_omitted_contact_and_invalid_gap_destroy_potential_control():
    phi = s.Integer(100)
    D = s.Rational(1, 2)
    g = s.S.One
    H = g * phi * phi / (2 * (D + 2))
    bad = model.V.subs({model.D: D, model.g: g, model.phi: phi, model.H: H})
    assert bad < 0
    # Even for a healthy heavy mass, the bare positive pole by itself leaves
    # an unbounded negative quartic along its heavy potential valley.
    D = s.Integer(10)
    H = g * phi * phi / (2 * (D + 2))
    without_contact = phi * phi / 2 + (D + 2) * H * H / 2 - g * H * phi * phi / 2
    assert without_contact < 0
    assert model.CONTACT < 0


def test_independent_literal_vertex_and_Feynman_exchange_signs():
    phi, H, g, C, channel, M2 = s.symbols("phi H g C channel M2")
    assert s.diff(g * H * phi**2 / 2, phi, phi, H) == g
    assert s.diff(C * phi**4 / 24, phi, 4) == C
    exchange = (s.I * g) ** 2 * s.I / (channel - M2) / s.I
    assert s.cancel(exchange - g * g / (M2 - channel)) == 0


def test_independent_complete_original_massive_tree_geometric_matching():
    S, T, U, D, g = s.symbols("s t u D g")
    C = -g * g * (3 / D - 2 / D**2)
    raw = C + g * g * sum(1 / (D + 2 - z) for z in (S, T, U))
    lam, gam = g * g / (2 * D**3), g * g / D**4
    target = (
        2 * lam * sum((z - 2) ** 2 for z in (S, T, U)) + 3 * gam * S * T * U - 8 * gam
    )
    remainder = gam * sum((z - 2) ** 4 / (D - (z - 2)) for z in (S, T, U))
    shell = {U: 4 - S - T}
    assert s.cancel((raw - target - remainder).subs(shell)) == 0
    assert (
        s.expand((sum((z - 2) ** 3 for z in (S, T, U)) - 3 * S * T * U + 8).subs(shell))
        == 0
    )
    assert s.cancel((raw - C - target - remainder).subs(shell) + C) == 0
    assert C != 0


@pytest.mark.parametrize(
    "D,S,x",
    [
        (100, 4, s.Rational(-1, 2)),
        (100, 6, s.Rational(1, 3)),
        (1000, 10, s.Rational(4, 5)),
        (1000000, 100, s.Rational(-3, 7)),
    ],
)
def test_independent_physical_angles_and_relative_error(D, S, x):
    D, S = s.Integer(D), s.Integer(S)
    g = s.Rational(1, 7)
    lam, gam = g * g / (2 * D**3), g * g / D**4
    T, U = -(S - 4) * (1 - x) / 2, -(S - 4) * (1 + x) / 2
    old = (
        2 * lam * ((S - 2) ** 2 + (T - 2) ** 2 + (U - 2) ** 2)
        + 3 * gam * S * T * U
        - 8 * gam
    )
    C = -g * g * (3 / D - 2 / D**2)
    new = C + g * g * (1 / (D + 2 - S) + 1 / (D + 2 - T) + 1 / (D + 2 - U))
    R = gam * sum((z - 2) ** 4 / (D - (z - 2)) for z in (S, T, U))
    r = (S - 2) / D
    assert 0 < R and old >= lam * (S - 2) ** 2 > 0
    assert new - old == R
    assert R / old <= 6 * r * r / (1 - r)
    assert R / old < s.Rational(1, 60)


def test_exact_actual_energy_window_and_first_elastic_ratio_constants():
    r = (matching.ENERGY**2 - 2) / model.GAP
    assert r < s.Rational(32, 625) < 1
    rational = 6 * s.Rational(32, 625) ** 2 / (1 - s.Rational(32, 625))
    assert rational == s.Rational(6144, 370625) < s.Rational(1, 60)
    assert s.Rational(61, 60) ** 2 - 1 == s.Rational(121, 3600)
    assert 4 < matching.ENERGY**2 < model.MASS2 < 10**198
    assert model.LAMBDA > 2 * model.GAMMA


def test_independent_complete_forward_remainder_and_higher_coefficients():
    D, g, v = s.symbols("D g v", positive=True)
    gam = g * g / D**4
    R = gam * (v**4 / (D - v) + 16 / (D + 2) + v**4 / (D + v))
    assert s.cancel(R - 16 * gam / (D + 2) - 2 * gam * D * v**4 / (D**2 - v**2)) == 0
    assert s.diff(R, v, 2).subs(v, 0) == 0
    assert s.cancel(s.diff(R, v, 4).subs(v, 0) / 24 - 2 * gam / D) == 0
    assert R.subs(v, 0) > 0
    assert model.G2 == s.Rational(1, 2**26)


@pytest.mark.parametrize("D,S", [(100, 6), (1000, 10), (1000000, 100)])
def test_independent_full_angular_integrals_and_first_cut_diagnostic(D, S):
    # Finite-precision diagnostic at separate modest parameter fixtures, not
    # the continuum proof or a change of the actual model parameters.
    with mp.workdps(80):
        D, S = mp.mpf(D), mp.mpf(S)
        g = mp.mpf(1) / 7
        g2 = g * g
        C = -g2 * (3 / D - 2 / D**2)
        k = (S - 4) / 2
        a = D + 2 + k
        d = C + g2 / (D + 2 - S)
        A = lambda x: d + 2 * g2 * a / (a * a - k * k * x * x)
        logarithm = mp.log((a + k) / (a - k))
        closed1 = 2 * d + 2 * g2 * logarithm / k
        closed2 = (
            2 * d * d
            + (4 * d * g2 / k + 2 * g2 * g2 / (a * k)) * logarithm
            + 4 * g2 * g2 / (a * a - k * k)
        )
        literal1 = mp.quad(A, [-1, 0, 1])
        literal2 = mp.quad(lambda x: A(x) ** 2, [-1, 0, 1])
        assert abs(literal1 / closed1 - 1) < mp.mpf("1e-60")
        assert abs(literal2 / closed2 - 1) < mp.mpf("1e-60")
        lam, gam = g2 / (2 * D**3), g2 / D**4

        def old(x):
            T, U = -(S - 4) * (1 - x) / 2, -(S - 4) * (1 + x) / 2
            return (
                2 * lam * ((S - 2) ** 2 + (T - 2) ** 2 + (U - 2) ** 2)
                + 3 * gam * S * T * U
                - 8 * gam
            )

        old_norm = mp.quad(lambda x: old(x) ** 2, [-1, 0, 1])
        assert 1 < literal2 / old_norm < (mp.mpf(61) / 60) ** 2
        fourth = mp.quad(lambda x: A(x) * (35 * x**4 - 30 * x * x + 3) / 8, [-1, 0, 1])
        assert fourth > 0


@pytest.mark.parametrize(
    "n,ell",
    [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2), (4, 2), (4, 3), (4, 4)],
)
def test_independent_positive_even_power_Legendre_coefficients(n, ell):
    x = s.Symbol("x")
    P = s.diff((x * x - 1) ** (2 * ell), x, 2 * ell) / (
        2 ** (2 * ell) * s.factorial(2 * ell)
    )
    coefficient = s.integrate(s.expand(x ** (2 * n) * P), (x, -1, 1))
    assert coefficient >= 0
    assert (coefficient > 0) == (n >= ell)


def test_complete_angular_primitives_and_removable_threshold():
    a, k, g2, d, x = cut.a, cut.k, cut.g2, cut.d, cut.x
    log = s.log((a + k * x) / (a - k * x))
    first = d * x + g2 * log / k
    second = (
        d * d * x
        + (2 * d * g2 / k + g2 * g2 / (a * k)) * log
        + 2 * g2 * g2 * x / (a * a - k * k * x * x)
    )
    assert s.cancel(s.diff(first, x) - cut.A) == 0
    assert s.cancel(s.diff(second, x) - cut.A * cut.A) == 0
    assert s.cancel(s.limit(cut.I1, k, 0) - 2 * d - 4 * g2 / a) == 0
    assert s.cancel(s.limit(cut.I2, k, 0) - 2 * (d + 2 * g2 / a) ** 2) == 0


@pytest.mark.parametrize(
    "energy_squared,angle",
    [(5, s.Rational(-2, 3)), (16, s.Rational(3, 5)), (10**196, s.Rational(1, 3))],
)
def test_actual_hierarchy_with_enough_precision_for_contact_cancellation(
    energy_squared, angle
):
    # 1000 digits resolve the original low-energy amplitude and its much
    # smaller remainder after the very large contact/exchange cancellation.
    # This is a diagnostic; exact algebra and inequalities certify the result.
    with mp.workdps(1000):
        lam = mp.power(10, -600)
        gam = 1024 * mp.power(10, -800)
        D = 2 * lam / gam
        g2 = mp.mpf(1) / (2**26)
        C = -g2 * (3 / D - 2 / D**2)
        S = mp.mpf(energy_squared)
        x = mp.mpf(int(angle.p)) / int(angle.q)
        T, U = -(S - 4) * (1 - x) / 2, -(S - 4) * (1 + x) / 2
        raw = C + g2 * (1 / (D + 2 - S) + 1 / (D + 2 - T) + 1 / (D + 2 - U))
        target = (
            2 * lam * ((S - 2) ** 2 + (T - 2) ** 2 + (U - 2) ** 2)
            + 3 * gam * S * T * U
            - 8 * gam
        )
        remainder = gam * sum((z - 2) ** 4 / (D - (z - 2)) for z in (S, T, U))
        assert remainder > 0 and target > 0
        assert abs((raw - target) / remainder - 1) < mp.mpf("1e-100")
        assert 0 < (raw - target) / target < mp.mpf(1) / 60


def test_named_window_has_no_heavy_production_and_tree_is_not_exact_unitary():
    assert matching.ENERGY**2 < model.MASS2 < 4 * model.MASS2
    # Every named interaction contains an even number of light fields.
    phi = model.phi
    assert s.expand(model.V - model.V.subs(phi, -phi)) == 0
    # At one physical angle a real positive tree supplies a nonzero elastic
    # optical RHS, but the real tree alone has no continuous absorptive part.
    value = matching.TREE.subs(
        {
            matching.S: 5,
            matching.T: -s.Rational(1, 2),
            model.D: 100,
            model.g: s.Rational(1, 7),
        }
    )
    assert value > 0 and s.im(value) == 0 and value * value > 0
    assert "not the full real loop amplitude" in cut.data()["first_absorptive_scope"]


def test_named_contact_does_not_assert_full_parent_or_constant_matching():
    constant = 16 * model.GAMMA / (model.GAP + 2)
    assert constant > 0
    assert model.data()["separately_named_model"] == "V2S-T1"
    assert "not the original Proca/DHOST affine parent" in model.data()["model_scope"]
    assert "not the full angular quantum error premise" in matching.data()["boundary"]
