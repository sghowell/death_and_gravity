"""Independent geometric, variational and signed-budget checks."""

import json

import mpmath as mp
import pytest
import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_canonical_bounce_gap import audit, equations, map, nullstress, verify
from p8_vacuum_clock_transparent_map import gate
from p8_vacuum_curved_dirac_stress import audit as previous
from p8_vacuum_curved_dirac_stress import stress


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_unsupported_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_direct_connection_curvature_and_physical_null_Einstein_sign():
    a, ad, add, N, Nd, Ndd = s.symbols("a adot addot N Ndot Nddot", positive=True)
    metric = s.diag(N * N, -a * a, -a * a, -a * a)
    inverse = metric.inv()
    dt = lambda f: (
        s.diff(f, a) * ad
        + s.diff(f, ad) * add
        + s.diff(f, N) * Nd
        + s.diff(f, Nd) * Ndd
    )
    partial = lambda f, index: dt(f) if index == 0 else s.Integer(0)
    connection = {}
    for r in range(4):
        for mu in range(4):
            for nu in range(4):
                connection[r, mu, nu] = s.simplify(
                    sum(
                        inverse[r, l]
                        * (
                            partial(metric[l, nu], mu)
                            + partial(metric[l, mu], nu)
                            - partial(metric[mu, nu], l)
                        )
                        / 2
                        for l in range(4)
                    )
                )
    Ric = s.zeros(4)
    for mu in range(4):
        for nu in range(4):
            Ric[mu, nu] = s.simplify(
                sum(
                    partial(connection[r, nu, mu], r)
                    - partial(connection[r, r, mu], nu)
                    + sum(
                        connection[r, r, l] * connection[l, nu, mu]
                        - connection[r, nu, l] * connection[l, r, mu]
                        for l in range(4)
                    )
                    for r in range(4)
                )
            )
    scalar = s.simplify(
        sum(inverse[i, j] * Ric[i, j] for i in range(4) for j in range(4))
    )
    Ein = Ric - metric * scalar / 2
    gauge = {N: 1, Nd: 0, Ndd: 0}
    assert s.simplify(scalar.subs(gauge) + 6 * (add / a + ad * ad / (a * a))) == 0
    assert s.simplify(Ein[0, 0].subs(gauge) - 3 * ad * ad / (a * a)) == 0
    null = s.simplify((Ein[0, 0] + Ein[1, 1] / (a * a)).subs(gauge))
    assert s.simplify(null + 2 * (add / a - ad * ad / (a * a))) == 0
    assert null.subs({a: 1, ad: 0, add: 4}) == -8


@pytest.mark.parametrize(
    "x", [-s.Rational(1, 2), -s.Rational(1, 3), 0, s.Rational(1, 3), s.Rational(1, 2)]
)
def test_entire_central_interval_algebraic_factor_and_signed_margin(x):
    Hd = 4 * (1 - x * x) / (1 + x * x) ** 2
    remainder = 4 * (1 - 4 * x * x) * (13 + 3 * x * x) / (25 * (1 + x * x) ** 2)
    assert Hd - s.Rational(48, 25) == remainder
    assert remainder >= 0
    bound = nullstress.data()["full_uniform_absolute_free_null_stress_upper"]
    residual = analytic.KAPPA * (2 * Hd + 1) - bound - analytic.KAPPA / 2
    assert (
        residual
        >= equations.margins(analytic.KAPPA / 2)[
            "central_interval_null_equation_residual_lower"
        ]
    )
    assert residual > 4 * analytic.KAPPA


@pytest.mark.parametrize("p,q", [(0, 0), (1, 0), (0, 1), (-2, 3), (7, -11)])
def test_positive_scalar_metric_and_gauge_sources_do_not_repair_null_sign(p, q):
    C = s.Matrix([[2, 1], [0, 3]])
    G = C.T * C
    velocity = s.Matrix([p, q])
    assert G.det() > 0 and G[0, 0] > 0
    extra = (velocity.T * G * velocity)[0]
    gauge = s.Rational(2, 3) * (p * p + q * q)
    assert extra >= 0 and gauge >= 0
    assert (
        9 * analytic.KAPPA
        + extra
        + gauge
        - nullstress.data()["actual_bounce_absolute_free_null_stress_upper"]
        > 8 * analytic.KAPPA
    )


def test_complete_correlated_null_budget_uses_no_potential_magnitude():
    old = stress.data()["actual_uniform_exact_rational_enclosures"]
    r = old["complete_remainder_components"]
    l = old["complete_local_two_derivative_components"]
    n = nullstress.data()
    bounce = (
        r["complete_subtracted_energy_remainder"]
        + r["complete_subtracted_pressure_remainder"]
    )
    uniform = (
        bounce
        + l["local_two_derivative_energy"]
        + l["local_two_derivative_pressure"]
        + old["finite_Euler_energy_upper"]
        + old["finite_Euler_pressure_upper"]
    )
    assert n["actual_bounce_absolute_free_null_stress_upper"] == bounce
    assert n["full_uniform_absolute_free_null_stress_upper"] == uniform
    assert 0 < bounce < 10**416 and 0 < uniform < 10**595
    assert uniform < old["complete_potential_upper"]
    assert (
        old["complete_absolute_curved_energy_upper"]
        + old["complete_absolute_curved_pressure_upper"]
        > 10**788
    )


def test_nonzero_bounce_components_cancel_only_in_their_correlated_null_sum():
    with mp.workdps(700):
        m = mp.mpf(10) ** 200
        d = 3 * mp.mpf(10) ** 197
        tau = mp.mpf(10) ** -100
        Q = 16 * mp.pi**2
        for sign in (-1, 1):
            Mdot = sign * d / tau
            # Independent full covariant two-derivative action at H=0,ell=0.
            kinetic = -2 * Mdot * Mdot / (3 * Q)
            Fddot = -4 * Mdot * Mdot / (6 * Q)
            rho = kinetic
            P = kinetic - 2 * Fddot
            assert rho < 0 and P > 0
            assert abs(rho + P) < mp.mpf("1e-90")
            assert abs(rho) > mp.mpf("1e590")
            # Without the Newton anchor an extra constant-mass Hdot term remains.
            without_reference = P - 2 * m * m * 4 / (3 * Q)
            assert abs(without_reference - P) > mp.mpf("1e397")


@pytest.mark.parametrize("r0,dx", [(1, 1), (-3, 2), (5, -7)])
def test_finite_clock_map_first_variation_not_just_its_value(r0, dx):
    with mp.workdps(100):
        q = s.lambdify(gate.X, gate.data()["gate"], "mpmath")
        h = mp.mpf("1e-8")
        actual = (q(1 + h * dx) * (r0 + h) - q(1 - h * dx) * (r0 - h)) / (2 * h)
        assert abs(actual) < mp.mpf("1e-44")
        wrong = ((1 - (1 + h * dx)) * (r0 + h) - (1 - (1 - h * dx)) * (r0 - h)) / (
            2 * h
        )
        assert abs(wrong + r0 * dx) < mp.mpf("1e-80")
        assert abs(wrong) > mp.mpf(".9")


def test_small_additional_budget_is_conditional_not_a_universal_restriction():
    k = analytic.KAPPA
    free = nullstress.data()
    critical = (
        s.Rational(121, 25) * k - free["full_uniform_absolute_free_null_stress_upper"]
    )
    assert (
        equations.margins(critical)["central_interval_null_equation_residual_lower"]
        == 0
    )
    assert (
        equations.margins(100 * k)["central_interval_null_equation_residual_lower"] < 0
    )
    n = s.Symbol("actual_free_null", real=True)
    repair = -9 * k - n
    assert 9 * k + n + repair == 0
    assert "not asserted absent" in equations.data()["additional_terms_boundary"]
    assert "No physical bound" in equations.data()["not_a_full_candidate_exclusion"]


def test_leading_nonminimal_operator_is_outside_the_diagnostic_not_excluded():
    # This cancels only one pointwise equation, not a controlled global solution.
    t, xi, k = s.symbols("t xi kappa", real=True)
    F = xi * k * t * t  # xi*Phi^2 on Phi=sqrt(kappa)t
    null_at_bounce = -2 * s.diff(F, t, 2).subs(t, 0)
    assert s.simplify(null_at_bounce + 4 * xi * k) == 0
    assert (9 * k + null_at_bounce).subs(xi, s.Rational(9, 4)) == 0
    effective_planck = k - 2 * F
    assert effective_planck.subs({xi: s.Rational(9, 4), t: s.sqrt(2) / 3}) == 0
    assert effective_planck.subs({xi: s.Rational(9, 4), t: 1}) == -s.Rational(7, 2) * k


def test_scope_ancestry_decimal_and_exact_serializer():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert "NOT_FULL_INTERACTING_PARENT_OR_DHOST_ROW" in audit.matching()[-1]["status"]
    assert "nearby metric" in nullstress.data()["scope"]
    assert "classical action/first-variation" in map.data()["scope"]
    assert "No parent field equation" in map.data()["classical_stress_consequence"]
    for mod in audit.MODULES:
        json.dumps(verify.serialize(verify.payload(mod.data())))
    assert all(
        type(v) is str for v in nullstress.data()["decimal_diagnostics_only"].values()
    )
    with pytest.raises(ValueError, match="Inexact or nonfinite"):
        verify.serialize(s.Float(1))
