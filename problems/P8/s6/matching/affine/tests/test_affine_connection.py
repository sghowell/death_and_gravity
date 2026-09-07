"""Literal tensor, full-rank, Schur, source-typo, and exclusion controls."""

from itertools import permutations, product

import pytest
import sympy as sp
from p8_affine import connection as c


def _zero(value):
    if isinstance(value, sp.MatrixBase):
        return all(entry == 0 for entry in value)
    return value == 0


def test_all_exact_residuals_and_dimensions():
    assert len(c.checks()) == 30
    assert all(_zero(value) for value in c.checks().values())
    assert c.literal()["hessian"].shape == (64, 64)
    assert c.quotient()["hessian"].shape == (60, 60)


@pytest.mark.parametrize("p,s,cc,quartic", [
    (sp.Rational(1, 2), 1, 0, 0),
    (sp.Rational(2, 3), 2, sp.Rational(-3, 7), sp.Rational(5, 11)),
    (sp.Rational(1, 4), 3, sp.Rational(7, 5), sp.Rational(-2, 9)),
])
def test_quadratic_density_from_independent_full_curvature_and_epsilon(p, s, cc, quartic):
    """No Ricci symmetry, no torsion restriction, and a dense exact distortion."""
    k = {(a, b, d): sp.Rational(((19*a+7*b+3*d) % 13)-6, 11)
         for a, b, d in product(range(4), repeat=3)}
    metric = (-1, 1, 1, 1)

    def epsilon(a, b, d, e):
        indices = (a, b, d, e)
        if len(set(indices)) != 4:
            return 0
        return (-1)**sum(indices[i] > indices[j] for i in range(4) for j in range(i+1, 4))

    def curvature(a, b, d, e):
        return sum(k[a, z, d]*k[z, b, e]-k[a, z, e]*k[z, b, d] for z in range(4))

    ricci = sum(metric[b]*curvature(a, b, a, b)
                for a, b in product(range(4), repeat=2))
    double_dual = sum(
        sp.Rational(1, 4)*epsilon(gamma, 0, a, b)*metric[gamma]
        *epsilon(gamma, 0, d, e)*metric[a]*curvature(a, b, d, e)
        for gamma, a, b, d, e in product(range(4), repeat=5))
    gal4 = sum(
        epsilon(0, i, j, ell)*metric[ell]*epsilon(0, m, n, ell)
        *s**4*k[0, m, i]*k[0, n, j]
        for i, j, ell, m, n in product(range(4), repeat=5))
    independent = p*ricci+cc*s**2*double_dual+quartic*gal4
    vector = sp.Matrix([k[key] for key in c.INDICES])
    matrix = c.literal()["hessian"].subs({c.P: p, c.S: s, c.C: cc, c.F4: quartic})
    assert (vector.T*matrix*vector)[0]/2 == independent


def test_cubic_and_quartic_epsilon_index_order_is_not_frobenius_square():
    """A nonsymmetric affine Hessian must use tr(F²), not sum F_ij²."""
    h = sp.Matrix([[1, 2, 0], [-3, 4, 5], [7, 0, -2]])
    epsilon = sp.LeviCivita
    direct = sum(epsilon(i, j, ell)*epsilon(m, n, ell)*h[i, m]*h[j, n]
                 for i, j, ell, m, n in product(range(3), repeat=5))
    assert direct == sp.trace(h)**2-sp.trace(h*h)
    assert direct != sp.trace(h)**2-sum(value**2 for value in h)


def test_linear_source_from_coefficient_derivative_before_rest_frame():
    """Differentiate the unrestricted gradient product, then impose v=(s,0,0,0)."""
    metric = (-1, 1, 1, 1)
    s = sp.Rational(3, 2)
    pp, px, cc, cp, cx, cubic, quartic = map(sp.Rational, (2, 3, 5, 7, 11, 13, 17))
    h = sp.Matrix([[2, -1, 3, 5], [-1, 7, 2, -4], [3, 2, 11, 6], [5, -4, 6, 13]])
    k = {(a, b, d): sp.Rational(((5*a+11*b+17*d) % 19)-9, 7)
         for a, b, d in product(range(4), repeat=3)}
    eps = sp.Symbol("eps", real=True)
    coefficient_derivative = {}
    for e, alpha, beta in product(range(4), repeat=3):
        v = sp.Matrix([s*int(i == 0)+eps*h[e, i] for i in range(4)])
        x = sum(metric[i]*v[i]**2 for i in range(4))
        value = cc+cp*eps*s*int(e == 0)+cx*(x+s**2)
        coefficient_derivative[e, alpha, beta] = sp.diff(value*v[alpha]*v[beta], eps).subs(eps, 0)
    dp = [pp*s*int(i == 0)-2*px*s*h[i, 0] for i in range(4)]
    linear = sum(metric[b]*(-dp[a]*k[a, b, b]+dp[b]*k[a, b, a])
                 for a, b in product(range(4), repeat=2))
    for gamma, alpha, a, b in permutations(range(4)):
        for beta, d, e in permutations(i for i in range(4) if i != gamma):
            weight = (sp.LeviCivita(gamma, alpha, a, b)*sp.LeviCivita(gamma, beta, d, e)
                      *metric[a]*metric[gamma]/4)
            linear += weight*(-coefficient_derivative[d, alpha, beta]*k[a, b, e]
                              +coefficient_derivative[e, alpha, beta]*k[a, b, d])
    spatial = sp.Matrix(3, 3, lambda i, j: h[i+1, j+1]-eps*s*k[0, j+1, i+1])
    galileon = 2*cubic*s**2*sp.trace(spatial)+quartic*s**2*(
        sp.trace(spatial)**2-sp.trace(spatial*spatial))
    linear += sp.diff(galileon, eps).subs(eps, 0)
    substitutions = {c.S: s, c.PP: pp, c.PX: px, c.C: cc, c.CP: cp,
                     c.CX: cx, c.F3: cubic, c.F4: quartic}
    substitutions.update({c.H[a, b]: h[a, b] for a, b in c.H_PAIRS})
    source = c.literal()["source"].subs(substitutions)
    assert (source.T*sp.Matrix([k[key] for key in c.INDICES]))[0] == linear


def test_four_gauge_columns_and_trace_slice_are_complementary():
    data = c.quadratic()
    embed = c.quotient()["embedding"]
    gauge = data["gauge"]
    assert gauge.T*gauge == 4*sp.eye(4)
    assert gauge.T*embed == sp.zeros(4, 60)
    assert sp.Matrix.hstack(embed, gauge).rank() == 64


def test_literal_source_vanishes_on_projective_directions_before_cd():
    gauge = c.quadratic()["gauge"]
    assert c.literal()["hessian"]*gauge == sp.zeros(64, 4)
    assert c.literal()["source"].T*gauge == sp.zeros(1, 4)


def test_full_quotient_block_factorization():
    q = c.quotient()
    p = c.P
    assert tuple(map(len, q["blocks"])) == (9, 9, 6, 6, 9, 6, 9, 3, 3)
    assert q["block_determinants"] == (
        8192*p**14, -32*p**6, p**4*(8*p**2-1), p**4*(8*p**2-1),
        8192*p**14, p**4*(8*p**2-1), 8192*p**14, -16*p**6, -16*p**6)
    assert q["determinant"] == -2**52*p**72*(8*p**2-1)**3
    for first, second in permutations(q["blocks"], 2):
        assert q["hessian"].extract(first, second) == sp.zeros(len(first), len(second))


def test_exceptional_rank_and_source_compatibility_do_not_imply_uniqueness():
    data = c.exceptional()
    assert data["quotient_nullity"] == 3
    assert data["quotient_rank"] == 57
    assert data["quotient_kernel"].rank() == 3
    assert _zero(data["kernel_residual"])
    assert _zero(data["source_compatibility_residual"])
    # The familiar sourced particular solution stays finite; the full inverse does not.
    assert all(not value.has(sp.zoo, sp.nan) for value in
               c.eliminate()["solution"].subs(c.P, data["p"]))
    assert c.quotient()["determinant"].subs(c.P, data["p"]) == 0


def test_actual_tube_margin_and_inverse_bound_include_closed_endpoints():
    w = sp.symbols("w", real=True)
    assert (8*w-1).subs(w, sp.Rational(9, 40)) == sp.Rational(4, 5)
    assert sp.Rational(9, 40) > sp.Rational(9, 20)**2
    assert sp.Rational(11, 40) < sp.Rational(11, 20)**2
    bound = c.inverse_bound()
    assert len(bound["row_bounds"]) == 60
    assert bound["norm_upper"] == sp.Rational(880109, 36000) < 25
    assert all(_zero(value) for value in bound["inverse_residuals"])
    for matrix in c.quotient()["block_matrices"]:
        inverse = matrix.subs(c.P, sp.Rational(1, 2)).inv()
        assert max(sum(abs(value) for value in inverse.row(i))
                   for i in range(inverse.rows)) <= bound["norm_upper"]


def test_auxiliary_quadratic_form_is_indefinite_without_a_kinetic_verdict():
    matrix = c.quotient()["block_matrices"][-1]
    plus, minus = sp.Matrix([1, 1, 0]), sp.Matrix([1, -1, 0])
    assert (plus.T*matrix*plus)[0] == -4*c.P**2
    assert (minus.T*matrix*minus)[0] == 4*c.P**2
    assert c.calibration()["algebraic_only"] is True
    assert c.calibration()["UV_or_kinetic_health_claim"] is False


def test_schur_density_equals_literal_stationary_value():
    data = c.eliminate()
    k = data["solution"]
    full = c.quadratic()
    difference = (k.T*full["hessian"]*k)[0]/2+(full["source"].T*k)[0]/2
    assert sp.factor(difference) == 0
    assert _zero(data["full_euler_residual"])
    assert c.coefficients()["density_reconstruction_residual"] == 0


def test_corrected_q2_is_direct_not_imported_and_printed_version_fails():
    result = c.coefficients()
    p, s = c.P, c.S
    t = c.PP+c.F3*s**2
    assert sp.factor(result["Q1"]-(4*p*t-2*c.PP-c.CP*s**2)) == 0
    assert sp.factor(result["Q2"]-(4*(p+3*s**2*c.PX)*t-2*c.PP-c.CP*s**2)/s**2) == 0
    defect = c.source_formula_comparison()["printed_Q2_defect"]
    assert defect == 12*c.PX*(s**2+1)*t/s**2
    assert defect.subs({c.PX: 1, c.PP: 1, c.F3: 0, s: 2}) == 15
    assert defect.subs(c.PX, 0) == 0  # why an X-independent check would miss it


def test_pure_palatini_direct_connection_fixes_lower_derivative_typo():
    result = c.palatini_control()
    assert _zero(result["full_euler_residual"])
    assert result["density_residual"] == 0
    assert result["Q1"] == 0
    assert result["Q2"] == 6*c.PP*c.PX/c.P
    assert result["printed_Q2_defect"] != 0
    assert c.quadratic()["gauge"].T*result["solution"] == sp.zeros(4, 1)


@pytest.mark.parametrize("p,s", [(sp.Rational(1, 2), 1),
                                  (sp.sqrt(sp.Rational(9, 40)), 2),
                                  (sp.sqrt(sp.Rational(11, 40)), sp.Rational(1, 3)),
                                  (sp.Rational(1, 4), 1)])
def test_exact_domain_acceptance_does_not_assume_positive_hessian(p, s):
    assert c.require_domain(p, s) == (p, sp.sympify(s))


@pytest.mark.parametrize("p,s", [(0, 1), (-1, 1), (1, 0), (1, -1),
                                  (sp.sqrt(2)/4, 1), (sp.oo, 1), (1, sp.nan),
                                  (sp.I, 1), (sp.Symbol("unknown"), 1),
                                  (sp.Float("0.5"), 1), ("True", 1), ([1], 1)])
def test_exact_domain_rejects_unproved_or_singular_parameters(p, s):
    with pytest.raises((TypeError, ValueError)):
        c.require_domain(p, s)


@pytest.mark.parametrize("p,s", [(True, 1), (1, False), (0.5, 1), (1, 1.0),
                                  ("1/2", 1), (1, "2")])
def test_exact_domain_rejects_inexact_boolean_and_string_inputs(p, s):
    with pytest.raises(TypeError):
        c.require_domain(p, s)


@pytest.mark.parametrize("indices", [(True, 0, 0), (0.0, 0, 0), (-1, 0, 0), (0, 4, 0)])
def test_index_domain_controls(indices):
    with pytest.raises((TypeError, ValueError)):
        c.index(*indices)
