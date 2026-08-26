"""S2 core tests: series arithmetic, exact polynomial systems, exact algebraic identities."""
from fractions import Fraction

import numpy as np
import pytest
from flint import arb, fmpq, fmpq_mpoly_ctx
from p4 import css
from p4.validated import recursion, sonic
from p4.validated.arbseries import Series, precision
from p4.validated.biseries import BiSeries
from p4.validated.systems import (
    centre_system,
    sonic_constraint_poly,
    sonic_constraint_propagation,
    sonic_system,
)


def _fq(x):
    fr = Fraction(float(x))
    return fmpq(fr.numerator, fr.denominator)


def test_series_and_biseries_arithmetic():
    with precision(128):
        dm = 3
        a = BiSeries.from_blocks([Series([1, 2], dm), Series([3, 1, 1], dm), Series([0.5], dm)], 6, dm)
        b = BiSeries.from_blocks([Series([2], dm), Series([1, -1], dm)], 6, dm)
        ab = a * b
        for n in range(4):
            man = Series([0], dm)
            for j in range(3):
                if 0 <= n - j < 2:
                    man = man + a[j] * b[n - j]
            assert all((ab[n][k] - man[k]).contains(arb(0)) for k in range(dm))
        assert [float(c) for c in a.deriv()[0].coeffs()] == [3.0, 1.0, 1.0]
        assert [float(c) for c in a.euler()[1].coeffs()] == [3.0, 1.0, 1.0]
        s = Series([4, 1, 0.5], 5)
        assert all(x.contains(arb(0)) for x in (s.inv() * s - 1).coeffs()[1:])
        assert all(x.contains(arb(0)) for x in (s.sqrt() * s.sqrt() - s).coeffs())


@pytest.mark.parametrize("seed", [0, 1])
def test_polynomial_systems_match_S1_rhs(seed):
    rng = np.random.default_rng(seed)
    S, C = sonic_system(), centre_system()
    for _ in range(10):
        y = np.array([1 + 2 * rng.random(), 0.3 + 3 * rng.random(), 0.05 + rng.random(), rng.uniform(-0.9, 0.9)])
        rhs = css.rhs_plain(0.0, y)
        args = [_fq(0)] + [_fq(z) for z in y]
        for r in range(4):
            lhs = sum(float(S.P[r][i](*args)) * rhs[i] for i in range(4))
            q = float(S.Q[r](*args))
            assert abs(lhs - q) < 1e-12 * (1 + abs(q))
        eps = rng.uniform(0.01, 0.3)
        Y = np.array([rng.uniform(0.5, 3), rng.uniform(0.1, 10), rng.uniform(-1, 1)])
        rhs = css.rhs_scaled3(np.log(eps), Y)
        args = [_fq(eps)] + [_fq(z) for z in Y]
        for r in range(3):
            lhs = sum(float(C.P[r][i](*args)) * rhs[i] for i in range(3))
            q = float(C.Q[r](*args))
            assert abs(lhs - q) < 1e-12 * (1 + abs(q))


def test_constraint_is_invariant_of_4d_flow():
    lam, Delta = sonic_constraint_propagation(sonic_system())      # raises if not exact
    assert lam.total_degree() == 10 and Delta.total_degree() == 7


def _hom_eval(poly, nums, den, ctx2):
    """p(numA/den, numN/den, numW/den, numV/den) * den^deg as a polynomial in (s, V0)."""
    deg = poly.total_degree()
    tot = ctx2.constant(0)
    for exps, coef in poly.terms():
        e = list(map(int, exps))[1:]                     # drop t
        term = ctx2.constant(0) + fmpq(coef)
        for k, ek in enumerate(e):
            term *= nums[k] ** ek
        term *= den ** (deg - sum(e))
        tot += term
    return tot


def test_sonic_closed_forms_exact_in_Q_V0_sqrt3():
    """Delta~(u0) = 0, l.Q_0 = 0 (KHA99 214-215) and C~(u0) = 0 (211) as identities in Q[V0, sqrt3]."""
    ctx2 = fmpq_mpoly_ctx.get(("s", "V0"))
    s, V0 = ctx2.gens()
    den = 8 * (1 - V0**2) * (1 - s * V0)
    numA = 2 * (7 + 2 * s * V0 - 3 * V0**2) * (1 - s * V0)
    numN = 8 * (1 - V0**2) * (s - V0)
    numW = (3 - 2 * s * V0 - 3 * V0**2) * (1 - s * V0)      # W0 = (3 - 2 sqrt3 V0 - 3 V0^2)/(8(1-V0^2))
    numV = V0 * den
    nums = [numA, numN, numW, numV]
    S = sonic_system()
    P3W, P3V, P4W, P4V = S.P[2][2], S.P[2][3], S.P[3][2], S.P[3][3]
    mod = s**2 - 3
    for poly in (P3W * P4V - P3V * P4W,                     # Delta~
                 P4W * S.Q[2] - P3W * S.Q[3],               # l . Q_0 with l = (P4W, -P3W)
                 sonic_constraint_poly(S.ctx)):
        h = _hom_eval(poly, nums, den, ctx2)
        _, r = divmod(h, mod)
        assert r.is_zero()


def test_level_matrix_formula_matches_finite_differences():
    ex = sonic.sonic_expansion("0.112439401388092", K=8)
    with precision(256):
        for n in (2, 5, 8):
            Mx, _ = recursion.extract_level_matrix(ex.sys, ex.eqs, ex.balls, n)
            Mf = recursion.level_matrix([ex.D], [ex.E], n)
            assert all((Mx[0][i, j] - Mf[0][i, j]).contains(arb(0)) for i in range(4) for j in range(4))
