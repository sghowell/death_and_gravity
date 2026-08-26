"""S2 core tests: series arithmetic, exact polynomial systems, exact algebraic identities."""
from fractions import Fraction

import numpy as np
import pytest
from flint import arb, fmpq, fmpq_mpoly_ctx
from p4 import css
from p4.validated import recursion, sonic
from p4.validated.arbseries import Series, precision
from p4.validated.biseries import BiSeries
from p4.validated.shootsys import remap
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


def _row(sys, r, ctx, dgens):
    """Residual row  sum_i P[r][i] dU_i - Q[r]  of a PolySystem as a polynomial in ``ctx``
    (whose first generators are those of the system, followed by the derivative generators ``dgens``)."""
    g = {k: k for k in range(sys.d + 1)}
    row = -remap(sys.Q[r], ctx, g)
    for i in range(sys.d):
        if not sys.P[r][i].is_zero():
            row += remap(sys.P[r][i], ctx, g) * dgens[i]
    return row


def test_centre_system_is_the_reduced_sonic_system_exactly():
    """Exact fmpq_mpoly identity: the centre rows (N, 3, 4) are the sonic rows 2-4 with the momentum
    constraint A = 1 + 2WT/S (T = 1 + V^2/3 + (4/3)NV, S = 1 - V^2) substituted and the change of
    variables (N, W, V) = (n/t, w t^2, v t), d/dx = theta, i.e. (N', W', V') = ((n' - n)/t,
    t^2(w' + 2w), t(v' + v)); after clearing the A-denominator by S and the negative powers of t by t
    the quotients are the units 1, t^3 S, t^2 S (remainders exactly zero)."""
    Ssys, Csys = sonic_system(), centre_system()
    cs = fmpq_mpoly_ctx.get(("t", "A", "N", "W", "V", "dA", "dN", "dW", "dV"))
    cc = fmpq_mpoly_ctx.get(("t", "n", "w", "v", "dn", "dw", "dv"))
    t, n, w, v, dn, dw, dv = cc.gens()
    S = 1 - v**2 * t**2
    Tc = 1 + v**2 * t**2 * fmpq(1, 3) + fmpq(4, 3) * n * v                # T after the substitution
    # image of each sonic generator as (numerator, t-weight, S-weight): X = num / (t^a S^b)
    img = {0: (t, 0, 0), 1: (S + 2 * w * t**2 * Tc, 0, 1), 2: (n, 1, 0), 3: (w * t**2, 0, 0),
           4: (v * t, 0, 0), 5: None, 6: (dn - n, 1, 0), 7: (t**2 * (dw + 2 * w), 0, 0), 8: (t * (dv + v), 0, 0)}

    def substitute(poly):
        """t^Amax S^Bmax * poly(images), a polynomial in cc (Amax, Bmax the maximal weights)."""
        terms = [(list(map(int, e)), fmpq(c)) for e, c in poly.terms()]
        assert all(e[0] == 0 and e[5] == 0 for e, _ in terms)              # autonomous; no A' in rows 2-4
        wts = [(sum(e[k] * img[k][1] for k in range(9) if e[k]), sum(e[k] * img[k][2] for k in range(9) if e[k]))
               for e, _ in terms]
        Amax, Bmax = max(a for a, _ in wts), max(b for _, b in wts)
        tot = cc.constant(0)
        for (e, c), (a, b) in zip(terms, wts):
            term = cc.constant(0) + c
            for k in range(9):
                if e[k]:
                    term *= img[k][0] ** e[k]
            tot += term * t ** (Amax - a) * S ** (Bmax - b)
        return tot

    units = [cc.constant(1), t**3 * S, t**2 * S]
    for r in (1, 2, 3):
        sub = substitute(_row(Ssys, r, cs, cs.gens()[5:]))
        crow = _row(Csys, r - 1, cc, (dn, dw, dv))
        q, rem = divmod(sub, crow)
        assert rem.is_zero() and q == units[r - 1], (r, q, rem)


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
