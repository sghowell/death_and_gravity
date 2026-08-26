"""Variational augmentation of a polynomial system and derivative Taylor models.

For P(t,u) theta u = Q(t,u) and a scalar parameter lambda entering through the
initial/low-order data, y := du/dlambda solves the linearised system.  Together
they form the *augmented* polynomial system in z = (u, y):
    P~ = [[P, 0], [sum_l y_l dP/du_l, P]],   Q~ = (Q, DQ . y),
(from d/dlambda of P theta u = Q), which is again of the form P~ theta z = Q~ with
exact fmpq_mpoly entries, so ``recursion`` and ``tailbound`` apply verbatim.
This is what makes the *tail* of du/dlambda certifiable (the plain tail
certificates of A1/A2 bound |u_n| but say nothing about |du_n/dlambda|).

Derivative Taylor models.  If u_n(lambda) is known as a degree-m Taylor
polynomial at the centre c (point run) and the (m+1)-th coefficient is
enclosed over the interval X = [c-w, c+w] (interval-base-point run, A1 sec. 3),
then for every lambda in X
    u_n'(lambda) = sum_{k=1}^{m} k u_{n,k} (lambda-c)^{k-1} + (m+1) u_{n,m+1}(xi) (lambda-c)^m,
so the derivative is enclosed by the shifted polynomial on [-w, w] plus
(m+1) sup_X |u_{n,m+1}| w^m  (Lagrange remainder of u' at degree m-1).
"""
from __future__ import annotations

from flint import arb

from . import recursion
from .arbseries import Series, to_arb
from .polysys import PolySystem, make_ctx
from .shootsys import remap


def augment(sys, names=None, skip=()):
    """Augmented PolySystem in (t, u_1..u_d, y_j for j not in ``skip``).

    Components in ``skip`` have y_j = du_j/dlambda identically 0 (e.g. T = e^x, which
    does not depend on the shooting parameter): their y-equations are dropped and
    y_j is set to 0 in the coupling terms.  Unknown order: (u_1..u_d, y_j, j kept)."""
    d = sys.d
    keep = [j for j in range(d) if j not in skip]
    dy = len(keep)
    old = tuple(sys.ctx.names())
    names = tuple(old) + tuple(f"d{old[1 + j]}" for j in keep) if names is None else names
    ctx = make_ctx(names)
    gmap = {k: k for k in range(d + 1)}
    gens = ctx.gens()
    ycol = {j: gens[d + 1 + k] for k, j in enumerate(keep)}          # y_j generator
    z = ctx.constant(0)
    dP, dQ = sys.dP(), sys.dQ()
    n = d + dy
    P = [[z] * n for _ in range(n)]
    Q = [z] * n
    for r in range(d):
        Q[r] = remap(sys.Q[r], ctx, gmap)
        for i in range(d):
            if not sys.P[r][i].is_zero():
                P[r][i] = remap(sys.P[r][i], ctx, gmap)
    for k, r in enumerate(keep):
        acc = z
        for l in keep:
            if not dQ[r][l].is_zero():
                acc = acc + remap(dQ[r][l], ctx, gmap) * ycol[l]
        Q[d + k] = acc
        for i in range(d):
            if sys.P[r][i].is_zero():
                continue
            if i in ycol:
                P[d + k][d + keep.index(i)] = P[r][i]
            acc = z
            for l in keep:
                if not dP[r][i][l].is_zero():
                    acc = acc + remap(dP[r][i][l], ctx, gmap) * ycol[l]
            P[d + k][i] = acc
    return PolySystem(ctx, P, Q, euler=sys.euler)


def augmented_level_equations(eqs, d, extra=None):
    """Level equations for the augmented system: the u-rows unchanged, the y-rows the same
    combinations of residual rows shifted by d, plus ``extra[i]`` (list of (coef, row, shift))
    appended to y-row i (used for the d/dlambda of parameter-dependent combinations)."""
    rows = [list(r) for r in eqs.rows]
    for i, row in enumerate(eqs.rows):
        new = [(c, r + d, s) for (c, r, s) in row]
        if extra and i in extra:
            new += list(extra[i])
        rows.append(new)
    return recursion.LevelEquations(rows)


# ----------------------------------------------------------------------------
# derivative Taylor models from (point run, interval run) delta-Series data
# ----------------------------------------------------------------------------
def poly_part(c, m):
    """Degree-<=m polynomial part of a delta-Series (or arb) as a Series with cap m+1."""
    if isinstance(c, Series):
        return Series(c.coeffs(m + 1), m + 1)
    return Series([c] + [arb(0)] * m, m + 1)


def derivative_poly(c, m):
    """Polynomial part of d/ddelta of the delta-Series c, degrees 0..m-1 (cap m)."""
    co = c.coeffs(m + 1)
    return Series([co[k] * k for k in range(1, m + 1)], m)


def derivative_balls(point_coefs, interval_coefs, m, w):
    """Enclosures of du_n/dlambda over |lambda-c| <= w for all n, i (list over n of d-lists).

    ``point_coefs[n][i]`` delta-Series at the centre (degrees >= m needed),
    ``interval_coefs[n][i]`` delta-Series at the interval base point (degree m+1 needed)."""
    dw = arb(0, w)
    wpow = to_arb(w) ** m
    out = []
    for n in range(len(point_coefs)):
        row = []
        for i in range(len(point_coefs[n])):
            top = interval_coefs[n][i][m + 1].abs_upper() * wpow * (m + 1)
            row.append(derivative_poly(point_coefs[n][i], m)(dw) + arb(0, top))
        out.append(row)
    return out


def derivative_point(point_coefs):
    """du_n/dlambda at the centre (the delta^1 coefficients)."""
    return [[c[1] for c in cn] for cn in point_coefs]


def augmented_coefs(balls_u, balls_y):
    return [list(bu) + list(by) for bu, by in zip(balls_u, balls_y)]


def check_level_residuals(sys_aug, eqs_aug, coefs_aug, n_from):
    """Orders n_from..K whose augmented level residuals do not contain 0 (must be empty)."""
    res = recursion.level_residuals(sys_aug, eqs_aug, coefs_aug, n_from=n_from)
    return [n_from + j for j, r in enumerate(res) if not all(recursion.contains_zero(v) for v in r)]


def certify_augmented(sys_aug, eqs_aug, coefs_aug, nu=None, **kw):
    """(D~, E~, certificate) for the augmented truncation ``coefs_aug`` (balls over the interval)."""
    from .tailbound import certify_tail
    D, E = recursion.structure_matrices(sys_aug, eqs_aug, coefs_aug, dm=1)
    cert = certify_tail(sys_aug, eqs_aug, coefs_aug, D, E, nu=nu, **kw)
    return D, E, cert
