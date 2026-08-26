"""The two CSS systems in cleared-denominator polynomial form (exact fmpq_mpoly).

Sonic point (t = x, theta = d/dx, unknowns u = (A, N, W, V), S = 1 - V^2,
Phi_A = S F_A = (1-A) S + 2 W (1 + V^2/3), F_N = -2 + A - 2W/3):
    row 1:  S A'                                   = A Phi_A
    row 2:  N'                                     = N F_N
    row 3:  3(1+NV) S W' + 4 W (N+V) V'             = -3 W S e
    row 4:  (4V + N + 3NV^2) S W' + 4W(1+V^2+2NV) V' = -W S f
with e, f the KHA row-3/row-4 source terms of ``p4.css.coeffs`` (S1, corrected
transcription).  Rows 3-4 are 3WS x (row 3) and WS x (row 4) of ``p4.css``.

Regular centre (t = eps = e^x, theta = t d/dt = d/dx, unknowns
u = (n, w, v) = (N e^x, W e^{-2x}, V e^{-x}), constraint-reduced,
A - 1 = 2 W T / S with T = 1 + V^2/3 + (4/3) N V,  S = 1 - v^2 t^2):
    row N:  S n'                                                    = (4/3) n w t^2 (1 + 2nv + v^2 t^2)
    row 3:  3 S (1+nv) w' + 4 w (n + v t^2) v'                       = Q3
    row 4:  (n + 4 v t^2 + 3 n v^2 t^2) S w' + 4 w t^2 (1 + v^2 t^2 + 2nv) v' = Q4
obtained from rows 3-4 above by the substitution (multiplied by 3wS resp. t w S)
with F_A replaced by G (momentum constraint) -- see notes/s2-validated-expansions.md.
"""
from __future__ import annotations

from flint import fmpq

from .polysys import PolySystem, make_ctx


def F(n, d=1):
    return fmpq(n, d)


def sonic_system():
    ctx = make_ctx(("t", "A", "N", "W", "V"))
    t, A, N, W, V = ctx.gens()
    z = ctx.constant(0)
    S = 1 - V**2
    PhiA = (1 - A) * S + 2 * W * (1 + V**2 * F(1, 3))       # S * F_A
    FN = -2 + A - W * F(2, 3)
    P = [[S, z, z, z],
         [z, 1 + z, z, z],
         [z, z, 3 * (1 + N * V) * S, 4 * W * (N + V)],
         [z, z, (4 * V + N + 3 * N * V**2) * S, 4 * W * (1 + V**2 + 2 * N * V)]]
    Q = [A * PhiA,
         N * FN,
         N * V * W * PhiA - 4 * N * V * W * S * FN - 6 * N * V * W * S - F(8, 3) * N * V * W**2,
         -W * N * (S * PhiA + 4 * (1 + V**2) * S * FN + 2 * (1 + 3 * V**2) * S)]
    return PolySystem(ctx, P, Q, euler=False)


def sonic_constraint_poly(ctx):
    """C~ = (A-1) S - 2 W (1 + V^2/3 + (4/3) N V); the momentum constraint is C~ = 0."""
    t, A, N, W, V = ctx.gens()
    return (A - 1) * (1 - V**2) - 2 * W * (1 + V**2 * F(1, 3) + F(4, 3) * N * V)


def sonic_constraint_propagation(sys):
    """Exact identity  S Delta~ dC~/dx = Lambda C~  along the 4D CSS flow.

    Delta~ = det of the fluid principal block of P, W' = PW/Delta~, V' = PV/Delta~ (Cramer),
    A' = Q1/S, N' = Q2.  Returns (Lambda, Delta~); raises if the division is not exact."""
    t, A, N, W, V = sys.ctx.gens()
    C = sonic_constraint_poly(sys.ctx)
    S = 1 - V**2
    P3W, P3V, P4W, P4V = sys.P[2][2], sys.P[2][3], sys.P[3][2], sys.P[3][3]
    Q1, Q2, Q3, Q4 = sys.Q
    Delta = P3W * P4V - P3V * P4W
    PW = Q3 * P4V - Q4 * P3V
    PV = P3W * Q4 - P4W * Q3
    num = (C.derivative(1) * Q1 * Delta + C.derivative(2) * Q2 * S * Delta
           + (C.derivative(3) * PW + C.derivative(4) * PV) * S)
    lam, rem = divmod(num, C)
    if not rem.is_zero():
        raise ArithmeticError("constraint is not an invariant of the 4D flow (remainder != 0)")
    return lam, Delta


def centre_system():
    ctx = make_ctx(("t", "n", "w", "v"))
    t, n, w, v = ctx.gens()
    z = ctx.constant(0)
    t2 = t**2
    S = 1 - v**2 * t2
    P = [[S, z, z],
         [z, 3 * S * (1 + n * v), 4 * w * (n + v * t2)],
         [z, (n + 4 * v * t2 + 3 * n * v**2 * t2) * S, 4 * w * t2 * (1 + v**2 * t2 + 2 * n * v)]]
    QN = F(4, 3) * n * w * t2 * (1 + 2 * n * v + v**2 * t2)
    Q3 = -(6 * w * S * (1 + n * v) + 4 * w * v * (n + v * t2)
           + n * v * w * (2 * S + F(8, 3) * w * t2 * (3 + 5 * n * v + 2 * v**2 * t2)))
    Q4 = -(2 * w * (n + 4 * v * t2 + 3 * n * v**2 * t2) * S
           + 4 * w * v * t2 * (1 + v**2 * t2 + 2 * n * v)
           + n * w * (S * (-2 + 2 * v**2 * t2 - F(8, 3) * n * v * w * t2)
                      + F(16, 3) * w * t2 * (1 + v**2 * t2) * (1 + 2 * n * v + v**2 * t2)))
    return PolySystem(ctx, P, [QN, Q3, Q4], euler=True)


def centre_Ahat_fraction(ctx):
    """A^ = (A-1) e^{-2x} = num/den with num = 2 w T, den = S (polynomials in (t, n, w, v))."""
    t, n, w, v = ctx.gens()
    T = 1 + v**2 * t**2 * F(1, 3) + F(4, 3) * n * v
    return 2 * w * T, 1 - v**2 * t**2
