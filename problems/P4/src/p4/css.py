"""Koike-Hara-Adachi CSS Einstein-Euler system at k = 1/3 (gamma = 4/3).

Polar-areal metric ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2, fluid
velocity V relative to the static observer, similarity variables
    s = -ln(-t),  x = ln(-r/t),  N = alpha/(a e^x),  A = a^2,  W = 4 pi r^2 a^2 rho.
Equations (KHA95 eq. 18, re-derived from scratch with sympy in S1; the digest's
row 3 has a transcription slip -- its last term must read 2NV(1+4W/(9(1-V^2)))):
    A_x/A = F_A := 1 - A + 2W(1+V^2/3)/(1-V^2)                      (Hamiltonian)
    N_x/N = F_N := -2 + A - 2W/3                                     (slicing)
    (A_s + A_x)/A = G := -(8/3) N V W/(1-V^2)                        (momentum)
    row3: W_s/W + 4V V_s/(3(1-V^2)) + a W_x + b V_x + e = 0
    row4: 4V W_s/W + 4(1+V^2) V_s/(1-V^2) + c W_x + d V_x + f = 0
with a, b, c, d, e, f given in ``coeffs``.  CSS: drop all s-derivatives.
The momentum constraint F_A = G (KHA99 eq. 211) is then an invariant of the flow.

Sonic point at x = 0 (gauge): det [[a,b],[c,d]] = 0, i.e. the fluid speed relative
to the x = const line, v_rel = (1+NV)/(N+V), equals 1/sqrt(3).
"""
from __future__ import annotations

import numpy as np

SQ3 = np.sqrt(3.0)
K_EOS = 1.0 / 3.0


# ----------------------------------------------------------------------------
# field equations (generic in the number type)
# ----------------------------------------------------------------------------
def coeffs(A, N, W, V, Am1=None):
    """Coefficient functions of the KHA system at gamma = 4/3.

    ``Am1`` = A - 1 may be supplied separately to avoid cancellation near the
    regular centre (A -> 1).  Returns a dict; all entries are expressions in
    (A, N, W, V) only.
    """
    if Am1 is None:
        Am1 = A - 1
    V2 = V * V
    omV2 = 1 - V2
    FA = -Am1 + 2 * W * (1 + V2 / 3) / omV2          # A_x / A
    FN = -1 + Am1 - 2 * W / 3                        # N_x / N   (= -2 + A - 2W/3)
    G = -(8.0 / 3.0) * N * V * W / omV2              # (A_s + A_x)/A
    a = (1 + N * V) / W
    b = 4 * (N + V) / (3 * omV2)
    c = (4 * V + N + 3 * N * V2) / W
    d = 4 * (1 + V2 + 2 * N * V) / omV2
    e = -N * V * FA / 3 + (4.0 / 3.0) * V * N * FN + 2 * N * V * (1 + 4 * W / (9 * omV2))
    f = N * omV2 * FA + 4 * (1 + V2) * N * FN + 2 * N * (1 + 3 * V2)
    return dict(FA=FA, FN=FN, G=G, a=a, b=b, c=c, d=d, e=e, f=f,
                sa=1 / W, sb=4 * V / (3 * omV2), sc=4 * V / W, sd=4 * (1 + V2) / omV2)


def fluid_residuals(A, N, W, V, Wx, Vx, Ws=0, Vs=0, Am1=None):
    """(row3, row4) of the time-dependent fluid equations."""
    q = coeffs(A, N, W, V, Am1)
    r3 = q["sa"] * Ws + q["sb"] * Vs + q["a"] * Wx + q["b"] * Vx + q["e"]
    r4 = q["sc"] * Ws + q["sd"] * Vs + q["c"] * Wx + q["d"] * Vx + q["f"]
    return r3, r4


def det_and_numerators(q):
    """Delta = ad - bc and the numerators P, Q with W_x = P/Delta, V_x = Q/Delta."""
    Delta = q["a"] * q["d"] - q["b"] * q["c"]
    P = q["b"] * q["f"] - q["d"] * q["e"]
    Q = q["c"] * q["e"] - q["a"] * q["f"]
    return Delta, P, Q


def constraint(y):
    """Momentum constraint F_A - G (vanishes on CSS solutions)."""
    A, N, W, V = y
    q = coeffs(A, N, W, V)
    return q["FA"] - q["G"]


def v_rel(N, V):
    """Fluid velocity relative to the x = const line (HM01's V_z up to sign)."""
    return (1 + N * V) / (N + V)


# ----------------------------------------------------------------------------
# ODE right-hand sides
# ----------------------------------------------------------------------------
def rhs_plain(x, y):
    """CSS ODE in the plain variables y = (A, N, W, V); autonomous in x."""
    A, N, W, V = y
    q = coeffs(A, N, W, V)
    Delta, P, Q = det_and_numerators(q)
    return np.array([A * q["FA"], N * q["FN"], P / Delta, Q / Delta])


def to_scaled(x, y):
    """(A, N, W, V) -> (Ah, Nh, Wh, Vh) with A = 1 + Ah e^{2x}, N = Nh e^{-x}, W = Wh e^{2x}, V = Vh e^x."""
    A, N, W, V = y
    return np.array([(A - 1) * np.exp(-2 * x), N * np.exp(x), W * np.exp(-2 * x), V * np.exp(-x)])


def from_scaled(x, yh):
    Ah, Nh, Wh, Vh = yh
    return np.array([1 + Ah * np.exp(2 * x), Nh * np.exp(-x), Wh * np.exp(2 * x), Vh * np.exp(x)])


def rhs_scaled(x, yh):
    """CSS ODE in the centre-scaled variables (all O(1) as x -> -inf on regular solutions)."""
    Ah, Nh, Wh, Vh = yh
    e_x, e_2x = np.exp(x), np.exp(2 * x)
    Am1 = Ah * e_2x
    A, N, W, V = 1 + Am1, Nh / e_x, Wh * e_2x, Vh * e_x
    q = coeffs(A, N, W, V, Am1=Am1)
    Delta, P, Q = det_and_numerators(q)
    # e^{-2x} A F_A written without cancellation: F_A = e^{2x} (-Ah + 2 Wh (1+V^2/3)/(1-V^2))
    dAh = A * (-Ah + 2 * Wh * (1 + V * V / 3) / (1 - V * V)) - 2 * Ah
    dNh = Nh * (Am1 - 2 * W / 3)                       # Nh (F_N + 1)
    dWh = (P / Delta) / e_2x - 2 * Wh
    dVh = (Q / Delta) / e_x - Vh
    return np.array([dAh, dNh, dWh, dVh])


def mass_function(x, yh):
    """m~ = (1 - 1/A) e^x = 2m/(-t) up to normalisation; -> 0 at a regular centre."""
    Ah = yh[0]
    Am1 = Ah * np.exp(2 * x)
    return Am1 / (1 + Am1) * np.exp(x)


# ----------------------------------------------------------------------------
# constraint-reduced (3D) system: A eliminated by the momentum constraint F_A = G
#   A - 1 = 2 W (1 + V^2/3 + (4/3) N V) / (1 - V^2).
# The constraint surface is invariant but repelling for the backward-x flow of the
# 4D system (the centre fixed point has one unstable direction *transverse* to it),
# so shooting toward the centre must be done on the reduced system.
# ----------------------------------------------------------------------------
def Am1_constraint(N, W, V):
    return 2 * W * (1 + V * V / 3 + (4.0 / 3.0) * N * V) / (1 - V * V)


def rhs_plain3(x, y3):
    """Reduced CSS ODE, y3 = (N, W, V)."""
    N, W, V = y3
    Am1 = Am1_constraint(N, W, V)
    q = coeffs(1 + Am1, N, W, V, Am1=Am1)
    Delta, P, Q = det_and_numerators(q)
    return np.array([N * q["FN"], P / Delta, Q / Delta])


def rhs_scaled3(x, yh3):
    """Reduced CSS ODE in centre-scaled variables yh3 = (Nh, Wh, Vh)."""
    Nh, Wh, Vh = yh3
    e_x, e_2x = np.exp(x), np.exp(2 * x)
    N, W, V = Nh / e_x, Wh * e_2x, Vh * e_x
    Am1 = Am1_constraint(N, W, V)
    q = coeffs(1 + Am1, N, W, V, Am1=Am1)
    Delta, P, Q = det_and_numerators(q)
    return np.array([Nh * (Am1 - 2 * W / 3), (P / Delta) / e_2x - 2 * Wh, (Q / Delta) / e_x - Vh])


def scaled3_to_full(x, yh3):
    """(Nh, Wh, Vh) -> scaled 4-vector (Ah, Nh, Wh, Vh) with Ah from the constraint."""
    Nh, Wh, Vh = yh3
    e_x = np.exp(x)
    Am1 = Am1_constraint(Nh / e_x, Wh * e_x**2, Vh * e_x)
    return np.array([Am1 / e_x**2, Nh, Wh, Vh])


def mass_function3(x, yh3):
    return mass_function(x, scaled3_to_full(x, yh3))


def centre_check(x, yh):
    """Regular-centre diagnostics (KHA99 eqs. 220-222): all -> 0 as x -> -inf."""
    Ah, Nh, Wh, Vh = yh
    return dict(mass=mass_function(x, yh), Ah_minus_2Wh3=Ah - 2 * Wh / 3, NV_plus_half=Nh * Vh + 0.5)


# ----------------------------------------------------------------------------
# sonic point (x = 0): zeroth- and first-order data
# ----------------------------------------------------------------------------
def sonic_data(V0):
    """Zeroth-order sonic data (KHA99 sec. IV.3.2 at gamma = 4/3), as (A0, N0, W0, V0)."""
    N0 = (SQ3 - V0) / (1 - SQ3 * V0)
    A0 = (7 + 2 * SQ3 * V0 - 3 * V0**2) / (4 * (1 - V0**2))
    W0 = 0.375 * (1 - 2 * V0 / SQ3 - V0**2) / (1 - V0**2)
    return np.array([A0, N0, W0, V0])


def first_order_quadratic(V0):
    """Coefficients (c2, c1, c0) of the quadratic c2 V1^2 + c1 V1 + c0 = 0 for the
    first-order sonic coefficient V1 (derived with sympy from the series equations)."""
    c2 = 12 * SQ3 * (1 - V0**2) ** 2
    c1 = 12 * V0 * (V0**2 - 1) * (3 * V0**3 - 2 * SQ3 * V0**2 - 5 * V0 + 2 * SQ3)
    c0 = V0 * (9 * V0**6 - 9 * SQ3 * V0**5 - 21 * V0**4 + 10 * SQ3 * V0**3
               + 39 * V0**2 + 3 * SQ3 * V0 - 27)
    return c2, c1, c0


def first_order_discriminant(V0):
    """Discriminant of the V1 quadratic: real branches iff > 0.
    = 432 V0 (1-V0^2)^4 (3V0^3 - 5 sqrt3 V0^2 + 3V0 + 3 sqrt3)."""
    return 432 * V0 * (1 - V0**2) ** 4 * (3 * V0**3 - 5 * SQ3 * V0**2 + 3 * V0 + 3 * SQ3)


def sonic_branches(V0, complex_ok=False):
    """First-order sonic data.  Returns (y0, [y1_branch0, y1_branch1]) with
    y1 = (A1, N1, W1, V1); branches sorted by V1 (ascending real part)."""
    y0 = sonic_data(V0)
    A0, N0, W0, _ = y0
    q = coeffs(A0, N0, W0, V0)
    A1, N1 = A0 * q["FA"], N0 * q["FN"]
    c2, c1, c0 = first_order_quadratic(V0)
    disc = c1 * c1 - 4 * c2 * c0
    if disc < 0 and not complex_ok:
        raise ValueError(f"no real analytic branch at V0={V0}: discriminant {disc}")
    sq = np.sqrt(disc + 0j) if disc < 0 else np.sqrt(disc)
    V1s = [(-c1 - sq) / (2 * c2), (-c1 + sq) / (2 * c2)]
    out = []
    for V1 in V1s:
        W1 = -(q["e"] + q["b"] * V1) / q["a"]        # rank-1 zeroth-order relation
        out.append(np.array([A1, N1, W1, V1]))
    return y0, out


def sonic_eigenvalues(V0):
    """The nonzero eigenvalues mu = grad(Delta).y1 of the desingularised flow at the
    sonic point for the two branches (saddle iff mu_0 mu_1 < 0)."""
    from .tps import TPS
    y0, branches = sonic_branches(V0, complex_ok=True)
    mus = []
    for y1 in branches:
        ser = [TPS(np.array([y0[i], y1[i]], dtype=complex)) for i in range(4)]
        Delta, _, _ = det_and_numerators(coeffs(*ser))
        mus.append(Delta.c[1])
    return mus
