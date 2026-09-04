"""Exact nonlinear lapse removal for the D-only M0 witness.

Dimensionless time u=t/tau. The Hamiltonian is divided by sqrt(hat_h)
and M^2/tau^2. sigma is the perturbation of the dimensionless trace of
the ADM momentum; shear2 is its traceless square, NOT a scalar-only ansatz.
rho=tau^2 R[hat_h]. Spatial momentum constraints are still imposed by shift.
"""

from functools import cache

import sympy as sp
from p8 import jets as j
from p8.matter import ia_completion
from p8.rational_candidates import X, reconstruct

u = sp.Symbol("u", real=True)
N = sp.Symbol("N", positive=True)
sigma, rho, shear2 = sp.symbols("sigma rho shear2", real=True)
d = 1 + u**2
H = 2*u/d


@cache
def functions():
    old = reconstruct("D")
    F = old["F"].subs(j.t, u)
    a3 = 4/d**3
    _, a4, a5 = ia_completion(-sp.Rational(1, 2), 0, 0, a3, X)
    omega = (N**-4 - 1)/(2*d**3)
    return {"F": F, "a3": a3, "a4": a4, "a5": a5,
            "omega": omega, "omega_u": sp.diff(omega, u),
            "Theta": old["Theta"].subs(j.t, u),
            "Lambda": old["Lambda"].subs(j.t, u),
            "delta": old["delta"].subs(j.t, u),
            "J": sp.factor((old["Sigma_total"] + 3*old["Theta"]**2).subs(j.t, u))}


@cache
def transformation_checks():
    """Full kinetic square and spatial-conformal IBP, not a quadratic jet."""
    f = functions()
    omega = f["omega"]
    c = f["a3"]/N**3
    dv = -3*c**2/4
    e = (-X**2*f["a4"]).subs(X, N**-2)
    v, k, kk, z, w = sp.symbols("V K KK z W", real=True)
    # hat K has trace k and squared norm kk. K_old = hat K+(c V/2+W) I.
    shift = c*v/2 + w
    trace_old = k + 3*shift
    square_old = kk + 2*k*shift + 3*shift**2
    before = (square_old-trace_old**2)/2 + c*trace_old*v + dv*v**2
    after = (kk-k**2)/2 - 2*w*k - 3*w**2
    delta_N = N*sp.diff(omega, N)
    # R[e^(2 omega)hat_h] = e^(-2 omega)(Rhat-4 Lap omega-2 |grad omega|^2).
    # IBP of -2 N exp(omega) Lap omega adds 2 N exp(omega)(1+delta_N)delta_N z^2.
    spatial = -delta_N**2*z**2 + e*z**2 + 2*(1+delta_N)*delta_N*z**2
    return {
        "Ia_A4": sp.cancel(f["a4"] - (-f["a3"] + X**2*f["a3"]**2/4)),
        "Ia_A5": sp.cancel(f["a5"] + X*f["a3"]**2),
        "lapse_velocity_cancellation": sp.cancel(sp.diff(omega, N) + c/(2*N**2)),
        "kinetic_square": sp.cancel(before-after),
        "spatial_gradient_IBP": sp.cancel(spatial),
        "background_identity": omega.subs(N, 1),
        "linear_chart_dictionary": sp.cancel(sp.diff(omega, N).subs(N, 1)-f["delta"]),
    }


@cache
def legendre_checks():
    """Independent six-component symmetric-matrix Legendre transform."""
    A, B, W, FF = sp.symbols("A B W FF", nonzero=True)
    ks = sp.symbols("kxx kyy kzz kxy kxz kyz", real=True)
    ps = sp.symbols("pxx pyy pzz pxy pxz pyz", real=True)
    km = sp.Matrix([[ks[0], ks[3], ks[4]], [ks[3], ks[1], ks[5]], [ks[4], ks[5], ks[2]]])
    pm = sp.Matrix([[ps[0], ps[3], ps[4]], [ps[3], ps[1], ps[5]], [ps[4], ps[5], ps[2]]])
    lag = A*((sp.trace(km*km)-sp.trace(km)**2)/2 - 2*W*sp.trace(km) - 3*W**2) + FF + B*rho/2
    # p=1/2 dL/dK for diagonal components, 1/4 for off-diagonal coordinates.
    p_from_k = [sp.diff(lag, key)/(2 if i < 3 else 4) for i, key in enumerate(ks)]
    sol = sp.solve([p_from_k[i]-ps[i] for i in range(6)], ks)
    ham = (2*sp.trace(pm*km)-lag).subs(sol)
    expected = 2*(sp.trace(pm*pm)-sp.trace(pm)**2/2)/A - 2*W*sp.trace(pm) - FF - B*rho/2
    return {"six_component_Legendre": sp.cancel(ham-expected)}


@cache
def hamiltonian():
    f = functions()
    omega, wu = f["omega"], f["omega_u"]
    # p^i_j=-H delta^i_j+delta p^i_j, with trace(delta p)=sigma.
    Q = -3*H**2/2 + H*sigma - sigma**2/6 + shear2
    return (2*N*sp.exp(-3*omega)*Q - 2*wu*(-3*H+sigma)
            - N*sp.exp(3*omega)*f["F"].subs(X, N**-2)
            - N*sp.exp(omega)*rho/2)


@cache
def lapse_jets():
    """Exact N derivatives at N=1, through fourth order, factored by invariants."""
    h = hamiltonian()
    pieces = {"A": h.subs({sigma: 0, rho: 0, shear2: 0}),
              "B": sp.diff(h, sigma).subs(sigma, 0),
              "C": sp.diff(h, rho),
              "D": sp.diff(h, sigma, 2)/2,
              "E": sp.diff(h, shear2)}
    out = {}
    for key, expr in pieces.items():
        out[key] = tuple(sp.factor(sp.diff(expr, N, k).subs(N, 1)) for k in range(5))
    return out


@cache
def audit_checks():
    f, jets = functions(), lapse_jets()
    # These three bridges come from the nonlinear Hamiltonian, not gamma.py.
    return {**transformation_checks(), **legendre_checks(),
            "background_lapse": jets["A"][1],
            "nonlinear_to_J": sp.cancel(jets["A"][2] + 2*f["J"]),
            "nonlinear_to_Theta": sp.cancel(jets["B"][1] - 2*f["Theta"]),
            "nonlinear_to_Lambda": sp.cancel(jets["C"][1] + f["Lambda"]/2)}
