"""Nonlinear CD/M1 spatial chart and invariant canonical Hamiltonian.

The chart leaves lapse, shift and the physical matter metric unchanged as
physical objects: h_ij = T**(-1/2) hat_h_ij is a change of variables, not a
claim of minimal matter coupling to hat_h. The lapse boundary primitive is
local in spacetime. Only its finite N jets are needed at this checkpoint.
Spatial momentum constraints have NOT yet been eliminated.
"""

import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8.matter import ia_completion

from . import njets as nj

P8 = Path(__file__).resolve().parents[4]
WITNESS = P8 / "certificates" / "witness-CD_matter.json"
u = sp.Symbol("u", real=True)
N = sp.Symbol("N", positive=True)
X = sp.Symbol("X", positive=True)
sigma, rho, eta, shear2, z = sp.symbols("sigma rho eta shear2 z", real=True)
VARIABLES = (sigma, rho, eta, shear2, z)
WEIGHTS = (1, 1, 1, 2, 2)
d = 1 + u**2
H = 4*u/d
velocity = 1/(10*d**6)


@cache
def functions():
    """Read the pinned prior witness as input, not as a new reconstruction.

    verify.py checks both the prior source manifest and the witness digest.
    All new chart/Legendre/jet identities below are derived independently.
    """
    raw = json.loads(WITNESS.read_text())
    local = {"t": u, "X": X, "Matrix": sp.Matrix}
    result = {key: sp.sympify(value, locals=local)
              for key, value in raw["covariant_functions"].items()}
    result.update({key: sp.sympify(value, locals=local)
                   for key, value in raw["principal"].items()})
    result["delta"] = 1/(2*d**3)
    result["w"] = velocity*(3*result["delta"]-1)
    return result


@cache
def chart():
    f = functions()
    T = (-2*f["F2"]).subs(X, N**-2)
    omega = -sp.log(T)/4
    c = -1/(N*d**3)
    L0 = -sp.diff(T, u)/N
    W = -sp.diff(T, u)/(4*N*T)
    # Psi(u,N) = integral_1^N psi_N(u,n) dn, Psi(u,1)=0.
    psi_N = -3*sp.diff(T, u)/(2*N**4*d**3*T**sp.Rational(7, 4))
    return {"T": T, "omega": omega, "c": c, "L0": L0, "W": W,
            "psi_N": psi_N, "A": T**sp.Rational(1, 4)}


@cache
def transformation_checks():
    f, ch = functions(), chart()
    T, c, L0 = ch["T"], ch["c"], ch["L0"]
    _, a4, a5 = ia_completion(f["F2"], sp.diff(f["F2"], X), 0, f["A3"], X)
    dv = (X*(f["A3"]+a4+X*a5)).subs(X, N**-2)
    acc = (4*X*sp.diff(f["F2"], X)-X**2*a4).subs(X, N**-2)
    V, k, kk, W = sp.symbols("V k kk W", real=True)
    shift = c*V/(2*T)+W
    trace_old = k+3*shift
    square_old = kk+2*k*shift+3*shift**2
    before = T*(square_old-trace_old**2)/2+L0*trace_old+c*trace_old*V+dv*V**2
    after = T*(kk-k**2)/2+(-2*T*W+L0)*k-3*T*W**2+3*L0*W+3*L0*c*V/(2*T)
    delta_N = N*sp.diff(ch["omega"], N)
    # Curvature IBP differentiates N*e^omega*T, not only N*e^omega.
    spatial = acc+T*delta_N**2+2*(T+N*sp.diff(T, N))*delta_N
    return {
        "Ia_kinetic_degeneracy": sp.cancel(dv+3*c**2/(4*T)),
        "lapse_velocity_chart": sp.cancel(sp.diff(ch["omega"], N)+c/(2*T*N**2)),
        "kinetic_square_with_linear_V": sp.cancel(before-after),
        "spatial_gradient_IBP": sp.cancel(spatial),
        "boundary_primitive_derivative": sp.simplify(ch["psi_N"]+T**sp.Rational(-3, 4)*3*L0*c/(2*T*N**2)),
        "background_identity": sp.simplify(ch["omega"].subs(N, 1)),
        "linear_chart_dictionary": sp.cancel(sp.diff(ch["omega"], N).subs(N, 1)-f["delta"]),
        "reduced_linear_K_before_boundary": sp.simplify(T**sp.Rational(-3, 4)*(-2*T*ch["W"]+L0)
                                                         +sp.diff(T, u)/(2*N*T**sp.Rational(3, 4))),
        "reduced_potential_before_boundary": sp.cancel(-3*T*ch["W"]**2+3*L0*ch["W"]
                                                        -9*sp.diff(T, u)**2/(16*N**2*T)),
    }


@cache
def legendre_checks():
    """Vary all six symmetric metric velocities and the matter velocity."""
    A, P, U, B, E = sp.symbols("A P U B E", nonzero=True)
    ks = sp.symbols("kxx kyy kzz kxy kxz kyz", real=True)
    ps = sp.symbols("pxx pyy pzz pxy pxz pyz", real=True)
    v, r = sp.symbols("v r", real=True)
    km = sp.Matrix([[ks[0], ks[3], ks[4]], [ks[3], ks[1], ks[5]], [ks[4], ks[5], ks[2]]])
    pm = sp.Matrix([[ps[0], ps[3], ps[4]], [ps[3], ps[1], ps[5]], [ps[4], ps[5], ps[2]]])
    lag = N*(A*(sp.trace(km*km)-sp.trace(km)**2)/2+P*sp.trace(km)+U+B*rho/2)+E*v**2/(2*N)
    p_from_k = [sp.diff(lag, key)/(N*(2 if i < 3 else 4)) for i, key in enumerate(ks)]
    solution = sp.solve([p_from_k[i]-ps[i] for i in range(6)]+[sp.diff(lag, v)-r], (*ks, v))
    ham = (2*N*sp.trace(pm*km)+r*v-lag).subs(solution)
    expected = (2*N*(sp.trace(pm*pm)-sp.trace(pm)**2/2)/A+N*P*sp.trace(pm)/A
                -3*N*P**2/(4*A)-N*U-N*B*rho/2+N*r**2/(2*E))
    # The local divergence identity also includes the shift boundary.
    psi, psi_u, psi_N, Dn = sp.symbols("psi psi_u psi_N Dn")
    divergence = psi_u+psi_N*Dn+N*psi*sp.trace(km)
    boundary_remainder = psi_N*Dn-divergence
    return {"seven_velocity_Legendre": sp.cancel(ham-expected),
            "canonical_boundary_IBP": sp.expand(boundary_remainder+psi_u+N*psi*sp.trace(km))}


@cache
def generic_lapse_jets(*, omit_boundary=False):
    """Derive N jets symbolically before inserting large rational time functions.

    F=f0+f1*X+f2*X^2. R=1/d^3, Rd=dR/du, Rdd=d^2R/du^2.
    Derivatives, not Taylor coefficients, are returned. A fourth-order jet
    of the boundary primitive is exact for this purpose, not a truncation
    of any spacetime integral or a single Fourier mode.
    """
    R, Rd, Rdd, HH, l, f0, f1, f2 = sp.symbols("R Rd Rdd HH l f0 f1 f2", real=True)
    shift_X = nj.add(nj.Npower(-2), nj.constant(-1))
    T = nj.add(nj.constant(1), nj.scale(shift_X, R))
    Tu = nj.scale(shift_X, Rd)
    powT = lambda exponent: nj.power(T, sp.Rational(exponent))
    primitive_derivative = nj.scale(nj.mul(Tu, nj.Npower(-4), powT(sp.Rational(-7, 4))), -3*R/2)
    psi = (sp.Integer(0), *(sp.expand(primitive_derivative[k-1]/k) for k in range(1, 5)))
    if omit_boundary:
        psi = nj.constant(0)
    psi_u = tuple(sp.expand(sp.diff(value, R)*Rd+sp.diff(value, Rd)*Rdd) for value in psi)
    P = nj.add(nj.scale(nj.mul(Tu, nj.Npower(-1), powT(sp.Rational(-3, 4))), -sp.Rational(1, 2)),
               nj.scale(psi, -1))
    FF = nj.add(nj.constant(f0), nj.scale(nj.Npower(-2), f1), nj.scale(nj.Npower(-4), f2))
    potential = nj.add(FF, nj.scale(nj.mul(Tu, Tu, nj.Npower(-2), powT(-1)), sp.Rational(9, 16)))
    U = nj.add(nj.mul(powT(sp.Rational(-3, 4)), potential), nj.scale(nj.mul(psi_u, nj.Npower(-1)), -1))
    Na = nj.mul(nj.Npower(1), powT(sp.Rational(-1, 4)))
    Nm = nj.mul(nj.Npower(1), powT(sp.Rational(3, 4)))
    pieces = {
        "A": nj.add(nj.scale(Na, -3*HH**2), nj.scale(nj.mul(Na, P), -3*HH),
                    nj.scale(nj.mul(Na, P, P), -sp.Rational(3, 4)), nj.scale(nj.mul(nj.Npower(1), U), -1),
                    nj.scale(Nm, l**2/2)),
        "B": nj.add(nj.scale(Na, 2*HH), nj.mul(Na, P)),
        "C": nj.scale(Nm, -sp.Rational(1, 2)),
        "L": nj.scale(Nm, l),
        "D": nj.scale(Na, -sp.Rational(1, 3)),
        "E": nj.scale(Na, 2),
        "M": nj.scale(Nm, sp.Rational(1, 2)),
        "Z": nj.scale(Na, sp.Rational(1, 2)),
    }
    result = {key: nj.derivative_coefficients(row) for key, row in pieces.items()}
    return {"symbols": (R, Rd, Rdd, HH, l, f0, f1, f2), "jets": result,
            "primitive_jet": sum(value*(N-1)**k for k, value in enumerate(psi)),
            "primitive_u_jet": sum(value*(N-1)**k for k, value in enumerate(psi_u))}


@cache
def lapse_jets(*, omit_boundary=False):
    f = functions()
    generic = generic_lapse_jets(omit_boundary=omit_boundary)
    poly = sp.Poly(f["F"], X)
    values = (1/d**3, sp.diff(1/d**3, u), sp.diff(1/d**3, u, 2), H, velocity,
              poly.nth(0), poly.nth(1), poly.nth(2))
    mapping = dict(zip(generic["symbols"], values))
    return {key: tuple(sp.factor(expr.subs(mapping)) for expr in row)
            for key, row in generic["jets"].items()}


@cache
def audit_checks():
    f, jets = functions(), lapse_jets()
    return {**transformation_checks(), **legendre_checks(), **nj.identities(),
            "background_lapse": jets["A"][1],
            "nonlinear_to_J_including_matter": sp.cancel(jets["A"][2]+2*f["J"]),
            "nonlinear_to_Theta": sp.cancel(jets["B"][1]-2*f["Theta"]),
            "nonlinear_to_Lambda": sp.cancel(jets["C"][1]+f["Lambda"]/2),
            "nonlinear_to_matter_mixing": sp.cancel(jets["L"][1]+f["w"]),
            "matter_equation": sp.cancel(sp.diff(d**6*velocity, u))}
