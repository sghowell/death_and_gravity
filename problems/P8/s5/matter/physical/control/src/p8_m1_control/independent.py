"""Independent Hamiltonian and original-cosmic-time normalization checks.

This module does not import ``model`` or ``oscillator``.  It starts with the
pinned physical Hamiltonian, and differentiates its original time functions
before evaluating a point.  Point fixtures are not all-time bounds.

Conventions: H2=p^T A p/2+p^T B Q+Q^T C Q/2 has measure a^3.  If
T^T T=A^-1, Y=a^(3/2) T Q and pi=a^(3/2) T^-T p.  The final momentum is
P=pi+S Y, with L=|Ydot+Omega Y|^2/2-Y^T W Y/2.  No antisymmetric connection
is discarded, and positivity of A alone does not imply positivity of W.
"""

from functools import cache

import sympy as sp
from p8_m1 import nonlinear
from p8_m1_physical import quadratic


def _clean(matrix):
    return matrix.applyfunc(sp.cancel)


def _scalar(matrix):
    return matrix[0, 0]


def _symmetric(matrix):
    return (matrix+matrix.T)/2


def oscillator_from_jets(A, B, C, Bdot, T, Tdot, Tddot, H, Hdot):
    """Whiten a physical Hamiltonian using exact zeroth/first/second jets.

    All entries and derivatives are in one physical cosmic-time system of
    units.  In particular, q in these inputs is k_com^2/a^2, not ell^2 times
    that quantity.  The returned kinetic residual checks the supplied T.
    """
    inverse = _clean(T.inv())
    inverse_dot = -inverse*Tdot*inverse
    F = 3*H*sp.eye(2)/2+Tdot*inverse
    Bhat = _clean(T*B*inverse+F)
    Bhatdot = _clean(Tdot*B*inverse+T*Bdot*inverse+T*B*inverse_dot
                    +3*Hdot*sp.eye(2)/2+Tddot*inverse+Tdot*inverse_dot)
    S, Sdot = _symmetric(Bhat), _symmetric(Bhatdot)
    Omega = _clean(-(Bhat-Bhat.T)/2)
    Chat = _clean(inverse.T*C*inverse)
    W = _clean(Chat-S*S+S*Omega-Omega*S-Sdot)
    return {"kinetic_residual": _clean(T*A*T.T-sp.eye(2)),
            "whitened_cross": Bhat, "whitened_potential": Chat,
            "momentum_boundary": _clean(S), "connection": Omega,
            "potential": W, "volume_connection": _clean(F)}


@cache
def generic_checks():
    """Exact two-species identities with arbitrary coupled matrix entries."""
    a, b = sp.symbols("a b", positive=True)
    c, a1, b1, c1, nu = sp.symbols("c a1 b1 c1 nu", real=True)
    T = sp.Matrix([[a, 0], [c, b]])
    Td = sp.Matrix([[a1, 0], [c1, b1]])
    inverse = T.inv()
    A = inverse*inverse.T
    B = sp.Matrix(2, 2, sp.symbols("b00 b01 b10 b11", real=True))
    c00, c01, c11 = sp.symbols("c00 c01 c11", real=True)
    C = sp.Matrix([[c00, c01], [c01, c11]])
    Y = sp.Matrix(sp.symbols("Y0 Y1", real=True))
    pi = sp.Matrix(sp.symbols("pi0 pi1", real=True))
    velocity = sp.Matrix(sp.symbols("V0 V1", real=True))
    F = nu*sp.eye(2)+Td*inverse
    Bhat, Chat = T*B*inverse+F, inverse.T*C*inverse
    old_p, old_Q = T.T*pi, inverse*Y
    transformed = (_scalar(old_p.T*A*old_p)/2+_scalar(old_p.T*B*old_Q)
                   +_scalar(old_Q.T*C*old_Q)/2+_scalar(pi.T*F*Y))
    normalized = (_scalar(pi.T*pi)/2+_scalar(pi.T*Bhat*Y)
                  +_scalar(Y.T*Chat*Y)/2)
    checks = {"Hamiltonian_whitening": sp.cancel(transformed-normalized)}
    checks.update({f"kinetic_whitening_{i}{j}": sp.cancel((T*A*T.T-sp.eye(2))[i, j])
                   for i in range(2) for j in range(2)})

    # Use independent arbitrary entries in the normalized frame.  This proves
    # the boundary/connection algebra without using a cosmological identity.
    s0, s1, s2, ds0, ds1, ds2, omega = sp.symbols(
        "s0 s1 s2 ds0 ds1 ds2 omega", real=True)
    S = sp.Matrix([[s0, s1], [s1, s2]])
    Sdot = sp.Matrix([[ds0, ds1], [ds1, ds2]])
    Omega = sp.Matrix([[0, omega], [-omega, 0]])
    Bhat = S-Omega
    W = C-S*S+S*Omega-Omega*S-Sdot
    P = sp.Matrix(sp.symbols("P0 P1", real=True))
    shifted_pi = P-S*Y
    before = (_scalar(shifted_pi.T*shifted_pi)/2
              +_scalar(shifted_pi.T*Bhat*Y)+_scalar(Y.T*C*Y)/2
              -_scalar(Y.T*Sdot*Y)/2)
    after = _scalar(P.T*P)/2-_scalar(P.T*Omega*Y)+_scalar(Y.T*W*Y)/2
    checks["symmetric_canonical_boundary"] = sp.expand(before-after)
    original_lagrangian = (_scalar((velocity-Bhat*Y).T*(velocity-Bhat*Y))/2
                           -_scalar(Y.T*C*Y)/2)
    covariant_lagrangian = (_scalar((velocity+Omega*Y).T*(velocity+Omega*Y))/2
                            -_scalar(Y.T*W*Y)/2)
    boundary_dot = _scalar(velocity.T*S*Y)+_scalar(Y.T*Sdot*Y)/2
    checks["symmetric_Lagrangian_boundary"] = sp.expand(
        original_lagrangian-covariant_lagrangian+boundary_dot)
    checks.update({f"canonical_velocity_{i}": sp.expand(sp.diff(after, P[i])
                   -(P-Omega*Y)[i]) for i in range(2)})
    checks.update({f"canonical_momentum_{i}": sp.expand(-sp.diff(after, Y[i])
                   -(-W*Y-Omega*P)[i]) for i in range(2)})

    dw0, dw1, dw2 = sp.symbols("dw0 dw1 dw2", real=True)
    Wdot = sp.Matrix([[dw0, dw1], [dw1, dw2]])
    energy_rate = (_scalar(P.T*(-C*Y-Omega*P))
                   +_scalar(Y.T*C*(P-Omega*Y))+_scalar(Y.T*Wdot*Y)/2)
    covariant_rate = _scalar(Y.T*(Wdot+Omega*C-C*Omega)*Y)/2
    checks["covariant_energy_identity"] = sp.expand(energy_rate-covariant_rate)
    covariant = Wdot+Omega*C-C*Omega
    checks.update({f"covariant_derivative_symmetric_{i}{j}": sp.expand(
        (covariant-covariant.T)[i, j]) for i in range(2) for j in range(2)})
    return checks


def generic_negative_controls():
    """Nonzero omission residuals; they are not instability diagnoses."""
    y0, y1, p0, p1, nu, sdot, omega, w0, w1 = sp.symbols(
        "y0 y1 p0 p1 nu sdot omega w0 w1", real=True)
    return {
        "missing_volume_generator": nu*(p0*y0+p1*y1),
        "missing_symmetric_boundary_derivative": sdot*y0*y1,
        "reversed_covariant_commutator": 2*omega*(w1-w0)*y0*y1,
        "discarded_connection_in_velocity": sp.Matrix([-omega*y1, omega*y0]),
    }


@cache
def cosmic_background_jets(u_value, kcom_squared):
    """Differentiate the original witness in u=t/tau, setting tau=1.

    kcom_squared is held constant.  a=(1+u^2)^2, ell=sqrt(1+u^2), so the
    physical and local dimensionless momenta have different time drifts.
    """
    u_value, kcom_squared = sp.sympify(u_value), sp.sympify(kcom_squared)
    u, d, raw = nonlinear.u, nonlinear.d, nonlinear.functions()
    formulas = {quadratic.H: nonlinear.H, quadratic.l: nonlinear.velocity,
                quadratic.theta: raw["Theta"], quadratic.lam: raw["Lambda"],
                quadratic.w: raw["w"], quadratic.J: raw["J"],
                quadratic.q: kcom_squared/d**4}
    jets = tuple({symbol: sp.cancel(sp.diff(value, u, order).subs(u, u_value))
                  for symbol, value in formulas.items()} for order in range(3))
    ell = sp.sqrt(1+u_value**2)
    return {"jets": jets, "ell": ell, "x": u_value/ell, "y": 1/ell,
            "local_q": kcom_squared/ell**6, "physical_q": kcom_squared/ell**8,
            "scale_factor": ell**4, "volume_root": ell**6}


def _jet(expression, background, order=2):
    """Chain-rule differentiation before point evaluation, without expansion."""
    values, first, second = background["jets"]
    symbols = tuple(key for key in values if expression.has(key))
    at = lambda value: sp.cancel(value.subs(values, simultaneous=True))
    result = [at(expression)]
    if order >= 1:
        result.append(sp.cancel(sum(at(sp.diff(expression, key))*first[key]
                                    for key in symbols)))
    if order >= 2:
        result.append(sp.cancel(sum(at(sp.diff(expression, key))*second[key]
                                    for key in symbols)
                                +sum(at(sp.diff(expression, key, other))*first[key]*first[other]
                                     for key in symbols for other in symbols)))
    return tuple(result)


def _physical_factor_jets(chart, background):
    J, theta, lam, w, q = quadratic.J, quadratic.theta, quadratic.lam, quadratic.w, quadratic.q
    if chart == "unitary":
        factors, delta = (2*J/theta**2, sp.Integer(1), w/theta), 0
    elif chart == "gamma":
        R, D = lam**2-J/q, lam**2-(J+w**2/2)/q
        factors, delta = (2*J/R, R/D, -lam*w/R), 1
    else:
        raise ValueError("Use unitary or gamma scalar chart")
    return tuple(_jet(value, background) for value in factors), delta


def _triangular_jets(factors):
    (d1, d1d, d1dd), (d2, d2d, d2dd), (chi, chid, chidd) = factors
    a, b = sp.sqrt(d1), sp.sqrt(d2)
    ad, bd = d1d/(2*a), d2d/(2*b)
    add, bdd = d1dd/(2*a)-d1d**2/(4*a**3), d2dd/(2*b)-d2d**2/(4*b**3)
    return (sp.Matrix([[a, 0], [chi*b, b]]),
            sp.Matrix([[ad, 0], [chid*b+chi*bd, bd]]),
            sp.Matrix([[add, 0], [chidd*b+2*chid*bd+chi*bdd, bdd]]))


@cache
def cosmic_time_normalization(chart, u_value, kcom_squared):
    """Independent exact point fixture, including physical normalization jets."""
    background = cosmic_background_jets(u_value, kcom_squared)
    base = quadratic.symbolic(chart)
    values, first, _ = background["jets"]
    if values[quadratic.J].is_nonpositive is True or values[quadratic.q].is_nonpositive is True:
        raise ValueError("Require positive J and physical q")
    if chart == "unitary" and values[quadratic.theta] == 0:
        raise ValueError("The unitary velocity chart is singular at theta=0")
    if chart == "gamma":
        domain = sp.cancel(values[quadratic.lam]**2
                           -(values[quadratic.J]+values[quadratic.w]**2/2)/values[quadratic.q])
        if domain.is_nonpositive is True:
            raise ValueError("The gamma velocity chart requires q*lambda^2>J+w^2/2")
    at = lambda matrix: matrix.subs(values, simultaneous=True).applyfunc(sp.cancel)
    A, B = at(base["A"]), at(base["B"])
    C = at(sp.hessian(base["density"], quadratic.Q))
    Bdot = base["B"].applyfunc(lambda value: _jet(value, background, 1)[1])
    factor_jets, delta = _physical_factor_jets(chart, background)
    if any(value[0].is_nonpositive is True for value in factor_jets[:2]):
        raise ValueError("The requested point does not have positive scalar kinetic factors")
    T, Tdot, Tddot = _triangular_jets(factor_jets)
    result = oscillator_from_jets(A, B, C, Bdot, T, Tdot, Tddot,
                                  values[quadratic.H], first[quadratic.H])
    return {**background, **result, "A": A, "B": B, "C": C,
            "alpha": at(base["alpha"]), "T": T, "Tdot": Tdot, "Tddot": Tddot,
            "factor_jets": factor_jets, "delta": delta}


def _bounce_q(q):
    q = sp.sympify(q)
    if q.is_number and q <= 6:
        raise ValueError("The positive gamma bounce velocity chart requires q>6")
    return q


@cache
def bounce_from_cosmic_time(q):
    """Exact bounce W from Hamiltonian time jets, without the main derivation.

    First time jets, B and H vanish at the bounce, so W=Chat-Sdot.  The
    off-diagonal is stored as sqrt(d2/d1) times the reported rational factor.
    This special-time simplification is not an all-time oscillator formula.
    """
    q = _bounce_q(q)
    background = cosmic_background_jets(sp.Integer(0), q)
    base = quadratic.symbolic("gamma")
    values, first, _ = background["jets"]
    factors, _ = _physical_factor_jets("gamma", background)
    (d1, _, d1dd), (d2, _, d2dd), (chi, _, chidd) = factors
    C = sp.hessian(base["density"], quadratic.Q).subs(values, simultaneous=True)
    Bd = base["B"].applyfunc(lambda value: _jet(value, background, 1)[1])
    b11, b22 = Bd[0, 0]-chi*Bd[0, 1], chi*Bd[0, 1]+Bd[1, 1]
    b12_factor = d1*Bd[0, 1]/d2
    b21_factor = chi*Bd[0, 0]+Bd[1, 0]-chi**2*Bd[0, 1]-chi*Bd[1, 1]+chidd
    W11 = ((C[0, 0]-2*chi*C[0, 1]+chi**2*C[1, 1])/d1
           -b11-3*first[quadratic.H]/2-d1dd/(2*d1))
    W22 = C[1, 1]/d2-b22-3*first[quadratic.H]/2-d2dd/(2*d2)
    W12 = (C[0, 1]-chi*C[1, 1])/d2-(b12_factor+b21_factor)/2
    checks = {"background_H": values[quadratic.H],
              "background_B": base["B"].subs(values, simultaneous=True),
              "first_factor_jets": sp.Matrix([value[1] for value in factors])}
    return {"d1": d1, "d2": d2, "chi": chi, "ratio": sp.cancel(d2/d1),
            "potential11": sp.cancel(W11), "potential22": sp.cancel(W22),
            "potential12_factor": sp.cancel(W12), "connection_factor": sp.Integer(0),
            "connection_dot_factor": sp.cancel((b21_factor-b12_factor)/2),
            "vanishing_checks": checks}


def bounce_expected(q):
    """Closed known-answer formulas for the positive gamma chart at u=0.

    q=8 is a useful indefinite-W control, not a dynamical-instability proof.
    These formulas retain finite-q mixing and are not principal K=G data.
    """
    q = _bounce_q(q)
    d1 = 2398*q/(200*q-1199)
    d2 = (200*q-1199)/(200*(q-6))
    chi = 20*q/(200*q-1199)
    ratio = sp.cancel(d2/d1)
    W11 = (239800*q**2-5994805*q-2875202)/(1199*(200*q-1199))
    W22 = (200*q**3-3603*q**2+21620*q-43338)/((q-6)*(200*q-1199))
    # Earlier Hamiltonian notes factored sqrt(d1/d2); use the main convention.
    W12 = -80*q*(199*q+5103)/(200*q-1199)**2
    return {"d1": d1, "d2": d2, "chi": chi, "ratio": ratio,
            "potential11": W11, "potential22": W22, "potential12_factor": W12,
            "connection_factor": sp.Integer(0),
            "connection_dot_factor": -80*q*(1000*q-12297)/(200*q-1199)**2,
            "determinant": (1199*q**3-37198*q**2+166106*q+86676)/(1199*(q-6))}


def unitary_connection_expected(u_value):
    """Physical Omega12 for positive lower-Cholesky T, away from theta=0."""
    u_value = sp.sympify(u_value)
    if u_value == 0:
        raise ValueError("The unitary velocity chart is singular at the bounce")
    raw, u = nonlinear.functions(), nonlinear.u
    theta, lam, J, H, matter = (value.subs(u, u_value) for value in
                              (raw["Theta"], raw["Lambda"], raw["J"],
                               nonlinear.H, nonlinear.velocity))
    return -sp.sign(theta)*sp.sqrt(2)*matter*(theta-H*lam)/sp.sqrt(J)
