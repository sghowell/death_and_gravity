"""Both physical scalar principal charts, including the moving boundary."""
from functools import cache

import sympy as sp
from p8_affine_aligned import dynamics as parent

from . import action

old = parent.old
delta_J = sp.Symbol("delta_J", real=True)
pb, ps = sp.symbols("canonical_P_b canonical_P_s", real=True)
theta_dot, lam_dot = sp.symbols("theta_dot lambda_dot", real=True)


def clean(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.factor))


@cache
def scalar():
    prior = parent.scalar()
    L = prior["after_shift"]+delta_J*old.n**2
    lapse = prior["lapse"]
    reduced = L.subs(old.temporal, prior["temporal"]).subs(old.n, lapse)
    K = clean(sp.hessian(reduced, (old.vd, old.sd, old.sigmad))/2)
    before = (-2*old.q*old.shift*old.vd+old.pm*old.sd+old.pi*old.sigmad-L).subs(
        {old.vd: old.theta*old.n-old.ell*old.matter/2,
         old.sd: old.pm-old.w*old.n,
         old.sigmad: old.temporal+old.pi/(parent.ZETA*old.q)}, simultaneous=True)
    after_temporal = before.subs(old.temporal, old.pi)
    R = parent.first_order()["R"]
    Jnew = old.J+delta_J
    hamiltonian = parent.first_order()["H"].subs(old.J, Jnew)
    return {"L": L, "K": K, "lapse": lapse, "R": R,
            "regular_lapse": -R/Jnew, "H": hamiltonian,
            "shift_constraint": sp.factor(sp.diff(L, old.shift)/(2*old.q)
                                          -old.theta*old.n+old.vd+old.ell*old.matter/2),
            "temporal_constraint_unchanged": sp.factor(sp.diff(L-prior["after_shift"], old.temporal)),
            "unitary_kinetic_change": clean(K-prior["kinetic"]-sp.diag(delta_J/old.theta**2, 0, 0)),
            "regular_lapse_equation": sp.factor(sp.diff(after_temporal, old.n).subs(old.n, -R/Jnew)),
            "regular_Hamiltonian": sp.factor(after_temporal.subs(old.n, -R/Jnew)-hamiltonian),
            "no_Theta_inverse": not sp.denom(sp.together(hamiltonian)).has(old.theta)}


@cache
def gamma():
    light_H = scalar()["H"].subs({old.pi: 0, old.sigma: 0})
    H = old.a**3*light_H.subs(
        {old.v: pb/(2*old.a**3*old.q), old.pm: ps/old.a**3}, simultaneous=True)
    H -= old.H*old.shift*pb
    coords, momenta = (old.shift, old.matter), (pb, ps)
    A = sp.hessian(H, momenta)
    B = sp.Matrix(2, 2, lambda i, j: sp.diff(H, momenta[i], coords[j]))
    C = sp.hessian(H, coords)
    alpha = clean(A.inv())
    beta = clean(-alpha*B)
    gamma_matrix = clean(B.T*alpha*B-C)

    def principal(matrix, power):
        return clean(matrix.applyfunc(lambda value: sp.limit(
            value/(old.a**3*old.q**power), old.q, sp.oo)))

    K = principal(alpha/2, 0)
    beta_q = principal(beta, 1)
    gamma_q = principal(gamma_matrix, 1)
    # (a^3 q)'=H a^3 q at fixed comoving momentum.
    beta_q_dot = beta_q.diff(old.theta)*theta_dot+beta_q.diff(old.lam)*lam_dot
    G = clean((old.H*beta_q+beta_q_dot-gamma_q)/2)
    J0 = old.J+delta_J+old.w**2/2
    D = old.q*old.lam**2-J0
    expected_K = sp.Matrix([[old.q*J0/D, -old.q*old.lam*old.w/(2*D)],
                           [-old.q*old.lam*old.w/(2*D), (old.q*old.lam**2-old.J-delta_J)/(2*D)]])
    return {"H": H, "A": clean(A), "alpha": alpha, "beta": beta, "gamma": gamma_matrix,
            "K_principal": K, "G_principal": G,
            "beta_leading": beta_q, "gamma_leading": gamma_q,
            "finite_q_kinetic": clean(alpha/(2*old.a**3)),
            "finite_q_kinetic_identity": clean(alpha/(2*old.a**3)-expected_K),
            "beta_leading_symmetric": clean(beta_q-beta_q.T),
            "gamma_has_no_q_squared_principal": principal(gamma_matrix, 2),
            "principal_time_boundary_retained": clean(G-sp.Matrix([
                [(-old.lam*theta_dot+old.H*old.lam*old.theta+old.theta*lam_dot-old.theta**2)/old.lam**2,
                 old.ell/2], [old.ell/2, sp.Rational(1, 2)]]))}


@cache
def cone_checks():
    bg = old.background()
    speed = sp.Symbol("squared_physical_speed", real=True)
    dJ = action.deformation()["delta_J"]
    out = {"actual_matter_mixing_identity": sp.factor(bg["w"]+bg["ell"]*bg["lam"]),
           "actual_moving_gradient_identity": sp.factor(
               bg["J"]+bg["w"]**2/2-bg["theta"]*(bg["H"]*bg["lam"]+sp.diff(bg["lam"], old.u))
               +bg["lam"]*sp.diff(bg["theta"], old.u)+bg["theta"]**2)}
    for chart, denominator, cross_sign in (("unitary", old.theta, 1), ("gamma", old.lam, -1)):
        K0 = sp.Matrix([[(old.J+old.w**2/2)/denominator**2, cross_sign*old.w/(2*denominator)],
                        [cross_sign*old.w/(2*denominator), sp.Rational(1, 2)]])
        K = K0+sp.diag(delta_J/denominator**2, 0)
        out[chart+"_characteristic"] = sp.factor(
            (K0-speed*K).det()/K.det()-(1-speed)*(old.J/(old.J+delta_J)-speed))
        out[chart+"_kinetic_determinant"] = sp.factor(K.det()-(old.J+delta_J)/(2*denominator**2))
        out[chart+"_gradient_determinant"] = sp.factor(K0.det()-old.J/(2*denominator**2))
        velocity = sp.Matrix(sp.symbols("clock_velocity matter_velocity", real=True))
        out[chart+"_positive_kinetic_squares"] = sp.factor(
            (velocity.T*K*velocity)[0]
            -(velocity[1]+cross_sign*old.w*velocity[0]/denominator)**2/2
            -(old.J+delta_J)*velocity[0]**2/denominator**2)
    actual = {symbol: bg[name] for symbol, name in ((old.J, "J"), (old.w, "w"),
              (old.theta, "theta"), (old.lam, "lam"), (old.H, "H"), (old.ell, "ell"))}
    actual |= {delta_J: dJ, theta_dot: sp.diff(bg["theta"], old.u), lam_dot: sp.diff(bg["lam"], old.u)}
    kg, gg = gamma()["K_principal"], gamma()["G_principal"]
    out["actual_regular_gamma_cone_margin"] = clean(
        (kg-gg-sp.diag(delta_J/old.lam**2, 0)).subs(actual, simultaneous=True))
    return out


@cache
def center():
    bg = old.background()
    actual = {symbol: bg[name].subs(old.u, 0) for symbol, name in
              ((old.J, "J"), (old.w, "w"), (old.lam, "lam"))}
    actual[delta_J] = 4*action.epsilon
    canonical_point = dict(actual)
    canonical_point.update({old.a: 1, old.H: 0, old.theta: 0, old.ell: sp.Rational(1, 10)})
    # Every other first background jet vanishes at the center, including
    # a', q', J', delta_J', lambda', ell' and w'. H'=4, Theta'=3.
    beta_prime = 4*gamma()["beta"].diff(old.H)+3*gamma()["beta"].diff(old.theta)
    center_gradient = clean(((beta_prime-gamma()["gamma"])/(2*old.q)).subs(canonical_point, simultaneous=True))
    principal_gradient = clean(center_gradient.applyfunc(lambda value: sp.limit(value, old.q, sp.oo)))
    return {"finite_q_kinetic": clean(gamma()["finite_q_kinetic"].subs(actual, simultaneous=True)),
            "independently_differentiated_principal_gradient": principal_gradient,
            "center_gradient_identity": clean(principal_gradient-sp.Matrix([[6, sp.Rational(1, 20)],
                                                                             [sp.Rational(1, 20), sp.Rational(1, 2)]])),
            "velocity_chart_threshold": 6+16*action.epsilon,
            "principal_squared_speeds": (sp.Integer(1), sp.Rational(1199, 800)/(sp.Rational(1199, 800)+4*action.epsilon)),
            "regular_phase_Hamiltonian_does_not_require_this_velocity_threshold": True}


@cache
def checks():
    s, g = scalar(), gamma()
    out = {name: s[name] for name in ("shift_constraint", "temporal_constraint_unchanged",
           "unitary_kinetic_change", "regular_lapse_equation", "regular_Hamiltonian")}
    out.update({name: g[name] for name in ("finite_q_kinetic_identity", "beta_leading_symmetric",
               "gamma_has_no_q_squared_principal", "principal_time_boundary_retained")})
    out.update(cone_checks())
    out["independent_center_time_jet_gradient"] = center()["center_gradient_identity"]
    return out
