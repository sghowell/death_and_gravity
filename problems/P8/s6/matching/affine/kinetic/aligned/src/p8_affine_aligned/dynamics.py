"""Source-aligned quadratic action with the physical matter metric retained."""
from functools import cache

import sympy as sp
from p8_affine_retuned import dynamics as previous
from p8_affine_retuned import geometry

from . import alignment, modes

old = previous.old
u, q, ZETA = old.u, old.q, modes.ZETA


@cache
def scalar():
    t, n, sigma, sigmad = old.temporal, old.n, old.sigma, old.sigmad
    dd = previous.DD
    # In the old temporal variable, B_N=-d changes the electric term.
    literal = (old.action()["base"]+(t+dd*n)**2/2-q*sigma**2/2
               +ZETA*q*(sigmad-t-dd*n)**2/2)
    shifted = literal.subs(t, t-dd*n)
    proca = t**2/2-q*sigma**2/2+ZETA*q*(sigmad-t)**2/2
    expected = old.action()["base"]+proca
    temporal = ZETA*q*sigmad/(1+ZETA*q)
    C = ZETA*q/(2*(1+ZETA*q))
    lapse = (old.vd+old.ell*old.matter/2)/old.theta
    reduced = (old.action()["base"]+C*sigmad**2-q*sigma**2/2).subs(n, lapse)
    kinetic = geometry.clean(sp.hessian(reduced, (old.vd, old.sd, sigmad))/2)
    target = sp.Matrix([[(old.J+old.w**2/2)/old.theta**2, old.w/(2*old.theta), 0],
                        [old.w/(2*old.theta), sp.Rational(1, 2), 0], [0, 0, C]])
    return {"literal_before_temporal_shift": literal, "after_shift": expected, "proca": proca,
            "temporal": temporal, "lapse": lapse, "C": C, "reduced": reduced,
            "kinetic": kinetic, "expected_kinetic": target,
            "triangular_temporal_shift": sp.factor(shifted-expected),
            "temporal_Euler": sp.factor(sp.diff(proca, t).subs(t, temporal)),
            "temporal_action": sp.factor(proca.subs(t, temporal)-C*sigmad**2+q*sigma**2/2),
            "zero_curl_auxiliary_limit": sp.factor(literal.subs({ZETA: 0, t: -dd*n, sigma: 0})-old.action()["base"]),
            "joint_auxiliary_determinant": sp.factor(sp.hessian(literal, (n, old.shift, t)).det()
                                                      +4*q**2*old.theta**2*(1+ZETA*q)),
            "physical_kinetic_split": geometry.clean(kinetic-target)}


@cache
def first_order():
    data = scalar()
    before = (-2*q*old.shift*old.vd+old.pm*old.sd+old.pi*old.sigmad-data["after_shift"]).subs(
        {old.vd: old.theta*old.n-old.ell*old.matter/2,
         old.sd: old.pm-old.w*old.n,
         old.sigmad: old.temporal+old.pi/(ZETA*q)}, simultaneous=True)
    after = before.subs(old.temporal, old.pi)
    R = (q*(old.theta*old.shift+old.lam*old.v)+old.w*old.pm/2
         -3*old.ell*old.matter*old.theta/2)
    lapse = -R/old.J
    target = previous.first_order()["H"].subs({previous.DD: 0, previous.ZETA: ZETA})
    return {"H": target, "R": R, "lapse": lapse,
            "vector_temporal_Euler": sp.factor(sp.diff(before, old.temporal).subs(old.temporal, old.pi)),
            "lapse_Euler": sp.factor(sp.diff(after, old.n).subs(old.n, lapse)),
            "full_Hamiltonian": sp.factor(after.subs(old.n, lapse)-target),
            "no_Theta_inverse": not sp.denom(sp.together(target)).has(old.theta)}


@cache
def cones():
    bg = old.background()
    K0 = sp.Matrix([[(bg["J"]+bg["w"]**2/2)/bg["theta"]**2, bg["w"]/(2*bg["theta"])],
                    [bg["w"]/(2*bg["theta"]), sp.Rational(1, 2)]])
    K = sp.diag(1, 1, sp.Rational(1, 2))
    K[:2, :2] = K0
    G = sp.diag(1, 1, sp.Rational(1, 2))
    G[:2, :2] = previous.principal()["G0"]
    speed = sp.Symbol("c_squared", real=True)
    return {"K": K, "G": G, "speed_squared": speed,
            "matter_cone_residual": geometry.clean(K-G),
            "characteristic_residual": sp.factor((G-speed*K).det()/K.det()-(1-speed)**3),
            "finite_q_massive_phase_speed_not_used": True}


@cache
def checks():
    data, canonical, cone = scalar(), first_order(), cones()
    out = {"all_component_alignment_input": alignment.rolling()["first_variation"],
           "nonlinear_mass_expansion_input": alignment.quadratic_mass()["residual"],
           "physical_matter_cone": cone["matter_cone_residual"],
           "three_scalar_luminal_characteristics": cone["characteristic_residual"]}
    for name in ("triangular_temporal_shift", "temporal_Euler", "temporal_action",
                 "zero_curl_auxiliary_limit", "joint_auxiliary_determinant", "physical_kinetic_split"):
        out[name] = data[name]
    for name in ("vector_temporal_Euler", "lapse_Euler", "full_Hamiltonian"):
        out["first_order_"+name] = canonical[name]
    return out
