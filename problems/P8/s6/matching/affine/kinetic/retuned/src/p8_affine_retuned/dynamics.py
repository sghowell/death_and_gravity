"""All actual scalar constraints for the separately mass-retuned trace.

Principal signs, the primary null and an isolated mass parameter do not
establish a nonlinear constraint theorem or a controlled heavy-field EFT.
"""
from functools import cache

import sympy as sp
from p8_affine_kinetic import scalar as old
from p8_affine_ricci import dynamics as ricci
from p8_affine_ricci import source as original_source

from . import geometry

u, q = old.u, old.q
ZETA = sp.Symbol("zeta_retuned", positive=True)
DD = sp.Symbol("trace_lapse_source", real=True)


@cache
def source():
    original = original_source.stationary_variation()
    actual = geometry.clean(geometry.update()["N_full"]*original["actual"])
    expected = sp.Matrix([-15*original_source.H*original_source.n/(4*original_source.h), 0, 0, 0])
    bg = old.background()
    dd = 15*bg["H"]/(4*bg["h"])
    ratio = sp.factor(dd/bg["theta"])
    return {"actual": actual, "expected": expected, "d": sp.factor(dd), "f": ratio,
            "background": geometry.clean(geometry.update()["N_full"]*original["background"]),
            "first_variation_residual": geometry.clean(actual-expected)}


@cache
def action():
    t, n, sigma, sigmad = old.temporal, old.n, old.sigma, old.sigmad
    extra = (t+DD*n)**2/2-q*sigma**2/2+ZETA*q*(sigmad-t)**2/2
    before = old.action()["base"]+extra
    tsol = (ZETA*q*sigmad-DD*n)/(1+ZETA*q)
    coefficient = ZETA*q/(2*(1+ZETA*q))
    after = coefficient*(sigmad+DD*n)**2-q*sigma**2/2
    lapse = (old.vd+old.ell*old.matter/2)/old.theta
    reduced = (old.action()["base"]+after).subs(n, lapse)
    velocities = sp.Matrix([old.vd, old.sd, sigmad])
    kinetic = geometry.clean(sp.hessian(reduced, tuple(velocities))/2)
    squares = ((old.sd+old.w*old.vd/old.theta)**2/2
               +old.J*old.vd**2/old.theta**2+coefficient*(sigmad+DD*old.vd/old.theta)**2)
    auxiliaries = (n, old.shift, t)
    return {"before": before, "extra": extra, "temporal": tsol, "C": coefficient,
            "after_temporal": old.action()["base"]+after, "lapse": lapse,
            "reduced": reduced, "kinetic": kinetic,
            "temporal_Euler": sp.factor(sp.diff(extra, t).subs(t, tsol)),
            "temporal_action": sp.factor(extra.subs(t, tsol)-after),
            "shift_residual": sp.factor(sp.diff(before, old.shift)/(2*q)
                                        -old.theta*n+old.vd+old.ell*old.matter/2),
            "kinetic_squares": sp.factor((velocities.T*kinetic*velocities)[0]-squares),
            "kinetic_determinant": sp.factor(kinetic.det()-coefficient*old.J/(2*old.theta**2)),
            "all_auxiliary_determinant": sp.factor(sp.hessian(before, auxiliaries).det()
                                                   +4*q**2*old.theta**2*(1+ZETA*q))}


@cache
def first_order():
    data = action()
    first = (-2*q*old.shift*old.vd+old.pm*old.sd+old.pi*old.sigmad-data["before"]).subs(
        {old.vd: old.theta*old.n-old.ell*old.matter/2,
         old.sd: old.pm-old.w*old.n,
         old.sigmad: old.temporal+old.pi/(ZETA*q)}, simultaneous=True)
    temporal = old.pi-DD*old.n
    after_temporal = sp.factor(first.subs(old.temporal, temporal))
    R = (q*(old.theta*old.shift+old.lam*old.v)+old.w*old.pm/2
         -3*old.ell*old.matter*old.theta/2+DD*old.pi/2)
    constant = (old.pm**2/2+q*old.ell*old.shift*old.matter+q*old.matter**2/2
                -q*old.v**2-3*old.ell**2*old.matter**2/4
                +(1+1/(ZETA*q))*old.pi**2/2+q*old.sigma**2/2)
    lapse = -R/old.J
    target = constant+R**2/old.J
    return {"before": first, "temporal": temporal, "lapse": lapse, "H": target,
            "R": R, "constant": constant,
            "temporal_Euler": sp.factor(sp.diff(first, old.temporal).subs(old.temporal, temporal)),
            "before_lapse_identity": sp.factor(after_temporal-constant+old.J*old.n**2+2*R*old.n),
            "lapse_Euler": sp.factor(sp.diff(after_temporal, old.n).subs(old.n, lapse)),
            "Hamiltonian_identity": sp.factor(after_temporal.subs(old.n, lapse)-target)}


@cache
def principal():
    """Exact time-dependent point transformation rho=sigma+f(u)*v.

The old spatial boundary is already included in G0. No time boundary
of a coefficient-dependent new term is silently discarded.
"""
    bg = old.background()
    f, dd = source()["f"], source()["d"]
    C = action()["C"]
    G0 = geometry.clean(ricci.gradients()["matrix"].subs(ricci.COUPLING, 0))
    driving = sp.Matrix([sp.diff(f, u), -f*bg["ell"]/2])
    rho, rhod = sp.symbols("rho rho_dot", real=True)
    x = sp.Matrix([old.v, old.matter, rho])
    velocity = sp.Matrix([old.vd, old.sd, rhod])
    mass = sp.Matrix([-f, 0, 1])
    G = sp.zeros(3)
    G[:2, :2] = G0-C*driving*driving.T/q
    G += mass*mass.T/2
    hvec = sp.Matrix([*driving, 0])
    K = sp.zeros(3)
    K[:2, :2] = sp.Matrix([
        [(bg["J"]+bg["w"]**2/2)/bg["theta"]**2, bg["w"]/(2*bg["theta"])],
        [bg["w"]/(2*bg["theta"]), sp.Rational(1, 2)]])
    K[2, 2] = C
    linear = sp.Matrix([0, 0, -2*C*(driving.T*x[:2, :])[0]])
    original_extra = C*(old.sigmad+dd*(old.vd+bg["ell"]*old.matter/2)/bg["theta"])**2-q*old.sigma**2/2
    transformed = original_extra.subs({old.sigma: rho-f*old.v,
        old.sigmad: rhod-f*old.vd-sp.diff(f, u)*old.v}, simultaneous=True)
    expected_extra = C*(rhod-(hvec.T*x)[0])**2-q*(rho-f*old.v)**2/2
    lagrangian = (velocity.T*K*velocity)[0]+(velocity.T*linear)[0]-q*(x.T*G*x)[0]
    expected = ((velocity[:2, :].T*K[:2, :2]*velocity[:2, :])[0]
                -q*(x[:2, :].T*G0*x[:2, :])[0]+expected_extra)
    ratio = sp.factor((driving.T*G0.inv()*driving)[0]/2)
    numerator, denominator = sp.fraction(ratio)
    witness = sp.expand(1000*denominator-numerator)
    return {"G0": G0, "G": G, "K": K, "driving": driving, "mass": mass,
            "linear_velocity_coefficient": linear, "x": x, "velocity": velocity,
            "rho": rho, "rho_dot": rhod, "L": lagrangian,
            "gradient_relative_ratio": ratio, "bound_numerator": numerator,
            "bound_denominator": denominator, "bound_witness": witness,
            "point_transformation_residual": sp.factor(transformed-expected_extra),
            "principal_action_residual": sp.factor(lagrangian-expected)}


@cache
def center():
    bg = old.background()
    pb, ps, pt = sp.symbols("P_b P_s P_trace", real=True)
    canonical = old.a**3*first_order()["H"].subs(
        {old.v: pb/(2*old.a**3*q), old.pm: ps/old.a**3, old.pi: pt/old.a**3}, simultaneous=True)
    canonical -= old.H*old.shift*pb
    values = {symbol: bg[name] for symbol, name in
              ((old.a, "a"), (old.J, "J"), (old.theta, "theta"), (old.lam, "lam"),
               (old.ell, "ell"), (old.w, "w"), (old.H, "H"))}
    values[DD] = source()["d"]
    actual = sp.factor(canonical.subs(values, simultaneous=True).subs(u, 0))
    momentum = geometry.clean(sp.hessian(actual, (pb, ps, pt)))
    kinetic = geometry.clean(momentum.inv()/2)
    expected = sp.zeros(3)
    expected[:2, :2] = ricci.center()["kinetic"]
    expected[2, 2] = action()["C"]
    return {"H": actual, "momentum_hessian": momentum, "kinetic": kinetic,
            "expected": expected, "velocity_domain": "q>6 and zeta>0"}


@cache
def checks():
    src, data, canonical, energy = source(), action(), first_order(), principal()
    out = {"all_component_background_trace": src["background"],
           "all_component_rolling_source": src["first_variation_residual"],
           "smooth_source_over_Theta": sp.factor(src["f"]-15/(4*(1+u**2)**3-1)),
           "point_transformation": energy["point_transformation_residual"],
           "full_principal_action": energy["principal_action_residual"],
           "regular_center_kinetic": geometry.clean(center()["kinetic"]-center()["expected"])}
    for name in ("temporal_Euler", "temporal_action", "shift_residual", "kinetic_squares",
                 "kinetic_determinant", "all_auxiliary_determinant"):
        out["scalar_"+name] = data[name]
    for name in ("temporal_Euler", "before_lapse_identity", "lapse_Euler", "Hamiltonian_identity"):
        out["first_order_"+name] = canonical[name]
    return out
