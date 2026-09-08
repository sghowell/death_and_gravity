"""Literal homogeneous action, physical stress source and regular reduction."""
from functools import cache

import sympy as sp
from p8_aligned_margin import dynamics as margin

old = margin.old
rho, pressure = sp.symbols("normalized_vector_energy normalized_vector_pressure", real=True)
Fv, Fn = sp.symbols("curvature_force lapse_force", real=True)
p, pm, shifted, Y = sp.symbols("homogeneous_p homogeneous_pm shifted_p scaled_shifted_p", real=True)
delta = sp.Symbol("physical_metric_lapse_shift", real=True)
Jnew = old.J+margin.delta_J


def clean(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.factor))


@cache
def action():
    # Restrict the literal local density before any nonzero-q shift division.
    L = margin.scalar()["L"].subs({old.q: 0, old.temporal: 0})+Fv*old.v+Fn*old.n
    velocities = {old.vd: old.theta*old.n-(p+3*old.ell*old.matter)/6,
                  old.sd: pm-old.w*old.n}
    before = sp.factor((p*old.vd+pm*old.sd-L).subs(velocities, simultaneous=True))
    numerator = old.theta*(p+3*old.ell*old.matter)-old.w*pm-Fn
    lapse = numerator/(2*Jnew)
    target = -(p+3*old.ell*old.matter)**2/12+pm**2/2+numerator**2/(4*Jnew)-Fv*old.v
    return {"L": L, "before_lapse": before, "lapse": lapse, "H": target,
            "metric_scalar_momentum": sp.factor(sp.diff(L, old.vd).subs(velocities)-p),
            "matter_scalar_momentum": sp.factor(sp.diff(L, old.sd).subs(velocities)-pm),
            "lapse_equation": sp.factor(sp.diff(before, old.n).subs(old.n, lapse)),
            "reduced_Hamiltonian": sp.factor(before.subs(old.n, lapse)-target),
            "no_Theta_or_spatial_momentum_inverse": not sp.denom(sp.together(target)).has(old.theta, old.q)}


@cache
def reduction():
    data = action()
    vd, sd = sp.diff(data["H"], p), sp.diff(data["H"], pm)
    pd = -sp.diff(data["H"], old.v)-3*old.H*p
    pmd = -sp.diff(data["H"], old.matter)-3*old.H*pm
    conserved = sp.factor(pmd+3*old.ell*vd+3*old.H*pm)
    Sd = pd-9*old.H*old.ell*old.matter+3*old.ell*sd
    sub = {p: shifted-3*old.ell*old.matter, pm: -3*old.ell*old.v}
    # Zero physical matter momentum perturbation is preserved. The exact
    # invariant is a^3*(pm+3ell*v), not pm separately.
    vdot = sp.factor(vd.subs(sub, simultaneous=True))
    Sdot = sp.factor(Sd.subs(sub, simultaneous=True))
    scaled_vdot = vdot.subs(shifted, Y/(4*old.a**3))
    scaled_Ydot = 4*old.a**3*(Sdot+3*old.H*shifted).subs(shifted, Y/(4*old.a**3))
    vector = sp.Matrix([scaled_vdot, scaled_Ydot])
    A = old.theta**2/(2*Jnew)-sp.Rational(1, 6)
    B = 3*old.theta*old.w*old.ell/(2*Jnew)
    C = 9*old.ell**2*(1+old.w**2/(2*Jnew))
    matrix = sp.Matrix([[B, A/(4*old.a**3)], [-4*old.a**3*C, -B]])
    force = sp.Matrix([-old.theta*Fn/(2*Jnew),
                       4*old.a**3*(Fv+3*old.ell*old.w*Fn/(2*Jnew))])
    n = (old.theta*Y/(4*old.a**3)+3*old.w*old.ell*old.v-Fn)/(2*Jnew)
    return {"A": A, "B": B, "C": C, "matrix": clean(matrix), "force": clean(force),
            "lapse": n, "matter_dot": -3*old.ell*old.v-old.w*n,
            "preserved_physical_matter_momentum": conserved,
            "scaled_regular_system": clean(vector-matrix*sp.Matrix([old.v, Y])-force),
            "trace_free_scaled_matrix": sp.factor(sp.trace(matrix)),
            "lapse_reconstruction": sp.factor(data["lapse"].subs(sub, simultaneous=True).subs(
                shifted, Y/(4*old.a**3))-n)}


@cache
def source():
    n, v = old.n, old.v
    # 1/2 T^{mu nu} delta g_{mu nu}, with actual physical curvature
    # perturbation v+delta*n and background T^{00}=rho,T^{ii}=p/a^2.
    contraction = -rho*n+3*pressure*(v+delta*n)
    target = 3*pressure*v+(3*delta*pressure-rho)*n
    return {"curvature_force": 3*pressure, "lapse_force": 3*delta*pressure-rho,
            "physical_stress_variation": sp.factor(contraction-target),
            "lapse_force_includes_physical_spatial_metric_variation":
                sp.factor(sp.diff(contraction, n)-(3*delta*pressure-rho))}


@cache
def checks():
    a, r, s = action(), reduction(), source()
    out = {name: a[name] for name in ("metric_scalar_momentum", "matter_scalar_momentum",
           "lapse_equation", "reduced_Hamiltonian")}
    out.update({name: r[name] for name in ("preserved_physical_matter_momentum",
               "scaled_regular_system", "trace_free_scaled_matrix", "lapse_reconstruction")})
    out.update({name: s[name] for name in ("physical_stress_variation",
               "lapse_force_includes_physical_spatial_metric_variation")})
    return out
