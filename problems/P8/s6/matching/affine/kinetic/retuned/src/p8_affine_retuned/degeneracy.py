"""Full trace source and preservation of the original primary lapse null.

Primary degeneracy/velocity rank is not a complete nonlinear Dirac analysis.
The physical free-matter metric is not replaced by the chart metric.
"""
from functools import cache

import sympy as sp
from p8_affine import connection as old

from . import geometry

h = sp.Symbol("h", positive=True)
H0, T, KHAT = sp.symbols("H0 spatial_H_trace K_hat", real=True)


@cache
def source():
    p, s = old.P, old.S
    delta = 4*p**2-1
    actual = geometry.source_centering()["Tstar"]
    trace = sum(old.H[i, i] for i in range(1, 4))
    temporal = (6*old.F3*p*s**3+3*delta*old.PP*s/(2*p)
                -3*delta*old.PX*s*old.H[0, 0]/p-delta*trace/s)
    expected = sp.Matrix([temporal, 0, 0, 0])
    # Tstar_mu = phi_mu*(alpha + c*Box(phi) + d*vHv); b*H.v is zero.
    alpha = 6*old.F3*p*s**2+3*delta*old.PP/(2*p)
    c = -delta/s**2
    d = -delta/s**4-3*delta*old.PX/(p*s**2)
    gradient = sp.Matrix([s, 0, 0, 0])
    rebuilt = gradient*(alpha+c*old.literal()["trace_H"]+d*old.literal()["vHv"])
    return {"actual": actual, "expected": expected,
            "alpha": alpha, "c": c, "d": d,
            "generic_component_residual": geometry.clean(actual-expected),
            "covariant_reconstruction": geometry.clean(actual-rebuilt)}


@cache
def trace_kinetic():
    p, s = old.P, old.S
    x = -s**2
    px = -1/(8*h*p)
    cx = 2*(1-4*p)*px/x-2*(p-2*p**2)/x**2
    values = {old.PX: px, old.CX: cx, old.PP: 0, old.CP: 0, old.F3: 0}
    values.update({item: 0 for item in old.H_SYMBOLS})
    values.update({old.H[0, 0]: H0, old.H[1, 1]: T/3,
                   old.H[2, 2]: T/3, old.H[3, 3]: T/3})
    # Include the actual fR ADM time boundary, f=2p², f_X=-1/(2h).
    f, fx = 2*p**2, -1/(2*h)
    action = sp.factor(old.eliminate()["reduced_density"].subs(values)
                       -2*f*T**2/(3*s**2)-4*fx*H0*T)
    square = -3*H0*s**2+8*T*h*p**2
    target = -square**2/(48*h**2*p**2*s**2)
    matrix = geometry.clean(sp.hessian(action, (H0, T)))
    null = sp.Matrix([1, 3*s**2/(8*h*p**2)])
    tstar = sp.factor(source()["actual"][0].subs(values))
    return {"action": action, "expected": target, "Hessian": matrix, "null": null,
            "square": square, "Tstar_kinetic": tstar,
            "Tstar_expected": -(4*p**2-1)*square/(8*h*p**2*s),
            "null_action_residual": geometry.clean(matrix*null),
            "null_source_residual": sp.factor((sp.Matrix([[sp.diff(tstar, H0), sp.diff(tstar, T)]])*null)[0])}


@cache
def chart():
    """h_physical=(4p²)^(-1/2)*h_hat removes the lapse velocity exactly."""
    p, s = old.P, old.S
    omega_x, omega_phi = 1/(16*h*p**2), -old.PP/(2*p)
    trace = -s*KHAT-3*s**2*omega_phi+6*s**2*omega_x*H0
    values = {old.H[0, 0]: H0, old.H[1, 1]: trace/3,
              old.H[2, 2]: trace/3, old.H[3, 3]: trace/3,
              old.PX: -1/(8*h*p)}
    tstar = sp.factor(source()["actual"][0].subs(values))
    expected = (4*p**2-1)*KHAT+6*p*old.F3*s**3
    ratio = sp.factor(3*(4*p**2-1)**2/(8*geometry.update()["gamma_t"]*p**2))
    return {"Tstar_normal": tstar, "expected": expected,
            "omega_X": omega_x, "omega_phi": omega_phi,
            "trace_mass_relative_kinetic_addition": ratio,
            "uniform_relative_upper": sp.Rational(19, 1080),
            "nonlinear_secondary_constraint_rank_claim": False}


@cache
def velocity_rank():
    """Six metric, three vector and one matter velocity before constraints.

The single negative metric-trace direction is unreduced gravitational
inertia, not a physical ghost verdict. Lapse and vector temporal have
no time derivative in this chart; secondary constraint rank is separate.
"""
    p = old.P
    f, delta = 2*p**2, 4*p**2-1
    gt = geometry.update()["gamma_t"]
    zeta = sp.Symbol("positive_curl", positive=True)
    trace, sx, sy, xy, xz, yz = sp.symbols("Ktrace Ksx Ksy Kxy Kxz Kyz", real=True)
    vectors = sp.symbols("E1 E2 E3", real=True)
    matter = sp.Symbol("chi_velocity", real=True)
    K = sp.Matrix([[trace/3+sx, xy, xz], [xy, trace/3+sy, yz], [xz, yz, trace/3-sx-sy]])
    L = (f*(sp.trace(K*K)-sp.trace(K)**2)+delta**2*trace**2/(2*gt)
         +zeta*sum(value**2 for value in vectors)/2+matter**2/2)
    variables = (trace, sx, sy, xy, xz, yz, *vectors, matter)
    matrix = geometry.clean(sp.hessian(L, variables))
    expected = sp.zeros(10)
    expected[0, 0] = -4*f/3+delta**2/gt
    expected[1:3, 1:3] = 2*f*sp.Matrix([[2, 1], [1, 2]])
    for i in (3, 4, 5):
        expected[i, i] = 4*f
    for i in (6, 7, 8):
        expected[i, i] = zeta
    expected[9, 9] = 1
    return {"Hessian": matrix, "expected": expected, "zeta": zeta,
            "trace_pivot": sp.factor(matrix[0, 0]),
            "trace_pivot_expected": -4*f*(1-chart()["trace_mass_relative_kinetic_addition"])/3,
            "positive_block_determinant": sp.factor(matrix[1:, 1:].det()),
            "determinant_expected": 768*f**5*zeta**3,
            "primary_velocity_rank": 10, "physical_inertia_not_inferred": True}


@cache
def checks():
    data, kinetic, ch, rank = source(), trace_kinetic(), chart(), velocity_rank()
    return {"all_four_trace_source_components": data["generic_component_residual"],
            "full_covariant_source_reconstruction": data["covariant_reconstruction"],
            "full_original_trace_kinetic_square": sp.factor(kinetic["action"]-kinetic["expected"]),
            "original_primary_null_vector": kinetic["null_action_residual"],
            "source_along_original_primary_null": kinetic["null_source_residual"],
            "source_uses_original_kinetic_square": sp.factor(kinetic["Tstar_kinetic"]-kinetic["Tstar_expected"]),
            "full_spatial_chart_source": sp.factor(ch["Tstar_normal"]-ch["expected"]),
            "no_lapse_velocity_in_chart_source": sp.diff(ch["Tstar_normal"], H0),
            "full_ten_velocity_block": geometry.clean(rank["Hessian"]-rank["expected"]),
            "trace_pivot_in_terms_of_uniform_bound": sp.factor(rank["trace_pivot"]-rank["trace_pivot_expected"]),
            "positive_nine_velocity_block": sp.factor(rank["positive_block_determinant"]-rank["determinant_expected"])}
