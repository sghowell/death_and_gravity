"""Full clock-lapse Taylor remainder from the actual lower-coefficient ODE."""
from functools import cache

import sympy as sp
from p8_affine_aligned import source_bounds

r, h = sp.Symbol("r", real=True), sp.Symbol("h", positive=True)
H, k, f1, RQ, RQx = sp.symbols("H k f1 R_Q R_Qx", real=True)


@cache
def ode():
    data = source_bounds.ode_bridge()
    x = source_bounds.dictionary.x
    actual_h = source_bounds.dictionary.lift()["h"]
    slope = sp.factor(sp.diff(data["forcing"], x).subs(x, -1))
    y = x+1
    difference = sp.factor(data["forcing"]-slope*y)
    return {"slope": slope, "forcing_difference": difference,
            "actual_forcing_remainder_identity": sp.factor(difference-slope*y**2/(actual_h-y)),
            "slope_H_bridge": sp.factor(slope+9*source_bounds.alignment.parent.old.background()["H"]/(8*actual_h**2))}


@cache
def expansion():
    N = 1+r
    s = 1/N
    y = 1-s**2
    Q = f1*y**2/2+RQ
    full = -y*k/h-3*H*y*(1-s)/h+sp.Rational(3, 2)*s*Q
    coordinate = N*full
    quadratic = -2*r*k/h-6*H*r**2/h+3*f1*r**2
    normal_remainder = (r**2*(3+2*r)*k/(h*N**2)
                        +3*H*r**3*(5+6*r+2*r**2)/(h*N**3)
                        -sp.Rational(3, 4)*f1*r**3*(16+39*r+40*r**2+20*r**3+4*r**4)/N**5
                        +sp.Rational(3, 2)*RQ/N)
    coordinate_remainder = (r**2*k/(h*N)+3*H*r**3*(3+2*r)/(h*N**2)
                            -sp.Rational(3, 4)*f1*r**3*(12+23*r+16*r**2+4*r**3)/N**4
                            +sp.Rational(3, 2)*RQ)
    # R_Qx is dR_Q/dx. Spatial phi derivatives vanish in the clock chart.
    derivative = sp.diff(coordinate_remainder, r)+sp.Rational(3, 2)*RQx*sp.diff(y, r)
    derivative_expected = (y*k/h+3*H*r**2*(9+11*r+4*r**2)/(h*N**3)
                           -sp.Rational(3, 2)*f1*r**2*(18+40*r+40*r**2+20*r**3+4*r**4)/N**5
                           +3*RQx/N**3)
    return {"N": N, "y": y, "normal": full, "coordinate": coordinate, "quadratic": quadratic,
            "normal_remainder": normal_remainder, "coordinate_remainder": coordinate_remainder,
            "coordinate_r_remainder": derivative_expected, "coordinate_k_remainder": r**2/(h*N),
            "normal_full_Taylor_identity": sp.factor(full-quadratic-normal_remainder),
            "coordinate_full_Taylor_identity": sp.factor(coordinate-quadratic-coordinate_remainder),
            "coordinate_r_derivative_identity": sp.factor(derivative-derivative_expected),
            "coordinate_k_derivative_identity": sp.factor(sp.diff(coordinate_remainder, k)-r**2/(h*N)),
            "physical_quadratic_source_bridge": sp.factor(quadratic.subs(f1, -9*H/(8*h**2))
                                                          +2*r*k/h+(6*H/h+27*H/(8*h**2))*r**2)}


@cache
def checks():
    data, equation = expansion(), ode()
    out = {name: data[name] for name in ("normal_full_Taylor_identity", "coordinate_full_Taylor_identity",
                                         "coordinate_r_derivative_identity", "coordinate_k_derivative_identity",
                                         "physical_quadratic_source_bridge")}
    out.update({name: equation[name] for name in ("actual_forcing_remainder_identity", "slope_H_bridge")})
    return out


@cache
def preparation():
    t, a3, k0, k1, h1 = sp.symbols("t a3 k0 k1 h1", real=True)
    data = expansion()
    # r and its first two time jets vanish. General smooth coefficient
    # variations cannot enter the cubic time coefficient multiplying r.
    substitution = {r: a3*t**3/6, k: k0+k1*t, h: h+h1*t, RQ: 0}
    series = sp.series(data["normal"].subs(substitution, simultaneous=True), t, 0, 4).removeO()
    expected = -a3*k0*t**3/(3*h)
    return {"prepared_full_source_jet_through_three": sp.expand(series-expected),
            "third_source_jet": -2*a3*k0/h,
            "time_variable": t, "lapse_third": a3, "trace_initial": k0}
