"""An explicit higher-operator boundary ambiguity, not a physical invariant."""
from functools import cache

import sympy as sp


@cache
def curvature_boundary():
    """Divergence of W(phi) R X grad(phi), retaining every complement term.

    X=grad(phi)^2, T=Box(phi), V=phi^mu phi_munu phi^nu, and
    Rv=grad(phi).grad(R). Coefficient extraction alone is not invariant
    under a boundary redistribution involving the displayed complement.
    """
    W, W1, R, X, T, V, Rv = sp.symbols("W W1 R X T V Rv", real=True)
    return {
        "symbols": (W, W1, R, X, T, V, Rv),
        "divergence": W1*R*X**2 + W*X*Rv + W*R*X*T + 2*W*R*V,
        "Delta_F2": W1*X**2,
        "Delta_A1": sp.S.Zero,
        "Delta_Xi_in_named_projection": -4*W1*X,
        "complement": W*X*Rv + W*R*X*T + 2*W*R*V,
        "full_action_changed": False,
        "Xi_projection_is_full_higher_order_invariant": False,
    }


@cache
def tensor_boundary():
    """Literal one-dimensional Euler test with theta=t and R=-qdot^2/2.

    The physical metric is diag(1,-exp(q),-exp(-q),-1). Its determinant
    is constant; X=1, Box(theta)=0, and V=0. With W=t the remaining two
    terms are an exact total derivative. Dropping the W-proportional
    term at t=0 before variation would leave a false second-order term.
    """
    time = sp.Symbol("t", real=True)
    q = sp.Function("q")(time)
    R = -sp.diff(q, time)**2/2
    apparent = R
    complement = time*sp.diff(R, time)
    density = apparent + complement

    def euler(expr):
        return sp.simplify(sp.diff(expr, q)
                           - sp.diff(sp.diff(expr, sp.diff(q, time)), time)
                           + sp.diff(sp.diff(expr, sp.diff(q, time, 2)), time, 2))

    return {
        "time": time, "q": q, "R": R,
        "density": density,
        "total_derivative": sp.diff(time*R, time),
        "full_euler": euler(density),
        "apparent_euler": euler(apparent),
        "complement_euler": euler(complement),
        "complement_coefficient_at_center": sp.S.Zero,
    }


@cache
def actual_center_scaling():
    """S6.20 center jets in the prescribed theta=T clock, tau=1 units.

    a=MF2*r/beta1 is the relative-metric correction coefficient, twice
    the own-f inverse parameter in S6.30. At the symmetric center a'=r'=0.
    The displayed W=MF2*r^2*a*a'*s' has W=0 but W' generally nonzero.
    These are finite-parameter identities, not an action at c=2.
    """
    c, MF2 = sp.symbols("c MF2", positive=True)
    delta = c - 2
    r = sp.Integer(2)
    a = c*delta/16
    a2 = (13*c**2 - 22*c + 8)/8
    t = -4*(c + 2)/c
    W1 = sp.factor(MF2*r**2*a*a2*t)
    return {
        "symbols": (c, MF2), "relative_correction_a": a,
        "relative_correction_a_theta2": a2,
        "s_theta": t, "W_center": sp.S.Zero, "W_theta_center": W1,
        "W_theta_over_delta_limit": -8*MF2,
        "normalized_Xi_shift_at_X1": sp.factor(-4*W1/(5*MF2)),
        "center_projection_not_physical_remainder": True,
    }


@cache
def checks():
    b = curvature_boundary()
    W, W1, R, X, T, V, Rv = b["symbols"]
    e = tensor_boundary()
    a = actual_center_scaling()
    c, MF2 = a["symbols"]
    time = sp.Symbol("z", real=True)
    scale = sp.Function("a")(time)
    phi = sp.Function("phi")(time)
    scalar = sp.Function("R")(time)
    weight = sp.Function("W")
    p = sp.diff(phi, time)
    lhs = sp.diff(scale**3*weight(phi)*scalar*p**3, time)/scale**3
    rhs = b["divergence"].subs({
        W: weight(phi), W1: sp.diff(weight(phi), phi), R: scalar, X: p**2,
        T: sp.diff(p, time) + 3*sp.diff(scale, time)*p/scale,
        V: p**2*sp.diff(p, time), Rv: p*sp.diff(scalar, time)})
    return {
        "literal_FLRW_covariant_divergence": sp.simplify(lhs - rhs),
        "named_projection_and_complement_exhaust_divergence": sp.expand(
            b["divergence"] - R*b["Delta_F2"] - b["complement"]),
        "projected_Xi_shift": sp.expand(b["Delta_Xi_in_named_projection"]
                                        + 2*sp.diff(b["Delta_F2"], X)),
        "tensor_density_is_boundary": sp.simplify(e["density"] - e["total_derivative"]),
        "full_tensor_Euler_of_boundary_vanishes": e["full_euler"],
        "vanishing_center_coefficient_cancels_apparent_Euler": sp.simplify(
            e["apparent_euler"] + e["complement_euler"]),
        "center_projection_alone_is_not_variationally_zero": sp.simplify(
            e["apparent_euler"] - sp.diff(e["q"], e["time"], 2)),
        "actual_center_W_derivative": sp.factor(a["W_theta_center"]
            + MF2*(c - 2)*(c + 2)*(13*c**2 - 22*c + 8)/8),
        "actual_boundary_shift_can_start_at_delta": sp.limit(
            a["W_theta_center"]/(c - 2), c, 2, dir="+") + 8*MF2,
    }


def calibration():
    b = curvature_boundary()
    e = tensor_boundary()
    return {
        "covariant_boundary": {key: value for key, value in b.items() if key != "symbols"},
        "literal_tensor_Euler_control": {
            "full_euler": e["full_euler"], "apparent_euler": e["apparent_euler"],
            "complement_euler": e["complement_euler"],
            "complement_coefficient_at_center": e["complement_coefficient_at_center"]},
        "actual_center": {key: value for key, value in actual_center_scaling().items()
                          if key != "symbols"},
        "full_S6_action_or_controlled_remainder_computed": False,
    }
