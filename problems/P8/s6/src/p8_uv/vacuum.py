"""Independently expanded massive vacuum scalar amplitude and matching control."""

from functools import cache

import sympy as sp

s, transfer, w, v = sp.symbols("s t_transfer w v", real=True)
mass, heavy_mass = sp.symbols("m M_h", positive=True)
coupling, mediator = sp.symbols("lambda g", real=True)


def labelled_contact(dot_products, coefficient=coupling):
    """Expand lambda*X_psi^2 with partial psi=i sum_j k_j e_j."""
    if sp.shape(dot_products) != (4, 4) or dot_products != dot_products.T:
        raise ValueError("Expected a symmetric four-leg dot-product matrix")
    labels = sp.symbols("e0:4")
    kinetic = -sum(dot_products[i, j]*labels[i]*labels[j]
                   for i in range(4) for j in range(4))
    return sp.expand(coefficient*kinetic**2).coeff(labels[0]).coeff(labels[1]).coeff(labels[2]).coeff(labels[3])


@cache
def contact():
    ds, dt, dw = [(z-2*mass**2)/2 for z in (s, transfer, w)]
    dots = sp.Matrix([[mass**2, ds, dt, dw], [ds, mass**2, dw, dt],
                      [dt, dw, mass**2, ds], [dw, dt, ds, mass**2]])
    return sp.factor(labelled_contact(dots))


def forward_coefficient(amplitude):
    """b2=(1/2)d_v^2 A at t=0, s=w=2m^2, with on-shell crossing."""
    crossed = amplitude.subs({s: 2*mass**2+v-transfer/2,
                              w: 2*mass**2-v-transfer/2}, simultaneous=True)
    return sp.factor(sp.diff(crossed, v, 2).subs({v: 0, transfer: 0})/2)


@cache
def heavy_exchange():
    return mediator**2/4*sum((z-2*mass**2)**2/(heavy_mass**2-z) for z in (s, transfer, w))


@cache
def derive():
    amplitude = contact()
    expected = 2*coupling*sum((z-2*mass**2)**2 for z in (s, transfer, w))
    exchange = heavy_exchange()
    b2_heavy = forward_coefficient(exchange)
    lam_eff = mediator**2/(8*heavy_mass**2)
    h, kinetic = sp.symbols("h Y", real=True)
    h_lagrangian = -heavy_mass**2*h**2/2+mediator*h*kinetic/2
    solution = sp.solve(sp.diff(h_lagrangian, h), h)[0]
    residuals = {"labelled_vs_pairing": sp.expand(amplitude-expected),
                 "s_t_crossing": sp.expand(amplitude-amplitude.xreplace({s: transfer, transfer: s})),
                 "s_w_crossing": sp.expand(amplitude-amplitude.xreplace({s: w, w: s})),
                 "forward_b2": forward_coefficient(amplitude)-4*coupling,
                 "heavy_stationary_matching": sp.expand(h_lagrangian.subs(h, solution)-lam_eff*kinetic**2),
                 "heavy_exact_b2": sp.cancel(b2_heavy-mediator**2/(2*(heavy_mass**2-2*mass**2))),
                 "heavy_leading_matching": sp.limit(heavy_mass**2*b2_heavy, heavy_mass, sp.oo)-mediator**2/2}
    if any(value != 0 for value in residuals.values()):
        raise ValueError("Vacuum amplitude or heavy-field control failed")
    return {"lagrangian": "Y/2-m^2*psi^2/2+lambda*Y^2; Y=(partial psi)^2",
            "amplitude": str(amplitude), "on_shell_sum": "s+t_transfer+w=4*m^2",
            "b2_definition": "(1/2) partial_v^2 A(v,0) at v=0; v=s+t_transfer/2-2*m^2",
            "b2": str(forward_coefficient(amplitude)),
            "heavy_field_tree_regression": {"amplitude": str(exchange), "b2_exact": str(b2_heavy),
                                            "lambda_leading": str(lam_eff), "domain": "M_h>2m>0",
                                            "warning": "derivative-coupled massive mediator, not a UV completion"},
            "exact_residuals": {key: "0" for key in residuals},
            "diagnostic_members": {"lambda>0": "POSITIVE_LEADING_SCALAR_COEFFICIENT",
                                   "lambda<0": "NEGATIVE_LEADING_SCALAR_COEFFICIENT",
                                   "lambda=0": "TREE_SATURATION_NOT_A_STANDALONE_EXCLUSION"},
            "scope": "Vacuum scalar decoupling patch only; no finite-gravity or bounce positivity verdict"}
