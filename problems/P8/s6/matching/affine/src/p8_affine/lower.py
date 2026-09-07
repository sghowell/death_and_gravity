"""Lower-order integration-by-parts and a smooth clock-tube dictionary.

The corrected source formula is an explicit hypothesis of this module,
to be discharged by the separately derived literal connection action.
The arbitrary solution of the linear ODE is fixed by q(u,-1)=0.
"""
from functools import cache

import sympy as sp

from . import dictionary as d


def source_lower(p, px, p_phi, f_phi, cubic, variable, delta, scalar=0,
                 *, printed_Q2=False):
    mismatch = p_phi-cubic*variable
    bracket = p-3*px if printed_Q2 else p-3*variable*px
    return {"P": scalar+3*variable*mismatch**2/delta,
            "Q1": -2*f_phi+4*p*mismatch/delta,
            "Q2": 2*f_phi/variable-4*bracket*mismatch/(variable*delta)}


@cache
def ode():
    values = d.lift()
    h, w = values["h"], values["w"]
    variable = d.x
    coefficient = 1/(2*variable)+3/(4*w)
    forcing = -3*values["f_phi"]/(2*w)
    factor = sp.sqrt(-variable)/w**sp.Rational(3, 4)
    z = sp.Symbol("z", negative=True)
    # Integral on the fixed compact x interval; no integration at x=0.
    integrand = (factor*forcing).subs(variable, z)
    q = sp.Integral(integrand, (z, -1, variable))/factor
    return {"h": h, "w": w, "coefficient": coefficient, "forcing": forcing,
            "integrating_factor": factor, "q": q, "integrand": integrand,
            "basepoint": sp.Integer(-1), "integration_variable": z}


@cache
def identities():
    values, equation = d.lift(), ode()
    x, p, px, f_phi = d.x, values["p"], values["px"], values["f_phi"]
    q, qphi, target_F = sp.symbols("q qphi target_F", real=True)
    cubic = -(q+f_phi)/(4*p*x)
    parent_scalar = target_F+x*qphi-3*x*(q+2*f_phi)**2/(16*p**2)
    reduced = source_lower(p, px, values["p_phi"], f_phi, cubic, x, 1, parent_scalar)
    required_qx = equation["forcing"]-equation["coefficient"]*q
    factor = equation["integrating_factor"]
    # Differentiating a definite integral uses the fundamental theorem,
    # not a CAS guess for a special-function antiderivative.
    return {
        "integrating_factor_derivative": sp.simplify(sp.diff(factor, x)/factor-equation["coefficient"]),
        "ODE_from_general_p": sp.simplify(equation["coefficient"]-(1/(2*x)-3*px/(2*p))),
        "forcing_from_general_p": sp.simplify(equation["forcing"]-3*px*f_phi/p),
        "cubic_solves_Q1": sp.simplify(reduced["Q1"]-q),
        "ODE_solves_Q2": sp.simplify(reduced["Q2"]-2*required_qx),
        "scalar_after_boundary": sp.simplify(reduced["P"]-x*qphi-target_F),
        "clock_forcing_zero": sp.simplify(equation["forcing"].subs(x, -1)),
        "clock_p_phi_zero": sp.simplify(values["p_phi"].subs(x, -1)),
    }


def divergence_identities():
    """Differentiate a vector density in a normal frame, all ten H jets."""
    signs = (-1, 1, 1, 1)
    v = sp.symbols("v0:4", real=True)
    jets = {(a, b): sp.Symbol(f"H{a}{b}", real=True) for a in range(4) for b in range(a, 4)}
    hessian = sp.Matrix(4, 4, lambda a, b: jets[min(a, b), max(a, b)])
    q, qx, qphi = sp.symbols("q qx qphi", real=True)
    norm = sum(signs[a]*v[a]**2 for a in range(4))
    norm_derivative = [2*sum(signs[b]*v[b]*hessian[a, b] for b in range(4)) for a in range(4)]
    dq = [qphi*v[a]+qx*norm_derivative[a] for a in range(4)]
    divergence = sum(signs[a]*(dq[a]*v[a]+q*hessian[a, a]) for a in range(4))
    box = sum(signs[a]*hessian[a, a] for a in range(4))
    vhv = sum(signs[a]*signs[b]*v[a]*v[b]*hessian[a, b] for a in range(4) for b in range(4))
    return {"ten_Hessian_jet_vector_divergence": sp.expand(divergence-q*box-2*qx*vhv-qphi*norm),
            "integration_by_parts_scalar_sign": sp.expand(q*box+2*qx*vhv-divergence+qphi*norm)}


def palatini_control():
    """Compare the printed formula with the independently known pure-p control.

    A literal four-index action replay is separate. Here the direct
    square (3/(2p)) (p_phi v+2p_x H.v)^2 fixes the cross coefficient.
    """
    p, px, pphi, x = sp.symbols("p px pphi x", nonzero=True)
    corrected = source_lower(p, px, pphi, pphi, 0, x, 2*p)
    printed = source_lower(p, px, pphi, pphi, 0, x, 2*p, printed_Q2=True)
    cross = 6*pphi*px/p
    return {"corrected_cross_residual": sp.factor(corrected["Q2"]-cross),
            "printed_cross_residual": sp.factor(printed["Q2"]-cross),
            "constant_X_normalization_can_hide_typo": sp.factor((printed["Q2"]-cross).subs(x, 1)),
            "printed_at_nonunit": sp.factor((printed["Q2"]-cross).subs({p: 2, px: 3, pphi: 5, x: -2}))}
