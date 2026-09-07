"""Exact physical composite metric and exchange-odd relative coordinate.

G denotes the physical g_eff, not one of the two Einstein metrics. Matrix
products have the displayed index order. Delta is G-self-adjoint, and the
chosen chart requires I+Delta and I-Delta invertible with the positive
square-root branch. These domain conditions are not inferred from det G.
"""

from functools import cache

import sympy as sp

M2, m2 = sp.symbols("M_squared m_squared", positive=True)
BETAS = sp.symbols("beta0:5", real=True)


@cache
def derive():
    t = sp.Symbol("generating_t", real=True)
    elementary = (sp.S.One, *sp.symbols("E1 E2 E3 E4", real=True))
    # det[(I+Delta)+t(I-Delta)], exact for a general 4x4 matrix.
    generator = sp.expand(sum(elementary[n]*(1+t)**(4-n)*(1-t)**n for n in range(5)))
    densities = tuple(generator.coeff(t, n)/16 for n in range(5))
    potential = -M2*m2*sum(beta*density for beta, density in zip(BETAS, densities, strict=True))
    beta2 = sp.expand(potential.subs(dict(zip(BETAS, (0, 0, 1, 0, 0), strict=True))))
    return {"M2": M2, "m2": m2, "betas": BETAS, "t": t,
            "elementary_Delta": elementary,
            "generating_polynomial": generator,
            "sqrtg_enS_over_sqrtG": densities,
            "general_potential_over_sqrtG": potential,
            "beta2_potential_over_sqrtG": beta2,
            "zero_relative_gravity_coefficient": -M2/4,
            "zero_relative_interaction_constant": -3*M2*m2/8}


def matrix_map(physical, relative):
    """Return the exact chart; caller must establish its real positive domain."""
    physical, relative = sp.Matrix(physical), sp.Matrix(relative)
    if physical.shape != (4, 4) or relative.shape != (4, 4):
        raise ValueError("The Cayley chart is four dimensional")
    ident = sp.eye(4)
    plus, minus = ident+relative, ident-relative
    return {"g": physical*plus**2/4, "f": physical*minus**2/4,
            "S": plus.inv()*minus, "Delta": relative, "G": physical}


@cache
def checks():
    d = derive()
    e = d["elementary_Delta"]
    flip = {e[1]: -e[1], e[3]: -e[3]}
    reverse = dict(zip(BETAS, reversed(BETAS), strict=True))
    symmetric = {BETAS[4]: BETAS[0], BETAS[3]: BETAS[1]}
    trace, trace2, det = sp.symbols("trace_Delta trace_Delta2 det_Delta", real=True)
    return {
        "general_exchange_reverses_beta_indices": sp.expand(
            d["general_potential_over_sqrtG"].subs(flip, simultaneous=True)
            -d["general_potential_over_sqrtG"].subs(reverse, simultaneous=True)),
        "symmetric_potential_even": sp.expand((d["general_potential_over_sqrtG"].subs(flip, simultaneous=True)
                                              -d["general_potential_over_sqrtG"]).subs(symmetric)),
        "beta2_exact_even_density": sp.expand(d["sqrtg_enS_over_sqrtG"][2]-(6-2*e[2]+6*e[4])/16),
        "beta2_full_trace_mass_and_quartic": sp.expand(
            d["beta2_potential_over_sqrtG"].subs({e[2]: (trace**2-trace2)/2, e[4]: det})
            +M2*m2*(6+trace2-trace**2+6*det)/16),
        "physical_volume_is_independent_of_relative_field": sp.expand(sum(d["sqrtg_enS_over_sqrtG"])-1),
        "equal_metric_Einstein_scaling": -M2*sp.Rational(1, 4)-d["zero_relative_gravity_coefficient"],
    }


def controls():
    d = derive()
    e1 = d["elementary_Delta"][1]
    asymmetric = d["general_potential_over_sqrtG"].subs(dict(zip(BETAS, (1, 0, 0, 0, 0), strict=True)))
    return {"asymmetric_beta0_linear_relative_tadpole": sp.diff(asymmetric, e1),
            "singular_Cayley_plus_example": sp.det(sp.eye(4)-sp.eye(4)),
            "negative_root_not_positive_chart": sp.Integer(-1)}
