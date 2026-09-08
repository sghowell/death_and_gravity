"""Physical high-frequency characteristic cone after all rolling constraints.

The original two-scalar cone is exactly luminal. A nonzero temporal
source mixes it with a healthy longitudinal vector and widens the cone.
This is a principal limit, not the massive phase speed at finite q.
"""
from functools import cache

import sympy as sp

from . import dynamics as d
from . import geometry


@cache
def characteristic():
    data = d.principal()
    K = geometry.clean(data["K"].applyfunc(lambda value: sp.limit(value, d.q, sp.oo)))
    G = geometry.clean(data["G"].applyfunc(lambda value: sp.limit(value, d.q, sp.oo)))
    speed = sp.Symbol("speed_squared", real=True)
    bg = d.old.background()
    kappa = sp.factor(d.source()["d"]**2/(2*bg["J"]))
    polynomial = sp.factor((G-speed*K).det()/K.det())
    expected = (1-speed)*((1-speed)**2-kappa*speed)
    minor = geometry.clean((G-K).extract((0, 2), (0, 2)))
    return {"K": K, "G": G, "speed_squared": speed, "kappa": kappa,
            "polynomial": polynomial, "expected": expected, "matter_cone_minor": minor,
            "speed_plus": 1+kappa/2+sp.sqrt(kappa**2+4*kappa)/2,
            "speed_minus": 1+kappa/2-sp.sqrt(kappa**2+4*kappa)/2,
            "new_curl_absent_from_principal_speeds": not kappa.has(d.ZETA),
            "zero_curl_is_a_separate_rank_limit": True,
            "violation_below_an_EFT_cutoff_claim": False}


@cache
def general_isotropic_mass():
    """The same front-cone mechanism for any positive retained gamma.

No general mass-retuning/open-tube theorem is inferred from this
quadratic lemma; the literal certified action keeps gamma=1.
"""
    gamma, J, theta = sp.symbols("gamma J theta", positive=True)
    w, f, speed = sp.symbols("w f speed_squared", real=True)
    K = sp.Matrix([[(J+w**2/2)/theta**2, w/(2*theta), 0],
                   [w/(2*theta), sp.Rational(1, 2), 0],
                   [0, 0, 1/(2*gamma)]])
    mass = sp.Matrix([-f, 0, 1])
    G = sp.zeros(3)
    G[:2, :2] = K[:2, :2]
    G += mass*mass.T/(2*gamma)
    kappa = f**2*theta**2/(2*gamma*J)
    expected = (1-speed)*((1-speed)**2-kappa*speed)
    return {"residual": sp.factor((G-speed*K).det()/K.det()-expected),
            "minor_residual": sp.factor((G-K).extract((0, 2), (0, 2)).det()+f**2/(4*gamma**2))}


@cache
def checks():
    data, general = characteristic(), general_isotropic_mass()
    kappa = data["kappa"]
    return {"original_two_scalar_cone_exactly_luminal": geometry.clean(d.principal()["G0"]-data["K"][:2, :2]),
            "complete_three_scalar_characteristic_polynomial": sp.factor(data["polynomial"]-data["expected"]),
            "negative_matter_cone_minor": sp.factor(data["matter_cone_minor"].det()+d.source()["f"]**2/4),
            "front_speed_product": sp.factor(data["speed_plus"]*data["speed_minus"]-1),
            "front_speed_sum": sp.factor(data["speed_plus"]+data["speed_minus"]-2-kappa),
            "general_positive_retained_mass_characteristic": general["residual"],
            "general_positive_retained_mass_minor": general["minor_residual"]}
