"""Exact inputs to the reciprocal-scale integrating-factor argument.

The all-time theorem is the written argument, not a finite point scan.
No G1 reconstruction, scalar-health verdict or EFT cutoff is imported.
"""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_variable_beta import model, tensors


def affine_tail_bound(alpha, a0, z0, momentum=1):
    """Exact future null affine bound for a fixed conserved momentum.

    The caller must supply values from a solution satisfying the written
    theorem. Positive numbers alone do not certify a background solution.
    """
    values = (alpha, a0, z0, momentum)
    if any(isinstance(value, bool) or not isinstance(value, (int, Fraction)) for value in values):
        raise TypeError("Require positive exact rational coefficients and momentum magnitude")
    if any(value <= 0 for value in values):
        raise ValueError("The tail starts at u0>0 with a0,z0,alpha and momentum strictly positive")
    alpha, a0, z0, momentum = map(Fraction, values)
    return alpha/(a0*z0*momentum)


@cache
def derive():
    a, alpha, h, z, z0, p = sp.symbols("a alpha h z z0 P", positive=True)
    hp = sp.Symbol("h_prime", real=True)
    y = alpha/a**2
    c = h/z
    zp = (c-y)*p/(2*y**3)

    def derivative(value):
        return a*h*sp.diff(value, a)+hp*sp.diff(value, h)+zp*sp.diff(value, z)

    D = h/y-z
    B = (hp+2*h*h)/y
    F = p/(2*y*y*z)
    primitive = alpha/(z0*a)
    source = model.derive()
    shift = tensors.derive()
    u = sp.Symbol("u", real=True)
    scale = (1+u*u)**2
    cd_y = 2/scale**2
    cd_h = sp.diff(scale, u)/scale
    residuals = {
        "literal_f_null_source": sp.factor(source["rho_f"]+source["pressure_f"]
            -(source["c"]-source["y"])*source["P"]/(source["c"]*source["y"]**3)),
        "literal_relative_shift_sign_prefactor": sp.factor(shift["literal_shift_coefficient"]
                                                          -shift["shift_coefficient"]),
        "reciprocal_scale_derivative": sp.factor(derivative(alpha/a)+h*alpha/a),
        "lapse_order_identity": sp.factor(c-y-y*D/z),
        "integrating_factor_equation": sp.factor(derivative(D)+F*D-B),
        "affine_primitive_derivative": sp.factor(derivative(primitive)+alpha*h/(z0*a)),
        "affine_bound_remainder": sp.factor(alpha*c/a+derivative(primitive)
                                             -alpha*h*(z0-z)/(a*z*z0)),
        "CD_geometric_sign_polynomial": sp.factor(sp.diff(cd_h/cd_y, u)
                                                   -2*(1+u*u)**2*(1+7*u*u)),
    }
    return {"residuals": residuals, "D": D, "B": B, "F": F,
            "affine_tail_bound": primitive,
            "variables": (a, alpha, h, hp, z, z0, p),
            "scope": "Specified reciprocal two-metric domain; not a physical-g no-go or cutoff"}


def checks():
    residuals = derive()["residuals"]
    if any(value != 0 for value in residuals.values()):
        raise ValueError("An exact reciprocal geometric input failed")
    return residuals
