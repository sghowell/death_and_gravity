"""Exact upper-p-cut density and independent two-by-two Gram identities."""
from functools import cache

import sympy as sp
from p8_vector_metric_dispersion import spectral

from . import block

h, z = block.h, block.z
ell = sp.Symbol("cut_atanh_sqrt_z", real=True)


@cache
def data():
    item = block.data()
    cs = tuple(spectral.actual(value) for value in spectral.coefficients())
    polynomial = cs[0]+cs[1]*(z+sp.Rational(1,3))+cs[2]*(z**2+z/3+sp.Rational(1,5))
    D = block.clean(item["local_matrix"]+polynomial/4-sp.sqrt(z)*item["spectral_matrix"]*ell/4)
    W = block.clean(sp.sqrt(z)*item["spectral_matrix"]/8)
    X = sp.factor(D.det()-sp.pi**2*W.det())
    Y = sp.factor(D[0,0]*W[1,1]+D[1,1]*W[0,0]-2*D[0,1]*W[0,1])
    adjD, adjW = block.clean(D.adjugate()), block.clean(W.adjugate())
    return {"real_cut_matrix": D, "positive_cut_weight": W,
            "substitution": {str(ell): "atanh(sqrt(momentum_fraction)) for 0<z<1"},
            "determinant_real_part": X, "minus_determinant_imaginary_part_over_pi": Y,
            "determinant_modulus_squared": X**2+sp.pi**2*Y**2,
            "density_numerator": X*(-adjW)+Y*adjD,
            "positive_Gram_numerator": adjD*W*adjD+sp.pi**2*W.det()*adjW}


@cache
def generic_checks():
    a,b,c,d,e,f = sp.symbols("cut_a cut_b cut_c cut_d cut_e cut_f", real=True)
    D,W = sp.Matrix([[a,b],[b,c]]),sp.Matrix([[d,e],[e,f]])
    B = D-sp.I*sp.pi*W
    X,Y = D.det()-sp.pi**2*W.det(),a*f+c*d-2*b*e
    numerator = Y*D.adjugate()-X*W.adjugate()
    gram = D.adjugate()*W*D.adjugate()+sp.pi**2*W.det()*W.adjugate()
    out = {"cut_determinant_real_part": sp.factor(sp.re(B.det())-X),
           "cut_determinant_imaginary_sign": sp.factor(sp.im(B.det())+sp.pi*Y),
           "adjugate_is_linear_in_two_dimensions": block.clean(B.adjugate()-D.adjugate()+sp.I*sp.pi*W.adjugate()),
           "inverse_positive_spectral_Gram_numerator": block.clean(numerator-gram),
           "inverse_density_imaginary_part": block.clean(
               ((D.adjugate()-sp.I*sp.pi*W.adjugate())*(X+sp.I*sp.pi*Y)).applyfunc(sp.im)/sp.pi-numerator),
           "inverse_cut_Gram_cross_terms_cancel": block.clean(
               D.adjugate()*W*W.adjugate()-W.adjugate()*W*D.adjugate())}
    return out


@cache
def checks():
    item, parent = data(), block.data()
    W,D = item["positive_cut_weight"],item["real_cut_matrix"]
    out = generic_checks().copy()
    out["strict_cut_weight_determinant"] = sp.factor(W.det()-8*z**3*(1-z)**2/(6561*h**2))
    out["cut_real_threshold_matches_direct_endpoint"] = block.clean(D.subs({z:0,ell:0})-parent["threshold_matrix"])
    # Polynomial division on the two banks: d=z-i0, Im I_n=+pi*z^(n-1/2)/2.
    cs = tuple(spectral.actual(value) for value in spectral.coefficients())
    jump = sum((value*z**(j+sp.Rational(1,2))/8 for j,value in enumerate(cs)),sp.zeros(2))
    out["three_radial_moments_have_the_same_cut_weight"] = block.clean(jump-W)
    out["cut_Gram_fixture_strictly_positive"] = sp.factor(
        W.det().subs({z:sp.Rational(1,2),h:1})-sp.Rational(1,26244))
    return out
