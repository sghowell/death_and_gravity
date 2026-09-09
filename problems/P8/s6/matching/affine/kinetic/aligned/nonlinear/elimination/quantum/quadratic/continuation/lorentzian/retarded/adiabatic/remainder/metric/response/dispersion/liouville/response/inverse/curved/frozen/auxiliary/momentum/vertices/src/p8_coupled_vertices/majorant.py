"""Continuous-time positive majorants for all unnormalized physical H3/H4."""
from fractions import Fraction
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_momentum import constraints
from p8_m1_tree.majorant import Series as _Series

from . import coefficients


@cache
def coefficient_bounds():
    rows={}
    maxima=[sp.Integer(0)]*5
    for powers,values in coefficients.jets().items():
        bounds=[]
        proofs=[]
        for order,value in enumerate(values):
            numerator,denominator=sp.fraction(sp.cancel(value))
            num,den=sp.Poly(numerator,model.u),sp.Poly(denominator,model.u)
            if den.nth(0)<=0 or any(coefficient<0 or degree[0]%2 for degree,coefficient in den.terms()):
                raise ValueError("The actual lapse-jet denominator lacks this positive-even certificate")
            upper=sum(abs(coefficient)*sp.Rational(1,2)**degree[0] for degree,coefficient in num.terms())
            exact=upper/den.nth(0)
            cauchy=10**4*sp.factorial(order)*100**order
            bound=min(exact,cauchy)
            bounds.append(bound)
            proofs.append({"numerator_absolute_half_interval_sum":upper,
                           "denominator_positive_even_constant_lower":den.nth(0),
                           "independent_S6_75_Cauchy_bound":cauchy})
            maxima[order]=max(maxima[order],bound)
        rows[powers]={"derivative_bounds":tuple(bounds),"proof_data":tuple(proofs)}
    return {"rows":rows,"common_N_derivative_bounds":tuple(maxima),
            "time_interval":(sp.Rational(-1,2),sp.Rational(1,2))}


def geometry(field,K):
    perturbation=3*field
    metric=1+perturbation
    inverse=1+sum(3**(n-1)*perturbation**n for n in range(1,5))
    determinant_remainder=(6*metric**3).from_degree(1)
    volume=determinant_remainder.analytic(sp.Rational(1,2))
    inverse_volume=determinant_remainder.analytic(sp.Rational(-1,2))
    connection=Fraction(9,2)*K*inverse*perturbation
    curvature=9*inverse*(6*K*connection+18*connection**2)
    return {"metric":metric,"inverse":inverse,"volume":volume,
            "inverse_volume":inverse_volume,"connection":connection,"curvature":curvature}


def bound(seed,derivative_bound,inverse_transfer_bound,*,include_vector=True):
    """Every individual linear field/momentum component <=seed per labelled leg.

    K bounds each spatial derivative on every used Fourier mask.
    inverse_transfer_bound bounds 1/|k_mask| on all nonempty proper masks.
    This is an unnormalized vertex estimate, not a normalized free-mode bound.
    """
    args=tuple(constraints.exact(value) for value in (seed,derivative_bound,inverse_transfer_bound))
    if any(value<=0 for value in args) or type(include_vector) is not bool:
        raise ValueError("Require positive exact phase, derivative and inverse-transfer bounds and native flag")
    seed,K,invk=map(Fraction,args)
    f=_Series.at(1,4*seed)
    geo=geometry(f,K)
    # |H|<=2 and |ell|<=1/10 on I. Each free momentum entry
    # is <=2+(2+1/6+1+2)*f < 2+6*f.
    pi=2+6*f
    matter=Fraction(1,10)+f
    corrections=[]
    for degree in range(1,4):
        source=3*K*pi+9*geo["connection"]*pi
        current=K*matter*f
        if include_vector:
            # Pi^j F_ij contributes <=6*K*f^2 and W_i div Pi <=3*K*f^2.
            current+=9*K*f*f
        source+=Fraction(3,2)*geo["inverse"]*current
        # Full tensor lift inverse: sqrt(2)/|k| times ||source||_2
        # <=sqrt(6)*entry_bound/|k| <3*entry_bound/|k|.
        correction=3*invk*source.coefficients[degree]
        pi+=_Series.at(degree,correction)
        corrections.append(correction)
    mixed=3*pi*geo["metric"]*geo["inverse_volume"]
    dp=(2*mixed).from_degree(1)
    dc=(matter*geo["inverse_volume"]).from_degree(1)
    divergence=3*K*f*geo["inverse_volume"]
    shear=(12*mixed*mixed).from_degree(2)
    zeta=Fraction(1,10**6)
    electric=(9/zeta)*geo["metric"]*f*f*geo["inverse_volume"]**2
    magnetic=324*zeta*K*K*geo["inverse"]**2*f*f
    vector=9*geo["inverse"]*f*f
    gradient=9*K*K*geo["inverse"]*f*f
    invariants=(dp,dc,divergence,shear,electric,magnetic,vector,gradient,geo["curvature"])
    linear=sum(invariants,_Series.at(0,0))
    quadratic=dp*dp+dp*divergence+divergence*divergence+dc*dc
    B=tuple(Fraction(value) for value in coefficient_bounds()["common_N_derivative_bounds"])
    L=tuple(value*linear for value in B[:4])
    Q=tuple(value*quadratic for value in B[:3])
    n1=20*L[1]
    force2=(B[3]/2)*n1*n1+L[2]*n1+Q[1]
    h=B[0]+L[0]+Q[0]+10*L[1]*L[1]
    h+=(B[3]/6)*n1**3+Fraction(1,2)*L[2]*n1*n1+Q[1]*n1
    h+=(B[4]/24)*n1**4+Fraction(1,6)*L[3]*n1**3+Fraction(1,2)*Q[2]*n1*n1+10*force2*force2
    # Keep a majorant even for moving-boundary terms whose integrated
    # higher coefficients vanish by flat TT orthogonality.
    h=geo["volume"]*h+36*pi*geo["metric"]+96*f+120*f*f+Fraction(1,10)*f
    return {"phase_component_seed":seed,"derivative_bound":K,"inverse_transfer_bound":invk,
            "four_labelled_leg_degree_seed":f,"geometry":geo,
            "momentum":pi,"York_correction_majorants":tuple(corrections),
            "invariants":invariants,"actual_N_derivative_bounds":B,
            "n1_majorant":n1,"second_lapse_force_majorant":force2,
            "physical_Hamiltonian_degree_majorant":h,
            "cubic_labelled_kernel_upper":h.coefficients[3],
            "quartic_labelled_kernel_upper":h.coefficients[4],
            "normalization":"Raw canonical phase components in fixed local units; no M*tau factors or free mode functions included."}


def serialize_bound(data):
    def convert(value):
        if isinstance(value,_Series):
            return tuple(map(str,value.coefficients))
        if isinstance(value,Fraction):
            return str(value)
        if isinstance(value,dict):
            return {str(key):convert(v) for key,v in value.items()}
        if isinstance(value,(tuple,list)):
            return tuple(convert(v) for v in value)
        return value
    return convert(data)
