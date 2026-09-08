"""Finite subtracted mass-insertion loop on a specified complex momentum ball."""
from functools import cache

import sympy as sp
from p8_affine_retuned import bounds as exact_parent

from . import pole


@cache
def dimensional_logarithms():
    """Actual radial Laurent coefficients, not just their pole residues.

    Integrals are divided by 1/(16*pi^2); ell=log(Delta/mu_bar^2).
    The d-dimensional angular coefficients are analytic at epsilon=0,
    so their epsilon corrections multiply poles into local polynomials.
    """
    epsilon = sp.Symbol("epsilon_DR", positive=True)
    Delta = sp.Symbol("Delta", positive=True)
    ell = sp.Symbol("log_Delta_over_mu_squared", real=True)
    data, identities = {}, {}
    expected = ((1, -ell), (-2*Delta, Delta*(2*ell-1)),
                (3*Delta**2, Delta**2*(-3*ell+2)))
    for order, (residue, finite) in enumerate(expected):
        ratio = sp.prod(j-epsilon for j in range(2, order+2))
        integral = (Delta**order*sp.exp(epsilon*(sp.EulerGamma-ell))
                    *ratio*sp.gamma(epsilon-order))
        actual_residue = sp.limit(epsilon*integral, epsilon, 0)
        actual_finite = sp.limit(integral-residue/epsilon, epsilon, 0)
        data[order] = {"pole": actual_residue, "finite": actual_finite}
        identities["dimensional_radial_pole_"+str(order)] = sp.expand(actual_residue-residue)
        identities["dimensional_radial_finite_"+str(order)] = sp.expand(actual_finite-finite)
        identities["finite_log_coefficient_is_minus_pole_"+str(order)] = sp.expand(sp.diff(actual_finite, ell)+residue)
    return {"integrals": data, "checks": identities}


@cache
def logarithmic_form():
    t = pole.x*(1-pole.x)
    # F0,F2,F4 are homogeneous of external-momentum degrees 0,2,4.
    F0 = pole.m2**2*(pole.TA**2+2*pole.TA2)/8
    F2 = pole.m2*(pole.PA2/2+t*(pole.P2*(pole.TA**2-2*pole.TA2)/4+pole.TA*pole.PA))
    F4 = t**2*(pole.P2**2*(pole.TA**2+2*pole.TA2)/8+2*pole.P2*pole.PA2
               +pole.P2*pole.TA*pole.PA+pole.PA**2)-t*pole.P2*pole.PA2/2
    z = sp.Symbol("momentum_scaling_squared", real=True)
    v = t*pole.P2/pole.m2
    original = (F0+z*F2+z**2*F4)*sp.log(1+z*v)
    local = sp.series(original, z, 0, 3).removeO()
    subtracted = (F0*(sp.log(1+z*v)-z*v+z**2*v**2/2)
                  +z*F2*(sp.log(1+z*v)-z*v)+z**2*F4*sp.log(1+z*v))
    return {"F0": F0, "F2": F2, "F4": F4, "subtracted": subtracted,
            "local_Taylor_polynomial": local,
            "homogeneous_parameter_integrand": sp.expand(F0+F2+F4-pole.parameter_integral()["integrand"]),
            "exact_through_four_derivatives_subtraction": sp.expand(original-local-subtracted)}


@cache
def constants():
    t = pole.x*(1-pole.x)
    f0 = sp.Rational(904, 2187)
    f2 = sp.Rational(8, 81)+sp.Rational(2032, 2187)*t
    f4_absolute = sp.Rational(3640, 2187)*t**2+sp.Rational(8, 81)*t
    actual = sp.integrate(f0*t**3/3+f2*t**2/2+f4_absolute*t, (pole.x, 0, 1))
    return {"f0_envelope": f0, "f2_envelope": f2, "f4_absolute_envelope": f4_absolute,
            "integrated_absolute_remainder_constant": sp.factor(actual)}


@cache
def proof_checks():
    c = constants()
    upper4 = sp.Rational(32*120, 98415)
    lower4 = sp.Rational(32*101, 98415)*sp.Rational(64, 125)**2
    return {"fourth_derivative_pole_compact_lower_above_one_over_125": bool(lower4 > sp.Rational(1, 125)),
            "fourth_derivative_pole_upper_below_one_over_25": bool(upper4 < sp.Rational(1, 25)),
            "mixed_fourth_coefficient_between_twice_extreme_diagonals": 2*101 < 220 < 2*120,
            "complex_mass_direction_norm": bool(sp.Rational(4, 9) > sp.Rational(28, 81) > 0),
            "actual_f2_direction_envelope": sp.Rational(592+1440, 2187) == sp.Rational(2032, 2187),
            "actual_f4_direction_envelope": sp.Rational(904+864+1440+432, 2187) == sp.Rational(3640, 2187),
            "integrated_remainder_constant": c["integrated_absolute_remainder_constant"] == sp.Rational(4852, 229635),
            "integrated_remainder_constant_below_one_over_40": bool(c["integrated_absolute_remainder_constant"] < sp.Rational(1, 40)),
            "unit_complex_momentum_ball_log_margin": 1-sp.Rational(1, 4) == sp.Rational(3, 4),
            "unit_ball_remainder_constant_below_one_over_30": sp.Rational(1, 40)/sp.Rational(3, 4) == sp.Rational(1, 30),
            "loop_prefactor_remainder_denominator": 64*30 == 1920,
            "rational_pi_lower_control": 1920*9 == 17280,
            "subtraction_not_unknown_higher_operator_or_full_quantum_matching": True}


def scale_bound(planck_time_product, reference_mass_time_product, momentum_squared_time_squared=1):
    L, R, Q = [exact_parent.exact(value, name) for value, name in
               ((planck_time_product, "M_tau"), (reference_mass_time_product, "m0_tau"),
                (momentum_squared_time_squared, "Hermitian_momentum_squared_tau_squared"))]
    if L.is_positive is not True or R.is_positive is not True or Q.is_nonnegative is not True or Q > R**2:
        raise ValueError("Require positive scales and 0<=Hermitian momentum squared<=m0 squared")
    return {"M_tau": L, "m0_tau": R, "Hermitian_momentum_squared_tau_squared": Q,
            "momentum_ball_ratio": Q/R**2,
            "subtracted_loop_kernel_over_reference_density": Q**3/(17280*L**2*R**2),
            "full_curved_or_quantum_EFT_remainder": False}


@cache
def checks():
    d = logarithmic_form()
    return {**dimensional_logarithms()["checks"],
            **{name: d[name] for name in ("homogeneous_parameter_integrand", "exact_through_four_derivatives_subtraction")}}
