"""Limiting causal source kernels and exact finite-error inequalities."""

from fractions import Fraction as Q
from functools import cache
from math import factorial

import sympy as sp
from p8_composite_asymptotics import limiting


@cache
def checks():
    old = limiting.general_jets()
    v, a = old["v"], old["a"]
    z = sp.Symbol("z", positive=True)
    q = sp.Symbol("q", positive=True)

    def derivative(value):
        return (sp.diff(value, z)*v*z+sp.diff(value, v)*old["vprime"]
                +sp.diff(value, a)*old["aprime"])

    factor = z/sp.sqrt(v)
    normalization = sp.cancel(derivative(derivative(factor))/factor)
    u = sp.Symbol("u", real=True)
    pulse = u**2*(1-u)**2
    return {
        "physical_relative_normalization": sp.simplify(normalization-old["N_relative"]),
        "physical_relative_full_potential": sp.simplify(
            (z*z*v*q+z*z*(v-1))/factor**2-normalization
            -(v*v*q-old["heavy_growth_coefficient"])),
        "physical_relative_source_weight": sp.simplify(z/factor-sp.sqrt(v)),
        "physical_relative_output_weight": sp.simplify(z/factor-sp.sqrt(v)),
        "source_pulse_integral": sp.integrate(pulse, (u, 0, 1))-sp.Rational(1, 30),
        "source_first_moment": sp.integrate((1-u)*pulse, (u, 0, 1))-sp.Rational(1, 60),
        "pulse_zero_endpoints": pulse.subs(u, 0)+pulse.subs(u, 1),
        "pulse_zero_endpoint_derivatives": sp.diff(pulse, u).subs(u, 0)+sp.diff(pulse, u).subs(u, 1),
        "positive_pulse_maximum": sp.Rational(1, 16)-pulse-(u-sp.Rational(1, 2))**2*(sp.Rational(1, 2)-(u-sp.Rational(1, 2))**2),
    }


def rational_checks():
    # exp(1)>sum_0^3 1/n!=8/3. The tail after n=5 is bounded by
    # (1/720)/(1-1/7)=7/4320, proving exp(1)<11/4.
    exp_lower = sum((Q(1, factorial(n)) for n in range(4)), Q(0))
    exp_upper = sum((Q(1, factorial(n)) for n in range(6)), Q(0))+Q(7, 4320)
    first_moment = Q(1, 60)
    error = Q(1, 600)
    kinetic_margin = Q(1, 4)-Q(140, 19)*Q(1, 64)
    light_factor = 1-Q(1, 64)/6
    values = {
        "exp_lower_partial_sum": exp_lower,
        "exp_upper_geometric_tail": exp_upper,
        "exp_upper_below_11_over_4": Q(11, 4)-exp_upper,
        "z1_upper": 7*Q(11, 4)+5/Q(8, 3)-11,
        "z_upper_below_11": 11-(7*Q(11, 4)+5/Q(8, 3)-11),
        "relative_kernel_coefficient_lower": kinetic_margin,
        "relative_kernel_above_one_eighth": kinetic_margin-Q(1, 8),
        "light_kernel_linear_lower_factor": light_factor,
        "pulse_integral": Q(1, 30),
        "pulse_linear_moment": first_moment,
        "finite_full_minus_locked_lower": first_moment-2*error,
        "finite_locked_positive_lower": light_factor*first_moment-error,
        "finite_full_positive_lower": (1+light_factor)*first_moment-error,
        # (2/3)R_eps-L_eps >= I/3-(2/3+1)delta = I/6.
        "one_third_relative_mismatch_margin": first_moment/3-Q(5, 3)*error,
    }
    if not (exp_lower == Q(8, 3) and kinetic_margin == Q(41, 304)
            and values["finite_full_minus_locked_lower"] == Q(1, 75)
            and values["one_third_relative_mismatch_margin"] == Q(1, 360)):
        raise ValueError("A pulse, band or finite-response margin failed")
    if not all(value > 0 for value in values.values()):
        raise ValueError("A positive comparison margin failed")
    return {key: str(value) for key, value in values.items()}


def source_checks():
    """Literal action/clock/source units, one real TT polarization eij*eij=1.

Define the probe convention by delta S_probe=(1/2) integral sqrt|g|
T^{ij} delta g_ij. With T^{ij}=Pi eij/Ae² and delta g_ij=-Ae² Gamma eij,
the source is -Ae³ Pi Gamma/2 in T time. J_u is its literal u-action
coefficient; Jbar=J_u/m belongs to the physical action divided by m.
"""
    m, planck, ae, eps, sigma = sp.symbols("m M A_e epsilon sigma", positive=True)
    pulse, gamma, g, f, z = sp.symbols("p Gamma g f z", real=True)
    jbar = planck**2*eps*sigma*pulse/4
    ju = m*jbar
    pi = -2*m*ju/ae**3
    normalization = planck**2*eps**2*sigma**2/8
    return {
        "physical_to_u_source": sp.cancel(-ae**3*pi*gamma/(2*m)-ju*gamma),
        "normalized_u_source": sp.cancel(ju/m-jbar),
        "physical_probe_stress_units": sp.cancel(pi+planck**2*m**2*eps*sigma*pulse/(2*ae**3)),
        "rescaled_source_action": sp.cancel(jbar*(eps*sigma*(z*g+f)/(1+eps*z))/normalization
                                             -2*pulse*(z*g+f)/(1+eps*z)),
        "absolute_metric_response_scale": sp.cancel(eps*sigma*(z*g+f)/(1+eps*z)
                                                      -(sigma*g+eps*sigma*f/(eps*z))/(1+1/(eps*z))),
    }
