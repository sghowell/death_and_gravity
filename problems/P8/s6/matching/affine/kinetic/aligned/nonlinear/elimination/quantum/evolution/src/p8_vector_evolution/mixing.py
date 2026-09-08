"""Phase-separated mixing and differentiated physical quadratic forms."""
from functools import cache

import sympy as sp
from p8_vector_state import comparison, wkb

from . import majorants

C = majorants.RESIDUAL_CONSTANT
C1 = majorants.RESIDUAL_DERIVATIVE_CONSTANT
MIXING_CONSTANT = 4*C1+138*C


@cache
def identities():
    r, theta, eta = sp.symbols("residual_rate theta diagonal_phase", real=True)
    original = sp.Matrix([[-sp.I*r, -sp.I*r*sp.exp(2*sp.I*theta)],
                          [sp.I*r*sp.exp(-2*sp.I*theta), sp.I*r]])
    T = sp.diag(sp.exp(-sp.I*eta), sp.exp(sp.I*eta))
    actual = T.inv()*original*T-sp.diag(-sp.I*r, sp.I*r)
    expected = sp.Matrix([[0, -sp.I*r*sp.exp(2*sp.I*(theta+eta))],
                          [sp.I*r*sp.exp(-2*sp.I*(theta+eta)), 0]])
    rd, aa, ad, psi, psid = sp.symbols("r_dot a a_dot phase_rate phase_rate_dot")
    phase = sp.Symbol("total_phase")
    g, gd = r*aa/(2*psi), (rd*aa+r*ad)/(2*psi)-r*aa*psid/(2*psi**2)
    exponential = sp.exp(-2*sp.I*phase)
    ibp = -gd*exponential+2*sp.I*psi*g*exponential+gd*exponential-sp.I*r*aa*exponential
    return {"diagonal_phase_removed_exactly": (actual-expected).applyfunc(sp.simplify),
            "oscillatory_integration_by_parts_endpoint_and_bulk_identity": sp.simplify(ibp)}


@cache
def physical_derivative_identities():
    vr, vi, pr, pi = sp.symbols("v_real v_imag p_real p_imag", real=True)
    A, B, Ad, Bd, omega, lam, H, d, rho = sp.symbols("A B A_dot B_dot omega lambda H d residual", real=True)
    v2, p2, cross = vr**2+vi**2, pr**2+pi**2, vr*pr+vi*pi
    numerator = A*p2+B*omega**2*v2
    exact = (Ad*p2+Bd*omega**2*v2+2*lam*B*omega**2*v2-3*H*numerator
             +sp.diff(numerator, vr)*(pr+d*vr)+sp.diff(numerator, vi)*(pi+d*vi)
             +sp.diff(numerator, pr)*(-d*pr-omega**2*vr)
             +sp.diff(numerator, pi)*(-d*pi-omega**2*vi))
    target = ((Ad-(2*d+3*H)*A)*p2+omega**2*(Bd+(2*lam+2*d-3*H)*B)*v2
              +2*omega**2*(B-A)*cross)
    reference = exact+sp.diff(numerator, pr)*rho*vr+sp.diff(numerator, pi)*rho*vi
    return {"physical_readout_time_derivative_before_momentum_integration": sp.expand(exact-target),
            "reference_residual_extra_time_derivative": sp.expand(reference-target-2*A*rho*cross)}


@cache
def product_identities():
    ar, ai, br, bi, fr, fi, pr, pi = sp.symbols("ar ai br bi fr fi pr pi", real=True)
    A, B, f, p = ar+sp.I*ai, br+sp.I*bi, fr+sp.I*fi, pr+sp.I*pi
    v, dv = A*f+B*sp.conjugate(f), A*p+B*sp.conjugate(p)
    na, nb = ar**2+ai**2, br**2+bi**2
    norm = lambda value: sp.expand_complex(value*sp.conjugate(value))
    return {"phase_independent_squared_mode_difference": sp.expand(norm(v)-norm(f)
                -2*nb*norm(f)-2*sp.re(A*sp.conjugate(B)*f**2)-(na-nb-1)*norm(f)),
            "phase_independent_squared_derivative_difference": sp.expand(norm(dv)-norm(p)
                -2*nb*norm(p)-2*sp.re(A*sp.conjugate(B)*p**2)-(na-nb-1)*norm(p)),
            "phase_independent_mixed_readout_difference": sp.expand(sp.re(dv*sp.conjugate(v))-sp.re(p*sp.conjugate(f))
                -2*nb*sp.re(p*sp.conjugate(f))-2*sp.re(A*sp.conjugate(B)*p*f)
                -(na-nb-1)*sp.re(p*sp.conjugate(f)))}


@cache
def proof_checks():
    m, amax = wkb.MASS_TIME_MIN, comparison.AMAX
    out = {}
    for kind in ("transverse", "longitudinal"):
        d = majorants.bounds(kind)
        out[kind+"_residual_value_bound"] = bool(d["residual_bracket_upper"] <= C)
        out[kind+"_residual_derivative_bound"] = bool(d["residual_time_derivative_integer_upper"] <= C1)
        out[kind+"_reference_tail_derivative_bound"] = bool(d["reference_tail_density_derivative_integer_upper"] <= majorants.TAIL_DERIVATIVE_CONSTANT)
    out.update({
        "phase_corrected_frequency_lower_quarter": bool(C/m**6 < sp.Rational(1, 4)),
        "phase_corrected_frequency_derivative_upper_seven": bool((C1+4*C)/m**6 < 1),
        "coefficient_norm_exponential_below_two": bool(2*C/m**5 < sp.Rational(1, 2)),
        "integration_by_parts_small_product_absorbed": bool(2*C**2/m**5 < C),
        "mixing_constant_from_endpoint_bulk_and_feedback": MIXING_CONSTANT == 2*(2*C1+69*C),
        "mixing_below_one_for_every_momentum": bool(MIXING_CONSTANT/m**6 < 1),
        "kinetic_weight_derivative_below_two": bool(sp.Rational(4, 3)+sp.Rational(4, 9) < 2),
        "potential_weight_derivative_below_two": bool(sp.Rational(28, 27)+sp.Rational(28, 81) < 2),
        "pressure_weight_derivatives_below_two": bool(sp.Rational(2, 3) < 2),
        "momentum_fraction_derivative_envelope": 2*sp.Integer(2)/4 == 1,
        "c1_derivative_envelope": sp.Integer(4)+1 == 5,
        "differentiated_kinetic_weight_envelope": sp.Integer(2)+2*3+3*2 == 14,
        "differentiated_potential_weight_envelope": sp.Integer(2)+(2*2+2*3+3*2)*sp.Rational(3, 2) == 26,
        "reference_mixed_product_below_two": bool(sp.Integer(3) < 4),
        "differentiated_radial_tail_integrable": bool(sp.Integer(2)-4 < -1),
        "reference_squared_physical_derivative_below_three_omega": bool(sp.Rational(9, 4)+25/m**2 < 3),
        "differentiated_quadratic_form_constant": bool(204/m+24 < 25),
        "derivative_integral_integer_constant": bool((75*amax**2+6)*amax**3/24 < 40),
        "mixing_constant_dominates_residual": bool(MIXING_CONSTANT >= C)})
    return out


def checks():
    return {**identities(), **physical_derivative_identities(), **product_identities()}
