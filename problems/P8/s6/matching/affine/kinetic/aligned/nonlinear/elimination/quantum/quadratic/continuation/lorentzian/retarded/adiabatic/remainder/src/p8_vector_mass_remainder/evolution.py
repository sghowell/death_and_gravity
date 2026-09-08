"""Exact mixing-variation algebra and continuous momentum integral bounds."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact
from p8_vector_mass_adiabatic import bounds as local_bounds
from p8_vector_regularity import preparation
from p8_vector_state import comparison, wkb

from . import envelopes, tail

AMAX = comparison.AMAX


@cache
def constants():
    by_sector = {}
    for sector in ("T", "L"):
        data = envelopes.reference(sector)
        C = data["base_residual_over_inverse_frequency_eighth_upper"]
        dC = data["delta_residual_over_inverse_frequency_eighth_upper"]
        B = data["delta_W_over_frequency_upper"]
        m = wkb.MASS_TIME_MIN
        E = B+data["delta_c_upper"]/m
        K = AMAX*((12+4*E)/m+20*AMAX*B)/2
        by_sector[sector] = {
            "delta_mixing_ninth_inverse_frequency_coefficient": sp.ceiling(16*(dC+2*C*B)),
            "delta_mixing_eighth_inverse_frequency_coefficient": sp.ceiling(16*C*B*AMAX),
            "delta_reference_bilinear_over_nu_squared_upper": sp.ceiling(K),
            "varied_reference_tail_over_inverse_frequency_fifth_upper":
                tail.reference_tail(sector)["varied_reference_tail_over_inverse_frequency_fifth_upper"]}
    common = {key: max(row[key] for row in by_sector.values()) for key in by_sector["T"]}
    return {"by_sector": by_sector, "common": common,
            "frozen_initial_and_evolved_mixing": preparation.constants()["common"]}


def bound(planck_time_product, mass_time_product):
    L, m = exact(planck_time_product, "M_tau"), exact(mass_time_product, "m0_tau")
    if not bool(L > 0) or not bool(m >= 1000):
        raise ValueError("Require M*tau>0 and m0*tau>=1000")
    data, old = constants()["common"], constants()["frozen_initial_and_evolved_mixing"]
    beta6 = old["initial_mixing_envelopes"][6]+old["evolution_mixing_envelope"]/m**4
    I4, I7, I8 = AMAX**3/(24*m), AMAX**3/(135*m**4), AMAX**3/(192*m**5)
    D9 = data["delta_mixing_ninth_inverse_frequency_coefficient"]
    D8 = data["delta_mixing_eighth_inverse_frequency_coefficient"]
    K = data["delta_reference_bilinear_over_nu_squared_upper"]
    T = data["varied_reference_tail_over_inverse_frequency_fifth_upper"]
    parts = {"exact_minus_reference_variation": 3*(30*AMAX*(D9*I8+D8*I7)+6*beta6*K*I4)/L**2,
             "varied_reference_minus_adiabatic": 3*T/(54*m**2*L**2)}
    local = local_bounds.operator_bound(L, m)["C4_to_C0_upper_bound"]
    return {"components": parts, "nonlocal_subtracted_C10_to_C0_upper_bound": sum(parts.values()),
            "finite_local_component_upper_bound": local,
            "complete_mass_response_C10_to_C0_upper_bound": sum(parts.values())+local,
            "source_norm": "Maximum of time derivatives zero through ten on [-1/2,1/2], with source zero on an initial neighborhood",
            "scope": "Homogeneous prepared mass-source response; not a coupled feedback contraction or a spatial response bound"}


@cache
def algebra_checks():
    rho, W, theta = sp.symbols("residual frequency phase", real=True)
    K = sp.I*rho/(2*W)*sp.Matrix([[-1, -sp.exp(2*sp.I*theta)], [sp.exp(-2*sp.I*theta), 1]])
    eta = sp.diag(1, -1)
    out = {"Bogoliubov_Wronskian_preserving_transport": (sp.conjugate(K).T*eta+eta*K).applyfunc(sp.simplify)}
    aa, bb, f, p = sp.symbols("mixing_A mixing_B reference_f reference_p")
    A, B, omega2 = sp.symbols("readout_A readout_B omega_squared", real=True)
    v, pp = aa*f+bb*sp.conjugate(f), aa*p+bb*sp.conjugate(p)
    exact_Q = (A*pp*sp.conjugate(pp)+B*omega2*v*sp.conjugate(v))/2
    Q = (A*p*sp.conjugate(p)+B*omega2*f*sp.conjugate(f))/2
    Z = (A*p**2+B*omega2*f**2)/2
    target = (aa*sp.conjugate(aa)+bb*sp.conjugate(bb))*Q+aa*sp.conjugate(bb)*Z+sp.conjugate(aa)*bb*sp.conjugate(Z)
    out["exact_mode_readout_mixing_decomposition"] = sp.expand(exact_Q-target)
    m = sp.Symbol("positive_mass", positive=True)
    for power, expected in ((4, 1/(8*sp.pi*m)), (7, 1/(15*sp.pi**2*m**4)), (8, 1/(64*sp.pi*m**5))):
        actual = m**(3-power)*sp.gamma(sp.Rational(power-3, 2))/(8*sp.pi**sp.Rational(3, 2)*sp.gamma(sp.Rational(power, 2)))
        out["radial_mixing_integral_"+str(power)] = sp.simplify(actual-expected)
    actual = m**-2*sp.gamma(1)/(8*sp.pi**sp.Rational(3, 2)*sp.gamma(sp.Rational(5, 2)))
    out["radial_reference_tail_integral_five"] = sp.simplify(actual-1/(6*sp.pi**2*m**2))
    phase, parameter = sp.symbols("phase variation_parameter", real=True)
    squeezed = 2*sp.Rational(5, 4)*sp.Rational(3, 4)*sp.cos(2*(phase+parameter))
    out["nonzero_initial_mixing_phase_control"] = sp.diff(squeezed, parameter).subs({phase: sp.pi/4, parameter: 0})+sp.Rational(15, 4)
    return out
