"""Same-event singular logarithm controlled by distinguished emissions."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantitative_soft_cutoff.conditioning import parameters

from . import source


def tail_remaining_log_upper(index, resolution, infrared):
    a, x, eta = parameters(index, resolution, infrared)
    return a * eta * (2 - s.log(x))


def logratio_mark_upper(index, resolution, infrared):
    a, _, eta = parameters(index, resolution, infrared)
    return 11200 * a * eta * (1 - s.log(eta)) / source.KAPPA


def same_event_log_upper(index, resolution, infrared):
    a, _, eta = parameters(index, resolution, infrared)
    return 14600 * a * eta * (1 - s.log(eta)) / source.KAPPA


@cache
def data():
    a, x, h, w, v, y, eta, L = s.symbols("a x h w v y eta L", positive=True)
    J = h + (1 - h) * s.log(1 - h)
    checks = {
        "endpoint_log_antiderivative": s.diff(J, h) + s.log(1 - h),
        "endpoint_log_zero": J.subs(h, 0),
        "endpoint_log_full": s.limit(J, h, 1, dir="-") - 1,
        "log_subadditivity_positive_cross_term": (1 + w / y) * (1 + v / y)
        - 1
        - (w + v) / y
        - w * v / y**2,
        "harmonic_recurrence_for_full_tail": s.expand_func(
            s.polygamma(0, a + 2) - s.polygamma(0, a + 1) - 1 / (a + 1)
        ),
        "distinguished_emission_factor": (
            s.Symbol("n", integer=True, positive=True)
            / s.factorial(s.Symbol("n", integer=True, positive=True))
            - 1 / s.factorial(s.Symbol("n", integer=True, positive=True) - 1)
        ),
        "small_tail_integrated_log": s.integrate(1 + s.log(x / w), (w, 0, eta))
        - eta * (2 + s.log(x / eta)),
        "small_tail_budget": s.Integer(1400 * 4 + 1400 * 2 * 2 - 11200),
        "large_tail_budget": s.Rational(1400 * 4 * 3, 2) - 8400,
        "same_event_total_budget": s.Integer(2 * 1700 + 11200 - 14600),
        "same_event_rounding_slack": (
            14600 * (1 + L) - (1700 * (L + 2) + 11200 * (1 + L))
        )
        - 1700 * L,
    }
    gates = {
        "Campbell_applied_before_conditioning": True,
        "retains_exact_CDF_ratio_one_minus_w_over_x_power_a": True,
        "zero_tail_endpoint_defined_by_limit": True,
        "harmonic_psi_a_plus_one_bounded_by_one_for_a_le_one": True,
        "small_eta_and_large_eta_regimes_both_covered": True,
        "infinite_log_sum_bounded_by_total_tail_over_positive_remaining_energy": True,
        "remaining_energy_positive_almost_surely_on_cut_for_a_positive": True,
        "no_conditioned_Poisson_independence_assumption": True,
        "large_tail_constant_dominated": bool(s.Integer(8400) < 11200),
    }
    return {
        "checks": {name: s.simplify(value) for name, value in checks.items()},
        "gates": gates,
        "whole_exact_tail_log_identity": "E[T_eta*abs(ln(x-R))|R<=x]=a*integral_0^eta(1-w/x)^a[-ln(x-w)+psi(a+1)+EulerGamma]dw <=a*eta*(ln(1/x)+2). At a=0 the tail and mark vanish.",
        "whole_same_event_logratio": "For eta<=x/4, one distinguished emission and the split remaining energy q<=x/2 or q>=x/2 give E[ln(1+T_eta/(x-R))|cut]<=a*eta/x*[4+2a*(2+ln(x/eta))]. The complementary eta>x/4 regime uses the exact R-weighted remaining-log moment.",
        "whole_same_event_log_bound": "The unchanged current moduli 1700*T/kappa and 1400*R_eta/kappa give a total same-event log error <=14600*a*eta*(1+ln(1/eta))/kappa.",
    }
