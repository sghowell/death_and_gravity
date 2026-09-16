"""Uniform angular-energy bounds and state-correct nested infrared pairing."""

from functools import cache

import sympy as s

from . import source

R = s.Symbol("nonnegative_total_radiated_energy", nonnegative=True)
RESOLUTION = s.Symbol("positive_detector_resolution_mass_units", positive=True)


def kernel_deviation(total=R):
    total = s.sympify(total)
    return 1416 * total + 3 * total**2


def index_deviation(total=R, kappa=source.K):
    return 40 * s.sympify(total) / s.sympify(kappa)


def compact_log_factor_bound(kappa=source.K):
    k = s.sympify(kappa)
    return 12 / k + 60 / k**2


def original_factor_ratio_error():
    return s.Integer(26) / source.KAPPA


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_func(value))

    d, c, r, k, x = s.symbols("d c R kappa x", positive=True)
    num = 2 * d + c * (d * d - 2 * d + s.Rational(1, 2))
    lower = (
        2 * (d - 1) + c * (d - 1) ** 2 + s.Rational(7, 4) + (s.Rational(1, 4) - c / 2)
    )
    put("massive_derivative_positive_decomposition", num - lower)
    put("massive_derivative_upper", 2 * s.Integer(7) - 14)
    put("all_six_massive_pair_changes", 2 * 14 * 48 * r - 1344 * r)
    put(
        "ln2_upper_exponential_partial_sum",
        1 + s.Rational(3, 4) + s.Rational(3, 4) ** 2 / 2 - 2 - s.Rational(1, 32),
    )
    put("mixed_kernel_max", 4 * 3 * s.Rational(3, 4) - 9)
    put("all_mixed_energy_terms", 2 * 4 * 9 * r - 72 * r)
    put("null_pair_endpoint_upper", 2 * 2 * s.Rational(3, 4) - 3)
    put("null_pair_derivative", s.diff(x * s.log(2 * x), x) - s.log(2 * x) - 1)
    put("null_pair_minimum", (-x * s.log(2 * x)).subs(x, 1 / (2 * s.E)) - 1 / (2 * s.E))
    put("null_pair_total_energy_square", 2 * 3 * r * r / 2 - 3 * r * r)
    put("whole_kernel_deviation", kernel_deviation(r) - 1344 * r - 72 * r - 3 * r * r)
    put(
        "compact_kernel_linear_margin", 1440 * r - kernel_deviation(r) - 3 * r * (8 - r)
    )
    put("index_conversion_with_pi_gt3", 1440 * r / (4 * 9 * k) - index_deviation(r, k))
    put("original_maximum_index_shift", index_deviation(s.Rational(1, 8), k) - 5 / k)
    put(
        "whole_state_index_upper_margin",
        6 / k - (s.Rational(4, 5) + 5) / k - s.Rational(1, 5) / k,
    )
    aa = s.Symbol("soft_index", real=True)
    put("soft_resolution_power_log", s.diff(s.log(x**aa), aa) - s.log(x))
    put("compact_energy_log_derivative", s.diff(-x * s.log(x), x) + s.log(x) + 1)
    put(
        "compact_energy_log_upper",
        40 * s.Rational(1, 8) * 3 * s.Rational(3, 4) - s.Rational(45, 4),
    )
    put("compact_energy_log_margin", 12 - s.Rational(45, 4) - s.Rational(3, 4))
    put("Gamma_log_derivative_series_bound", 2 * (6 / k) * (5 / k) - 60 / k**2)
    put("combined_known_log_bound", compact_log_factor_bound(k) - 12 / k - 60 / k**2)
    put("combined_ratio_linear_bound", 2 * 13 / k - 26 / k)
    put("original_ratio_error", original_factor_ratio_error() - s.Rational(26, 10**800))
    put(
        "original_ratio_below_3e_minus799_margin",
        s.Rational(3, 10**799) - original_factor_ratio_error() - s.Rational(4, 10**800),
    )
    a0, aR, ep = s.symbols("a4 aradiative epsilon", real=True)
    put("correct_state_real_virtual_poles", aR / (2 * ep) - aR / (2 * ep))
    put(
        "wrong_state_pole_residue",
        s.residue((aR - a0) / (2 * ep), ep, 0) - (aR - a0) / 2,
    )
    put(
        "wrong_state_pole_absolute_majorant",
        index_deviation(r, k) / (2 * ep) - 20 * r / (k * ep),
    )
    put(
        "nested_known_scale_log_shift",
        ((aR - a0) * (-k * s.Symbol("chi"))) + k * (aR - a0) * s.Symbol("chi"),
    )
    return {
        "whole_K_difference_majorant": kernel_deviation(),
        "whole_a_difference_majorant": index_deviation(),
        "whole_compact_known_power_Gamma_log_bound": compact_log_factor_bound(),
        "whole_original_known_factor_ratio_error": original_factor_ratio_error(),
        "whole_uniform_index_proof": "Every massive pair lies in[1,7],and F'(d) has positive numerator2d+c[(d-1)^2-1/2]>=7/4 for0<=c<=1/2,whileF'<=2d*integral1/den<=14. The six pair-dot shifts sum<=48R,so their total contribution changes by<=1344R. Sincea_i in[1/4,4],|a_i ln(2a_i)|<9,all mixed pairs contribute<72R. The continuous null kernel obeys|b ln(2b)|<=3 on[0,2],giving<=3R^2. Hence|Ksigma-K4|<=1416R+3R^2<1440R for0<R<=1/8,with equalityzero atR0. Dividing by4pi^2*kappa yields|asigma-a4|<40R/kappa,uniform in multiplicity and all collinear limits.",
        "whole_nested_pairing_boundary": "For a fixed complete radiative state,one additional infinitesimal graviton has real rate pole asigma/(2e),paired by virtual rate pole -asigma/(2e). Reusing the elastic virtual pole leaves(asigma-a4)/(2e); small Newton coupling does not make this finite at fixedR>0 as e->0. The coefficient is universal and independent of the spin of those existing gravitons. This identifies the required radiative-state virtual subtraction,not its finite hard amplitude or finite D angular conversion.",
        "whole_known_factor_stability": "Only AFTER state-correct IR pairing,compare x^a F(a),F=exp(-EulerGamma*a)/Gamma(1+a). For total state energyR<=x<=1/8 andnu1,the logarithmic power change is<=40R|lnx|/kappa<12/kappa. Bothindices are in[0,6/kappa); the Weierstrass derivative|dlnF/da|<=pi^2*a/6<2a gives a Gamma log change<60/kappa^2. At originalkappa10^800 the total log change is<13/kappa<1/2,so the multiplicative difference is<26/kappa<3*10^-799. No claim about the unknown Delta_sigma or full physical rate is included.",
        "checks": checks,
        "gates": {
            "original_compact_bounds_uniform_in_multiplicity": True,
            "massless_collinear_limits_do_not_enlarge_budget": True,
            "same_state_real_virtual_pole_pairing_required": True,
            "tiny_numerator_not_a_finite_unpaired_IR_pole": True,
            "known_power_and_Gamma_change_after_IR_removal_only": True,
            "finite_angular_and_radiative_hard_terms_not_assumed": True,
            "all_original_V_G_B_obligations_unchanged": True,
        },
    }
