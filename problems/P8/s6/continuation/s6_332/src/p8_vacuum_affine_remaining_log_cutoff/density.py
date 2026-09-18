"""Common-space signed-density L1 cutoff, not configuration-law TV."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantitative_soft_cutoff.conditioning import parameters

from . import source


def conversion_density_cutoff_upper(index, resolution, infrared):
    a, _, eta = parameters(index, resolution, infrared)
    return 42000 * a * eta * (1 - s.log(eta)) / source.KAPPA


def common_space_cutoff_upper(index, resolution, infrared):
    a, _, eta = parameters(index, resolution, infrared)
    return 63000 * a * eta * (1 - s.log(eta)) ** 2 / source.KAPPA


def original_common_space_cutoff_upper(resolution, infrared):
    _, _, eta = parameters(0, resolution, infrared)
    return 51000 * eta * (1 - s.log(eta)) ** 2 / source.KAPPA**2


@cache
def data():
    m, n, r, eta, t, u = s.symbols("m n r eta t u", real=True)
    z = s.Symbol("z", positive=True)
    f = lambda v: v * (1 - s.log(v))
    checks = {
        "same_event_density_normalization": m * (1 - r) - n - ((m - n) - r * m),
        "tail_entropy_subadditivity": s.expand(
            f(t)
            + f(u)
            - f(t + u)
            - ((t + u) * s.log(t + u) - t * s.log(t) - u * s.log(u))
        ),
        "distinguished_tail_entropy_integral": s.integrate(1 - s.log(z), (z, 0, eta))
        - eta * (2 - s.log(eta)),
        "conversion_density_budget": s.Integer(2 * 11000 + 20000 - 42000),
        "whole_density_budget": s.Integer(42000 + 14600 + 6300 - 62900),
        "original_parameter_budget": s.Rational(63000 * 4, 5) - 50400,
        "exact_zero_index_density": common_space_cutoff_upper(
            0, s.Rational(1, 8), s.Rational(1, 64)
        ),
    }
    gates = {
        "common_underlying_probability_space_is_specified": True,
        "both_normalization_error_terms_retained": True,
        "positive_endpoint_equalities_null_for_nonzero_index": True,
        "zero_index_mark_vanishes_exactly": True,
        "strict_general_roundup": bool(s.Integer(62900) < 63000),
        "strict_original_roundup": bool(s.Integer(50400) < 51000),
        "TV_of_signed_density_not_of_configuration_pushforward_laws": True,
        "quantitative_limiting_insertion_not_joint_regulator_rate": True,
        "no_full_state_evaluation_outside_its_physical_cut": True,
    }
    return {
        "checks": {name: s.simplify(value) for name, value in checks.items()},
        "gates": gates,
        "whole_common_space_comparison": "On the full fixed-Born Poisson space let A={R<=x}, A_eta={R_eta<=x}, p=P(A), p_eta=P(A_eta). Compare f=1_A*Z(sigma;x-R)/p and f_eta=1_A_eta*Z(sigma_eta;x-R_eta)/p_eta, with Z=deltaDelta+deltaa*ln(remaining energy). Conditional probabilities are never replaced by 1.",
        "whole_density_split": "For any mark, E|f_eta-f| <=E[|M_eta-M||A]+r*E[|M_eta||A]+E[|M_eta|*1_(A_eta minus A)|A_eta], r=1-p/p_eta<=a*eta/x. The finite-conversion contribution is <=42000*a*eta*L_eta/kappa.",
        "whole_quantitative_statement": "For0<=a<=1 and0<eta<=x<=1/8, the common-space signed-density L1 distance is <=63000*a*eta*L_eta^2/kappa, hence <51000*eta*L_eta^2/kappa^2 at original a<4/(5kappa). It implies conditioned-mean convergence. Total variation means the total variation norm, without a factor1/2.",
        "whole_negative_boundary": "For a>0 the retained emission configuration has finite count and the full configuration infinite count almost surely, even after the positive-probability energy cut. Those configuration laws are singular. This is not a claim about every angular-energy pushforward, which can coalesce atoms.",
    }
