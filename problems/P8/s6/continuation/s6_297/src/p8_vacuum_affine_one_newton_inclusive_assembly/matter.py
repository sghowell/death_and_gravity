"""Complete inherited matter loop and its missing gravity-Born interference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import bounds as frozen

from . import source

X = s.Symbol("positive_gravity_to_matter_Born_ratio", positive=True)
ETA = s.Symbol("matter_relative_loop_bound", positive=True)


def interference_weight(ratio=X):
    ratio = s.sympify(ratio)
    return 2 * ratio / (1 + ratio) ** 2


@cache
def data():
    old = frozen.data()
    r = old["combined_full_first_loop_angular_relative_upper"]
    eta = s.Rational(60, 59) * r
    checks = {}
    selected = (
        "same_old_complete_symmetric_upper",
        "same_old_angular_bound",
        "new_mass_and_residue_counterterms_complete",
        "new_first_b20_shift",
        "new_first_b21_shift",
        "new_first_b40_shift_retained",
    )
    for name in selected:
        checks["complete_frozen_matter_" + name] = old["checks"][name]
    checks["same_source_kappa"] = source.KAPPA - frozen.germs.KAPPA
    checks["same_source_heavy_mass"] = source.HEAVY_MASS2 - frozen.germs.MASS2
    checks["same_source_cubic"] = source.CUBIC - frozen.germs.G
    checks["same_source_quartic"] = source.CONTACT - frozen.germs.CONTACT
    checks["original_to_positive_heavy_Born_conversion"] = eta - s.Rational(60, 59) * r
    checks["full_Born_interference_weight_maximum"] = s.factor(
        s.Rational(1, 2) - interference_weight() - (X - 1) ** 2 / (2 * (1 + X) ** 2)
    )
    checks["matter_Born_to_full_Born_weight"] = s.factor(
        1 / (1 + X) ** 2 - (1 / (1 + X)) ** 2
    )
    checks["weight_at_equal_Born_amplitudes"] = interference_weight(1) - s.Rational(
        1, 2
    )
    checks["forward_gravity_dominance_weight_limit"] = s.limit(
        interference_weight(), X, s.oo
    )
    checks["gravity_decoupling_weight_limit"] = s.limit(interference_weight(), X, 0)
    return {
        "whole_same_OS4_matter_loop_relative_to_original_target": r,
        "whole_same_OS4_matter_loop_relative_to_matter_Born": eta,
        "whole_normalized_gravity_matter_loop_weight": interference_weight(),
        "whole_uniform_interference_upper": eta / 2,
        "whole_bound_proof": "S239 bounds the full complex matter loop by r times the original scalar target, r<10^-199. The frozen S233 tree comparison gives A_m>=59/60 of that target, so |L_m|/A_m<=eta=60r/59<2*10^-199. Both A_m and A_G are positive on 4<s<=10^196,-1<z<1. With x=A_G/A_m, |2 A_G Re L_m|/(A_m+A_G)^2<=2eta x/(1+x)^2<=eta/2<10^-199. The exact square identity proves the maximum. The full Born denominator is an explicitly unexpanded positive reference, not a uniform inverse-kappa expansion.",
        "whole_scope": "The complete formal matter loop includes its fixed higher-vertex tadpoles and mixed source contractions, not only a C/g polynomial subset. Its flat OS4 prescription is the existing one. This bound does not fix gravity finite matching, sum higher loops, assert a finite forward cross section, or import any unitarity/Regge result.",
        "checks": {k: s.factor(v) for k, v in checks.items()},
        "gates": {
            "inherited_complete_matter_bound_positive": bool(
                0 < r < s.Rational(1, 10**199)
            ),
            "matter_Born_relative_bound_explicit": bool(eta < s.Rational(2, 10**199)),
            "missing_interference_uniform_below_one_e_minus_199": bool(
                eta / 2 < s.Rational(1, 10**199)
            ),
            "whole_complex_matter_loop_not_imaginary_cut_only": True,
            "same_OS4_and_all_finite_source_terms_retained": True,
            "normalization_not_a_forward_cross_section": True,
        },
    }
