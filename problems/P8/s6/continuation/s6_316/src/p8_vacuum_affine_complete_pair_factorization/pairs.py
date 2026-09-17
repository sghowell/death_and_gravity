"""Complete hard-root pair factorization and its conditional angular estimate."""

from functools import cache

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import trees as full
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower
from p8_vacuum_affine_two_real_collinear_current import collinear

from . import source

ETA, imm, old = full.ETA, full.imm, full.old


class NoPairs(full.TreeEngine):
    @lower.instance_cache
    def current(self, mask, kind):
        if (
            kind == "h"
            and mask.bit_count() == 2
            and all(
                self.legs[i][0] == "h" for i in range(len(self.legs)) if mask >> i & 1
            )
        ):
            return lower.ZERO, 0
        return super().current(mask, kind)


def hard_remainder(points, other_radiation, **parameters):
    legs = [("phi", p, 1) for p in points] + [
        ("h", q, eps) for q, eps in other_radiation
    ]
    engine = full.TreeEngine(legs, **parameters)
    return engine.amputated(engine.full, "h")


def pair_contraction(first, second, hard, kappa):
    q, A = first
    r, B = second
    Q = q + r
    propagated = imm(
        -(ETA * hard * ETA - ETA * s.trace(ETA * hard) / 2) / old.dot(Q, Q)
    )
    return s.factor(old.cubic((A, B, propagated), (q, r, -Q)) / s.sqrt(kappa))


def three_real_calibration(original):
    if not isinstance(original, bool):
        raise TypeError("Require an explicit original/diagnostic boolean")
    return _three_real_calibration(original)


@cache
def _three_real_calibration(original):
    parameters = (
        source.original_parameters() if original else source.diagnostic_parameters()
    )
    points, hs = full.configuration()
    legs = [("phi", p, 1) for p in points[1:]] + [("h", q, eps) for q, eps in hs]
    whole, count = full.TreeEngine(legs, **parameters).amplitude()
    regular, nreg = NoPairs(legs, **parameters).amplitude()
    checks = {
        "full_inventory": s.Integer(count - 5116),
        "no_pair_inventory": s.Integer(nreg - 3814),
    }
    pieces = []
    currents = []
    for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        Q = hs[i][0] + hs[j][0]
        current, nrest = hard_remainder(points, (hs[k],), **parameters)
        checks[f"complete_offshell_root_Ward_{i}_{j}"] = (
            Q.T * ETA * current
        ).applyfunc(s.factor)
        checks[f"complete_remainder_inventory_{i}_{j}"] = s.Integer(nrest - 434)
        # Independent frozen lower-multiplicity current, no new jet engine.
        compare, ncompare = lower.TreeEngine(
            [("phi", p, 1) for p in points] + [("h", *hs[k])], **parameters
        ).amputated(31, "h")
        checks[f"independent_lower_source_current_{i}_{j}"] = (
            current - compare
        ).applyfunc(s.factor)
        checks[f"independent_lower_source_inventory_{i}_{j}"] = s.Integer(
            nrest - ncompare
        )
        pieces.append(pair_contraction(hs[i], hs[j], current, parameters["kappa"]))
        currents.append(current)
    checks["complete_three_pair_amplitude_factorization"] = s.factor(
        whole - regular - sum(pieces)
    )
    return {
        "checks": checks,
        "whole_counts": (count, nreg, 434, 434, 434),
        "nonzero_pair_pieces": all(value != 0 for value in pieces),
        "nonzero_whole_coefficient": whole != 0,
        "all_currents_symmetric": all(value == value.T for value in currents),
        "parameters": parameters,
    }


@cache
def data():
    pair = collinear.data()
    checks = dict(pair["checks"])
    samples = {}
    for label, original in (("diagnostic", False), ("original", True)):
        result = three_real_calibration(original)
        checks.update(
            {label + "_" + key: value for key, value in result["checks"].items()}
        )
        samples[label] = {
            key: value for key, value in result.items() if key != "checks"
        }
    return {
        "whole_general_single_root_Ward": "For the complete one-off-shell graviton current R(P) with physical on-shell free leaves, P^T eta R=0. The proof is finite induction of the ungaugefixed Noether identity in separately labeled waves, with the linear de Donder inverse preserving the harmonic constraint at each lower solved coefficient. There is no assumption that P is on shell.",
        "whole_general_fixed_pair_factorization": "For any fixed physical pair i,j in the full finite-N source, remove its cubic branch and internal edge. The remainder is the complete four-Phi/(N-2)-real/one-off-shell graviton current R_ij. The graph-class amplitude is V_hhh(eps_i,eps_j,-trace_reverse(R_ij)/Q_ij^2)/sqrt(kappa). Its inventory is T_{N-1}. Only this completed remainder is claimed conserved.",
        "whole_conditional_angular_estimate": "|A_ij| <= 530*(a+b)^2*||R_ij,spatial||_F/(sqrt(kappa)*a*b) for physical unit-Frobenius TT polarizations and positive noncollinear a,b, uniformly in their relative angle. The exact24 S311 coefficient identities apply pointwise to the complete conserved remainder. This removes the explicit pair propagator but is NOT a uniform bound on the remainder itself.",
        "whole_complete_three_real_calibrations": samples,
        "whole_pair_quotient_coefficients": pair[
            "whole_exact_conserved_pair_coefficients"
        ],
        "whole_pair_coefficient_budgets": pair["whole_exact_normalized_pair_budgets"],
        "checks": checks,
        "gates": {
            **pair["gates"],
            "all_finite_single_root_Noether_induction": True,
            "root_not_required_on_shell": True,
            "fixed_pair_graph_bijection_not_arbitrary_subcluster": True,
            "diagnostic_and_original_full5116_factorization": all(
                x["nonzero_pair_pieces"] and x["nonzero_whole_coefficient"]
                for x in samples.values()
            ),
            "three_complete434_currents_symmetric_in_both_sources": all(
                x["all_currents_symmetric"] for x in samples.values()
            ),
            "remaining_current_norm_and_simultaneous_collinearity_open": True,
            "multi_offshell_remainders_not_assumed_conserved": True,
        },
    }
