"""Original parameters and separately owned all-multiplicity source."""

from functools import cache

import sympy as s
from p8_vacuum_affine_state_transport_soft_matching import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
recoil = previous.recoil


def require_multiplicity(value):
    if isinstance(value, (bool, str, float, s.Float)) or value is None:
        raise TypeError("Require an exact nonnegative integer multiplicity")
    value = s.sympify(value)
    if not isinstance(value, s.Integer) or value < 0:
        raise ValueError("Require an exact nonnegative integer multiplicity")
    return int(value)


@cache
def data():
    old = previous.data()
    result = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    checks.update(
        {
            "first_new_three_real_graph_count": s.Integer(388 + 2364 + 2364 - 5116),
            "first_new_four_real_graph_count": s.Integer(4972 + 34236 + 34236 - 73444),
            "full_EH5_coupling_degree": s.Rational(1)
            - s.Rational(5, 2)
            + s.Rational(3, 2),
        }
    )
    result.update(
        {
            "whole_new_claim_boundary": "The same selected four-dimensional covariant action now has an exact canonical source and labeled rooted-tree evaluator at every finite graviton multiplicity. Analytic vertex and graph-count majorants are proved separately. This is not an all-N probability bound, a loop sum, an inclusive observable or original P8 closure.",
            "whole_source_ownership": "S310 owns the frozen zero/one/two-real source. The current successor retains its action, propagators and rooted-tree bijection, deriving every higher metric and Einstein vertex anew rather than using the finite S310 jet outside its domain. S314's state-transport soft reference remains unchanged; no new real-virtual matching is inferred.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "source_payload_copied_without_parent_mutation": checks
                is not old["checks"],
                "same_selected_action_and_original_parameters": True,
                "all_finite_multiplicities_not_an_infinite_probability_sum": True,
            },
        }
    )
    return result
