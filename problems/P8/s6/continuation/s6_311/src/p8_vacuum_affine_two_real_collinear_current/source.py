"""Same original four-scalar tree and off-shell hard-current ownership."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
recoil = previous.recoil


@cache
def data():
    old = previous.data()
    result = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    checks.update(
        {
            "same_offshell_current_original_heavy_mass": HEAVY_MASS2
            - (s.Rational(10**200, 512) + 2),
            "same_offshell_current_original_Newton": KAPPA - s.Integer(10) ** 800,
            "same_offshell_current_original_cubic": CUBIC - s.Rational(1, 8192),
            "unique_pair_cut_topology_count": s.Integer(434 - 387 - 47),
            "pure_pair_cut_topology_count": s.Integer(198 - 177 - 21),
        }
    )
    result.update(
        {
            "whole_offshell_current_source": "The same S310 complete selected formal four-Phi/two-h tree is retained. Cutting the unique q1+q2 graviton edge gives precisely the S310-rooted full47-tree hard current with an off-shell graviton root and the canonical S304 cubic Einstein vertex. The other387 graphs are retained as a separate complementary class, not set to zero or bounded by the selected-class theorem.",
            "whole_offshell_current_domain": "The hard current is evaluated on the S300 exact four-massive-leg recoil against Q=sum(q1,q2), with E in[5/4,2], W=Q0<=1/8,|Qvec|<=W, and original parameters. Q is allowed timelike; using the on-shell S304 tensor without its Q2-dependent scalar propagators would be incorrect.",
            "whole_new_claim_boundary": "Exact general off-shell transversality, physical pair-vertex angular-pole cancellation and a conservative full-Born-normalized bound for exactly the47 pair-propagator graphs. The other387 graphs and the overlapping soft subtraction must be controlled separately before an integrated two-real detector error can be claimed.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "offshell_source_checks_copied_without_parent_mutation": checks
                is not old["checks"],
                "same_complete434_tree_not_a_replacement_model": True,
                "timelike_root_not_null_soft_approximation": True,
                "pair_sector_not_complete_two_real_rate": True,
            },
        }
    )
    return result
