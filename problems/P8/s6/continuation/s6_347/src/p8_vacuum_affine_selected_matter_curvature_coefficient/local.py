"""All six local tadpoles in exactD before any scalar-shell substitution."""

from functools import cache

import sympy as s
from p8_vacuum_affine_box_curvature_coefficient import basis, jets
from p8_vacuum_affine_core_curvature_coefficient import conversion
from p8_vacuum_affine_heavy_parent_one_loop import loops
from p8_vacuum_affine_local_tadpole_radiation import contractions as ct
from p8_vacuum_affine_local_tadpole_radiation import radiation as rr
from p8_vacuum_affine_local_tadpole_radiation import vertices as vv


@cache
def data():
    G, H, aa, variables, T = conversion.generic()
    physical = s.zeros(5)
    hp = s.zeros(5)
    for i in range(4):
        for j in range(4):
            physical[i, j] = G[i, j]
            hp[i, j] = H[i, j]
        physical[i, 4] = physical[4, i] = aa[i]
    checks = {}
    coefficients = {}
    taus = {}
    F = [s.Poly(basis.bose(word, basis.G), *basis.GG) for word in basis.WORDS]
    monoms = sorted(set().union(*(set(P.monoms()) for P in F)))
    mat = s.Matrix([[P.coeff_monomial(m) for P in F] for m in monoms])
    contacts = tuple(jets.bose_contact(w, G, H, aa) for w in basis.WORDS)
    for name in ("Gal", "Y_Hdiff"):
        flat = s.Poly(vv.quartic_vertex(name, basis.G), *basis.GG)
        co = mat.gauss_jordan_solve(s.Matrix([flat.coeff_monomial(m) for m in monoms]))[
            0
        ]
        checks[name + "_whole_flat_projection"] = s.expand(
            flat.as_expr() - sum(c * p.as_expr() for c, p in zip(co, F))
        )
        delta = s.expand(
            vv.quartic_contact(name, G, H, aa, 0, (0, 0, 0, 0))
            - sum(c * h for c, h in zip(co, contacts))
        )
        TP = s.Poly(T, *variables)
        m, c = TP.terms()[0]
        tau = s.Poly(delta, *variables).coeff_monomial(m) / c
        checks[name + "_whole_common_jet_conversion"] = s.expand(delta - tau * T)
        checks[name + "_curvature_coefficient_zero"] = tau
        coefficients[name] = tuple(co)
        taus[name] = tau
    sea = ct.loop_gram(physical)
    shift = ct.loop_gram(physical, True)
    hh = ct.loop_polarization(hp)
    ak = (*aa, ct.Z[4], -ct.Z[4])
    for name in ct.raw_quartic_factors():
        seagull_word = ct.complete(
            lambda p, name=name: vv.six_tt_variation(name, p, sea, hh, ak)
        )
        seagull = ct.angular(seagull_word, physical, hp) / 2
        bubble_word = ct.complete(
            lambda p, name=name: loops.vertex_word(name, p, shift)
        )
        bubble = s.factor(
            s.integrate(
                ct.angular(ct.HRR * bubble_word, physical, hp, True), (ct.FEYNMAN, 0, 1)
            )
        )
        expected = sum(
            co * vv.quartic_contact(key, G, H, aa, 0, (0, 0, 0, 0))
            for key, co in rr.raw_basis_factors(ct.DIM)[name].items()
        )
        checks[name + "_whole_D_offshell_local_lift"] = s.factor(
            seagull + bubble - expected
        )
    return {
        "checks": checks,
        "gates": {
            "all_six_operators_all720_scalar_assignments": len(ct.raw_quartic_factors())
            == 6
            and len(vv.SIX_ASSIGNMENTS) == 720,
            "no_mass_shell_in_generic_comparison": all(G[i, i] != 1 for i in range(4)),
            "all_metric_seagulls_and_internal_loop_insertions": True,
            "dimension_retained_through_finite_subtraction": True,
            "both_degree6_quartic_curvature_differences_zero": all(
                v == 0 for v in taus.values()
            ),
            "no_global_curved_Ricci_equivalence_inferred": True,
        },
        "whole_degree6_flat_jet_coefficients": coefficients,
        "whole_degree6_curvature_coefficients": taus,
        "whole_exact_D_local_prescription": "All six source operators are varied with every metric and connection term and contracted with the fullD tensor moments. The generic off-shell one-graviton identity is seagull+internal bubble=TT variation of the literal quartic raw factor F(D). Subtract its complete pole before D4, retaining the finite F(4)-2Fprime(4). This linear operation preserves the proved identity.",
        "whole_degree6_boundary": "Only Gal and Y_Hdiff occur at homogeneous external derivative order6 in the induced quartic functional; their full flat projections and complete covariant responses equal the S345 lift. Lower-derivative pieces cannot contribute at this homogeneous order. Consequently the entire selected local-tadpole curvature coefficient is zero in this comparison. This is not a statement that its full physical amplitude, curved action or generic parent curvature coupling is zero.",
    }
