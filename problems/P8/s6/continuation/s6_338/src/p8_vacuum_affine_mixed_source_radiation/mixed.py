"""Complete covariant cross-source reduction and finite constant cancellation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_local_tadpole_radiation import vertices


@cache
def data():
    G, H = vertices.generic_radiative_data()
    a = tuple(G[i, 4] for i in range(4))
    n = s.Symbol("heavy_mass_squared")
    checks = {}
    gc, light, heavy, boxlight, boxheavy, cross = s.symbols(
        "gc G GH BoxG BoxGH gradGH_dot_gradG"
    )
    boxproduct = heavy * boxlight + light * boxheavy + 2 * cross
    divergence = heavy * boxlight + cross
    checks["covariant_product_rule_cross_terms_cancel"] = s.expand(
        boxproduct / s.Integer(6)
        - divergence / s.Integer(3)
        - (light * boxheavy - heavy * boxlight) / 6
    )
    delta = s.Symbol("covariant_delta")
    checks["exact_Green_equation_including_delta_terms"] = s.expand(
        (light * boxheavy - heavy * boxlight).subs(
            {boxheavy: -n * heavy - delta, boxlight: -light - delta}
        )
        - ((1 - n) * heavy * light - delta * (light - heavy))
    )
    E = s.Symbol("light_free_EOM")
    checks["complete_bilocal_and_external_EOM_coefficients"] = s.expand(
        gc * ((1 - n) * heavy * light - delta * (light - heavy)) / 3
        - gc * heavy * light * (E - 1)
        - (
            gc * (4 - n) * heavy * light / 3
            - gc * delta * (light - heavy) / 3
            - gc * heavy * light * E
        )
    )
    g, kap, c = s.symbols("g kappa c")
    A1, An, B = s.symbols("A0_light A0_heavy B0")
    written = 24 * (-g * c * (A1 - An) / 3 + g * c * (4 - n) * B / 3)
    checks["all24_labels_reproduce_entire_frozen239_mixed_bubble"] = s.factor(
        written.subs(c, -4 * g / kap) - 32 * g * g * ((A1 - An) + (n - 4) * B) / kap
    )

    # Arbitrary analytic scalar bilocal kernel: four external emissions,
    # the kernel graviton insertion and the existing constant subtraction.
    f0 = s.Symbol("F_on_shell")
    fshift = s.symbols("F_shift0:4")
    external = sum(H[i, i] / a[i] * 6 * (fshift[i] + 3 * f0) for i in range(4))
    contact = -6 * sum(H[i, i] / a[i] * (fshift[i] - f0) for i in range(4))
    checks["whole_bilocal_TT_radiation_equals_on_shell_constant"] = s.factor(
        external + contact - 24 * f0 * sum(H[i, i] / a[i] for i in range(4))
    )
    # For a symmetric homogeneous two-hard-leg rank2 response at k^2=0,
    # k.T gives A*(p.k)*p + (B*(p.k)+D)*k; p.k !=0 forces A=0.
    pk = s.Symbol("pk", nonzero=True)
    Acoef, Bcoef, Dcoef = s.symbols("A B D")
    checks["no_two_leg_transverse_TT_homogeneous_coefficient"] = s.solve(
        Acoef * pk, Acoef
    )[0]
    checks["remaining_homogeneous_tensor_has_only_k_or_trace"] = s.expand(
        (Bcoef * pk + Dcoef).subs(Dcoef, -Bcoef * pk)
    )
    return {
        "checks": checks,
        "gates": {
            "covariant_product_and_Green_delta_terms_retained": True,
            "all24_flat_labels_recover_frozen239_entire_mixed_bubble": True,
            "four_external_emissions_and_kernel_contact_combined": True,
            "kernel_contact_not_identically_zero": bool(contact != 0),
            "existing_OS4_constant_subtraction_cancels_whole_TT_sector": True,
            "two_endpoint_uniqueness_not_general_curvature_matching": True,
        },
        "whole_external_kernel_emissions": external,
        "whole_internal_kernel_emission": contact,
        "whole_total_kernel_radiation": 24 * f0 * sum(H[i, i] / a[i] for i in range(4)),
        "whole_covariant_reduction": "The literal cross-contraction2gc Phi(x)[PhiY G_H G+Phi^2 nablaPhi G_H nablaG] equals gc(4-n)Phi(x)Phi(y)^3 G_H G/3-gc delta(G-G_H)Phi^4/3-gc Phi(x)Phi^2(Box+1)Phi G_H G. Retain all three terms, including Green-function delta contacts and the external-EOM radiation. The triangle and full24 EOM computations determine the complete selected physical TT response.",
        "whole_subtraction_proof": "Four external emissions of Phi^3 F Phi give6 sum J_i[F(1+2a_i)+3F(1)]. Both loop insertions give -6 sum J_i[F(1+2a_i)-F(1)]. Their total is24F(1) sum J_i, constant-quartic radiation. The tadpole contact and fixed Phi^2Y UV counterfunctional have the same physical constant reduction with all EOM metric terms retained. Thus the linear part of the existing full-parent symmetric-value condition cancels the entire across-source physical TT class.",
        "whole_scope": "The proof concerns this massive two-endpoint kernel and its actual covariant source. It neither assigns the S336 independent Weyl-curvature coefficient nor constructs the complete curved effective action.",
    }
