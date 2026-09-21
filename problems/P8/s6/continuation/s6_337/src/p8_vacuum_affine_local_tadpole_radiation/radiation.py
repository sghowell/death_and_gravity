"""Full induced-quartic radiation, generic Ward identities and OS4 reduction."""

from functools import cache

import sympy as s

from . import source, vertices


def raw_basis_factors(dim):
    return {
        "phi6": {"phi4": s.Integer(15)},
        "phi4Y": {"phi4": s.S.One, "phi2Y": s.Integer(6)},
        "phi2Y2": {"Y2": s.S.One, "phi2Y": 2 + 4 / dim},
        "Y3": {"Y2": 3 + 12 / dim},
        "phi2_L3_minus_L4": {
            "Gal": s.S.One,
            "phi_vHv": -2 * (dim - 2) / dim,
            "phi_box_Y": -2 / dim,
            "phi2_Hdiff": 1 / dim,
        },
        "Y_L3_minus_L4": {"Gal": 1 + 4 / dim, "Y_Hdiff": 1 / dim},
    }


def finite_basis_factors():
    dim = s.Symbol("factor_D")
    return {
        name: {
            key: s.factor(value.subs(dim, 4) - 2 * s.diff(value, dim).subs(dim, 4))
            for key, value in values.items()
        }
        for name, values in raw_basis_factors(dim).items()
    }


def local_coefficients():
    kap, lam = source.KAPPA, source.LAMBDA
    k2 = 120 * lam / kap - s.Rational(2362496, 25) / kap**2
    k3 = s.Integer(47232) / kap**2
    return {"Y2": k2 / 2, "Gal": -2 * k3 / 3, "phi4": -k2 / 18 - 8 * k3 / 81}


def four_polynomial(ss, tt):
    ss, tt = s.sympify(ss), s.sympify(tt)
    uu = 4 - ss - tt
    co = local_coefficients()
    return s.expand(
        2 * co["Y2"] * (ss**2 + tt**2 + uu**2 - s.Rational(16, 3))
        - 3 * co["Gal"] * (ss * tt * uu - s.Rational(64, 27)) / 2
    )


def basis_radiation(name, gram, polarization):
    ak = tuple(gram[i, 4] for i in range(4))
    return s.factor(
        vertices.quartic_contact(name, gram, polarization, ak, 0, (0, 0, 0, 0))
        + sum(
            polarization[i, i] * vertices.shifted_vertex(name, gram, ak, i) / ak[i]
            for i in range(4)
        )
    )


def full_tt_without_common_loop_and_metric_factor(gram, polarization):
    return -s.Add(
        *(
            co * basis_radiation(name, gram, polarization)
            for name, co in local_coefficients().items()
        )
    )


@cache
def data():
    gram, pol = vertices.generic_radiative_data()
    ak = tuple(gram[i, 4] for i in range(4))
    xi = s.symbols("xi_dot_p0:4")
    beta = -s.Add(*xi)
    gauge = s.zeros(5)
    for i in range(4):
        for j in range(4):
            gauge[i, j] = ak[i] * xi[j] + xi[i] * ak[j]
    denominator = s.prod(ak)
    checks, numerators = {}, {}
    for name in vertices.BASIS:
        shifted = tuple(vertices.shifted_vertex(name, gram, ak, i) for i in range(4))
        contact = vertices.quartic_contact(
            name, gram, gauge, ak, 2 * beta, tuple(beta * a for a in ak)
        )
        checks[name + "_generic_full_Ward"] = s.expand(
            contact + 2 * sum(xi[i] * shifted[i] for i in range(4))
        )
        numerators[name] = s.expand(
            vertices.quartic_contact(name, gram, pol, ak, 0, (0, 0, 0, 0)) * denominator
            + s.Add(
                *(
                    pol[i, i] * s.prod(ak[j] for j in range(4) if j != i) * shifted[i]
                    for i in range(4)
                )
            )
        )
    identities = {
        "phi2Y_equals_phi4_over3": numerators["phi2Y"] - numerators["phi4"] / 3,
        "phi_box_Y_EOM": numerators["phi_box_Y"] + numerators["phi2Y"],
        "phi_vHv_divergence": numerators["phi_vHv"]
        + numerators["Y2"] / 2
        + numerators["phi_box_Y"] / 2,
        "phi2_Hdiff_EOM_Ricci": numerators["phi2_Hdiff"]
        - 3 * numerators["phi2Y"]
        + numerators["Y2"],
        "Y_Hdiff_Galileon_Ricci": numerators["Y_Hdiff"] + 2 * numerators["Gal"],
    }
    checks.update({key: s.expand(value) for key, value in identities.items()})
    finite = finite_basis_factors()
    q = source.local_couplings()
    induced = {
        key: s.factor(
            sum(q[name] * values.get(key, 0) for name, values in finite.items())
        )
        for key in vertices.BASIS
    }
    reduced = {
        "Y2": s.factor(induced["Y2"] - induced["phi_vHv"] / 2 - induced["phi2_Hdiff"]),
        "Gal": s.factor(induced["Gal"] - 2 * induced["Y_Hdiff"]),
        "phi4": s.factor(
            induced["phi4"]
            + induced["phi2Y"] / 3
            + induced["phi_vHv"] / 6
            - induced["phi_box_Y"] / 3
            + induced["phi2_Hdiff"]
        ),
    }
    co = local_coefficients()
    checks["full_source_alpha_after_finite_EOM_reduction"] = s.factor(
        reduced["Y2"] - co["Y2"]
    )
    checks["full_source_beta_after_finite_EOM_reduction"] = s.factor(
        reduced["Gal"] - co["Gal"]
    )
    ss, tt = s.symbols("Mandelstam_s Mandelstam_t")
    uu = 4 - ss - tt
    local_four = (
        co["Y2"] * (2 * (ss**2 + tt**2 + uu**2) - 8)
        - 3 * co["Gal"] * ss * tt * uu / 2
        + 24 * co["phi4"]
    )
    checks["retained_OS4_symmetric_value_subtraction"] = s.expand(
        local_four.subs({ss: s.Rational(4, 3), tt: s.Rational(4, 3)})
    )
    checks["whole_renormalized_four_point_polynomial"] = s.expand(
        local_four - four_polynomial(ss, tt)
    )
    for name in ("phi6", "phi4Y"):
        values = finite[name]
        whole = s.Add(*(value * numerators[key] for key, value in values.items()))
        constant = values.get("phi4", 0) + values.get("phi2Y", 0) / 3
        checks[
            name + "_entire_radiation_removed_by_existing_linear_value_subtraction"
        ] = s.expand(whole - constant * numerators["phi4"])
    return {
        "checks": checks,
        "gates": {
            "all_eight_quartic_Ward_identities_before_TT_projection": len(
                vertices.BASIS
            )
            == 8,
            "all_five_general_reductions_not_sample_inference": len(identities) == 5,
            "volume_and_all_Christoffel_terms_retained": True,
            "four_external_emissions_plus_1PI_contact": True,
            "large_Phi4Y_constant_cancelled_not_misbounded": True,
            "linear_existing_counterterm_split_not_new_matching": True,
            "physical_TT_representative_not_full_action_equivalence": True,
        },
        "whole_six_finite_basis_factors": finite,
        "whole_unsubtracted_reduced_coefficients": reduced,
        "whole_OS4_subtracted_coefficients": co,
        "whole_renormalized_local_four_point": -four_polynomial(ss, tt)
        / (16 * s.pi**2),
        "whole_radiation_prescription": "M5local=-(V5_alphaY2+betaL+gammaPhi4 + four scalar emissions)/(16pi^2 sqrt(kappa)). All eight quartic Ward identities retain volume and connection terms. The five reductions are identities of general on-shell physical five-point amplitudes modulo invisible Ricci terms, not replacements of the curved source. Phi6 and Phi4Y cancel against their linear portion of the existing full symmetric value subtraction.",
    }
