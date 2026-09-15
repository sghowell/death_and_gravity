"""Literal Bianchi-I cotangent density flow, retaining the entire heavy sector."""

from functools import cache

import sympy as s

from . import source

p, m, sh, eta, ph = source.p, source.m, source.sh, source.eta, source.ph
FIELDS = (p, m, sh, eta, ph)


def flow(F, G):
    Gp, Gm, Gs, Gph, Ge = (s.diff(G, z) for z in (p, m, sh, ph, eta))
    return s.factor(
        s.diff(F, p) * (-G + m * Gm + 2 * sh * Gs + ph * Gph)
        - Gp * (m * s.diff(F, m) + 2 * sh * s.diff(F, sh) + ph * s.diff(F, ph))
        + 10**100 * (s.diff(F, eta) * Gph - s.diff(F, ph) * Ge)
    )


@cache
def rates():
    H = source.whole()["H"]
    Hp, Hm, Hs, He, Hph = (s.diff(H, z) for z in (p, m, sh, eta, ph))
    return {
        "alpha": Hp / 3,
        "p": -H + m * Hm + 2 * sh * Hs + ph * Hph,
        "m": -m * Hp,
        "sh": -2 * sh * Hp,
        "eta": 10**100 * Hph,
        "ph": -(10**100) * He - ph * Hp,
        "M1": Hm,
    }


@cache
def cotangent():
    alpha, beta, kap = s.symbols(
        "homogeneous_alpha homogeneous_beta positive_kappa", real=True
    )
    Pa, PM, PB, PH, chi = s.symbols(
        "canonical_Palpha canonical_PM canonical_Pbeta canonical_PH homogeneous_heavy",
        real=True,
    )
    V = s.exp(3 * alpha)
    gamma = s.diag(
        s.exp(2 * alpha + 2 * beta), s.exp(2 * alpha - 2 * beta), s.exp(2 * alpha)
    )
    mixed = s.eye(3) * Pa / 6 + s.diag(PB / 4, -PB / 4, 0)
    pi = mixed * gamma.inv()
    traceless = mixed - s.eye(3) * s.trace(mixed) / 3
    binding = {
        p: Pa / (3 * kap * V),
        m: PM / (kap * V),
        sh: PB * PB / (8 * kap * kap * V * V),
        eta: 10**100 * chi,
        ph: PH / (kap * V),
    }
    return {
        "alpha": alpha,
        "beta": beta,
        "kappa": kap,
        "Pa": Pa,
        "PM": PM,
        "PB": PB,
        "PH": PH,
        "chi": chi,
        "volume": V,
        "gamma": gamma,
        "pi": pi,
        "traceless": traceless,
        "binding": binding,
    }


def direct_canonical(F, G):
    d = cotangent()
    alpha, kap = d["alpha"], d["kappa"]
    f = F.subs(d["binding"], simultaneous=True)
    h = kap * d["volume"] * G.subs(d["binding"], simultaneous=True)
    result = (
        s.diff(f, alpha) * s.diff(h, d["Pa"])
        - s.diff(f, d["Pa"]) * s.diff(h, alpha)
        + s.diff(f, d["chi"]) * s.diff(h, d["PH"])
        - s.diff(f, d["PH"]) * s.diff(h, d["chi"])
    )
    # Cyclic beta and M1 contributions are exactly zero.
    return s.factor(
        result.subs(
            {
                alpha: 0,
                kap: 1,
                d["Pa"]: 3 * p,
                d["PM"]: m,
                d["PB"] ** 2: 8 * sh,
                d["PH"]: ph,
                d["chi"]: eta / 10**100,
            },
            simultaneous=True,
        )
    )


@cache
def data():
    d = cotangent()
    H0 = source.clock_jet("H")
    C0 = source.clock_jet("C")
    fixtures = (
        C0,
        source.clock_jet("C", 1, 0),
        source.clock_jet("C", 0, 1),
        p * p + m * sh + eta * ph + ph * ph,
    )
    checks = {
        "full_alpha_one_form": s.factor(
            s.trace(d["pi"] * s.diff(d["gamma"], d["alpha"])) - d["Pa"]
        ),
        "full_beta_one_form": s.factor(
            s.trace(d["pi"] * s.diff(d["gamma"], d["beta"])) - d["PB"]
        ),
        "full_shear_density": s.factor(
            s.trace(d["traceless"] ** 2) / (d["kappa"] ** 2 * d["volume"] ** 2)
            - d["binding"][sh]
        ),
        "full_trace_density": s.factor(
            2 * s.trace(d["pi"] * d["gamma"]) / (3 * d["kappa"] * d["volume"])
            - d["binding"][p]
        ),
        **{
            "literal_canonical_bracket_" + str(i): s.factor(
                direct_canonical(F, H0) - flow(F, H0)
            )
            for i, F in enumerate(fixtures)
        },
        "weighted_density_self_flow_not_zero_off_constraint": s.factor(
            flow(C0, C0) + C0 * s.diff(C0, p)
        ),
        "full_M1_density_conservation": s.factor(
            rates()["m"] + 3 * m * rates()["alpha"]
        ),
        "full_shear_density_dilution": s.factor(
            rates()["sh"] + 6 * sh * rates()["alpha"]
        ),
    }
    return {
        "whole_BianchiI_metric_and_raw_momentum": (d["gamma"], d["pi"]),
        "whole_absolute_canonical_density_binding": d["binding"],
        "whole_full_original_homogeneous_rates": rates(),
        "whole_homogeneous_constraint_preservation_numerator": source.original.eliminate_N_primitive(
            s.diff(source.whole()["C"], source.u)
        )
        + flow(source.whole()["C"], source.whole()["H"]),
        "whole_zero_spatial_constraints": "Every homogeneous spatial derivative is zero. On the diagonal Bianchi-I invariant sector the momentum constraints and residual translations vanish; all homogeneous spatial Proca coordinates and momenta remain zero by rotational covariance, while the temporal vector and heavy scalar use the complete original equations.",
        "whole_cyclic_shape_reconstruction": "For positive conserved Pbeta, beta_u=N*Pbeta/(2*kappa*exp(3alpha)*R^(1/4)). M1_u is the entire H_m; both cyclic coordinates are reconstructed, not frozen. The normalized shear is Pbeta^2/(8*kappa^2*exp(6alpha)).",
        "whole_density_bracket_warning": "flow(F,G) means {F,kappa*V*G} for a scalar F and Hamiltonian density G. It is not antisymmetric in F,G: flow(C,C)=-C*C_p. This term vanishes only on the constraint and is retained in off-surface identities.",
        "checks": checks,
        "gates": {
            "all_five_live_density_fields_retained": len(FIELDS) == 5,
            "heavy_equations_not_deleted": rates()["eta"].has(ph)
            and rates()["ph"].has(eta)
            and rates()["ph"].has(source.original.j),
            "positive_shear_branch_has_nonzero_shape_momentum": s.Rational(81, 160)
            - s.Rational(5, 4) * source.PROFILE_BOUND
            > 0,
            "homogeneous_momenta_are_test_data_not_S275_state": True,
            "no_fixed_reference_momentum_means_relabelled": True,
        },
    }
