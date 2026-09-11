"""Complete canonical source, regular null limits and quantitative coefficient decay."""

from functools import cache

import sympy as s
from p8_vacuum_analytic_affine_parent import source as original


@cache
def data():
    K, z = s.symbols("kappa zeta", positive=True)
    Y, r, ry, rp, Q, Z, H, p = s.symbols(
        "Y r r_Y r_Phi Box_Phi Z_Phi Href Phi_mu", real=True
    )
    # Read the literal frozen source packet, rather than defining a replacement source.
    prior = original.data()["shifted_source_P8"]
    table = {str(x): x for x in prior.free_symbols}
    sub = {
        table["R"]: 1 + r / K,
        table["R_X"]: ry,
        table["R_u"]: rp / s.sqrt(K),
        table["positive_X"]: Y / K,
        table["Box_P8_u"]: Q / s.sqrt(K),
        table["uHu_P8"]: Z / K ** s.Rational(3, 2),
        table["H_clock"]: H,
    }
    actual = p / s.sqrt(K) * prior.subs(sub, simultaneous=True)
    RK = 1 + r / K
    expected = (
        p * r / K * (Q / Y - Z / Y**2)
        - 3 * p * r * H / K ** s.Rational(3, 2)
        + (3 * p * r / (K * K * RK) * (rp / 4 + ry * Z / (2 * Y)))
    )
    a, ay, ap = s.symbols("regular_a a_Y a_Phi", real=True)
    regular = {r: Y * Y * a, ry: 2 * Y * a + Y * Y * ay, rp: Y * Y * ap}
    V = p * a * (Y * Q - Z)
    B = -3 * p * Y * Y * a * H
    C = 3 * p * Y * Y * a * (Y * Y * ap / 4 + (2 * a + Y * ay) * Z / 2)
    regular_R = 1 + Y * Y * a / K
    U = V + B / s.sqrt(K) + C / (K * regular_R)
    J = U / s.sqrt(K * z)
    contact = U**2 / (2 * K)
    eps = s.Symbol("inverse_sqrt_kappa", positive=True)
    checks = {
        "literal_complete_source_canonical_chain": s.factor(actual - expected),
        "regular_all_Y_source_decomposition": s.factor(
            expected.subs(regular, simultaneous=True) - U / K
        ),
        "canonical_source_normalization": s.factor(s.sqrt(K / z) * (U / K) - J),
        "retained_source_squared_contact": s.factor(K * (U / K) ** 2 / 2 - contact),
        "nonzero_null_gradient_source": s.factor((U / K).subs(Y, 0) + p * a * Z / K),
        "nonzero_null_gradient_J": s.factor(J.subs(Y, 0) + p * a * Z / s.sqrt(K * z)),
        "canonical_source_vanishes": s.limit(J.subs(K, 1 / eps**2), eps, 0),
        "contact_vanishes": s.limit(contact.subs(K, 1 / eps**2), eps, 0),
        "canonical_source_leading_coefficient": s.simplify(
            s.limit((J / eps).subs(K, 1 / eps**2), eps, 0) - V / s.sqrt(z)
        ),
        "source_contact_leading_coefficient": s.simplify(
            s.limit((contact / eps**2).subs(K, 1 / eps**2), eps, 0) - V * V / 2
        ),
    }
    return {
        "literal_complete_canonical_source": expected,
        "regular_source_numerator": {"V_mu": V, "B_mu": B, "C_mu": C},
        "regular_U_mu": U,
        "canonical_J_mu": J,
        "complete_contact_per_component": contact,
        "norm_bound": "On a fixed compact canonical jet set let v,b,c bound Euclidean norms of V,B,C. Ustar=v+b/sqrt(kappa0)+2c/kappa0 bounds ||U|| for all kappa>=kappa0, since R>1/2. Thus ||J||<=Ustar/sqrt(kappa*zeta), and |kappa S^2/2|<=Ustar^2/(2kappa) at eta.",
        "full_fixed_vector_action": "-F(A)^2/4+A^2/(2zeta)-A.J+kappa S^2/2; A=sqrt(kappa*zeta) W, mass1000 fixed. The source-free vector is a free spectator in the classical decoupled limit.",
        "not_inferred": "No regulator-removed determinant, uniform on-shell Green norm, interacting-state limit or full contour theorem is inferred.",
        "checks": checks,
    }


@cache
def coefficients():
    K, K0 = s.symbols("kappa kappa0", positive=True)
    r, Y, a, ay = s.symbols("r Y regular_a a_Y", real=True)
    ry = Y * (2 * a + Y * ay)
    d4 = -7 * ry**2 / (4 * (K + r))
    d5 = ry**2 / ((K + r) * Y)
    # Strict bound R>1/2 gives these explicit, removable coefficients.
    b4 = 7 * ry**2 / (2 * K)
    b5 = 2 * Y * (2 * a + Y * ay) ** 2 / K
    return {
        "dependent_A4": d4,
        "dependent_A5": s.cancel(d5),
        "absolute_A4_bound": b4,
        "absolute_A5_bound": "2*abs(Y)*(2a+Y*a_Y)^2/kappa",
        "generic_finite_jet_bound": "For any fixed compact set strictly inside the canonical strip and any finite derivative order j, analytic f,r,a and the uniformly nonzero R denominator give finite C_j with ||Delta A4||_Cj+||Delta A5||_Cj<=C_j/kappa. The analogous source and contact norms are O(kappa^-1/2) and O(kappa^-1). Constants depend on the full compact jet set and j, not on kappa.",
        "checks": {
            "A4_exact_ratio_to_signed_majorant": s.factor(d4 + b4 / (2 * (1 + r / K))),
            "A5_exact_ratio_to_signed_majorant": s.factor(d5 - b5 / (2 * (1 + r / K))),
            "A5_quotient_analytic_at_null": s.cancel(
                d5 - Y * (2 * a + Y * ay) ** 2 / (K + r)
            ),
            "A5_null_limit": s.cancel(d5).subs(Y, 0),
            "A4_null_limit": d4.subs(Y, 0),
            "reciprocal_fixed_domain_factor": s.factor(
                1 / (K + r) - 1 / (K * (1 + r / K))
            ),
            "lower_domain_denominator": s.factor(K - K0 / 2 - K / 2 - (K - K0) / 2),
        },
    }
