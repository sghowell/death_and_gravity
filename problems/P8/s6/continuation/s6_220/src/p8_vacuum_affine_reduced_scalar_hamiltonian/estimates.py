"""Explicit infrared-safe scalar norm estimates and finite-band classical comparison."""

from functools import cache

import sympy as s
from p8_vacuum_affine_full_spatial_remainder import assembly
from p8_vacuum_affine_prepared_ward_reconstruction import norms as old_norms
from p8_vacuum_affine_quantum_retuning import profile

from . import scalar

DETECTOR = s.Integer(20)
SOURCE = s.Integer(10) ** 20
LOCAL = s.Integer(10) ** 35
CURRENT = s.Integer(10) ** 117
RECONSTRUCTED_TAIL = s.Integer(10) ** 76
PHASE_MAP = s.Integer(400)
GENERATOR = s.Integer(1000)


@cache
def projection_data():
    h13 = old_norms.leibniz_bound(13, 7)
    c12 = old_norms.leibniz_bound(12, 28)
    raw = 4 * (6 + 4 * h13 + 2 * c12)
    y = s.Symbol("transfer_squared", nonnegative=True)
    source_ward = 56 + 16 + 168 + 24 + 2 * (8 + 24)
    detector_ward = 4 * (46 + 2 * 22)
    clock = 60
    local = (source_ward + detector_ward + clock) * s.Integer(10) ** 30
    full = DETECTOR * SOURCE * assembly.CURRENT + LOCAL
    tail = DETECTOR * SOURCE * assembly.TAIL
    checks = {
        "new_source_norm_exact_addition": raw
        - old_norms.data()["projection_unrounded"]["source"]
        - 12,
        "source_Ward_explicit_sum": source_ward - 328,
        "detector_Ward_explicit_sum": detector_ward - 360,
        "both_Ward_and_clock_contacts_sum": local - 748 * s.Integer(10) ** 30,
        "complete_scalar_response_unrounded": full - 4 * s.Integer(10) ** 116 - LOCAL,
        "reconstructed_spatial_regulator_unrounded": tail - 6 * s.Integer(10) ** 75,
        "real_detector_projector_weight_margin": s.expand(16 * (1 + y) - (14 + 2 * y))
        - s.expand(2 + 14 * y),
        "source_two_spatial_derivatives": s.expand((1 + y) ** 8 - (1 + y) ** 6 * y * y)
        - s.expand((1 + y) ** 6 * (1 + 2 * y)),
    }
    return {
        "scalar_input": "s=(n,zeta,b), zeta=v+delta n, b=div beta. The nonzero-P scalar Fourier shift is beta=-iP b/|P|^2. These coordinates are not an assertion that beta itself lies in an unweighted infrared L2 space.",
        "norms": "V03[s]^2=integral_I integral_R3 (1+|P|^2)^3 |s_hat|^2; U138[s]^2=sum_r0..13 integral_I integral_R3 (1+|P|^2)^8 |partial_t^r s_hat|^2. Components use the Euclidean norm. The unitary Fourier transform and original source/detector germs are fixed.",
        "primitive_bounds": "On the unit slab both endpoint primitives have L2 norm<=1, and their j-th derivatives are input derivatives of order j-1 for j>=1. With eta=I n, c=I b-|P|^2 I(a^-2 eta), ||c||<= (1+|P|^2)||s||. These formulas are used after exact Ward cancellation, never by bounding beta itself.",
        "projection": "M[Qsyn_D]<=20 V03[D], Z136[Qsyn_G]<=1e20 U138[G]. For the detector, 2sqrt3<4, |H|<=2 and a^-2<=1 give14+2|P|^2<=16(1+|P|^2). The source uses the original Cauchy H_j<=7j!4^j and (a^-2)_j<=28j!4^j, with sqrt14<4; its exact new raw coefficient is the old coefficient plus12.",
        "local_bound": "With A=-a^3 rho/2 and B=-a P/2, original first stress jets<1e30 imply |A|<=4e30, |A'|<=28e30, |B|<=1e30, |B'|<=3e30. The explicit ordered source contributes328e30, synchronous detector360e30 and additional clock chart60e30. Use |tr Qsyn|<=22(1+|P|^2)||s|| and its time derivative<=46(1+|P|^2) times the first source-jet norm. All products have at most two factors(1+|P|^2), controlled by V03 U138. No detector time derivative is used.",
        "complete_clock_scalar_comparison": "The full S219 reference Gaussian ADM response plus the extra clock-map one-current contact extends from smooth sources supported away from P=0 to the stated scalar Sobolev spaces, with |Rscalar(D,G)|<1e117 V03[D]U138[G]. Density and the bounded projector give an infrared extension, not a continuous direction-independent value at P=0.",
        "regulator_boundary": "Replacing only the synchronous spatial response by S219's original homogeneous-anchored approximation, and keeping exact covariant Ward/clock terms, gives error<1e76 V03 U138/K for K>=2000. This reconstructed approximation is not identified with the literal sharp-band lapse/shift response; no finite-band Ward identity is asserted.",
        "source_unrounded": raw,
        "local_unrounded": local,
        "response_unrounded": full,
        "reconstructed_tail_unrounded": tail,
        "checks": checks,
        "gates": {
            "source_rounding": raw < SOURCE,
            "detector_rounding": 16 < DETECTOR,
            "full_local_Ward_and_clock_rounding": local < LOCAL,
            "full_scalar_current_rounding": full < CURRENT,
            "reconstructed_spatial_tail_rounding": tail < RECONSTRUCTED_TAIL,
            "infrared_extension_not_unweighted_shift_norm": True,
            "extra_clock_contact_not_deleted": True,
            "exact_Ward_only_after_original_covariant_limit": True,
            "derivative_loss_and_different_domains_explicit": True,
        },
    }


def phase_matrices():
    mapping = s.Matrix(
        [
            scalar.lapse_numerator() / (2 * (scalar.J0 + scalar.dJ)),
            scalar.v
            + scalar.delta * scalar.lapse_numerator() / (2 * (scalar.J0 + scalar.dJ)),
            scalar.pv / 2,
        ]
    ).jacobian(scalar.Z)
    actual = mapping.subs(scalar.actual_coefficients(), simultaneous=True)
    return actual.subs(scalar.q, 0), actual.diff(scalar.q) / scalar.a**2


@cache
def phase_data():
    q = s.Symbol("q_nonnegative", nonnegative=True)
    eps = profile.EPS
    lower = s.Rational(1215, 800) * s.Rational(4, 5) ** 18 - 4 * eps
    lrow = 2 * q + s.Rational(273, 100)
    lam = 1 + q
    # Symmetric base Hessian row-sum upper bounds; the rank-one term is separate.
    rows = (
        2 * q + s.Rational(9, 100),
        q + s.Rational(13, 200),
        s.Rational(1, 20),
        s.S.One,
    )
    margins = [s.Poly(3 * lam - r, q) for r in rows]
    checks = {
        "lapse_linear_row_margin": s.expand(3 * lam - lrow) - (q + s.Rational(27, 100)),
        "rank_one_Hessian_upper": s.Integer(50) * 9 - 450,
        "full_Hessian_plus_weighted_damping_margin": s.expand(
            GENERATOR * lam**2 - (450 * lam**2 + 3 * lam + 6)
        )
        - s.expand(541 + 1097 * q + 550 * q * q),
        "phase_map_safe_row_sum": s.Integer(200)
        + 101
        + s.Rational(1, 2)
        - s.Rational(603, 2),
        "phase_detector_response_factor": PHASE_MAP * CURRENT
        - 4 * s.Integer(10) ** 119,
        "retuning_T_stress_bound": s.Rational(1)
        + 3 * s.Rational(1, 2)
        - s.Rational(5, 2),
    }
    return {
        "coefficient_bounds": "On I, |H|<=2, |Theta|<=2, delta<=1/2, |E|<=1, ell<=1/10 and |w|<=1/10. Fixed |A|<epsilon and |Tcorr|<(5/2)epsilon<1/100, while Jnew>1/100. Thus the lapse numerator coefficient row sum is<=2q+273/100<3(1+q), 1/(2Jnew)<50 and |n|<=200(1+|P|^2)|Z|.",
        "phase_input_maps": "The actual coefficient-sector map s=L0(t)Z+|P|^2 L2(t)Z has V03[s_D]<=400 V05[Z_D]. Set C13=sum_r0..13 sum_j0..r binom(r,j) sup_I(sum_entries abs(L0^(j))+sum_entries abs(L2^(j))). Then U138[s_G]<=C13 U13,10[Z_G]. This follows by the finite jet convolution bound. C13 is finite by smoothness and the positive Jnew pivot; no unproved numerical13-jet bound on the fixed reference stress is asserted.",
        "phase_pullback": "The linear comparison pullback obeys |Rphase|<1e120 C13 V05 U13,10. V05 uses five spatial derivatives, U13,10 thirteen time and ten spatial derivatives. The clock second-chart contact is included. This is not a same-space inverse estimate, and no canonical1/kappa suppression is transferred to reduced source variance or to an unspecified norm.",
        "classical_generator": "For Z=(v,sigma,pv,ps), K=Jcanonical Hess(Hred)-3H diag(0,0,1,1). The base symmetric Hessian norm is<=3(1+q); the rank-one lapse term is<=450(1+q)^2. Therefore ||K(t,P)||2<1000(1+|P|^2)^2, including the crossing and current QG1 coefficient retuning.",
        "compact_momentum_propagator": "For each mathematical Lambda<infinity and |P|<=Lambda, the classical coefficient-sector IVP has its unique fundamental solution with norm<=exp(1000(1+Lambda^2)^2|t-s|). The zero-initial-data retarded solution has L2-time operator bound at most the same exponential on this unit slab. This follows from matrix Gronwall/Volterra iteration and the exact regular Hamiltonian, not a frozen eigenvalue argument.",
        "compact_band_scope": "Lambda labels compact-momentum estimates, not a chosen physical EFT cutoff or a uniform continuum inverse. The exponential grows quartically in Lambda. No stability, quantum constraint inversion or full nonlinear background follows; derivative loss prevents using the weak Gaussian bound as a direct Neumann contraction.",
        "new_J_lower": lower,
        "linear_lapse_row_upper": lrow,
        "base_Hessian_row_bounds": rows,
        "base_Hessian_row_positive_margins": [r.as_expr() for r in margins],
        "checks": checks,
        "gates": {
            "actual_retuned_lapse_lower": lower > s.Rational(1, 100),
            "actual_fixed_A_small": eps < s.Rational(1, 100),
            "actual_fixed_T_small": s.Rational(5, 2) * eps < s.Rational(1, 100),
            "all_base_Hessian_row_margins_positive": all(
                all(c >= 0 for c in r.all_coeffs()) and r.eval(0) > 0 for r in margins
            ),
            "phase_map_rounding": s.Rational(603, 2) < PHASE_MAP,
            "phase_response_rounding": PHASE_MAP * CURRENT < s.Integer(10) ** 120,
            "source_C13_finite_not_numerically_assumed": True,
            "retuned_coefficient_sector_not_full_quantum_constraints": True,
            "finite_band_not_physical_cutoff": True,
            "no_continuum_inverse_or_stability_from_Gronwall": True,
        },
    }
