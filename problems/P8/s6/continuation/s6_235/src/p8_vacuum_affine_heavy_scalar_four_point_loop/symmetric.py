"""Explicit symmetric on-shell first-loop contact, with sign and classical margin."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import model

from . import amplitude, zero_jet

NAME = "V2S-T1-OS4"
D, g = model.D, amplitude.g
z, xi, eta = amplitude.z, amplitude.xi, amplitude.eta
S0 = s.Rational(4, 3)
N = D + 2
CFIX = -g * g * (3 / D - 2 / D**2)
AFIX = CFIX + g * g / (N - S0)
TREE = CFIX + 3 * g * g / (N - S0)
LOWER = s.Rational(2, 3) + D * z
QC = (1 - z) ** 2 + N * z - S0 * (1 - z) ** 2 * xi * (1 - xi)
QD = QC - S0 * z * z * eta * (1 - eta)
LOG_BOUND = s.Integer(462)
B0, C0, D0 = s.symbols("B_symmetric C_symmetric D_symmetric", real=True)
LOOP = 3 * (AFIX**2 * B0 / 2 + 2 * AFIX * g * g * C0 + 2 * g**4 * D0) / (16 * s.pi**2)
CONTACT = -LOOP
VALLEY = model.VALLEY_QUARTIC.subs(model.g, g)
ABS_UPPER = 14793 * model.G2**2 / (192 * model.GAP**2)
VALLEY_RATIO_UPPER = 14793 * model.G2 * (model.GAP + 2) / (768 * (model.GAP - 1))
LOOP_TREE_RATIO_LOWER = 5985 * model.G2 * (3 * model.GAP + 2) / 2560


@cache
def data():
    A, B, y = s.symbols("positive_A positive_B centered_angle", positive=True)
    P1 = s.atan(y * s.sqrt(B / A)) / s.sqrt(A * B)
    P2 = y / (2 * A * (A + B * y * y)) + s.atan(y * s.sqrt(B / A)) / (
        2 * A * s.sqrt(A * B)
    )
    J1 = 2 * s.atan(s.sqrt(B / A) / 2) / s.sqrt(A * B)
    J2 = 1 / (2 * A * (A + B / 4)) + s.atan(s.sqrt(B / A) / 2) / (A * s.sqrt(A * B))
    r = s.Rational
    primitive = s.log(LOWER) / D
    box_primitive = (s.log(LOWER) + r(2, 3) / LOWER) / D**2
    point = amplitude.complete_loop().subs(
        {amplitude.svar: S0, amplitude.tvar: S0, amplitude.uvar: S0}
    )
    point = point.subs(
        {
            amplitude.n: N,
            amplitude.C: CFIX,
            amplitude.Bfun(S0): B0,
            amplitude.Cfun(S0): C0,
            amplitude.Dfun(S0, S0): D0,
        }
    )
    n = amplitude.n
    jet_B = zero_jet.B21
    checks = {
        "complete_on_shell_symmetric_constraint": 3 * S0 - 4,
        "full_tree_symmetric_exact_positive_value": s.cancel(
            TREE - 4 * g * g / (D * D * (3 * D + 2))
        ),
        "negative_fixed_channel_value": s.cancel(
            AFIX + 2 * g * g * (3 * D * D - 2) / (D * D * (3 * D + 2))
        ),
        "fixed_channel_strict_magnitude_lower": s.cancel(
            -AFIX / g**2
            - r(19, 10) / D
            - (3 * D * D - 38 * D - 40) / (10 * D * D * (3 * D + 2))
        ),
        "triangle_all_angle_lower_polynomial": s.expand(
            QC
            - LOWER
            - S0 * (1 - z) ** 2 * (xi - r(1, 2)) ** 2
            - r(2, 3) * z
            - r(2, 3) * z * z
        ),
        "box_all_angle_lower_polynomial": s.expand(
            QD
            - LOWER
            - S0 * (1 - z) ** 2 * (xi - r(1, 2)) ** 2
            - S0 * z * z * (eta - r(1, 2)) ** 2
            - r(2, 3) * z
            - r(1, 3) * z * z
        ),
        "triangle_upper_integral_primitive": s.cancel(s.diff(primitive, z) - 1 / LOWER),
        "box_upper_integral_primitive": s.cancel(
            s.diff(box_primitive, z) - z / LOWER**2
        ),
        "complete_symmetric_loop_from_six_boxes": s.cancel(point - LOOP),
        "new_finite_contact_cancels_only_symmetric_value": s.cancel(LOOP + CONTACT),
        "symmetric_upper_bracket_constant": r(9, 4) + 8 * LOG_BOUND - r(14793, 4),
        "symmetric_negative_bracket_constant": r(9, 4)
        + 2 * LOG_BOUND
        - 2 * r(19, 10) * 375
        + r(1995, 4),
        "exact_valley_margin_ratio": s.cancel(
            (14793 * g**4 / (192 * D * D)) / (24 * VALLEY)
            - 14793 * g * g * (D + 2) / (768 * (D - 1))
        ),
        "exact_large_loop_to_tree_lower_ratio": s.cancel(
            (5985 * g**4 / (640 * D * D)) / TREE - 5985 * g * g * (3 * D + 2) / 2560
        ),
        "off_shell_B21_lower_formula": s.cancel(
            jet_B - n * (s.log(n) - 1 + 1 / n) / (n - 1) ** 2
        ),
        "actual_Csym_lower_rational_margin": s.Rational(39300, 101)
        - 375
        - s.Rational(1425, 101),
        "first_complete_angle_primitive": s.simplify(
            s.diff(P1, y) - 1 / (A + B * y * y)
        ),
        "second_complete_angle_primitive": s.simplify(
            s.diff(P2, y) - 1 / (A + B * y * y) ** 2
        ),
        "first_angle_integral_endpoints": s.simplify(
            P1.subs(y, r(1, 2)) - P1.subs(y, -r(1, 2)) - J1
        ),
        "second_angle_integral_endpoints": s.simplify(
            P2.subs(y, r(1, 2)) - P2.subs(y, -r(1, 2)) - J2
        ),
        "first_angle_removable_endpoint": s.limit(J1, B, 0) - 1 / A,
        "second_angle_removable_endpoint": s.limit(J2, B, 0) - 1 / A**2,
    }
    return {
        "new_prescription_name": NAME,
        "explicit_new_finite_condition": "Extend the S234 base scheme by deltaC_fin=-A1_base(s0,s0,s0), s0=4/3. This on-shell crossing-symmetric point is subthreshold, not real2-to2 scattering kinematics. The condition fixes the four-point VALUE through first loop order; it does not fix higher derivatives or all-angle matching. No frozen S233/S234 prescription is edited.",
        "complete_symmetric_tree": TREE,
        "complete_symmetric_first_loop": LOOP,
        "defined_finite_contact": CONTACT,
        "actual_parameter_integral_domain": {
            "Q_C": QC,
            "Q_D": QD,
            "positive_common_lower": LOWER,
        },
        "uniform_subthreshold_integral_bounds": "All angle parameters lie in[0,1]. Both Q_C,Q_D>=2/3+D z. Hence0<B_symmetric<1/2, 0<C_symmetric<ln(1+3D/2)/D, 0<D_symmetric<ln(1+3D/2)/D^2. The actual logarithm is below462. These bounds apply to complete convergent parameter integrals at the stated point.",
        "strict_negative_loop_proof": "The actual n=M_H^2>10^197 gives ln(n)>394. Since C_symmetric>B21(n)>393/n>375/D and abs(A_symmetric)>(19/10)g^2/D, while abs(A_symmetric)<3g^2/D, the complete loop bracket is below-(1995/4)g^4/D^2. Thus A1_base(s0,s0,s0)<-5985g^4/(64pi^2 D^2)<0.",
        "large_unadjusted_loop_to_tree": "The positive complete tree at the point is4g^2/[D^2(3D+2)]. The unadjusted S234 first-loop correction has magnitude MORE THAN10^190 times this finely cancelled tree value. This is a failure of unadjusted value matching, not a proof of strong coupling, a full parent exclusion or impossible quantum renormalization.",
        "actual_loop_tree_ratio_lower": LOOP_TREE_RATIO_LOWER,
        "complete_absolute_loop_upper": ABS_UPPER,
        "finite_contact_and_classical_quartic_margin": "The new finite contact is positive and less than10^-6 times24q, q=g^2(D-1)/[6D^2(D+2)]. A CLASSICAL CONTACT-ONLY COMPARISON replacing C by C+deltaC_fin keeps the full heavy-square potential nonnegative and coercive with remaining quartic above(1-10^-6)q. This is not a quantum effective-potential theorem or a claim that all finite source/mass counterterms preserve the old bare minimum.",
        "actual_contact_to_valley_ratio_upper": VALLEY_RATIO_UPPER,
        "one_dimensional_integrals": {"J1": J1, "J2": J2},
        "reduced_integral_and_enclosure": "At s0 set a=(1-z)^2+n z,b=s0(1-z)^2,A=a-b/4,B=b. Then C_symmetric=integral(1-z)J1 dz and Dbar(s0,0)=integral z(1-z)J2 dz. The EXACT D_symmetric is retained between Dbar(s0,0) and Dbar(s0,0)/(1-s0/(4n))^2, not replaced by its lower bound. Log coordinate log[1+(n-1)z] resolves the heavy endpoint layer for independent diagnostics.",
        "boundary": "One explicit finite first-loop matching condition and a classical contact margin do not bound the remaining physical-angle amplitude, match b20/b21 to full quantum accuracy, control omitted loops, construct an exact UV theory or supply original V/G/B/P8.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "actual_large_gap_and_channel_lower_margin": model.GAP > 200
            and 3 * model.GAP**2 - 38 * model.GAP - 40 > 0,
            "actual_mass_logarithm_lower_domain": model.MASS2 > 10**197,
            "actual_common_logarithm_upper_domain": 1 + 3 * model.GAP / 2 < 10**198,
            "actual_positive_finite_contact_fits_valley": 0
            < VALLEY_RATIO_UPPER
            < s.Rational(1, 10**6),
            "actual_unadjusted_loop_exceeds_tree_by_named_factor": LOOP_TREE_RATIO_LOWER
            > 10**190,
            "actual_heavy_angle_enclosure_denominator_positive": 0
            < S0 / (4 * model.MASS2)
            < 1,
            "new_finite_condition_not_original_prescription_change": True,
            "contact_only_comparison_not_full_quantum_potential": True,
        },
    }
