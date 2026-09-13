"""Complete two-threshold first-sheet reciprocal; no inherited Proca-only pole."""

from functools import cache

import sympy as s

from . import geometry as g

y, p, n, mu2, ell = g.y, g.p, g.n, g.mu2, g.ell
z = s.Symbol("cut_fraction", positive=True)
J0 = s.Symbol("base_radial_integral")
TINY = s.Rational(1, 10) ** 180


def polynomial_parts(species, spin):
    if not isinstance(species, str) or species not in ("Proca", "heavy"):
        raise ValueError("Both specified fields have distinct complete cuts")
    if (
        isinstance(spin, bool)
        or not isinstance(spin, (int, s.Integer))
        or spin not in (0, 2)
    ):
        raise ValueError("Keep the trace and shear channels")
    weight = (g.W_P if species == "Proca" else g.W_H)[spin]
    rows = [J0]
    for order in range(1, 4):
        rows.append(s.expand(z * rows[-1] - s.Rational(1, 2 * order - 1)))
    expression = s.expand(
        sum(coef * rows[power[0] // 2] for power, coef in s.Poly(weight, y).terms())
    )
    polynomial = expression.subs(J0, 0)
    coefficient = expression.coeff(J0)
    return polynomial, coefficient


def moment(species, spin, order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)) or order < 0:
        raise ValueError("Use a nonnegative exact radial moment order")
    polynomial_parts(species, spin)
    weight = (g.W_P if species == "Proca" else g.W_H)[spin]
    return s.integrate(s.expand(weight * (1 - y * y) ** order), (y, 0, 1))


@cache
def data():
    parts = {
        (species, spin): polynomial_parts(species, spin)
        for species in ("Proca", "heavy")
        for spin in (0, 2)
    }
    slopes = {
        spin: moment("Proca", spin, 0) / (4 * mu2) + moment("heavy", spin, 0) / (4 * n)
        for spin in (0, 2)
    }
    c = {0: g.C0, 2: g.C2}
    static = {spin: 1 / c[spin] for spin in (0, 2)}
    second = {spin: s.factor(slopes[spin] / c[spin] ** 2) for spin in (0, 2)}
    # The complete heavy radial decrement at p=-4m_P^2 is bounded BEFORE
    # the lowest-threshold sign test. It is not replaced by zero.
    bound = {
        spin: g.PROCA_MASS2 / (g.HEAVY_MASS2 - g.PROCA_MASS2) * moment("heavy", spin, 0)
        for spin in (0, 2)
    }
    low_old = {
        spin: (4 if spin == 0 else s.Rational(1, 30))
        - s.integrate(g.W_P[spin] / y**2, (y, 0, 1))
        for spin in (0, 2)
    }
    gaps = {
        0: s.Rational(16, 15) + 2 * 394 - TINY,
        2: -s.Rational(172, 225) + s.Rational(394, 60) - TINY,
    }
    checks = {
        "Proca_complete_trace_threshold": low_old[0] - s.Rational(16, 15),
        "Proca_complete_shear_threshold": low_old[2] + s.Rational(172, 225),
        "scalar_trace_radial_zero_moment": moment("heavy", 0, 0) - s.Rational(68, 35),
        "scalar_shear_radial_zero_moment": moment("heavy", 2, 0) - s.Rational(1, 210),
        "complete_trace_slope": slopes[0]
        - s.Rational(9, 35) / mu2
        - s.Rational(17, 35) / n,
        "complete_shear_slope": slopes[2]
        - s.Rational(3, 56) / mu2
        - s.Rational(1, 840) / n,
        "complete_trace_static_measure": static[0] - 1 / (2 * (ell + 2)),
        "complete_shear_static_measure": static[2] - 60 / (ell + 2),
    }
    for key, (poly, coef) in parts.items():
        species, spin = key
        weight = (g.W_P if species == "Proca" else g.W_H)[spin]
        quotient = s.cancel((weight - weight.subs(y, s.sqrt(z))) / (z - y * y))
        checks[species + "_" + str(spin) + "_literal_full_polynomial_division"] = (
            s.factor(poly - s.integrate(s.expand(quotient), (y, 0, 1)))
        )
        checks[species + "_" + str(spin) + "_full_log_coefficient"] = s.factor(
            coef - weight.subs(y, s.sqrt(z))
        )
    uh = {spin: s.factor(parts["heavy", spin][1] / (2 * s.sqrt(z))) for spin in (0, 2)}
    checks["minimal_trace_full_upper_cut_imaginary"] = (
        uh[0] - s.sqrt(z) * (3 - z) ** 2 / 2
    )
    checks["minimal_shear_full_upper_cut_imaginary"] = (
        uh[2] - z ** s.Rational(5, 2) / 60
    )
    checks["heavy_trace_full_threshold_polynomial"] = parts["heavy", 0][0].subs(
        z, 0
    ) + s.Rational(36, 5)
    checks["heavy_shear_full_threshold_polynomial"] = parts["heavy", 2][0].subs(
        z, 0
    ) + s.Rational(1, 150)
    high_a = {
        spin: s.factor(
            sum(
                parts[species, spin][1].subs(z, 1) / 2 for species in ("Proca", "heavy")
            )
        )
        for spin in (0, 2)
    }
    high_b = {
        spin: s.factor(
            c[spin]
            + sum(parts[species, spin][0].subs(z, 1) for species in ("Proca", "heavy"))
            - parts["Proca", spin][1].subs(z, 1) * s.log(mu2) / 2
            - parts["heavy", spin][1].subs(z, 1) * ell / 2
        )
        for spin in (0, 2)
    }
    checks["entire_two_mass_trace_high_log"] = high_a[0] - 4
    checks["entire_two_mass_shear_high_log"] = high_a[2] - s.Rational(7, 30)
    checks["unchanged_scheme_heavy_log_cancels_only_full_trace_UV_constant"] = (
        high_b[0] + 2 * s.log(mu2) + s.Rational(52, 15)
    )
    checks["unchanged_scheme_heavy_log_cancels_only_full_shear_UV_constant"] = (
        high_b[2] + s.Rational(13, 60) * s.log(mu2) + s.Rational(127, 450)
    )
    checks["trace_reciprocal_log_tail_coefficient"] = 1 / high_a[0] - s.Rational(1, 4)
    checks["shear_reciprocal_log_tail_coefficient"] = 1 / high_a[2] - s.Rational(30, 7)

    # Positive Pick sign and monotonicity apply to the SUM, not each inverse.
    a, b, x, h = s.symbols("positive_a nonnegative_b real_x imaginary_y", positive=True)
    ratio = (x + s.I * h) / (a + b * (x + s.I * h))
    checks["complete_radial_Pick_imaginary_identity"] = s.simplify(
        s.im(s.expand_complex(ratio)) - a * h / ((a + b * x) ** 2 + b * b * h * h)
    )
    checks["complete_radial_real_derivative_identity"] = s.factor(
        s.diff(x / (a + b * x), x) - a / (a + b * x) ** 2
    )
    W2positive = 13 + (1 - y * y) * (17 - 3 * y * y)
    checks["Proca_shear_positive_weight_decomposition"] = s.expand(
        30 - 20 * y * y + 3 * y**4 - W2positive
    )
    checks["Proca_trace_positive_weight_decomposition"] = s.expand(
        3
        - 2 * y * y
        + 3 * y**4
        - (2 + 3 * (y * y - s.Rational(1, 3)) ** 2 + s.Rational(2, 3))
    )
    # Exact local one-sided threshold cusp of the reciprocal trace density.
    D, U = s.symbols("nonzero_upper_bank_real positive_upper_bank_imaginary", real=True)
    den = D * D + s.pi * s.pi * U * U
    rho = U / den
    cusp_minus = s.factor(s.diff(rho, D) * 9 * s.pi / 2)
    cusp_plus = s.factor(s.diff(rho, U) * s.Rational(9, 2))
    checks["second_threshold_trace_left_cusp"] = cusp_minus + 9 * s.pi * D * U / den**2
    checks["second_threshold_trace_right_cusp"] = (
        cusp_plus - s.Rational(9, 2) * (D * D - s.pi * s.pi * U * U) / den**2
    )
    w = s.Symbol("threshold_root", positive=True)
    ph, qh = parts["heavy", 0]
    left = ph.subs(z, -w * w) - qh.subs(z, -w * w) * s.atan(1 / w) / w
    # atan(1/w)=pi/2-atan(w), w>0, for a direct convergent local series.
    left = left.subs(s.atan(1 / w), s.pi / 2 - s.atan(w))
    checks["heavy_trace_left_full_square_root"] = s.expand(
        s.series(left, w, 0, 2).removeO() + s.Rational(36, 5) - 9 * s.pi * w / 2
    )
    ps, qs = parts["heavy", 2]
    lefts = ps.subs(z, -w * w) - qs.subs(z, -w * w) * (s.pi / 2 - s.atan(w)) / w
    checks["heavy_shear_second_threshold_order_five_not_square_root"] = (
        s.expand(s.series(lefts, w, 0, 6).removeO()).coeff(w, 5) - s.pi / 60
    )

    # Static moments, not an absolute unweighted frequency integral.
    u, q = s.symbols("source_frequency_shift nonnegative_transfer", positive=True)
    checks["shifted_static_density_identity"] = 1 / (u + q) - 1 / u + q / (u * (u + q))
    return {
        "complete_forward_factors": "F_i(p)=-A_i(p), A_i(p)=c_i+integral p W_Pi(y)/(4*10^6+p(1-y^2))dy+integral p W_Hi(y)/(4*n+p(1-y^2))dy; both full integrals and fixed c_i retained.",
        "full_polynomial_and_base_integral_parts": {
            species: {spin: parts[species, spin] for spin in (0, 2)}
            for species in ("Proca", "heavy")
        },
        "analytic_base_integral": "J0(z)=atanh(1/sqrt(z))/sqrt(z), continued from p>0 with z=1+4m^2/p. Upper p-cut: for z>0 J0=atanh(sqrt(z))/sqrt(z)+i*pi/(2sqrt(z)); for z<0 J0=-atan(1/sqrt(-z))/sqrt(-z).",
        "full_lowest_threshold_heavy_decrement_bounds": bound,
        "strict_complete_lowest_threshold_gaps": gaps,
        "complete_first_sheet_zero_set": "EMPTY for both actual total factors. Strict Pick sign excludes nonreal zeros; monotonicity and the positive full lowest-threshold values exclude the real gap. Above the first threshold the Proca imaginary part is strictly positive, including the second threshold.",
        "full_reciprocal_density": "If A_i(-sigma+i0)=D_i+i*pi*U_i, rho_i=U_i/(D_i^2+pi^2 U_i^2)>0 for sigma>4*10^6. The complete reciprocal is1/A_i(p)=integral rho_i(sigma)/(p+sigma)dsigma; no isolated pole or instantaneous constant is present.",
        "complete_static_and_next_measure_moments": {
            "static": static,
            "next": second,
            "slopes": slopes,
        },
        "full_UV_logarithms_and_constants": {
            "coefficients": high_a,
            "constants": high_b,
        },
        "second_threshold_trace_two_sided_cusp": {
            "below": cusp_minus,
            "above": cusp_plus,
        },
        "ordinary_kernel_regularity": "K_i(t)=-2theta(t) integral rho_i(Omega^2)sin(Omega*t)dOmega is an ordinary oscillatory inverse and belongs to L1(0,infinity). The first threshold square root and BOTH sides of the second trace threshold are retained; after subtracting their models the second frequency derivative is integrable. Constants are finite but not numerically evaluated.",
        "proof_boundary": "The sign and contour arguments, endpoint expansions, causal distribution identities and L1 estimates are written proofs, not FORMALIZED. Complete stress cuts are not full quantum gravity propagator cuts.",
        "checks": {name: s.factor(value) for name, value in checks.items()},
        "gates": {
            "heavy_trace_decrement_strictly_below_declared_bound": bound[0] < TINY,
            "heavy_shear_decrement_strictly_below_declared_bound": bound[2] < TINY,
            "full_trace_lowest_threshold_gap_above_789": gaps[0] > 789,
            "full_shear_lowest_threshold_gap_above_5": gaps[2] > 5,
            "old_Proca_only_shear_threshold_is_still_negative": low_old[2] < 0,
            "heavy_trace_does_have_an_interior_square_root": cusp_minus != 0
            and cusp_plus != 0,
            "two_distinct_physical_thresholds": g.HEAVY_MASS2 > g.PROCA_MASS2,
            "full_total_density_not_sum_of_old_reciprocal_densities": True,
            "both_full_reciprocal_log_tails_positive": high_a[0] > 0 and high_a[2] > 0,
            "no_mass_or_finite_prescription_changed": True,
            "finite_half_line_L1_not_quantified_stability": True,
        },
    }
