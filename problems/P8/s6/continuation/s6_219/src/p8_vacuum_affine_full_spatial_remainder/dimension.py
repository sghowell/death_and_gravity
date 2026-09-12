"""Full trace geometries on the original complex-dimensional endpoint domain."""

from functools import cache

import sympy as s
from p8_vacuum_affine_corrected_spatial_current import subtraction
from p8_vacuum_affine_matched_spatial_current import dimension as continued
from p8_vacuum_affine_ordered_scalar_symbol import density, geometry
from p8_vacuum_affine_uniform_uv_remainder import domain

SIGMA = s.Rational(1, 4)
eta = s.Symbol("dimension_offset", real=True)


def rational_bound(value, variables, radii):
    shifted = s.cancel(value.subs(geometry.d, 3 + eta))
    numerator, denominator = s.fraction(shifted)
    vs, rs = (eta, *variables), (SIGMA, *radii)
    denpoly, numpoly = s.Poly(denominator, *vs), s.Poly(numerator, *vs)
    if not all(c.is_Rational for c in denpoly.coeffs() + numpoly.coeffs()):
        raise ValueError("Require exact rational geometry coefficients")
    bound = lambda poly: sum(
        abs(c) * s.prod(radius**power for radius, power in zip(rs, exponents))
        for exponents, c in poly.terms()
    )
    constant = abs(denpoly.coeff_monomial(1))
    defect = bound(denpoly) - constant
    floor = constant - defect
    if floor <= 0:
        raise ValueError("The stated geometry domain has no proved denominator floor")
    return {
        "upper": s.factor(bound(numpoly) / floor),
        "denominator_floor": floor,
        "reconstruction": s.factor(shifted - numpoly.as_expr() / denpoly.as_expr()),
    }


@cache
def normalized_geometries(channel):
    """Real momentum directions, before analytic transverse moments."""
    D, td, G, tg, norm = geometry.channel_directions(channel)
    extra = geometry.d - 4
    tauD, tauG = s.trace(D) + extra * td, s.trace(G) + extra * tg
    BD, BG = D - tauD * s.eye(4) / 2, G - tauG * s.eye(4) / 2
    bd, bg = td - tauD / 2, tg - tauG / 2
    c, h, u, v = s.symbols("c h u v", real=True)
    X, Y, Z = geometry.X, geometry.Y, geometry.Z
    n, w = s.Matrix([c * X, c * Y, c * Z, u]), s.Matrix([h * X, h * Y, h * Z, v])
    P, Q = s.eye(4) - n * n.T, s.eye(4) - w * w.T
    MD, MG = geometry.magnetic(n, w, D, tauD), geometry.magnetic(n, w, G, tauG)
    md, mg = -(n.T * D * w)[0] - n.dot(w) * bd, -(n.T * G * w)[0] - n.dot(w) * bg
    out = {
        "00": s.trace(P * BD * Q * BG) + extra * bd * bg,
        "01": s.trace(P * BD * Q * MG.T) + extra * bd * mg,
        "10": s.trace(P * MD * Q * BG) + extra * md * bg,
        "11": s.trace(P * MD * Q * MG.T) + extra * md * mg,
        "TL": (w.T * BD * P * BG * w)[0],
        "LT": (n.T * BD * Q * BG * n)[0],
        "LL": (n.T * BD * w)[0] * (n.T * BG * w)[0],
        "LC": tauG * (n.T * BD * w)[0],
        "CL": tauD * (n.T * BG * w)[0],
        "CC": tauD * tauG,
    }
    return (X, Y, Z, c, h, u, v), {
        key: s.cancel(value / norm) for key, value in out.items()
    }


@cache
def normalized_reconstruction():
    checks = {}
    R = s.sqrt(1 - 2 * geometry.u * geometry.y + geometry.y**2)
    powers = {
        "00": 0,
        "01": 1,
        "10": 1,
        "11": 2,
        "TL": 0,
        "LT": 0,
        "LL": 0,
        "LC": 1,
        "CL": 1,
        "CC": 2,
    }
    for channel in density.CHANNELS:
        variables, values = normalized_geometries(channel)
        _X, _Y, _Z, c, h, u, v = variables
        replacement = {c: 1, h: -1 / R, u: geometry.u, v: (geometry.y - geometry.u) / R}
        checks[channel] = s.Matrix(
            [
                s.factor(
                    geometry.azimuth(
                        s.cancel(
                            value.subs(replacement, simultaneous=True)
                            * R ** powers[key]
                        )
                    )
                    - geometry.contractions(channel)[key]
                )
                for key, value in values.items()
            ]
        )
    return checks


@cache
def geometry_bounds():
    out, checks = {}, {}
    for channel in density.CHANNELS:
        for key, value in geometry.contractions(channel).items():
            row = rational_bound(value, (geometry.y, geometry.u), (domain.EPS, 1))
            out[channel + "_far_" + key] = row["upper"]
            checks[channel + "_far_" + key] = row["reconstruction"]
        variables, values = normalized_geometries(channel)
        for key, value in values.items():
            row = rational_bound(value, variables, (1,) * len(variables))
            out[channel + "_near_" + key] = row["upper"]
            checks[channel + "_near_" + key] = row["reconstruction"]
    inverse = geometry.invariant_matrix().inv()
    for i in range(6):
        for j in range(6):
            row = rational_bound(inverse[i, j], (), ())
            out[f"inverse_{i}_{j}"] = row["upper"]
            checks[f"inverse_{i}_{j}"] = row["reconstruction"]
    return out, checks


def moment_rows():
    out = {}
    for total in range(5):
        for A in range(total + 1):
            for B in range(total - A + 1):
                C = total - A - B
                out[A, B, C] = (
                    s.rf(s.Rational(1, 2), A)
                    * s.rf(s.Rational(1, 2), B)
                    * s.rf(SIGMA / 2, C)
                    / s.rf(1 - SIGMA / 2, total)
                )
    return out


@cache
def data():
    bounds, exact = geometry_bounds()
    mode = continued.mode_data()
    moments = moment_rows()
    checks = {
        **{
            key + "_normalized_all_momentum_geometry_reconstruction": value
            for key, value in normalized_reconstruction().items()
        },
        "all_sixty_full_scalar_geometry_and_36_inverse_denominators": s.Matrix(
            list(exact.values())
        ),
        "all_35_grouped_transverse_moment_bounds": len(moments) - 35,
        "full_complex_dimension_canonical_rate_bound": 3
        * ((1 + SIGMA) / 2 + s.Rational(101, 100))
        - s.Rational(981, 200),
        "same_far_radial_integrable_exponent": -2 + SIGMA + s.Rational(7, 4),
        "same_angular_integrable_exponent": -SIGMA / 2 + s.Rational(1, 8),
        "same_near_external_degree": 5 + SIGMA - s.Rational(21, 4),
        "fixed_six_invariant_basis_determinant": s.factor(
            geometry.invariant_matrix().det()
            - geometry.d**3 * (geometry.d - 1) ** 3 / 2
        ),
    }
    return {
        "explicit_geometry_majorants": bounds,
        "grouped_transverse_moment_majorants": moments,
        "domain": "Keep |d-3|<=1/4 and both original complex time discs. The exact30 full scalar far geometries on |y|<=1/100,|u|<=1 and every36 inverse-basis entry have nonzero denominator floors from their rational coefficient majorants. This is coefficient-wise analytic continuation, not a positive noninteger-dimensional Hilbert norm.",
        "all_real_momenta_geometry": "Normalize the two real internal directions separately: n=(cX,cY,cZ,u), w=(hX,hY,hZ,v), where |c|,|h|,|u|,|v|<=1 and X,Y,Z share the transverse orientation. All30 normalized scalar geometries are polynomials divided only by safe dimension factors. Their absolute coefficient bounds are explicit. The grouped Z moments have numerator(eta/2)_C, bounded by(1/8)_C; every transverse denominator has floor(7/8)_degree. All35 moments through total degree4 are bounded by1, so the same polynomial bounds survive the analytic azimuth average. At either zero internal momentum use any real direction; mass keeps amplitudes bounded and the point has measure zero.",
        "full_constraint_amplitude": "The normalized full C term contains rk r_l pk pl/(2 a^2 omega_k omega_l). For real internal momenta on the complex clock discs each r/(a omega) is bounded by2, so it is at most2|pk pl|. The continued longitudinal canonical rate is below981/200<5, the original transverse rate is smaller, and all four W8 defect bounds keep the modes and inverse summed phases nonzero. Thus LC,CL,CC add no momentum degree and no singular small-leg denominator.",
        "both_retarded_branches": "Use the same-d analytic pair(Fplus-Fsharp)/(2i), not conjugation of d. S217's all-four-order bound gives each frequency's real-time phase defect<1/40 and both phases have modulus<20/19<2. The state correction and finite time remainder themselves remain in physical dimension as the original prescription requires.",
        "dominated_limit": "The full normalized pair amplitudes, all ten analytic geometries and inverse-basis entries have finite d-uniform majorants. Applying the same complete source-jet Cauchy recurrence gives far C_j U^(5-j) r^-4, U=m+|P|. Multiplying by the dimension factor yields radial exponent-7/4 and growth U^(4+1/4-j). The unexpanded low and near raw/Taylor/log rows have worst degree21/4<6. Angular weight has envelope(1-u^2)^(-1/8). Thus the full scalar Q dimension limit commutes with all original integrals and equals its physical remainder. No numerical physical constant is inferred from this qualitative uniform bound; the physical constants come from the full ten-field contraction proof.",
        "checks": checks,
        "gates": {
            "all_96_new_geometry_and_basis_bounds": len(bounds) == 96,
            "all_exact_rational_majorants_nonnegative": all(
                v >= 0 for v in bounds.values()
            ),
            "all_grouped_transverse_moments_controlled": all(
                v <= 1 for v in moments.values()
            ),
            "continued_full_canonical_rates_below_five": s.Rational(981, 200) < 5,
            "all_four_real_domain_W8_defects": all(
                row["all_real_relative_W8_defect"] < s.Rational(1, 1000)
                for row in mode["rows"].values()
            ),
            "all_four_far_domain_W8_defects": all(
                row["far_relative_scaled_W8_defect"] < s.Rational(1, 1000)
                for row in mode["rows"].values()
            ),
            "both_retarded_branch_phase_bounds": all(
                subtraction.dimension_data()["gates"].values()
            ),
            "six_spatial_derivative_weight_still_suffices": 5 + SIGMA < 6,
            "angular_and_far_radial_envelopes_integrable": -SIGMA / 2 > -1
            and -2 + SIGMA < -1,
        },
    }
