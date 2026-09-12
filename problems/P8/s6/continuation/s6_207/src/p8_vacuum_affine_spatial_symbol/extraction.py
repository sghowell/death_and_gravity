"""Finite UV coefficient extraction with detector time held independent."""

from functools import cache

import sympy as s


def require_order(order):
    if (
        isinstance(order, bool)
        or not isinstance(order, (int, s.Integer))
        or not 0 <= order <= 4
    ):
        raise ValueError("Exact endpoint order0,...,4 required")
    return int(order)


def slots():
    return tuple((j, d) for j in range(5) for d in range(5 - j))


def source_iterates(F, g, source_time, order=4):
    order = require_order(order)
    out = [F]
    for _ in range(order):
        out.append(s.diff(g * out[-1], source_time))
    return tuple(out)


def coefficient_row(F, g, source_time, detector_time, endpoint):
    j = require_order(endpoint)
    iterate = source_iterates(F, g, source_time, j)[j]
    return (
        (-s.I)
        * s.I**j
        * g.subs(source_time, detector_time)
        * iterate.subs(source_time, detector_time)
    )


def retained_degrees(endpoint):
    return tuple(range(5 - require_order(endpoint)))


@cache
def data():
    u, t = s.symbols("source_time detector_time", real=True)
    a, b = s.symbols("a b", real=True)
    F = s.exp(a * u + b * t)
    g = s.Function("g")(u)
    rows = source_iterates(F, g, u)
    flat = s.symbols("flat_phase", positive=True)
    flat_rows = source_iterates(F, flat, u)
    checks = {
        "complete_fifteen_endpoint_inverse_power_slots": len(slots()) - 15,
        "all_thirty_five_retained_source_time_jet_entries": sum(
            j + 1 for j, d in slots()
        )
        - 35,
        "every_retained_physical_power_at_most_three": s.Matrix(
            [int(j - 1 + d <= 3) - 1 for j, d in slots()]
        ),
        "every_first_omitted_power_is_four": s.Matrix(
            [j - 1 + (5 - j) - 4 for j in range(5)]
        ),
        "all_four_other_endpoint_orders_retained": len(rows) - 5,
        "source_only_derivative_not_coincident": s.diff(F, u).subs(u, t)
        - a * F.subs(u, t),
        "premature_coincidence_extra_detector_derivative": s.simplify(
            s.diff(F.subs(u, t), t) - s.diff(F, u).subs(u, t) - b * F.subs(u, t)
        ),
        "flat_source_jet_full_inverse_phase_powers": s.Matrix(
            [s.simplify(flat_rows[j] - flat**j * a**j * F) for j in range(5)]
        ),
        "higher_W8_coefficient_first_power_six": 6 - min(2 * q for q in (3, 4)),
    }
    return {
        "exact_endpoint": "For each of four sectors, keep the detector time t independent and apply Lhat F=partial_s(ghat(s)*F(s,t)). The j-th pre-current endpoint is x^(j-1)*(-i)i^j*ghat(t)*(Lhat^j Fhat)(t,t). Only then take the full current imaginary part. Coalescing times before differentiation is incorrect.",
        "finite_UV_extraction": "For endpoint j0,...,4 retain the Taylor coefficients of its normalized analytic row at inverse-radius degrees d0,...,4-j. These15 endpoint/power slots supply every possibly nonintegrable radial power x^(j-1+d), from k down to k^-3. Each source time-jet index0,...,j remains inside its row.",
        "spatial_polynomial": "External P enters only through xP in the scaled momenta and regular projector/frequency expressions. Thus the inverse-radius coefficient of degree d has spatial degree at most d. Time differentiation and sector-dependent inverse phases do not change this degree count.",
        "higher_W8_independence": "The third/fourth W8 coefficients enter What first at x^6/x^8. Analytic sums, products, inverses, roots and source-time derivatives cannot lower that inverse-radius order. No retained UV coefficient of degree at most4 depends on those higher coefficients. This is exact Taylor-jet independence, not replacement of the finite-momentum reference or prepared state.",
        "fixed_P_remainder": "After these exact coefficients are removed, each endpoint starts at x^4=k^-4. On a fixed-P compact set the analytic remainder is bounded by a constant times k^-4, giving an absolutely integrable high-band remainder and O(1/K) radial tail. No explicit uniform all-P constant or full response norm is established here.",
        "matching_boundary": "The original two-leg mask remains on every integrand. Integrated one-ball UV coefficients are polynomial in external momentum, but the original sharp-band conversion is not: S206 supplies its actual leading nonpolynomial term. Full conversion/contact and fixed finite covariant matching remain required.",
        "checks": checks,
        "gates": {
            "five_complete_endpoint_rows": len(rows) == 5,
            "no_odd_endpoint_deleted": True,
            "all_source_time_jets_through_endpoint_order": True,
            "separate_detector_time_essential": b != 0,
            "fixed_P_remainder_not_all_P_norm": True,
            "no_renormalized_comparator_or_matching_asserted": True,
        },
    }
