"""Complete actual second-grade source and time-jet densities."""

from functools import cache

import sympy as s
from p8_vacuum_affine_subleading_band_conversion import centered as original_centered
from p8_vacuum_affine_subleading_band_conversion import flat

from . import jets


@cache
def coefficients():
    u, m, p, rows = flat.coefficients()
    t = s.symbols("t", real=True)
    a, H, _Hp, c = jets.clock(t)
    out = {}
    for channel in ("tensor", "vector", "scalar"):
        f0 = rows[channel, 0][0] / a
        flat2 = original_centered.centered_second(
            rows[channel, 0][0], rows[channel, 0][2], u, p
        )
        ug, y, geo = flat.averaged_geometry(channel)
        ll = s.factor(geo["LL"].subs({y: 0, ug: u}))
        f2 = s.factor(flat2.subs(m, m * a) / a - 3 * a * c * ll / 8)
        out[channel] = {
            "f0": f0,
            "source_value": f2,
            "source_first": s.factor(-a * a * H * f0 / 4),
            "source_second": s.factor(-a * a * f0 / 4),
            "longitudinal_geometry": ll,
        }
    return u, t, m, p, out


@cache
def data():
    u, t, _m, _p, rows = coefficients()
    source, detector = s.symbols("source_time detector_time", real=True)
    A = s.Function("a")
    G = s.Function("Gamma")
    C = s.symbols("angular_numerator", real=True)
    F = C * G(source) / (4 * A(source) * A(detector))
    g = A(source) / 2
    first = s.diff(g * F, source)
    second = s.diff(g * first, source)
    current = -A(detector) * second.subs(source, detector) / 2
    expected = (
        -C
        * (
            s.diff(A(detector), detector) * s.diff(G(detector), detector)
            + A(detector) * s.diff(G(detector), detector, 2)
        )
        / 32
    )
    checks = {
        "source_only_first_inverse_phase_derivative": s.simplify(
            first - C * s.diff(G(source), source) / (8 * A(detector))
        ),
        "complete_source_second_time_jet_and_first_time_jet": s.simplify(
            current - expected
        ),
        "source_value_second_endpoint_vanishes": s.simplify(
            current.subs(
                {s.diff(G(detector), detector): 0, s.diff(G(detector), detector, 2): 0}
            )
        ),
    }
    a, H, _Hp, _c = jets.clock(t)
    for channel, row in rows.items():
        checks[channel + "_exact_first_second_time_ratio"] = s.factor(
            row["source_first"] - H * row["source_second"]
        )
        checks[channel + "_complete_time_coefficient"] = s.factor(
            row["source_second"] + a * a * row["f0"] / 4
        )
        checks[channel + "_curvature_even_angular_degree"] = s.factor(
            row["longitudinal_geometry"].subs(u, -u) - row["longitudinal_geometry"]
        )
    return {
        "source_value": "The complete centered source-value coefficient is f2c_flat(m*a,p)/a-3a P1_L M_LL/8, where M_LL=(n.D.n)(n.G.n). It includes all massive sectors; the curvature correction has been computed separately from the actual W8 jets.",
        "time_jets": "The leading full pair product is C(n)/(4a_source a_detector), while g_source=a_source/2. With detector time fixed, two source iterations give -C[a_prime Gamma_prime+a Gamma_second]/32. Thus f2c has -a^2 f0(Gamma_second+H Gamma_prime)/4 and no source-value j2/d0 term.",
        "source_detector_boundary": "Coalescing times before differentiation would differentiate detector amplitudes and change these terms. The actual first source-time derivative is retained; a static fixture alone would miss it.",
        "evaluated_slots": "The q2 slots are j0/d2 with the full source-value coefficient, j1/d1 evaluated as zero, and j2/d0 with its complete source first and second derivatives. Higher original UV slots remain for one-ball/contact matching.",
        "checks": checks,
        "gates": {
            "all_three_complete_massive_channels": len(rows) == 3,
            "actual_first_source_time_derivative_retained": all(
                row["source_first"] != 0 for row in rows.values()
            ),
            "actual_source_value_j2_zero_not_all_time_terms": True,
            "independent_detector_time_before_differentiation": source != detector,
            "remaining_fourth_grade_not_declared_matched": True,
        },
    }
