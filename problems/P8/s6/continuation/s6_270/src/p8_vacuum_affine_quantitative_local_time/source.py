"""Whole source, primitive and mixed complex-lapse/real-time jet bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth import background as current
from p8_vacuum_affine_heavy_scalar_parent import clock as heavy
from p8_vacuum_affine_heavy_source_filtration import source as heavy_source
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical as parent
from p8_vacuum_affine_nonlinear_lapse_branch import source as bounce
from p8_vacuum_affine_physical_background_vertices import parent as physical

u, N = current.u, parent.N
R = s.Function("whole_time_source_R")(u, N)
F = s.Function("whole_time_source_F")(u, N)
j = s.Function("whole_time_normalized_source")(u, N)
primitive = s.Function("whole_time_primitive")(u, N)
Hclock = s.Function("whole_original_clock_Hubble")(u)
COORDS = bounce.COORDS
p, G, dp, ph, eta, sh, el, ma, wm, gm, gh, curv = COORDS
mu = bounce.mu
TIME = s.Rational(1, 10**2000)
LAPSE_RADIUS = s.Rational(1, 10**245)
IMAGE = s.Rational(1, 10**260)
JET = s.Integer(10) ** 30
U = R ** -s.Rational(3, 4)
IN = 3 * R ** -s.Rational(7, 4) * s.diff(R, u) * s.diff(R, N) / (4 * N)
B = -U * s.diff(R, u) / (2 * N) - primitive
Fhat = U * (F + 9 * s.diff(R, u) ** 2 / (16 * R * N**2)) - s.diff(primitive, u) / N
lower = -3 * Hclock * (R - 1) / N
NORMAL = (
    -3 * (p - B - (R - 1) * G) ** 2 / (4 * R ** s.Rational(1, 4))
    - Fhat
    + G**2 / (2 * U)
    - lower * G
    + 2 * sh * R ** -s.Rational(1, 4)
    + ((s.Rational(1, 10) + dp) ** 2 + ph**2) / (2 * U)
    + R ** -s.Rational(1, 4) * (gm + gh) / 2
    + U * (mu * eta**2 / 2 - j * eta / 10**100)
    - R ** s.Rational(3, 4) * curv / 2
    + el * R ** s.Rational(1, 4) / 2
    + ma * R ** s.Rational(1, 4) / 4
    + R ** -s.Rational(1, 4) * wm / 2
)
HAMILTONIAN = N * NORMAL
TEMPORAL = (
    -3 * (R - 1) * (p - B - (R - 1) * G) / (2 * R ** s.Rational(1, 4)) + lower - G / U
)


@cache
def eliminate_N_primitive(expression):
    rules = {}
    for term in expression.atoms(s.Derivative):
        if term.expr != primitive:
            continue
        count = dict(term.variable_count)
        if count.get(N, 0):
            rules[term] = s.diff(IN, u, count.get(u, 0), N, count[N] - 1)
    return expression.xreplace(rules)


CONSTRAINT = eliminate_N_primitive(s.diff(HAMILTONIAN, N))


def parent_binding():
    return {
        parent.R: R,
        parent.B: B,
        parent.F: Fhat,
        parent.j: j,
        parent.Href: Hclock,
        parent.p: p,
        parent.G: G,
        parent.pm: s.Rational(1, 10) + dp,
        parent.ph: ph,
        parent.h: eta / 10**100,
        parent.n: mu * 10**200,
        parent.shear: sh,
        parent.electric: el * parent.zeta,
        parent.magnetic: ma / parent.zeta,
        parent.wmass: wm,
        parent.gm: gm,
        parent.gh: gh,
        parent.curvature: curv,
    }


@cache
def complex_source_bounds():
    packet = current.source_data()
    radius, distance = s.Rational(1, 10000), s.Rational(11, 10000)
    delta = distance * (2 + distance) / (1 - distance) ** 2
    switch = packet["full_switch_modulus_base"]
    heavy_packet = heavy.data()
    suppressed = s.factorial(8) * 2**8 / heavy.family.LOCALIZER**8
    entire_heavy = heavy_packet["common_complete_coefficient_prefactor"] * suppressed
    constant, complement = 10**301 * switch, 10**207 * switch
    scalar = (
        packet["complete_complex_tree_scalar_modulus_bound"]
        + constant
        + complement
        + entire_heavy
    )
    metric = 1 + delta / (1 - distance**2) ** 3 + 8 * switch + entire_heavy
    normalized_source = entire_heavy + 10**25 * switch / 10**400
    rows = {}
    for i in range(6):
        for k in range(6 - i):
            cauchy = s.factorial(i) * s.factorial(k) / radius ** (i + k)
            # Only N is complexified for the original smooth profile.
            profile = (
                8 * current.prior_initial.PROFILE_BOUND * s.factorial(k) / radius**k
            )
            rows[(i, k)] = (
                cauchy * metric,
                cauchy * scalar + profile,
                cauchy * normalized_source,
            )
    return {
        "packet": packet,
        "radius": radius,
        "distance": distance,
        "X_image_distance": delta,
        "entire_heavy": entire_heavy,
        "constant": constant,
        "complement": complement,
        "scalar": scalar,
        "metric": metric,
        "normalized_source": normalized_source,
        "rows": rows,
    }


@cache
def data():
    source = complex_source_bounds()
    full = parent.full_trace()
    mapped = full["whole_reduced_Hamiltonian"].subs(parent_binding(), simultaneous=True)
    temporal = full["whole_temporal_solution"].subs(parent_binding(), simultaneous=True)
    rules = {
        current.R: R,
        current.F: F,
        current.N: N,
        s.diff(current.R, current.u): s.diff(R, u),
        s.diff(current.R, current.N): s.diff(R, N),
        current.primitive: primitive,
        s.diff(current.primitive, current.u): s.diff(primitive, u),
    }
    germs = physical.source_germs()
    X = physical.X
    remainder = germs["entire_R_difference"]
    clock = 1 + (N**-2 - 1) / (1 + u * u) ** 3
    checks = {
        "literal_whole_original_time_Hamiltonian": s.factor(mapped - HAMILTONIAN),
        "literal_whole_original_time_temporal_root": s.factor(temporal - TEMPORAL),
        "literal_original_primitive_N_contact": s.factor(
            current.IN.xreplace(rules) - IN
        ),
        "literal_original_B_with_full_primitive": s.factor(
            current.B.xreplace(rules) - B
        ),
        "entire_source_R_clock_zero": s.simplify(remainder.subs(X, 1)),
        "entire_source_first_R_clock_jet_zero": s.simplify(
            s.diff(remainder, X).subs(X, 1)
        ),
        "actual_full_R_clock_value": s.factor(clock.subs(N, 1) - 1),
        "actual_full_R_first_N_clock_value": s.factor(
            s.diff(clock, N).subs(N, 1) + 2 / (1 + u * u) ** 3
        ),
        "actual_full_normalized_source_clock_value": heavy_source.physical_source(
            u, X
        ).subs(X, 1),
        "same_kappa": heavy_source.KAPPA - 10**800,
        "same_heavy_mass_ratio": parent.source.MASS2 / 10**200
        - bounce.parent.source.MASS2 / 10**200,
        "whole_twelve_invariant_count": s.Integer(len(COORDS) - 12),
    }
    return {
        "whole_original_time_Hamiltonian_in_complete_density_invariants": HAMILTONIAN,
        "whole_original_time_auxiliary_constraint": CONSTRAINT,
        "whole_original_time_temporal_solution": TEMPORAL,
        "whole_original_primitive_N": IN,
        "whole_primitive_boundary": "I(u,1)=0 at every real time. Pure time derivatives are actual N integrals of the corresponding IN derivatives. Only N derivatives are replaced by the exact IN identity. Iu is not set to zero off the clock or off the bounce.",
        "whole_actual_function_bindings": source["packet"][
            "entire_current_function_bindings"
        ],
        "whole_actual_fixed_profiles": source["packet"]["whole_fixed_profile_bindings"],
        "whole_same_full_normalized_heavy_source": source["packet"][
            "complete_normalized_heavy_source"
        ],
        "whole_evaluated_real_time_and_complex_lapse_radii": [TIME, LAPSE_RADIUS],
        "whole_complex_source_Cauchy_radii": [source["radius"], source["distance"]],
        "whole_complex_X_image_distance": source["X_image_distance"],
        "whole_all_21_mixed_source_jet_ceilings": {
            str(key): [JET, JET, 1] for key in source["rows"]
        },
        "whole_source_jet_proof": "Retain the entire original analytic tree, switches, all fixed vacuum constants and heavy localizer on the S256 joint neighborhood. Its scalar analytic modulus is<1e7, metric modulus<2 and full normalized source modulus<1e-100. Cauchy radius1e-4 gives all mixed total-order-five RF jets<1e30 and j jets<1. For the original smooth reference profiles use their REAL time C5 bound and N-only Cauchy; no holomorphic time assumption is made for the full F. N stays in the nonzero branch because full R(u,1)=1 and |RN|<1e30, so |R-1|<1e30*1e-245.",
        "checks": checks,
        "gates": {
            "all_original_source_gates": all(source["packet"]["gates"].values()),
            "actual_small_joint_neighborhood": LAPSE_RADIUS + source["radius"]
            < source["distance"]
            and TIME + source["radius"] < source["distance"],
            "whole_complex_X_image": source["X_image_distance"] < s.Rational(1, 400),
            "complete_heavy_suppression": source["entire_heavy"]
            < s.Rational(1, 10**2700),
            "all_actual_constant_switch_terms_retained": source["constant"]
            < s.Rational(1, 10**2200),
            "whole_analytic_scalar_modulus": source["scalar"] < 10**7,
            "whole_analytic_metric_modulus": source["metric"] < 2,
            "whole_normalized_source_modulus": source["normalized_source"]
            < s.Rational(1, 10**100),
            **{
                "whole_mixed_RFj_jet_" + str(key): bool(
                    row[0] < JET and row[1] < JET and row[2] < 1
                )
                for key, row in source["rows"].items()
            },
            "entire_complex_R_stays_near_one": JET * LAPSE_RADIUS < s.Rational(1, 100),
            "real_time_inside_original_C5_source_slab": TIME < current.TIME_LENGTH,
            "actual_mass_adapted_ratio": s.Rational(1, 10**4)
            < parent.source.MASS2 / 10**200
            < s.Rational(1, 100),
            "all_spatial_matter_vector_channels_retained": all(
                HAMILTONIAN.has(z) for z in COORDS
            ),
            "full_primitive_and_Hclock_retained": HAMILTONIAN.has(primitive, Hclock),
            "no_full_complex_time_or_original_quantum_matching_claim": True,
        },
    }
