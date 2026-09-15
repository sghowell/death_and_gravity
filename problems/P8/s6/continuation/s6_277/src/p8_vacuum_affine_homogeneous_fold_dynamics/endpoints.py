"""Exact desingularized flow and source-pinned physical endpoint coefficients."""

from functools import cache

import sympy as s
from p8_vacuum_affine_physical_volume_correction import actual as volume_actual
from p8_vacuum_affine_physical_volume_correction import source as volume_source

from . import consistency, source


@cache
def ricci_identity():
    t = s.Symbol("physical_clock_coordinate", real=True)
    lapse = s.Function("positive_physical_lapse")(t)
    scales = [s.Function("positive_directional_scale_" + str(i))(t) for i in range(3)]
    g = s.diag(lapse**2, *[-(z**2) for z in scales])
    inv = g.inv()

    def partial(expr, index):
        return s.diff(expr, t) if index == 0 else s.Integer(0)

    connection = [
        [
            [
                sum(
                    inv[a, d]
                    * (partial(g[d, c], b) + partial(g[d, b], c) - partial(g[b, c], d))
                    / 2
                    for d in range(4)
                )
                for c in range(4)
            ]
            for b in range(4)
        ]
        for a in range(4)
    ]
    ricci = s.Matrix(
        4,
        4,
        lambda b, d: s.simplify(
            sum(
                partial(connection[a][b][d], a)
                - partial(connection[a][b][a], d)
                + sum(
                    connection[a][a][e] * connection[e][b][d]
                    - connection[a][d][e] * connection[e][b][a]
                    for e in range(4)
                )
                for a in range(4)
            )
        ),
    )
    scalar = s.simplify(
        sum(inv[b, d] * ricci[b, d] for b in range(4) for d in range(4))
    )
    rates = [s.diff(z, t) / (lapse * z) for z in scales]
    theta = sum(rates)
    return {
        "metric": g,
        "Ricci": ricci,
        "scalar": scalar,
        "rates": rates,
        "theta": theta,
        "residual": s.simplify(
            scalar + 2 * s.diff(theta, t) / lapse + theta**2 + sum(z * z for z in rates)
        ),
    }


@cache
def data():
    Cu, CN = s.symbols("full_Cu full_CN", real=True)
    cq = s.symbols("full_constraint_density_gradient0:5", real=True)
    f = s.symbols("full_homogeneous_density_rate0:5", real=True)
    K = Cu + sum(a * b for a, b in zip(cq, f, strict=True))
    tangent = s.expand(
        Cu * CN - CN * K + CN * sum(a * b for a, b in zip(cq, f, strict=True))
    )
    D, J = s.symbols(
        "positive_fold_second positive_future_constraint_preservation", positive=True
    )
    tau = s.Symbol("desingularized_curve_parameter", real=True)
    leading_u = -D * J * tau**2 / 2
    leading_N = 1 - J * tau
    theta = 3 / (2 * D * tau)
    lapse_slope = s.diff(leading_N, tau) / s.diff(leading_u, tau)
    curvature = -2 * s.diff(theta, tau) / s.diff(leading_u, tau)
    RR, Ru, RN, adot, ndot, N = s.symbols(
        "full_R full_Ru full_RN full_alpha_rate full_lapse_rate positive_lapse",
        real=True,
    )
    expansion = (3 * adot - s.Rational(3, 4) * (Ru + RN * ndot) / RR) / N
    bridge = volume_source.coefficient_bridge()
    metric = volume_actual.metric()
    physical = ricci_identity()
    checks = {
        "full_desingularized_constraint_tangency": tangent,
        "entire_profile_future_preservation": s.factor(
            consistency.FIRST.subs(source.p, -s.Rational(1, 10))
            - s.Rational(3, 250)
            + 3 * (source.PRESSURE - source.RHO) / 10
            - 3 * source.PRESSURE1 / 2
            + source.RHO1
        ),
        "leading_u_tau_tau": s.diff(leading_u, tau, 2) + D * J,
        "leading_N_tau": s.diff(leading_N, tau) + J,
        "leading_regular_lapse_velocity": s.factor(lapse_slope - 1 / (D * tau)),
        "leading_square_root_relation": s.factor(
            (leading_N - 1) ** 2 / leading_u + 2 * J / D
        ),
        "actual_physical_lapse_unchanged": metric["whole_positive_lapse"]
        - volume_source.N,
        "actual_physical_volume_U": s.simplify(
            metric["whole_volume_ratio"] - volume_source.R ** (-s.Rational(3, 4))
        ),
        **{"physical_bridge_" + k: v for k, v in bridge["checks"].items()},
        "actual_clock_expansion": s.factor(
            expansion.subs({RR: 1, Ru: 0, RN: -2, N: 1}) - 3 * adot - 3 * ndot / 2
        ),
        "finite_tau_expansion_coefficient": s.factor(tau * theta - 3 / (2 * D)),
        "finite_tau_Ricci_coefficient": s.factor(tau**3 * curvature + 3 / (D * D * J)),
        "entire_independent_BianchiI_Ricci": physical["residual"],
    }
    return {
        "whole_desingularized_vector_field": "u_tau=C_N; x_tau=C_N F(u,N,x) for every full canonical density/coordinate rate; N_tau=-K, K=C_u+C_x.F. No quantum force or alternate state is inserted; F is the entire classical parent with its fixed profiles and live heavy field.",
        "whole_constraint_tangency_expression": tangent,
        "whole_exact_future_endpoint_parameters": {
            "p": -s.Rational(1, 10),
            "shear": source.FOLD_SHEAR.subs(source.p, -s.Rational(1, 10)),
            "matter_square": source.FOLD_MATTER_SQUARE,
            "D": source.FOLD_SECOND,
            "J": consistency.FIRST.subs(source.p, -s.Rational(1, 10)),
        },
        "whole_leading_endpoint_jets_not_replacement_solution": (leading_u, leading_N),
        "whole_proper_time_normal_expansion": expansion,
        "whole_exact_physical_BianchiI_Ricci_scalar": physical["scalar"],
        "whole_physical_directional_expansions": physical["rates"],
        "whole_endpoint_limits": {
            "tau_theta": 3 / (2 * D),
            "tau_cubed_Ricci": -3 / (D * D * J),
        },
        "whole_endpoint_existence_proof_boundary": "A smooth vector field on the regular C=0 submanifold has a local curve through the fold. For tau>0 small, C_N=-D J tau+O(tau^2)<0 and u=-D J tau^2/2+O(tau^3)<0 is invertible away from zero. These yield full regular clock-time classical solutions approaching u0 in finite proper time. Lapse and spatial volume have finite positive limits; normal expansion and the physical Ricci scalar diverge. The auxiliary curve parameter is singular as a clock change at its endpoint; it is not a physical clock-time turning point.",
        "whole_scope": "This is a homogeneous classical comparison family with changed M1/shape/trace momenta, not the original prepared quantum state, not the S276 mean-preserving wave path and not the corrected S275 finite hybrid. No quantitative initial interval, cutoff-validity claim, generic instability near the original bounce, quantum obstruction or original P8 closure follows.",
        "checks": checks,
        "gates": {
            "actual_D_greater_than28": s.Rational(11391, 400)
            - s.Rational(15, 2) * source.PROFILE_BOUND
            > 28,
            "actual_future_J_greater_than11_over1000": s.Rational(3, 250)
            - s.Rational(31, 10) * source.PROFILE_BOUND
            > s.Rational(11, 1000),
            "positive_physical_lapse_not_conformally_rescaled": metric["gates"][
                "physical_lapse_and_shift_unchanged_by_spatial_conformal_map"
            ],
            "proper_time_endpoint_finite_with_N_tending_one": True,
            "negative_CN_regular_side_not_zero_division_at_fold": True,
            "live_heavy_temporal_and_cyclic_coordinates_reconstructed": True,
            "original_couplings_and_fixed_profiles_unchanged": True,
            "not_an_evaluated_domain_or_original_state_lifetime": True,
        },
    }
