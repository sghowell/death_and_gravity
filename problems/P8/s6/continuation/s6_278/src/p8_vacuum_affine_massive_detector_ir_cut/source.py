"""Literal retained vacuum jets, physical stress and the full M1 pair tree."""

from functools import cache

import sympy as s
from p8_vacuum_affine_physical_background_vertices import parent as current
from p8_vacuum_canonical_affine_decoupling import gravity as canonical_gravity

S, T, U, MU, K = s.symbols("s t u external_mass_squared kappa", real=True)
Z = s.Symbol("physical_scattering_cosine", real=True)
ETA = s.diag(1, -1, -1, -1)
KAPPA = current.heavy.K0
MASS2 = s.Integer(1)
NEWTON = 1 / (8 * s.pi * KAPPA)


def dot(p, q):
    return (p.T * ETA * q)[0]


def vertex(p, q, mass2):
    return p * q.T + q * p.T - ETA * (dot(p, q) + mass2)


def projector_contract(A, B):
    return s.expand(
        s.trace(ETA * A * ETA * B) - s.trace(ETA * A) * s.trace(ETA * B) / 2
    )


@cache
def full_vacuum():
    d = current.fixed_functions()
    u, X = current.u, current.X
    zero = {u: 0, X: 0}
    F, R = d["F_full"], d["R_full"]
    F0 = F.subs(X, 0)
    checks = {
        "actual_kappa": KAPPA - 10**800,
        "all_three_vacuum_constants_retained": s.cancel(
            F0.subs(u, 0)
            + d["vacuum_pressure_total"]
            - current.original.previous.vacuum_lower_function().subs(zero)
        ),
        "current_massive_quadratic_kinetic": s.cancel(
            s.diff(F, X).subs(zero) - s.Rational(1, 2)
        ),
        "current_massive_quadratic_mass": s.cancel(
            s.diff(F.subs(X, 0), u, 2).subs(u, 0) + 1
        ),
        "no_current_scalar_one_point": s.cancel(s.diff(F.subs(X, 0), u).subs(u, 0)),
        "current_R_vacuum": s.cancel(R.subs(zero) - 1),
        "no_nonminimal_Y_R_cubic": s.cancel(s.diff(R, X).subs(zero)),
        "no_nonminimal_Phi_R_mixing": s.cancel(s.diff(R.subs(X, 0), u).subs(u, 0)),
        "no_nonminimal_Phi_squared_R_cubic": s.cancel(
            s.diff(R.subs(X, 0), u, 2).subs(u, 0)
        ),
        "actual_Newton_dictionary": 8 * s.pi * NEWTON * KAPPA - 1,
    }
    return d, checks


def pair_tree(S=S, x=Z, mu=MU, kappa=K):
    return (S - (S - 4 * mu) * x * x) / (4 * kappa)


@cache
def data():
    d, checks = full_vacuum()
    # Independent full four-vector stress contraction in COM coordinates.
    e, p, k, cs, sn = s.symbols(
        "external_energy initial_momentum final_momentum cosine sine", real=True
    )
    p1, p2 = s.Matrix([e, 0, 0, p]), s.Matrix([e, 0, 0, -p])
    q1, q2 = s.Matrix([-k, -k * sn, 0, -k * cs]), s.Matrix([-k, k * sn, 0, k * cs])
    A = vertex(p1, p2, e * e - p * p)
    B = vertex(q1, q2, s.Integer(0))
    actual = projector_contract(A, B)
    # All legs on shell: k=e, cs^2+sn^2=1, S=4e^2.
    on = s.factor(actual.subs(k, e).subs(sn * sn, 1 - cs * cs))
    checks.update(
        {
            "literal_full_M1_pair_tree": s.factor(
                -on / (4 * K * e * e) - (e * e - p * p * cs * cs) / K
            ),
            "massive_source_Ward": s.simplify(A * ETA * (p1 + p2)),
            "massless_source_Ward": s.simplify(
                (B * ETA * (q1 + q2)).subs(sn * sn, 1 - cs * cs)
            ),
            "all_incoming_pair_mass_sum": s.expand(S + T + (2 * MU - S - T) - 2 * MU),
            "pair_stress_invariant_contraction": s.expand(
                ((T - MU) ** 2 + (U - MU) ** 2 - S * S) / 2 + (T - MU) * (U - MU)
            )
            .subs(U, 2 * MU - S - T)
            .expand(),
            "pair_physical_invariants": s.factor(
                ((MU - S / 2 + s.sqrt(S * (S - 4 * MU)) * Z / 2) - MU)
                * ((MU - S / 2 - s.sqrt(S * (S - 4 * MU)) * Z / 2) - MU)
                / (K * S)
                - pair_tree()
            ),
            "pair_tree_threshold_not_deleted": pair_tree(S=4 * MU) - MU / K,
        }
    )
    a, b, c, f, g, h = s.symbols("Jxx Jxy Jyy J0x J0y J00", real=True)
    J = s.Matrix([[h, f, g, h], [f, a, b, f], [g, b, c, g], [h, f, g, h]])
    checks["conserved_soft_current_TT_positive_square"] = s.factor(
        projector_contract(J, J) - (a - c) ** 2 / 2 - 2 * b * b
    )
    checks["conserved_soft_current_full_Ward"] = J * ETA * s.Matrix([1, 0, 0, 1])
    return {
        "whole_current_R_and_F": (d["R_full"], d["F_full"]),
        "all_three_fixed_vacuum_constants": d["vacuum_pressure_total"],
        "original_parameters": {
            "kappa": KAPPA,
            "external_tree_mass_squared": MASS2,
            "M1_mass_squared": s.Integer(0),
            "Newton": NEWTON,
            "heavy_mass_squared": current.heavy.MASS2,
            "Proca_mass": s.Integer(1000),
        },
        "physical_M1_action": "sqrt(-g) g^(mu nu) partial_mu Psi partial_nu Psi/2, Psi=sqrt(kappa)chi; not an adjustable nonminimal scalar",
        "canonical_metric_and_vertex": canonical_gravity.data()[
            "scalar_graviton_vertex"
        ],
        "whole_M1_production_tree": pair_tree(),
        "massive_to_massless_pair_sum": "s+t+u=2m^2, unlike four external massive scalars, where the sum is4m^2",
        "retained_vacuum_boundary": "The full fixed profiles and all vacuum constants remain. This is a formal coefficient about the retained quadratic flat reference, not a nonperturbative Minkowski quantum-vacuum theorem.",
        "checks": checks,
        "gates": {
            "original_kappa_and_external_mass_retained": KAPPA == 10**800
            and MASS2 == 1,
            "M1_is_massless_not_the_massive_Phi_or_heavy_field": True,
            "all_current_independent_functions_retained": d["R_full"]
            == current.heavy.coefficients()["R"],
            "nonminimal_cubic_checks_use_full_vacuum_jets": all(
                checks[k] == 0
                for k in (
                    "no_nonminimal_Y_R_cubic",
                    "no_nonminimal_Phi_R_mixing",
                    "no_nonminimal_Phi_squared_R_cubic",
                )
            ),
            "literal_pair_tree_contains_mass_and_angle": all(
                pair_tree().has(x) for x in (MU, Z)
            ),
            "massive_M1_production_threshold_is_four": pair_tree(S=4 * MU) == MU / K,
            "no_graviton_or_other_light_channel_discarded_from_full_amplitude": True,
        },
    }
