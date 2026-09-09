"""Independent negative controls and the second-order conservation boundary."""

from functools import cache

import sympy as sp
from p8_coupled_momentum import constraints, quadratic

from . import bridge, center, spatial


@cache
def data():
    c = center.data()
    N = sp.Symbol("positive_physical_lapse", positive=True)
    dc = sp.Symbol("matter_density_deviation", real=True)
    gradient = sp.Symbol("nonnegative_spatial_gradient_squared", nonnegative=True)
    full_center = (sp.Rational(1, 10) + dc) ** 2 / (2 * N**3) + gradient / (2 * N)
    amplitude = sp.Symbol("small_tensor_momentum_amplitude", real=True)
    pure = {value: 0 for value in c["fields"]}
    pure[c["fields"][5]] = amplitude
    lapse = sp.factor(c["n2"].subs(pure, simultaneous=True))
    density = sp.factor(c["rho"].subs(pure, simultaneous=True))
    # Formal order accounting, not an assertion that the second variation
    # alone is conserved on the undeformed background.
    e = sp.Symbol("formal_perturbation_parameter")
    D = sp.symbols("connection0 connection1 connection2", commutative=False)
    T = sp.symbols("stress0 stress1 stress2", commutative=False)
    ward = sp.expand(
        sum(e**i * D[i] for i in range(3)) * sum(e**i * T[i] for i in range(3))
    ).coeff(e, 2)
    return {
        "actual_full_classical_center_density": full_center,
        "pure_tensor_second_order_lapse": lapse,
        "pure_tensor_quadratic_matter_density": density,
        "pure_tensor_implicit_branch_lapse_Hessian": -sp.Rational(243, 80),
        "second_order_divergence_terms": ward,
        "checks": {
            "pure_tensor_second_order_lapse_is_nonzero": lapse
            - sp.Rational(40, 81) * amplitude**2,
            "pure_tensor_quadratic_density_is_negative_correction": density
            + amplitude**2 / 135,
            "full_classical_center_density_has_two_nonnegative_terms": full_center
            - (sp.Rational(1, 10) + dc) ** 2 / (2 * N**3)
            - gradient / (2 * N),
            "second_order_Ward_accounting_retains_connection_variations": sp.expand(
                ward - D[0] * T[2] - D[1] * T[1] - D[2] * T[0]
            ),
        },
    }


def malformed_field(kind):
    ctx = constraints.context(((1, 0, 0), (-1, 0, 0)), symbolic_scale=True)
    f, t, tp, W, Pi = quadratic.fields(ctx, "curvature", "matter")
    if kind == "nonlinear":
        f["zeta"] += ctx.leg(0) * ctx.leg(1)
    elif kind == "background":
        f["chi"] += 1
    elif kind == "wrong_context":
        other = constraints.context(((1, 0, 0), (-1, 0, 0)), symbolic_scale=True)
        f["chi"] = other.leg(0)
    elif kind == "nontracefree":
        t[0][0] += ctx.leg(0)
    elif kind == "nonsymmetric":
        t[0][1] += ctx.leg(0)
    elif kind == "nontransverse":
        t[0][1] += ctx.leg(0)
        t[1][0] += ctx.leg(0)
    elif kind == "missing_scalar":
        del f["chi_p"]
    elif kind == "short_vector":
        W.pop()
    elif kind == "short_tensor":
        t.pop()
    elif kind == "nonjet":
        f["chi"] = sp.Integer(0)
    elif kind == "noncontext":
        ctx = None
    else:
        raise ValueError("Unknown malformed-input control")
    return spatial.construct(ctx, f, t, tp, W, Pi)


def bad_cases():
    result = []
    for value in (
        True,
        0.1,
        "0",
        sp.oo,
        sp.nan,
        sp.Symbol("u"),
        sp.Rational(-11, 1000),
        sp.Rational(11, 1000),
    ):
        result.append(("time_" + str(value), spatial.at_time, (value,)))
    for i in (0, 1):
        for value in (None, True, 0, "undeclared", "tensor3"):
            args = ["curvature", "matter"]
            args[i] = value
            result.append(
                ("channel_" + str(i) + "_" + str(value), spatial.pair, tuple(args))
            )
    for kind in (
        "nonlinear",
        "background",
        "wrong_context",
        "nontracefree",
        "nonsymmetric",
        "nontransverse",
        "missing_scalar",
        "short_vector",
        "short_tensor",
        "nonjet",
        "noncontext",
    ):
        result.append(("field_" + kind, malformed_field, (kind,)))
    return result


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("An invalid physical observable input was accepted: " + name)
    return len(bad_cases())


@cache
def summary():
    c = center.data()
    b = bridge.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "missing_linear_lapse_shortcut_matrix_nonzero": c["matrices"][
            "missing_if_lapse_truncated"
        ]
        != sp.zeros(14),
        "omitting_either_canonical_boundary_shift_changes_actual_map": all(
            b[key] != sp.zeros(4)
            for key in (
                "missing_matter_boundary_shift_matrix",
                "missing_metric_boundary_shift_matrix",
            )
        ),
        "negative_scalar_principal_configuration_test": c[
            "negative_principal_configuration_test"
        ],
        "negative_tensor_momentum_density_Hessian": c["matrices"]["rho"][5, 5],
        "negative_transverse_Proca_momentum_density_Hessian": c["matrices"]["rho"][
            12, 12
        ],
        "full_classical_positive_density_not_refuted": True,
        "no_arbitrary_amplitude_or_unbounded_energy_claim": True,
        "no_automatic_positive_square_QEI_or_full_Ward_transfer": True,
    }
