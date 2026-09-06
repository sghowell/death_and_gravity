"""Exact C1 bounds for the actual past-prepared nonlinear Wick response.

This module bounds R[U] = integral k (F[U] - 1 - F1[U]) dk, with
F = 2 k |v_k|^2 and a common original plane-wave start.  It does not bound
the singular linear logarithmic operator, prove a SEE fixed point, or
assert that every C1 potential defines a smooth Hadamard spacetime.

All numerical interfaces use exact rational inputs.  The exponential and
logarithm evaluations below are upper bounds, not floating-point values.
The functions in identities() have no numerical input and check formal
identities used in the written proof.
"""

from __future__ import annotations

import sympy as sp

Z_MAX = sp.Rational(1, 4)
PREPARED_U1 = sp.Rational(61013499, 8192)
DELTA_CALIBRATION = sp.Rational(1, 10**14)


def rational(value: object) -> sp.Rational:
    """Read an exact rational, rejecting bools, binary floats and infinities."""
    if isinstance(value, (bool, float, sp.Float)):
        raise TypeError("an exact finite rational is required")
    if not isinstance(value, (int, str, sp.Rational)):
        raise TypeError("an exact finite rational is required")
    try:
        result = sp.Rational(value)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("an exact finite rational is required") from exc
    if result.is_finite is not True:
        raise ValueError("an exact finite rational is required")
    return result


def _nonnegative(value: object, name: str) -> sp.Rational:
    result = rational(value)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _positive(value: object, name: str) -> sp.Rational:
    result = rational(value)
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _integer(value: object, name: str, minimum: int) -> int:
    result = rational(value)
    if result.q != 1 or result < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return int(result)


def exponential_upper(argument: object) -> sp.Rational:
    """Return 1/(1-x) >= exp(x), for an exact 0 <= x < 1."""
    argument = _nonnegative(argument, "argument")
    if argument >= 1:
        raise ValueError("the geometric exponential majorant requires x < 1")
    return 1 / (1 - argument)


def log_ratio_upper(ratio: object, terms: object = 8) -> sp.Rational:
    """Bound log(r) for r >= 1 by an atanh series plus geometric tail."""
    ratio = _positive(ratio, "ratio")
    if ratio < 1:
        raise ValueError("ratio must be at least one")
    terms = _integer(terms, "terms", 1)
    argument = (ratio - 1) / (ratio + 1)
    partial = sum(
        (2 * argument ** (2 * j + 1) / (2 * j + 1) for j in range(terms)),
        sp.S.Zero,
    )
    tail = 2 * argument ** (2 * terms + 1) / (
        (2 * terms + 1) * (1 - argument**2)
    )
    return partial + tail


def dyson_coefficients(order: object, duration: object, split: object) -> dict:
    """Return c_n with integral k |F_n| dk <= c_n B^n, for n >= 3.

    The low-frequency estimate uses |sin(k r)/k| <= duration.  The
    high-frequency estimate uses |sin(k r)/k| <= 1/k.  The quadratic
    term must instead be treated by its exact logarithmic kernel.
    """
    order = _integer(order, "order", 3)
    duration = _positive(duration, "duration")
    split = _positive(split, "split")
    factorial = sp.factorial(order)
    infrared = split**2 * (2 * duration**2) ** order / (2 * factorial)
    ultraviolet = (
        split**2
        * (2 * duration / split) ** order
        / (factorial * (order - 2))
    )
    return {
        "order": order,
        "duration": duration,
        "split": split,
        "infrared": infrared,
        "ultraviolet": ultraviolet,
        "total": infrared + ultraviolet,
    }


def derivative_order_bound(
    order: object, derivative_cap: object, duration: object, split: object
) -> dict:
    """Bound one order of R'[U]-R'[V], per ||U'-V'||_infinity."""
    data = dyson_coefficients(order, duration, split)
    derivative_cap = _nonnegative(derivative_cap, "derivative_cap")
    amplitude_cap = derivative_cap * data["duration"]
    multiplier = data["order"] ** 2 * amplitude_cap ** (data["order"] - 1)
    return {
        **data,
        "derivative_cap": derivative_cap,
        "amplitude_cap": amplitude_cap,
        "infrared_lipschitz": multiplier * data["infrared"],
        "ultraviolet_lipschitz": multiplier * data["ultraviolet"],
        "lipschitz": multiplier * data["total"],
    }


def shared_order_bound(
    order: object,
    derivative_cap: object,
    duration: object,
    future_length: object,
    split: object,
) -> dict:
    """Use the marked-factor L1 simplex estimate for a common prehistory."""
    data = dyson_coefficients(order, duration, split)
    derivative_cap = _nonnegative(derivative_cap, "derivative_cap")
    future_length = _positive(future_length, "future_length")
    if future_length > data["duration"]:
        raise ValueError("future_length cannot exceed the full history")
    amplitude_cap = derivative_cap * data["duration"]
    order = data["order"]
    multiplier = (
        order
        * amplitude_cap ** (order - 1)
        * future_length
        / data["duration"]
        * (1 + (order - 1) * future_length / (2 * data["duration"]))
    )
    return {
        **data,
        "derivative_cap": derivative_cap,
        "amplitude_cap": amplitude_cap,
        "future_length": future_length,
        "lipschitz": multiplier * data["total"],
    }


def _domain(derivative_cap: object, duration: object) -> dict:
    derivative_cap = _nonnegative(derivative_cap, "derivative_cap")
    duration = _positive(duration, "duration")
    amplitude_cap = derivative_cap * duration
    strength = amplitude_cap * duration**2
    if strength > Z_MAX:
        raise ValueError("the rational calibration requires M T^3 <= 1/4")
    return {
        "derivative_cap": derivative_cap,
        "duration": duration,
        "amplitude_cap": amplitude_cap,
        "strength": strength,
        "split": 1 / duration,
        "exponential_upper": exponential_upper(2 * strength),
    }


def full_history_bounds(derivative_cap: object, duration: object) -> dict:
    """Give C0/C1 Lipschitz constants per D=||U'-V'||, with U(0)=V(0)=0.

    Both potentials obey ||U'||, ||V'|| <= M on their whole history,
    whose actual duration is at most T. No continuation beyond the
    observed endpoint is required when T is a strict upper bound.
    The reported R-value constant is T times the derivative constant.
    These are constants for the nonlinear mode functional R alone.
    """
    data = _domain(derivative_cap, duration)
    strength = data["strength"]
    quadratic = 2 * strength
    higher = 18 * strength**2 * data["exponential_upper"]
    derivative = quadratic + higher
    return {
        **data,
        "quadratic_derivative": quadratic,
        "higher_derivative": higher,
        "derivative": derivative,
        "value": data["duration"] * derivative,
    }


def shared_history_bounds(
    derivative_cap: object,
    duration: object,
    future_length: object,
    log_terms: object = 8,
) -> dict:
    """Bounds for a common past and a future interval of length at most L.

    T bounds the complete history duration, not an unproved extension of
    the physical potential beyond its observed endpoint.
    The logarithmic majorant is rational.  The sharper written quadratic
    coefficient is B L^2 [5/4 + log(T/L)/2].  Only the future interval
    contributes to the value difference, hence its multiplier is L.
    """
    data = _domain(derivative_cap, duration)
    future_length = _positive(future_length, "future_length")
    if future_length > data["duration"]:
        raise ValueError("future_length cannot exceed the full history")
    log_upper = log_ratio_upper(data["duration"] / future_length, log_terms)
    quadratic = (
        data["amplitude_cap"]
        * future_length**2
        * (sp.Rational(5, 4) + log_upper / 2)
    )
    higher = (
        18
        * data["strength"] ** 2
        * future_length
        / data["duration"]
        * data["exponential_upper"]
    )
    derivative = quadratic + higher
    return {
        **data,
        "future_length": future_length,
        "log_ratio_upper": log_upper,
        "quadratic_derivative": quadratic,
        "higher_derivative": higher,
        "derivative": derivative,
        "value": future_length * derivative,
    }


def identities() -> dict:
    """Return exact zero identities; these supplement, not replace, the proof."""
    time, inner, outer, epsilon, frequency = sp.symbols(
        "t r s epsilon k", positive=True
    )
    duration, future, strength = sp.symbols("T L z", positive=True)
    coupling = sp.symbols("u", real=True)
    left, right = outer - inner, time - inner
    abel = sp.log((epsilon**2 + 4 * right**2) / (epsilon**2 + 4 * left**2)) / 2
    oscillator = 1 - coupling / (frequency**2 + coupling) * sp.sin(
        sp.sqrt(frequency**2 + coupling) * time
    ) ** 2
    oscillator_quadratic = sp.sin(frequency * time) ** 2 / frequency**4 - (
        time * sp.sin(2 * frequency * time) / (2 * frequency**3)
    )
    order = sp.symbols("n", integer=True, positive=True)
    polynomial = sp.Function("P")(time)
    phase = sp.exp(-sp.I * frequency * time)
    return {
        "abel_time_derivative": sp.factor(
            sp.diff(abel, time) - right / (right**2 + epsilon**2 / 4)
        ),
        "abel_outer_endpoint": sp.simplify(abel.subs(outer, time)),
        "quadratic_kernel_derivative": sp.simplify(
            sp.diff(sp.log(right / left), time) - 1 / right
        ),
        "quadratic_kernel_outer_endpoint": sp.simplify(
            sp.log(right / left).subs(outer, time)
        ),
        "oscillator_second_order": sp.trigsimp(
            sp.diff(oscillator, coupling, 2).subs(coupling, 0) / 2
            - oscillator_quadratic
        ),
        "oscillator_radial_ibp": sp.trigsimp(
            frequency * oscillator_quadratic
            + sp.diff(sp.sin(frequency * time) ** 2 / frequency**2, frequency) / 2
        ),
        "phase_derivative_transfer": sp.simplify(
            sp.diff(phase * polynomial, time)
            + sp.I * frequency * phase * polynomial
            - phase * sp.diff(polynomial, time)
        ),
        "factorial_majorant_gap": sp.factor(
            sp.Rational(3, 2) - order / ((order - 1) * (order - 2))
            - (order - 3) * (3 * order - 2) / (2 * (order - 1) * (order - 2))
        ),
        "localized_quadratic_integral": sp.simplify(
            sp.integrate(
                (future - inner) * sp.log(duration / inner), (inner, 0, future)
            )
            - future**2 * (sp.log(duration / future) / 2 + sp.Rational(3, 4))
        ),
        "quadratic_total_kernel_mass": sp.simplify(
            sp.integrate(time - inner, (inner, 0, time)) - time**2 / 2
        ),
        "higher_tail_coefficient": sp.simplify(
            sp.Rational(3, 2) / strength * sp.Rational(3, 2) * (2 * strength) ** 3
            - 18 * strength**2
        ),
    }


def calibration() -> dict:
    """A genuine derivative-ball envelope containing A7's prepared potential.

    Time and potential are dimensionless in eta_star units. T=3 bounds
    the active history beginning at the free slice eta_star. The same
    state has already propagated freely from its original y=1/2 slice;
    moving the free phase origin does not reset its covariance. No
    potential beyond the observed interval is implicitly assumed.
    The factor two leaves a margin around the pinned A7 derivative cap.
    This does not assert that a full SEE contraction preserves that ball.
    """
    derivative_cap = 2 * DELTA_CALIBRATION * PREPARED_U1
    full = full_history_bounds(derivative_cap, 3)
    shared = shared_history_bounds(derivative_cap, 3, sp.Rational(1, 4))
    return {
        "delta": DELTA_CALIBRATION,
        "prepared_U1": PREPARED_U1,
        "margin_factor": sp.Integer(2),
        "full": full,
        "shared": shared,
        "full_derivative_rounded_upper": sp.Rational(1, 10**8),
        "shared_derivative_rounded_upper": sp.Rational(1, 10**10),
        "full_rounding_gap": sp.Rational(1, 10**8) - full["derivative"],
        "shared_rounding_gap": sp.Rational(1, 10**10) - shared["derivative"],
    }
