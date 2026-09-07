"""Exact necessary cone-error budgets, never an omitted-operator estimate."""

import sympy as sp


def _exact_finite(value):
    if isinstance(value, bool):
        raise TypeError("Exact finite real parameters required, not bool")
    value = sp.sympify(value)
    if value.has(sp.Float) or value.is_real is not True or value.is_finite is not True:
        raise TypeError("Exact finite real parameters with declared domains required")
    return value


def speed_budget(a, q, actual_null_stress, tensor_kinetic, kinetic_error=0, gradient_error=0):
    """A pointwise remainder lower bound; inputs/errors are physical-frame.

    A,q,GTh>0, n_h>=0, E_G,E_F>=0. If Gnew remains positive and
    |delta G|<=E_G, |delta F|<=E_F, then Fnew-Gnew >= the returned margin.
    Positivity of this margin excludes a subluminal corrected cone. The
    caller must independently prove those operator/solution error bounds.
    """
    a, q, nh, gt, eg, ef = map(_exact_finite, (a, q, actual_null_stress, tensor_kinetic, kinetic_error, gradient_error))
    if any(value.is_positive is not True for value in (a, q, gt)):
        raise ValueError("A, q and G_T_h must be strictly positive")
    if any(value.is_nonnegative is not True for value in (nh, eg, ef)):
        raise ValueError("Actual null stress and error budgets must be nonnegative")
    delta = a*nh/(4*q)
    squared_gap = sp.expand(2*delta+delta**2)
    margin = sp.expand(gt*squared_gap-eg-ef)
    return {"speed_gap": delta, "squared_speed_gap": squared_gap,
            "uncorrected_F_minus_G": gt*squared_gap,
            "corrected_F_minus_G_lower_bound": margin,
            "corrected_kinetic_lower_bound": gt-eg}


def checks():
    a, q, nh, gt, eg, ef = sp.symbols("A q n_h G_T_h E_G E_F", positive=True)
    d = speed_budget(a, q, nh, gt, eg, ef)
    m, target, tau = sp.symbols("m_0 M_target tau", positive=True)
    frozen = d["speed_gap"].subs({q: target**2*m**2/4, nh: target**2/(100*tau**2)})
    return {"exact_squared_speed_excess": sp.expand(d["squared_speed_gap"]-(1+d["speed_gap"])**2+1),
            "physical_coefficient_remainder_budget": sp.expand(d["corrected_F_minus_G_lower_bound"]-gt*d["squared_speed_gap"]+eg+ef),
            "conditional_free_M1_vacuum_scale_dictionary": sp.cancel(frozen-a/(100*(m*tau)**2))}


def controls():
    exact = speed_budget(2, 1, 12, sp.Rational(1, 14))
    protected = speed_budget(2, 1, 12, sp.Rational(1, 14), sp.Rational(1, 56), sp.Rational(1, 56))
    return {"actual_example_speed_gap": exact["speed_gap"],
            "actual_example_squared_gap": exact["squared_speed_gap"],
            "actual_example_F_minus_G": exact["uncorrected_F_minus_G"],
            "protected_kinetic_margin": protected["corrected_kinetic_lower_bound"],
            "protected_cone_margin": protected["corrected_F_minus_G_lower_bound"],
            "nonrolling_has_no_strict_gap": speed_budget(1, 1, 0, 1)["speed_gap"]}
