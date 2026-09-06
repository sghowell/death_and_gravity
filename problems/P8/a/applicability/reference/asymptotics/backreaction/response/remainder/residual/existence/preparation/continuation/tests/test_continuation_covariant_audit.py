"""Independent trace, sheet, causal-moment and exact-scale audits.

The frozen block is not asserted to linearize the actual rolling state.
No comparator pole is promoted to an actual nonlinear instability.
"""

from fractions import Fraction as Q

import sympy as sp


def test_actual_scale_trace_and_anomaly_without_a_flat_background_assumption():
    from p8a_continuation.decomposition import full_trace_in_s

    x = sp.Symbol("x", real=True)
    a, wick = sp.Function("a", positive=True)(x), sp.Function("S")(x)
    delta = sp.Symbol("delta", positive=True)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    q = wick/a**2
    # Direct A.11 trace in its actual auxiliary q convention.
    original = sp.diff(q, x, 2)+2*h*sp.diff(q, x)-u/(60*delta)
    original += (u**2/4+h**2*(u+h**2)/30)/a**2
    expanded = (sp.diff(wick, x, 2)-2*h*sp.diff(wick, x)
                +2*(u+h**2)*wick-a**2*u/(60*delta)
                +u**2/4+h**2*(u+h**2)/30)
    assert sp.simplify(a**2*original-expanded) == 0
    assert sp.simplify(expanded-full_trace_in_s(a, wick, delta, x)) == 0


def test_double_primitive_and_Einstein_sign_in_the_integrated_block():
    from p8a_continuation.resolvent import transfer

    s, scale, delta = sp.symbols("s af delta", positive=True)
    gamma, finite = sp.symbols("gamma df", real=True)
    # Integrating the trace once acts on u, then u=I X supplies another
    # primitive. Dropping either primitive or changing the sign is detected.
    log_block = (sp.log(s)+gamma)/2+finite
    exact = log_block-scale**2/(60*delta*s**2)
    assert sp.factor(exact*transfer(s, gamma+2*finite, scale**2/(30*delta))-1) == 0
    wrong_sign = log_block+scale**2/(60*delta*s**2)
    wrong_order = log_block-scale**2/(60*delta*s)
    assert sp.factor(wrong_sign-exact) != 0
    assert sp.factor(wrong_order-exact) != 0


def test_principal_sheet_angle_count_and_exclusion_of_imaginary_roots():
    angle = sp.Symbol("theta", positive=True)
    # Positive c requires -theta/sin(2theta)>0, excluding the first
    # quadrant. On the second quadrant this angular map is one-to-one.
    angular_log = sp.log(-angle/sp.sin(2*angle))-2*angle/sp.tan(2*angle)
    numerator = (sp.sin(2*angle)-2*angle*sp.cos(2*angle))**2
    numerator += 4*angle**2*sp.sin(2*angle)**2
    assert sp.trigsimp(sp.diff(angular_log, angle)
                       -numerator/(angle*sp.sin(2*angle)**2)) == 0
    assert (-angle/sp.sin(2*angle)).subs(angle, sp.pi/4) < 0
    assert (-angle/sp.sin(2*angle)).subs(angle, 3*sp.pi/4) > 0
    radius, beta = sp.symbols("radius beta", positive=True)
    imaginary_axis_imaginary_part = sp.im((sp.I*radius)**2
                                         *(sp.log(radius)+sp.I*sp.pi/2+beta))
    assert imaginary_axis_imaginary_part == -sp.pi*radius**2/2
    # The written endpoint/monotonicity argument proves the complete count,
    # rather than these two angular fixtures being treated as enumeration.


def test_literal_residue_cut_and_positive_feedback_resummation():
    from p8a_continuation.resolvent import cut_density, residue, transfer

    s, pole, c, r = sp.symbols("s p c r", positive=True)
    beta = sp.Symbol("beta", real=True)
    denominator = s**2*(sp.log(s)+beta)-c
    literal_residue = 2*pole**2/sp.diff(denominator, s).subs(s, pole)
    assert sp.simplify(literal_residue-residue(pole, beta)) == 0
    upper = 2/(sp.log(r)+beta-c/r**2+sp.I*sp.pi)
    lower = sp.conjugate(upper)
    assert sp.simplify((lower-upper)/(2*sp.pi*sp.I)-cut_density(r, beta, c)) == 0
    j = 2/(sp.log(s)+beta)
    assert sp.factor(j/(1-c*j/(2*s**2))-transfer(s, beta, c)) == 0


def test_elementary_root_brackets_and_norm_growth_with_fraction_only():
    # Independent interval arithmetic for beta, log(10), log(2), and c.
    delta, af = Q(1, 10**14), Q(5, 2)
    c = af**2/(30*delta)
    beta_lo, beta_hi = -Q(53, 60), Q(11, 30)
    left_upper = 10**12*(6*Q(5, 2)+beta_hi)-c
    right_lower = 4*10**12*(Q(1, 2)+12+beta_lo)-c
    assert left_upper < 0 < right_lower
    w_upper = 2*(1+6*Q(5, 2)+beta_hi)
    assert w_upper < 33
    assert (Q(8, 3)**10-1)/17-Q(8, 3) > 1000
    assert (Q(8, 3)**100-1)/17-Q(8, 3) > 10**28
    assert 2/(right_lower/(4*10**12)) == Q(240, 769)


def test_growing_pole_moment_is_a_forcing_condition_not_homogeneous_initial_data():
    t, p = sp.symbols("t p", positive=True)
    chi = t**2*(1-t)**2
    balanced = sp.diff(chi, t)-p*chi
    assert sp.simplify(sp.exp(-p*t)*balanced-sp.diff(sp.exp(-p*t)*chi, t)) == 0
    assert chi.subs(t, 0) == chi.subs(t, 1) == 0
    assert sp.diff(chi, t).subs(t, 0) == sp.diff(chi, t).subs(t, 1) == 0
    assert sp.simplify(sp.integrate(sp.exp(-p*t)*balanced, (t, 0, 1))) == 0
    # Conversely, nonnegative nonzero compact forcing has positive moment.
    assert chi.subs(t, sp.Rational(1, 2)) > 0
    # This exact finite-interval control is not the unknown full G_f[X].


def test_deleting_positive_pole_does_not_preserve_causal_inverse():
    from p8a_continuation.resolvent import transfer

    s = sp.Symbol("s", positive=True)
    pole, beta = sp.S(2), sp.S.Zero
    c = pole**2*sp.log(pole)
    residue = 2*pole/(1+2*sp.log(pole))
    block = (sp.log(s)+beta-c/s**2)/2
    deleted = transfer(s, beta, c)-residue/(s-pole)
    defect = sp.factor(deleted*block-1)
    assert sp.simplify(defect+residue*block/(s-pole)) == 0
    assert defect.subs(s, 4) != 0


def test_named_proper_time_growth_scale_does_not_depend_on_freezing_scale():
    af, delta, t0, w = sp.symbols("af delta T0 W", positive=True)
    beta = sp.EulerGamma-sp.Rational(19, 30)-sp.log(af/2)
    c = af**2/(30*delta)
    z = 2*c*sp.exp(2*beta)
    assert sp.diff(sp.simplify(z), af) == 0
    growth2 = (2*c/w)/(af*t0)**2
    restored = 192*sp.pi**2/(2880*sp.pi**2*delta*t0**2*w)
    assert sp.factor(growth2-restored) == 0
