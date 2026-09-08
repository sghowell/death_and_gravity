"""Independent full-connection, coupled-constraint and unit interfaces."""
from functools import cache

import sympy as sp
from p8_affine import connection
from p8_affine_kinetic import scalar as old

from . import dynamics, geometry, rank_two


def _matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.factor))


@cache
def two_trace_action():
    """The regular eight-trace Schur reduction at quadratic rolling order."""
    data = geometry.schur()
    m, n = data["M_clock"], data["N"]
    w, d = data["inverse_N_transpose_clock"], data["D_clock"]
    lift = _matrix(w*d.inv())
    projector = _matrix(sp.eye(60)-lift*n)
    return {"M": m, "N": n, "W": w, "D": d, "lift": lift,
            "projector": projector, "complement_dimension": 52,
            "quadratic_rolling_only": True}


@cache
def scalar_embedding():
    """Both vector temporals and the original lapse/shift before elimination."""
    z, D = rank_two.Z, rank_two.D
    q, n = old.q, old.n
    temporal = sp.Matrix(sp.symbols("trace_t1 trace_t2", real=True))
    sigma = sp.Matrix(sp.symbols("trace_s1 trace_s2", real=True))
    velocity = sp.Matrix(sp.symbols("trace_sd1 trace_sd2", real=True))
    time_shift = sp.Matrix(sp.symbols("trace_d1 trace_d2", real=True))
    space_shift = sp.Matrix(sp.symbols("trace_e1 trace_e2", real=True))
    mass_time, mass_space = temporal+time_shift*n, sigma+space_shift*n
    extra = ((mass_time.T*D.inv()*mass_time)[0]
             -q*(mass_space.T*D.inv()*mass_space)[0]
             +q*((velocity-temporal).T*z*(velocity-temporal))[0])/2
    before = old.action()["base"]+extra
    temporal_hessian = D.inv()+q*z
    solved = _matrix(temporal_hessian.inv()*(q*z*velocity-D.inv()*time_shift*n))
    lapse = (old.vd+old.ell*old.matter/2)/old.theta
    # Differentiate the exact reduced quadratic form, not a frequency ansatz.
    after = before.subs(dict(zip(temporal, solved, strict=True)), simultaneous=True).subs(n, lapse)
    pure_vector = _matrix(sp.hessian(after, tuple(velocity))/2)
    expected = rank_two.longitudinal()["kinetic"].subs(rank_two.Q, q)
    auxiliary = (n, old.shift, *temporal)
    hessian = sp.hessian(before, auxiliary)
    return {"before": before, "temporal": temporal, "velocity": velocity,
            "temporal_solution": solved, "lapse": lapse,
            "pure_vector_kinetic": pure_vector,
            "unchanged_shift": sp.factor(sp.diff(before, old.shift)/(2*q)
                                            -old.theta*n+old.vd+old.ell*old.matter/2),
            "pure_vector_residual": _matrix(pure_vector-expected),
            "full_auxiliary_rank_residual": sp.factor(
                hessian.det()+4*q**2*old.theta**2*temporal_hessian.det())}


@cache
def checks():
    g, independent = geometry.traces(), dynamics.geometry()
    data, roll, own_roll = two_trace_action(), geometry.rolling(), dynamics.rolling()
    m, n, w, d = data["M"], data["N"], data["W"], data["D"]
    lift, projector = data["lift"], data["projector"]
    full_m = connection.quadratic()["hessian"].subs(connection.P, sp.Rational(1, 2))
    embedding = geometry.schur()["embedding"]
    results = {
        "independent_full_trace_maps": g["full"]-independent["map"],
        "independent_all_block_response": _matrix(d-independent["D_two"]),
        "literal_rank_two_field_space": _matrix(d-sp.kronecker_product(rank_two.D, geometry.ETA)),
        "eight_trace_quotient_Euler": _matrix(m*w-n.T),
        "eight_trace_full_64_Euler": _matrix(full_m*embedding*w-g["full"].T),
        "eight_trace_right_inverse": _matrix(n*lift-sp.eye(8)),
        "eight_trace_mass_action": _matrix(lift.T*m*lift-d.inv()),
        "complement_annihilation": _matrix(n*projector),
        "complement_idempotence": _matrix(projector*projector-projector),
        "complement_cross_Euler": _matrix(projector.T*m*lift),
        "complement_dimension": sp.trace(projector)-52,
        "rank_one_Schur_interface": sp.factor(geometry.rank_one()["gamma"]-dynamics.coefficients()["gamma"]),
    }
    mapping = {own_roll["h"]: geometry.H, own_roll["H"]: roll["symbols"]["Hubble"]}
    for name in ("alpha", "d", "e"):
        results["independent_rolling_"+name] = sp.factor(roll[name]-own_roll[name].subs(mapping))
    results["independent_rolling_curl"] = sp.factor(roll["curl_lapse_value"]-own_roll["curl_n"].subs(mapping))
    results["independent_full_Vstar"] = _matrix(geometry.schur()["stationary"][:4, :]-independent["Vstar"])
    results["independent_full_Ustar"] = _matrix(geometry.schur()["stationary"][4:, :]-independent["Ustar"])
    scalar = scalar_embedding()
    for name in ("unchanged_shift", "pure_vector_residual", "full_auxiliary_rank_residual"):
        results["physical_rank_two_"+name] = scalar[name]
    # Exact agreement with the already frozen one-trace physical constraint action.
    one = dynamics.regular()
    at_v = {dynamics.gamma: sp.Rational(3, 8), dynamics.zeta: old.zeta,
            dynamics.d: old.dd, dynamics.e: old.e}
    results["frozen_single_V_action"] = sp.factor(one["before"].subs(at_v)-old.action()["full_before_temporal"])
    physical = dynamics.units(3, 2, 5)
    geometric = geometry.units(1, 0, mass_squared=3, time_scale=2, kinetic_coefficient=5)
    results["independent_nonunit_curl"] = physical["normalized_zeta"]-geometric["lambda_normalized"]
    results["nonunit_curl_value"] = physical["normalized_zeta"]-sp.Rational(5, 12)
    results["nonunit_mass_value"] = physical["isolated_mass_squared_physical"]-sp.Rational(8, 5)
    results["nonunit_mass_conversion"] = physical["isolated_mass_squared_normalized"]-4*physical["isolated_mass_squared_physical"]
    return results


@cache
def proof_checks():
    # Written inertia and kernel proofs supply the general domains; these
    # exact premises are not replaced by a finite scan of coefficient signs.
    bg = old.background()
    numerator, denominator = sp.fraction(sp.factor(bg["J"]))
    polynomial = sp.Poly(numerator, old.u)
    bounds = {
        "eight_trace_D_invertible": geometry.schur()["D_clock"].det() != 0,
        "original_quotient_invertible_on_clock": two_trace_action()["M"].det() != 0,
        "projector_rank_52": sp.trace(two_trace_action()["projector"]) == 52,
        "J_positive_constant": polynomial.coeff_monomial(1) > 0,
        "J_only_even_nonnegative_terms": all(monomial[0] % 2 == 0 and coefficient >= 0
                                               for monomial, coefficient in polynomial.terms()),
        "J_denominator_positive_for_all_real_u": denominator.is_positive is True,
        "rank_two_D_has_both_strict_signs": rank_two.D.det() < 0,
        "nonnull_rank_one_and_null_charts_separated": dynamics.coefficients()["gamma"].subs({dynamics.A: 1, dynamics.B: 1}) < 0,
        "no_low_frequency_EFT_conclusion": True,
    }
    return {name: bool(value) for name, value in bounds.items()}


@cache
def calibration():
    return {"all_connection_equations": 64, "projective_quotient": 60,
            "two_trace_retained": 8, "regular_two_trace_complement": 52,
            "null_single_trace_kernel_radical": 4,
            "unchanged_M1_and_physical_metric": True,
            "M_squared": 3, "tau": 2, "zeta_physical": 5,
            "nonunit": dynamics.units(3, 2, 5),
            "regular_complement_claim_is_quadratic_rolling_only": True}
