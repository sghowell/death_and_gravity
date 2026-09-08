"""Causal covariance response under a moving positive mode clock."""
from functools import cache

import sympy as sp

rate, physical_w2 = sp.symbols("acoustic_rate physical_frequency_squared", positive=True)
rate_dot, physical_d = sp.symbols("acoustic_rate_dot physical_canonical_rate", real=True)
bk, D = sp.symbols("acoustic_canonical_rate acoustic_frequency_squared", real=True)
dbk, dU = sp.symbols("varied_acoustic_canonical_rate varied_acoustic_mass_potential", real=True)
r, xi = sp.symbols("varied_log_acoustic_rate prepared_acoustic_time_shift", real=True)


def matrix(half_rate, frequency_squared):
    return sp.ImmutableMatrix([[2*half_rate, 2, 0], [-frequency_squared, 0, 1],
                               [0, -2*frequency_squared, -2*half_rate]])


def covariance_map():
    return sp.diag(rate, 1, 1/rate)


def varied_generator():
    return sp.ImmutableMatrix([[2*dbk, 0, 0], [-dU, 0, 0], [0, -2*dU, -2*dbk]])


@cache
def checks():
    S = covariance_map()
    St = sp.diag(rate_dot, 0, -rate_dot/rate**2)
    target_b = (physical_d+rate_dot/(2*rate))/rate
    transformed = St*S.inv()/rate+S*matrix(physical_d, physical_w2)*S.inv()/rate
    out = {"exact_acoustic_covariance_generator": sp.ImmutableMatrix(
        (transformed-matrix(target_b, physical_w2/rate**2)).applyfunc(sp.factor))}
    e = sp.Symbol("parameter", real=True)
    out["varied_generator_replays_both_sources"] = sp.ImmutableMatrix(
        sp.diff(matrix(bk+e*dbk, D+e*dU), e).subs(e, 0)-varied_generator())
    xx, yy, zz = sp.symbols("mode_variance mode_cross mode_momentum_variance", real=True)
    vector = sp.ImmutableMatrix([xx, yy, zz])
    rhs = matrix(bk, D)*vector
    out["mode_covariance_determinant_preserved"] = sp.factor(
        rhs[0]*zz+xx*rhs[2]-2*yy*rhs[1])
    dx, dy, dz = sp.symbols("varied_variance varied_cross varied_momentum_variance", real=True)
    dvector = sp.ImmutableMatrix([dx, dy, dz])
    drhs = matrix(bk, D)*dvector+varied_generator()*vector
    out["linearized_pure_covariance_constraint_preserved"] = sp.factor(
        drhs[0]*zz+dx*rhs[2]+rhs[0]*dz+xx*drhs[2]-2*rhs[1]*dy-2*yy*drhs[1])
    # W=delta Y at fixed physical u; Z=delta Y at fixed acoustic sigma.
    M, Ms, change = (sp.MatrixSymbol(name, 3, 3) for name in ("M", "M_sigma", "delta_M_fixed_u"))
    Y, Z = (sp.MatrixSymbol(name, 3, 1) for name in ("Y", "delta_Y_fixed_sigma"))
    W = Z+xi*M*Y
    Wdot = rate*M*W+rate*r*M*Y+rate*change*Y
    Zdot = Wdot-rate*r*M*Y-xi*rate*Ms*Y-xi*rate*M*M*Y
    out["fixed_acoustic_time_tangent_equation"] = sp.ImmutableMatrix(sp.expand(
        Zdot/rate-M*Z-(change-xi*Ms)*Y).as_explicit())
    # The physical covariance perturbation has opposite diagonal rate contacts.
    dY0, dY1, dY2, Y0, Y1, Y2, Yp0, Yp1, Yp2 = sp.symbols(
        "dY0 dY1 dY2 Y0 Y1 Y2 Yp0 Yp1 Yp2", real=True)
    Yvec = sp.ImmutableMatrix([Y0, Y1, Y2])
    dYvec = sp.ImmutableMatrix([dY0, dY1, dY2])
    Yprime = sp.ImmutableMatrix([Yp0, Yp1, Yp2])
    source_map = sp.diag(rate*r, 0, -r/rate)
    X = S.inv()*Yvec
    dX = S.inv()*(dYvec+xi*Yprime)-r*sp.diag(1, 0, -1)*X
    out["physical_covariance_pullback_including_rate_and_history"] = sp.ImmutableMatrix(
        (S*dX+source_map*X-xi*Yprime-dYvec).applyfunc(sp.factor))
    return out


@cache
def green_checks():
    f, fc, g, gc = sp.symbols("mode_t conjugate_mode_t mode_s conjugate_mode_s", complex=True)
    G = sp.I*(f*gc-fc*g)
    # Formal conjugation swaps the independently named mode/conjugate columns.
    conjugate = {f: fc, fc: f, g: gc, gc: g, sp.I: -sp.I}
    Gbar = G.xreplace(conjugate)
    covariance = -G*(fc*g+f*gc)
    target = -sp.I*((f*gc)**2-(fc*g)**2)
    out = {"retarded_Green_is_real": sp.expand(G-Gbar),
           "mode_variance_retarded_kernel": sp.expand(covariance-target)}
    w, t = sp.symbols("positive_frequency positive_time", positive=True)
    flat = {f: sp.exp(-sp.I*w*t)/sp.sqrt(2*w), fc: sp.exp(sp.I*w*t)/sp.sqrt(2*w),
            g: 1/sp.sqrt(2*w), gc: 1/sp.sqrt(2*w)}
    out["flat_sine_Green_normalization"] = sp.simplify(sp.expand_complex(G.subs(flat))-sp.sin(w*t)/w)
    out["flat_mode_variance_response_sign"] = sp.simplify(
        sp.expand_complex(covariance.subs(flat))+sp.sin(2*w*t)/(2*w**2))
    out["retarded_Green_zero_on_diagonal"] = sp.expand(G.subs({g: f, gc: fc}))
    amplitude = sp.Symbol("positive_mode_amplitude", positive=True)
    cross_rate, half_rate = sp.symbols("real_mode_cross_rate acoustic_half_rate", real=True)
    p = (cross_rate-sp.I/(2*amplitude**2))*amplitude
    pc = (cross_rate+sp.I/(2*amplitude**2))*amplitude
    out["retarded_Green_unit_derivative_jump"] = sp.expand(
        sp.I*((p+half_rate*amplitude)*amplitude-(pc+half_rate*amplitude)*amplitude)-1)
    aa, ac, bb, bc = sp.symbols("mix_A mix_A_conjugate mix_B mix_B_conjugate", complex=True)
    mixed = {f: aa*f+bb*fc, fc: ac*fc+bc*f, g: aa*g+bb*gc, gc: ac*gc+bc*g}
    out["Green_Bogoliubov_Wronskian_invariance"] = sp.expand(
        G.subs(mixed, simultaneous=True)-(aa*ac-bb*bc)*G)
    fixture = {f: (1-sp.I)/2, fc: (1+sp.I)/2, g: 1/sp.sqrt(2), gc: 1/sp.sqrt(2)}
    squeezed = covariance.subs(mixed, simultaneous=True).subs({aa: sp.Rational(5, 4), ac: sp.Rational(5, 4),
                                                              bb: sp.Rational(3, 4), bc: sp.Rational(3, 4)})
    out["nonzero_initial_squeezing_changes_variance_kernel"] = sp.simplify(
        (squeezed-covariance).subs(fixture)+sp.Rational(3, 2))
    return out
