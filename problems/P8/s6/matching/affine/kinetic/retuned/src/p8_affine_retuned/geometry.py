"""Exact source-centered mass update on all sixty quotient directions.

The literal new term is -mu*g^(ab)(T-Tstar)_a(T-Tstar)_b/2,
T=V+U, mu=11/9. It is not a change to the frozen S6.37 action.
"""
from functools import cache

import sympy as sp
from p8_affine import connection as old
from p8_affine_traces import geometry as traces

MU = sp.Rational(11, 9)
P = old.P
ETA = traces.ETA


def clean(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.factor)) if isinstance(value, sp.MatrixBase) else sp.factor(value)


def expanded_rational(value):
    """Exact normalization without factoring a many-source-variable numerator."""
    numerator, denominator = sp.together(value).as_numer_denom()
    return sp.expand(numerator)/sp.factor(denominator)


@cache
def update():
    data = old.quotient()
    combine = sp.eye(4).row_join(sp.eye(4))
    full = combine*traces.traces()["full"]
    n = full*data["embedding"]
    prior = traces.schur()
    w = prior["inverse_N_transpose"]*combine.T
    d = clean(n*w)
    time = 3*(2*P**3-1)/P
    space = (8*P+5)/(8*P**2)
    added = MU*n.T*ETA*n
    m = clean(data["hessian"]+added)
    middle = clean((ETA/MU+d).inv())
    inverse_n = clean(w*(sp.eye(4)+MU*ETA*d).inv())
    response = clean(n*inverse_n)
    gamma_t = 1/(1/time+MU)
    gamma_s = 1/(MU-1/space)
    return {"embedding": data["embedding"], "N_full": clean(full), "N": clean(n),
            "M_old": data["hessian"], "M": m, "added_hessian": clean(added),
            "W_old": clean(w), "D_old": d, "old_time": time, "old_space": space,
            "old_expected": sp.diag(time, space, space, space), "middle": middle,
            "W": inverse_n, "D": response,
            "gamma_t": sp.factor(gamma_t), "gamma_s": sp.factor(gamma_s),
            "expected_D": clean(sp.diag(gamma_t, -gamma_s, -gamma_s, -gamma_s)),
            "determinant_ratio": sp.factor((1+MU*time)*(1-MU*space)**3)}


@cache
def source_centering():
    data = update()
    stationary = old.eliminate()["solution"]
    tstar = clean(data["N_full"]*stationary)
    source = old.quadratic()["source"]-MU*data["N_full"].T*ETA*tstar
    full_m = old.quadratic()["hessian"]+MU*data["N_full"].T*ETA*data["N_full"]
    return {"stationary": stationary, "Tstar": tstar,
            "source_new": sp.ImmutableMatrix(source.applyfunc(expanded_rational)),
            "M_full": clean(full_m), "counterterm": expanded_rational(MU*(tstar.T*ETA*tstar)[0]/2),
            "full_Euler_residual": clean(full_m*stationary+source),
            "full_lift_residual": clean(full_m*data["embedding"]*data["W"]-data["N_full"].T)}


@cache
def complement():
    data = update()
    lift = clean(data["W"]*data["D"].inv())
    projection = clean(sp.eye(60)-lift*data["N"])
    return {"lift": lift, "projector": projection,
            "trace_right_inverse": clean(data["N"]*lift-sp.eye(4)),
            "complement_annihilation": clean(data["N"]*projection),
            "complement_idempotence": clean(projection**2-projection),
            "complement_cross_Euler": clean(projection.T*data["M"]*lift),
            "retained_mass_action": clean(lift.T*data["M"]*lift-data["D"].inv()),
            "complement_dimension": sp.factor(sp.trace(projection)-56)}


@cache
def checks():
    data, source, split = update(), source_centering(), complement()
    out = {"literal_generic_trace_response": clean(data["D_old"]-data["old_expected"]),
           "mass_retuned_response": clean(data["D"]-data["expected_D"]),
           "quotient_Euler_lift": clean(data["M"]*data["W"]-data["N"].T),
           "all64_Euler_lift": source["full_lift_residual"],
           "all64_unchanged_stationary_Euler": source["full_Euler_residual"],
           "projective_columns": clean(data["N_full"]*old.quadratic()["gauge"]),
           "clock_mass_Lorentz_matrix": clean(data["D"].subs(P, sp.Rational(1, 2))-ETA),
           "clock_stationary_unchanged_mass_center": clean(source["Tstar"]
                -data["N_full"]*source["stationary"])}
    out.update({name: value for name, value in split.items() if name not in ("lift", "projector")})
    return out
