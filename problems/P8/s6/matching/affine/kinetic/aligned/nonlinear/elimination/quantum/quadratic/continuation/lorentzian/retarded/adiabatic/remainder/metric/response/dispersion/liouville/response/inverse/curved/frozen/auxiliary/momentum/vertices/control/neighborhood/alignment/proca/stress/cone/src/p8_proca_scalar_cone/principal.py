"""Euler-first central scalar principal matrix in the physical matter frame."""
from functools import cache

import sympy as sp
from p8_offclock_scalar import principal as old
from p8_offclock_scalar import spatial as old_spatial

from . import background, spatial

q=old.q
N,P=background.N,background.P
speed=sp.Symbol("physical_speed_squared",real=True)


def clean_matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(spatial.clean))


def limit_matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(lambda x:old.infinity(spatial.clean(x))))


@cache
def blocks():
    T=sp.zeros(4)
    T[0,2]=1/(2*q)
    T[1,0]=-2*q
    T[2,1]=T[3,3]=1
    H=sp.Matrix((T.T*spatial.phase_matrix()*T).subs(spatial.k,sp.sqrt(q)))
    # Actual moving gamma generator from q_hat'= -2 H_hat q_hat.
    H[0,2]-=spatial.eta*spatial.at_lapse(background.data()["hat_Hubble_time_derivative"])
    H[2,0]=H[0,2]
    H=clean_matrix(H)
    A,B,C=H[2:,2:],H[2:,:2],H[:2,:2]
    A0,B0,C0=(clean_matrix(value.subs(spatial.eta,0)) for value in (A,B,C))
    Ad,Bd=(clean_matrix(value.diff(spatial.eta).subs(spatial.eta,0)) for value in (A,B))
    alpha=clean_matrix(A0.inv())
    beta=clean_matrix(-alpha*B0)
    betad=clean_matrix(alpha*Ad*alpha*B0-alpha*Bd)
    gamma=clean_matrix(C0-B0.T*alpha*B0)
    return {"gamma_canonical_map":sp.ImmutableMatrix(T),
            "actual_canonical_first_time_Hessian":H,
            "momentum_Hessian":A0,"momentum_coordinate_block":B0,"coordinate_Hessian":C0,
            "actual_momentum_Hessian_time_derivative":Ad,
            "actual_momentum_coordinate_time_derivative":Bd,
            "velocity_kinetic":alpha,"velocity_mixed":beta,
            "actual_velocity_mixed_time_derivative":betad,"velocity_potential":gamma,
            "high_q_kinetic":limit_matrix(alpha),
            "high_q_mixed_over_q":limit_matrix(beta/q),
            "high_q_mixed_time_derivative_over_q":limit_matrix(betad/q),
            "high_q_potential_over_q":limit_matrix(gamma/q),
            "Euler_first_gradient":limit_matrix((gamma+betad)/q)}


@cache
def formula():
    d=background.data()
    h,P2=d["fixed_phase_lapse_Hessian"],d["matter_momentum_squared"]
    r=-1/(4*N**3*h)
    A=sp.Matrix([[r,-r*P],[-r*P,N**-sp.Rational(1,2)+r*P2]])
    K=sp.Matrix([[-4*N**3*h+sp.sqrt(N)*P2,sp.sqrt(N)*P],
                 [sp.sqrt(N)*P,sp.sqrt(N)]])
    gbb=sp.factor(4*N**sp.Rational(3,2)*d["trace_force_time_derivative"]/3)
    G=sp.Matrix([[gbb,N**sp.Rational(3,2)*P],[N**sp.Rational(3,2)*P,N**sp.Rational(3,2)]])
    schur=sp.factor(gbb-N**sp.Rational(3,2)*P2)
    c2=sp.factor(schur/(-4*N**4*h))
    return {"high_q_momentum_Hessian":sp.ImmutableMatrix(A),
            "coordinate_time_kinetic":sp.ImmutableMatrix(K),
            "hat_momentum_gradient":sp.ImmutableMatrix(G),
            "clock_gradient_Schur":schur,
            "clock_physical_speed_squared":c2,"matter_physical_speed_squared":sp.Integer(1),
            "physical_clock_cone_margin":sp.factor(1-c2),
            "actual_lapse":N,"physical_spatial_metric_over_hat_metric":N,
            "physical_squared_momentum":q/N,
            "physical_cone_comparison_matrix":sp.ImmutableMatrix(N*K-G)}


@cache
def checks():
    b,d,f=blocks(),background.data(),formula()
    convert=lambda M:sp.ImmutableMatrix(M.applyfunc(spatial.at_lapse))
    oldJ=sp.diag(sp.Matrix([[0,1],[-1,0]]),sp.Matrix([[0,1],[-1,0]]))
    newJ=sp.zeros(4)
    newJ[:2,2:]=sp.eye(2)
    newJ[2:,:2]=-sp.eye(2)
    T=b["gamma_canonical_map"]
    K,G=f["coordinate_time_kinetic"],f["hat_momentum_gradient"]
    residual=sp.factor((G-speed*N*K).det().subs(P**2,d["matter_momentum_squared"])
        -N**2*sp.factor(K.det().subs(P**2,d["matter_momentum_squared"]))*(1-speed)*(f["clock_physical_speed_squared"]-speed))
    out={"actual_scalar_gamma_map_symplectic":clean_matrix(T.T*oldJ*T-newJ),
         "actual_high_q_kinetic_from_full_spatial_time_jet":clean_matrix(b["high_q_kinetic"]-convert(K)),
         "actual_Euler_first_gradient_from_full_spatial_time_jet":clean_matrix(b["Euler_first_gradient"]-convert(G)),
         "no_missing_q_leading_mixed_velocity_term":b["high_q_mixed_over_q"],
         "all_leading_mixed_time_derivative_in_clock_entry":clean_matrix(
             b["high_q_mixed_time_derivative_over_q"]-sp.diag(convert(G)[0,0],0)),
         "actual_spatial_potential_principal_matrix":clean_matrix(
             b["high_q_potential_over_q"]-convert(G-sp.diag(G[0,0],0))),
         "exact_physical_characteristic_factorization":residual,
         "exact_clock_kinetic_Schur":sp.factor(K[0,0]-K[0,1]**2/K[1,1]+4*N**3*d["fixed_phase_lapse_Hessian"]).subs(P**2,d["matter_momentum_squared"]),
         "exact_luminal_matter_cone_difference_offdiagonal":sp.factor((N*K-G)[0,1]),
         "exact_luminal_matter_cone_difference_bottom":sp.factor((N*K-G)[1,1]),
         "clock_gradient_at_original_bounce":G[0,0].subs(N,1)-12,
         "clock_speed_at_original_bounce":f["clock_physical_speed_squared"].subs(N,1)-sp.Rational(749375,749377)}
    # Differentiate the force at fixed canonical phase data, not along
    # the family of solutions parametrized by N.
    a=d["trace_a"]
    wrong=sp.diff(N*(d["trace_momentum_time_derivative"]-d["actual_trace_linear_time_derivative"])/(2*a),N)
    discrepancy=sp.factor(4*N**sp.Rational(3,2)*(wrong-d["trace_force_time_derivative"])/3)
    out["wrong_constraint_family_derivative_negative_control"]=sp.factor(
        discrepancy+sp.sqrt(N)*(5624995*N**4-8625003)/250000)
    # Zero-order independent reduction agrees with the old scalar block.
    indexes=(0,1,2,3)
    actual=spatial.phase_matrix().subs(spatial.eta,0)
    previous=old_spatial.phase_matrix().extract(indexes,indexes)
    # The previously omitted moving trace boundary is now retained. At
    # the center it adds 3*p0_dot*v^2/2, hence a 3*p0_dot Hessian.
    correction=sp.zeros(4)
    correction[0,0]=3*spatial.at_lapse(d["trace_momentum_time_derivative"])
    out["previous_central_spatial_scalar_matrix_plus_actual_boundary"]=clean_matrix(actual-previous-correction)
    return out
