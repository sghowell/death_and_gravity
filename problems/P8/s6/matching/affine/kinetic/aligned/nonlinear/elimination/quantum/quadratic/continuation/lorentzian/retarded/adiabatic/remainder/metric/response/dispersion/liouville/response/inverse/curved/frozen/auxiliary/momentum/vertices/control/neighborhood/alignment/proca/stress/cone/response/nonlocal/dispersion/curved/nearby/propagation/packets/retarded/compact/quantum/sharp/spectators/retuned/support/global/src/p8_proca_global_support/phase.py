"""Regular global density-canonical scalar phase before any singular chart."""
from functools import cache

import sympy as sp

from . import model

o=model.old
pv,ps=sp.symbols("global_metric_density_momentum global_matter_density_momentum",real=True)


@cache
def data():
    H=model.margin.scalar()["H"].subs({o.pi:0,o.sigma:0})
    phase=sp.Matrix([o.v,o.matter,pv,ps])
    substitution={o.shift:-pv/(2*o.a**3*o.q),o.pm:ps/o.a**3}
    density=sp.factor(o.a**3*H.subs(substitution,simultaneous=True))
    Omega=sp.zeros(4)
    Omega[:2,2:]=sp.eye(2)
    Omega[2:,:2]=-sp.eye(2)
    M=(Omega*sp.hessian(density,tuple(phase))).applyfunc(sp.factor)
    coeffs=[]
    for j in range(3):
        coeffs.append(M.applyfunc(lambda v,j=j:sp.Poly(sp.cancel(v),o.q).nth(j)))
    C=sp.Matrix([[1,0,0,0],[0,1,0,0]])
    Cgamma=sp.Matrix([[0,0,-1/(2*o.a**3*o.q),0],[0,1,0,0]])
    maps={"unitary":(model.k*C).col_join(C*M),
        "gamma":(model.k*Cgamma).col_join(model.dt(Cgamma)+Cgamma*M)}
    source=sp.Symbol("normalized_global_probe",real=True)
    force=-Omega*sp.Matrix([sp.diff(o.a**3*source*o.matter,state) for state in phase])/source
    return {"regular_density_Hamiltonian":density,"regular_density_phase":phase,
        "constant_symplectic_form":Omega,"regular_density_generator":M,
        "complete_polynomial_q_coefficients":coeffs,
        "original_relational_observable_row":sp.Matrix([[0,1,0,0]]),
        "normalized_physical_probe_force":force,
        "chart_velocity_maps_before_q_to_comoving_substitution":maps,
        "original_retarded_multiplier":"G_hat(t,s,k)=a(s)^3*U_density(t,s,k)[1,3]",
        "original_global_action_source_density":"kappa*a(u)^3*J*chi",
        "q_is_bilinear_complex_spatial_square_over_actual_scale_squared":True}


@cache
def checks():
    d=data()
    M,Omega=d["regular_density_generator"],d["constant_symplectic_form"]
    original=model.margin.scalar()
    return {"global_density_generator_exact_polynomial_degree_two_in_q":
        (M-sum((o.q**j*v for j,v in enumerate(d["complete_polynomial_q_coefficients"])),sp.zeros(4))).applyfunc(sp.factor),
        "global_density_generator_preserves_equal_time_CCR":(M*Omega+Omega*M.T).applyfunc(sp.factor),
        "global_regular_lapse_elimination_replayed_before_charts":original["regular_lapse_equation"],
        "global_regular_Hamiltonian_replayed_before_charts":original["regular_Hamiltonian"],
        "global_scalar_shift_constraint_replayed":original["shift_constraint"],
        "global_matter_source_is_original_local_density_momentum_force":
            d["normalized_physical_probe_force"]-sp.Matrix([0,0,0,o.a**3])}
