"""Fresh canonical Laurent chart for the actual retuned global bounce."""
from functools import cache

import sympy as sp
from p8_proca_global_support import charts, model, phase

o=model.old
k=sp.Symbol("positive_global_canonical_frequency",positive=True)


def clean(value):
    return value.applyfunc(sp.factor) if isinstance(value,sp.MatrixBase) else sp.factor(value)


@cache
def data():
    original=phase.data()
    Omega=original["constant_symplectic_form"]
    Bq=charts.data()["gamma"]["beta_q_principal"]
    exchange=sp.Matrix([[0,0,-1/(2*o.a**3*o.q),0],
        [0,1,0,0],[2*o.a**3*o.q,0,0,0],[0,0,0,1]])
    shift=sp.eye(4)
    shift[2:,:2]=-o.a**3*o.q*Bq
    scale=sp.diag(k*o.a**sp.Rational(3,2),k*o.a**sp.Rational(3,2),
        o.a**sp.Rational(-3,2),o.a**sp.Rational(-3,2))
    E=clean(scale*shift*exchange)
    Ei=clean(E.inv())
    M=clean((model.dt(E)+E*original["regular_density_generator"])*Ei)
    M=clean(M.subs(o.q,k*k/o.a**2))
    coefficients={j:M.applyfunc(lambda v,j=j:sp.Poly(sp.cancel(k**3*v),k).nth(j+3))
        for j in range(-3,2)}
    J=coefficients[1]
    K=2*charts.data()["gamma"]["K_principal"]
    G=2*charts.data()["gamma"]["G_principal"]
    E=clean(E.subs(o.q,k*k/o.a**2))
    Ei=clean(Ei.subs(o.q,k*k/o.a**2))
    return {"density_to_packet":E,"packet_to_density":Ei,
        "actual_symplectic_density_to_packet":E/sp.sqrt(k),
        "constant_symplectic_form":Omega,"complete_canonical_generator":M,
        "complete_Laurent_coefficients":coefficients,
        "principal_kinetic":K,"principal_gradient_before_actual_substitution":G,
        "leading_generator":J,"leading_Hamiltonian":clean(-Omega*J),
        "gamma_leading_symmetric_momentum_shift":Bq}


@cache
def checks():
    d=data()
    M,J,Omega=d["complete_canonical_generator"],d["leading_generator"],d["constant_symplectic_form"]
    E,Ei=d["density_to_packet"],d["packet_to_density"]
    K=d["principal_kinetic"]
    actual=model.substitution()
    G0=sp.Matrix([[(2*o.J+o.w**2)/o.lam**2,-o.w/o.lam],
        [-o.w/o.lam,1]])
    out={
        "actual_global_density_packet_is_scaled_symplectic":clean(E*Omega*E.T-k*Omega),
        "actual_global_packet_inverse_is_exact":clean(E*Ei-sp.eye(4)),
        "all_global_finite_Laurent_orders_retained":clean(M-sum(
            (k**j*value for j,value in d["complete_Laurent_coefficients"].items()),sp.zeros(4))),
        "actual_global_leading_inverse_kinetic":clean(J[:2,2:]*K-sp.eye(2)),
        "actual_global_leading_diagonal_blocks_vanish":J[:2,:2].row_join(J[2:,2:]),
        "actual_global_leading_gradient_full_chain":clean(
            (J[2:,:2]+G0/o.a**2).subs(actual,simultaneous=True)),
    }
    for j,value in d["complete_Laurent_coefficients"].items():
        out["global_Laurent_Hamiltonian_order_"+str(j)]=clean(value*Omega+Omega*value.T)
    return out
