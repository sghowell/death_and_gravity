"""All-order matrix Riccati induction with positive frequency SUM denominators."""
from functools import cache

import sympy as sp

rho=sp.Symbol("inverse_positive_frequency",real=True)


def order(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Integer)):
        raise TypeError("Require a positive exact integer Riccati order")
    if value<1:
        raise ValueError("The leading graph is supplied separately")
    return int(value)


def clean(value):
    return value.applyfunc(sp.factor) if isinstance(value,sp.MatrixBase) else sp.factor(value)


@cache
def principal():
    mass,a,c=sp.symbols("positive_clock_kinetic positive_scale positive_clock_speed",positive=True)
    mixing=sp.Symbol("real_clock_matter_mixing",real=True)
    V=sp.Matrix([[1,0],[mixing,1]])
    Vi=V.inv()
    D=sp.diag(mass,1)
    omega=sp.diag(c/a,1/a)
    K=Vi.T*D*Vi
    B=Vi.T*D*omega*Vi
    G=Vi.T*D*omega**2*Vi
    R=-sp.I*B
    return {"clock_kinetic":mass,"scale":a,"clock_speed":c,"mixing":mixing,
        "configuration_eigenbasis":V,"positive_frequency_matrix":omega,
        "leading_kinetic":K,"leading_spatial_Hamiltonian":G,
        "positive_graph_imaginary_part":B,"leading_negative_frequency_graph":R,
        "sum_denominators":[2*c/a,(c+1)/a,2/a]}


def solve_sylvester(forcing,V,omega):
    """Solve R0 K^-1 X + X K^-1 R0 = forcing, R0=-i B."""
    if forcing.shape!=(2,2) or V.shape!=(2,2) or omega.shape!=(2,2):
        raise ValueError("Require the two-field matrix domain")
    transformed=V.T*forcing*V
    solved=sp.Matrix(2,2,lambda i,j:sp.I*transformed[i,j]/(omega[i,i]+omega[j,j]))
    return clean(V.inv().T*solved*V.inv())


def next_coefficient(n,previous,hamiltonian,derivative,V,omega):
    """H/k=sum rho^j hamiltonian[j], R=sum rho^j previous[j]."""
    n=order(n)
    if len(previous)!=n:
        raise ValueError("Require every preceding coefficient through n-1")
    if not hamiltonian or any(value.shape!=(4,4) for value in hamiltonian):
        raise ValueError("Require complete rank-four Hamiltonian coefficients")
    if any(value.shape!=(2,2) for value in previous):
        raise ValueError("Require symmetric two-field graph coefficients")
    R=sum((rho**j*v for j,v in enumerate(previous)),sp.zeros(2))
    H=sum((rho**j*v for j,v in enumerate(hamiltonian)),sp.zeros(4))
    defect=rho*sum((rho**j*derivative(v) for j,v in enumerate(previous)),sp.zeros(2))
    defect+=R*H[2:,2:]*R+R*H[2:,:2]+H[:2,2:]*R+H[:2,:2]
    forcing=defect.applyfunc(lambda v:sp.expand(v).coeff(rho,n))
    return solve_sylvester(-forcing,V,omega)


@cache
def checks():
    d=principal()
    K,B,R=d["leading_kinetic"],d["positive_graph_imaginary_part"],d["leading_negative_frequency_graph"]
    V,omega=d["configuration_eigenbasis"],d["positive_frequency_matrix"]
    x,y,z=sp.symbols("symmetric_forcing_x symmetric_forcing_y symmetric_forcing_z",real=True)
    forcing=sp.Matrix([[x,y],[y,z]])
    solved=solve_sylvester(forcing,V,omega)
    return {
        "negative_frequency_leading_graph_solves_actual_type_Riccati":clean(R*K.inv()*R+d["leading_spatial_Hamiltonian"]),
        "positive_graph_B_has_exact_positive_square_decomposition":clean(V.T*B*V-sp.diag(d["clock_kinetic"]*d["clock_speed"]/d["scale"],1/d["scale"])),
        "leading_configuration_evolution_has_negative_frequency_sign":clean(K.inv()*R*V+sp.I*V*omega),
        "frequency_SUM_Sylvester_inverse_is_exact":clean(R*K.inv()*solved+solved*K.inv()*R-forcing),
        "symmetric_forcing_gives_symmetric_Riccati_coefficient":clean(solved-solved.T),
        "positive_graph_determinant":sp.factor(B.det()-d["clock_kinetic"]*d["clock_speed"]/d["scale"]**2),
    }


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,None,0,-1,sp.Rational(1,2))
    rejected=0
    for value in bad:
        try:
            order(value)
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(bad):
        raise ValueError("An invalid Riccati order was accepted")
    return {"rejected_inputs":rejected}
