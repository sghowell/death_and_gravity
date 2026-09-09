"""Exact matrix energy and canonical initial-column identities."""
from functools import cache

import sympy as sp


@cache
def checks():
    H,volume=sp.symbols("H positive_volume",real=True)
    a,b,c=sp.symbols("alpha11 alpha12 alpha22",real=True)
    ad,bd,cd=sp.symbols("alpha11_dot alpha12_dot alpha22_dot",real=True)
    x,y,z=sp.symbols("potential11 potential12 potential22",real=True)
    xd,yd,zd=sp.symbols("potential11_dot potential12_dot potential22_dot",real=True)
    r,rd=sp.symbols("mixing mixing_dot",real=True)
    alpha=sp.Matrix([[a,b],[b,c]])
    alpha_dot=sp.Matrix([[ad,bd],[bd,cd]])
    V=sp.Matrix([[x,y],[y,z]])
    Vdot=sp.Matrix([[xd,yd],[yd,zd]])
    A=sp.Matrix([[0,r],[-r,0]])
    Adot=sp.Matrix([[0,rd],[-rd,0]])
    Q=sp.Matrix(sp.symbols("coordinate0:2",real=True))
    v=sp.Matrix(sp.symbols("velocity0:2",real=True))
    acceleration=-alpha.inv()*((alpha_dot+3*H*alpha+2*A)*v+(Adot+3*H*A+V)*Q)
    derivative=volume*((v.T*alpha*acceleration)[0]+(v.T*V*Q)[0]
                        +(v.T*alpha_dot*v)[0]/2+(Q.T*Vdot*Q)[0]/2)
    derivative+=3*H*volume*((v.T*alpha*v)[0]+(Q.T*V*Q)[0])/2
    expected=volume*(-(v.T*(alpha_dot+3*H*alpha)*v)[0]
                     +(Q.T*(Vdot+3*H*V)*Q)[0]-2*(v.T*(Adot+3*H*A)*Q)[0])/2
    out={"full_two_component_energy_identity_with_antisymmetric_mixing":sp.factor(derivative-expected)}
    q,g,gd,R,Rd=sp.symbols("q gradient gradient_dot remainder remainder_dot",real=True)
    out["fixed_comoving_energy_potential_decomposition"]=sp.expand(
        -2*H*q*g+q*gd+Rd+3*H*(q*g+R)-q*(gd+H*g)-Rd-3*H*R)
    b11,b12,b22=sp.symbols("positive_b11 b12 positive_b22",real=True,nonzero=True)
    k,V0=sp.symbols("positive_frequency positive_initial_volume",positive=True)
    S=sp.Matrix([[sp.Symbol("s11",real=True),sp.Symbol("s12",real=True)],
                 [sp.Symbol("s12",real=True),sp.Symbol("s22",real=True)]])
    B=sp.Matrix([[b11,b12],[0,b22]])
    Q0=B.inv()/sp.sqrt(2*k*V0)
    P0=-sp.I*sp.sqrt(k*V0/2)*B.T+V0*S*Q0
    simplify=lambda matrix:sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))
    out["canonical_initial_column_Wronskian"]=simplify(
        Q0.conjugate().T*P0-P0.conjugate().T*Q0+sp.I*sp.eye(2))
    out["canonical_initial_column_isotropy"]=simplify(Q0.T*P0-P0.T*Q0)
    out["canonical_initial_coordinate_momentum_commutator"]=simplify(
        Q0*P0.conjugate().T-Q0.conjugate()*P0.T-sp.I*sp.eye(2))
    out["canonical_initial_momentum_momentum_commutator"]=simplify(
        P0*P0.conjugate().T-P0.conjugate()*P0.T)
    beta=S+A
    velocity0=(B.T*B).inv()*(P0/V0-beta*Q0)
    target=-sp.I*sp.sqrt(k/(2*V0))*B.inv()-(B.T*B).inv()*A*Q0
    out["initial_velocity_retains_antisymmetric_mixing"]=simplify(velocity0-target)
    wrong_P=-sp.I*sp.eye(2)/sp.sqrt(2)+A/sp.sqrt(2)
    defect=sp.eye(2)/sp.sqrt(2)*wrong_P-wrong_P.T*sp.eye(2)/sp.sqrt(2)
    out["wrong_full_beta_initial_shift_has_nonzero_isotropy_defect"]=simplify(defect-A)
    return out


@cache
def initial_mixing_negative_control():
    return sp.Matrix([[0,1],[-1,0]])
