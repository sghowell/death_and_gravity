"""Uniformly diagonalizable real principal speed in either complete Euler chart."""
from functools import cache

import sympy as sp


@cache
def data():
    J,dJ,a=sp.symbols("positive_J positive_margin positive_global_scale",positive=True)
    den=sp.Symbol("nonzero_chart_pivot",nonzero=True,real=True)
    w=sp.Symbol("actual_matter_mixing",real=True)
    c=sp.sqrt(J/(J+dJ))
    rows={}
    for name,sign in (("unitary",1),("gamma",-1)):
        K0=sp.Matrix([[(J+w*w/2)/den**2,sign*w/(2*den)],[sign*w/(2*den),sp.Rational(1,2)]])
        K=K0+sp.diag(dJ/den**2,0)
        V=sp.Matrix([[1,0],[-sign*w/den,1]])
        omega=sp.diag(c/a,1/a)
        S=V.row_join(V).col_join((sp.I*V*omega).row_join(-sp.I*V*omega))
        leading=sp.zeros(4)
        leading[:2,2:]=sp.eye(2)
        leading[2:,:2]=-K.inv()*K0/a**2
        rows[name]={"positive_principal_K":K,"positive_principal_G":K0,
            "exact_velocity_mode_basis":S,"exact_velocity_mode_inverse":S.inv().applyfunc(sp.simplify),
            "principal_generator":leading,"real_principal_frequencies":sp.diag(c/a,1/a,-c/a,-1/a),
            "basis_determinant":sp.factor(S.det()),
            "no_inverse_clock_matter_frequency_gap":True}
    return {"positive_J":J,"positive_margin":dJ,"physical_scale":a,
        "clock_squared_speed":c*c,"charts":rows}


@cache
def checks():
    d=data()
    out={}
    real,imag=sp.symbols("real_scalar_frequency imaginary_scalar_frequency",real=True)
    for name,c in d["charts"].items():
        S,Si=c["exact_velocity_mode_basis"],c["exact_velocity_mode_inverse"]
        lam=c["real_principal_frequencies"]
        principal=sp.I*(real+sp.I*imag)*lam
        out[name+"_actual_principal_velocity_modes_diagonalize_exactly"]=(
            c["principal_generator"]*S-sp.I*S*lam).applyfunc(sp.simplify)
        out[name+"_mode_inverse_is_exactly_the_literal_inverse"]=(S*Si-sp.eye(4)).applyfunc(sp.simplify)
        out[name+"_basis_determinant_has_no_frequency_gap_factor"]=sp.simplify(
            c["basis_determinant"]+4*sp.sqrt(d["clock_squared_speed"])/d["physical_scale"]**2)
        out[name+"_real_frequency_growth_cancels_in_Hermitian_part"]=(
            (principal+principal.conjugate().T)/2+imag*lam).applyfunc(sp.simplify)
    return out
