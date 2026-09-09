"""Exact all-real-frequency canonical chart for the actual global system."""
from functools import cache

import sympy as sp
from p8_proca_global_hadamard import canonical
from p8_proca_global_support import charts, model, phase

o=model.old
mu=sp.Symbol("regular_real_frequency",positive=True)
clean=canonical.clean

@cache
def data():
    original=phase.data()
    Omega=original["constant_symplectic_form"]
    qreg=mu**2/o.a**2
    Bq=charts.data()["gamma"]["beta_q_principal"]
    exchange=sp.Matrix([[0,0,-1/(2*o.a**3*qreg),0],
        [0,1,0,0],[2*o.a**3*qreg,0,0,0],[0,0,0,1]])
    shift=sp.eye(4)
    shift[2:,:2]=-o.a**3*qreg*Bq
    scale=sp.diag(mu*o.a**sp.Rational(3,2),mu*o.a**sp.Rational(3,2),
        o.a**sp.Rational(-3,2),o.a**sp.Rational(-3,2))
    E=clean(scale*shift*exchange)
    Ei=clean(E.inv())
    M=clean((model.dt(E)+E*original["regular_density_generator"])*Ei)
    M=clean(M.subs(o.q,(mu**2-1)/o.a**2))
    coeffs={j:M.applyfunc(lambda v,j=j:sp.Poly(sp.cancel(mu**3*v),mu).nth(j+3))
        for j in range(-3,2)}
    assert clean(M-sum((mu**j*v for j,v in coeffs.items()),sp.zeros(4)))==sp.zeros(4)
    assert clean(E*Omega*E.T-mu*Omega)==sp.zeros(4)
    assert clean(E*Ei-sp.eye(4))==sp.eye(4)*0
    for j,v in coeffs.items():
        assert clean(v*Omega+Omega*v.T)==sp.zeros(4),j
    assert clean(coeffs[1]-canonical.data()["leading_generator"])==sp.zeros(4)
    actual={j:clean(v.subs(model.substitution(),simultaneous=True)) for j,v in coeffs.items()}
    return {"density_to_regular":E,"regular_to_density":Ei,"generator":M,
        "coefficients":coeffs,"actual_coefficients":actual}

if __name__=="__main__":
    d=data()
    print("REGULAR_EXACT_SCALED_SYMPLECTIC_LAURENT_AND_LEADING_IDENTITIES_PASS",flush=True)
    for j,v in d["actual_coefficients"].items():
        print("ORDER",j,flush=True)
        print(v,flush=True)
        print("AT_CENTER",v.subs(model.u,0),flush=True)
