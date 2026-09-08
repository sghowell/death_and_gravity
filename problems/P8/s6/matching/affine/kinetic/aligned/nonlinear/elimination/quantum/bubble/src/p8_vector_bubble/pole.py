"""Ordered Proca numerator, angular contractions and exact pole polynomial."""
from functools import cache

import sympy as sp

m2 = sp.Symbol("m0_squared", positive=True)
x = sp.Symbol("feynman_x", real=True)
P2, PA, PA2 = sp.symbols("k_squared k_A_k k_A_squared_k", real=True)
TA, TA2 = sp.symbols("trace_A trace_A_squared", real=True)


@cache
def numerator():
    entries = sp.symbols("A0:10", real=True)
    A = sp.zeros(4)
    counter = 0
    for i in range(4):
        for j in range(i, 4):
            A[i, j] = A[j, i] = entries[counter]
            counter += 1
    p, q = [sp.Matrix(sp.symbols(prefix+"0:4", real=True)) for prefix in ("p", "q")]
    Gp, Gq = sp.eye(4)+p*p.T/m2, sp.eye(4)+q*q.T/m2
    actual = m2**2*sp.trace(Gp*A*Gq*A)
    expected = (m2**2*sp.trace(A*A)+m2*((p.T*A*A*p)[0]+(q.T*A*A*q)[0])+(p.T*A*q)[0]**2)
    Kp = ((p.T*p)[0]+m2)*sp.eye(4)-p*p.T
    return {"A": A, "p": p, "q": q, "expected": expected,
            "full_longitudinal_Proca_inverse": (Kp*Gp-((p.T*p)[0]+m2)*sp.eye(4)).applyfunc(sp.expand),
            "complete_two_insertion_numerator": sp.expand(actual-expected)}


@cache
def parameter_integral():
    t = x*(1-x)
    Delta = m2+t*P2
    # Relative to 1/(16pi² epsilon_DR), the angular integrals give
    # I0=1, I(l²)=-2Delta, I(l^4)=3Delta² in d=4 at pole order.
    integrand = (m2**2*TA2-m2*Delta*TA2+m2*(x**2+(1-x)**2)*PA2
                 +Delta**2*(TA**2+2*TA2)/8
                 -(1-2*x)**2*Delta*PA2/2+t*Delta*TA*PA+t**2*PA**2)
    actual = sp.integrate(sp.expand(integrand), (x, 0, 1))
    constant = m2**2*(TA**2+2*TA2)/8
    second = m2*(P2*(TA**2-2*TA2)/24+PA2/2+TA*PA/6)
    fourth = P2**2*(TA**2+2*TA2)/240-P2*PA2/60+P2*TA*PA/30+PA**2/30
    return {"integrand": integrand, "actual": sp.factor(actual),
            "constant": constant, "second": second, "fourth": fourth,
            "full_Feynman_parameter_pole": sp.expand(actual-constant-second-fourth)}


@cache
def actual_mass_direction():
    from p8_aligned_quantum import kernel
    h = sp.Symbol("h", positive=True)
    p = kernel.geometry.P
    masses = kernel.masses()
    alpha, beta = [sp.factor(sp.diff(masses[name], p).subs(p, sp.Rational(1, 2))*(-1/(2*h)))
                   for name in ("a", "b")]
    time2, space2 = sp.symbols("k_time_squared k_space_squared", nonnegative=True)
    substitute = {TA: alpha+3*beta, TA2: alpha**2+3*beta**2,
                  P2: time2+space2, PA: alpha*time2+beta*space2, PA2: alpha**2*time2+beta**2*space2}
    data = parameter_integral()
    out = {"h": h, "alpha": alpha, "beta": beta, "time2": time2, "space2": space2}
    out.update({name: sp.factor(data[name].subs(substitute)) for name in ("constant", "second", "fourth")})
    return out


@cache
def radial_residues():
    epsilon = sp.Symbol("epsilon_DR", positive=True)
    return {"radial_Gamma_residue_"+str(order):
            sp.limit(epsilon*sp.gamma(epsilon-order), epsilon, 0)-sp.Rational((-1)**order, sp.factorial(order))
            for order in range(3)}


@cache
def checks():
    data = numerator()
    out = {name: data[name] for name in ("full_longitudinal_Proca_inverse", "complete_two_insertion_numerator")}
    out["full_Feynman_parameter_pole"] = parameter_integral()["full_Feynman_parameter_pole"]
    out.update(radial_residues())
    return out
