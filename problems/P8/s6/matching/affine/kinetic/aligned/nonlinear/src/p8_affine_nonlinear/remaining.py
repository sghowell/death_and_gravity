"""The five shear velocities and the free physical-metric matter velocity."""
from functools import cache

import sympy as sp


@cache
def shear():
    # At a point choose an orthonormal basis for symmetric traceless tensors.
    # These are velocity-conjugate P_A=2*pi_TF,A/sqrt_hat, not pi_TF,A.
    U, minus_B = sp.symbols("U minus_B", positive=True)
    v, P = [sp.Matrix(sp.symbols(prefix+"0:5", real=True)) for prefix in ("shear", "P_shear")]
    L = U*minus_B*(v.T*v)[0]
    solution = P/(2*U*minus_B)
    substitution = dict(zip(v, solution, strict=True))
    H = (P.T*P)[0]/(4*U*minus_B)
    return {"U": U, "minus_B": minus_B, "L": L, "v": v, "P": P, "H": H,
            "five_shear_momentum_equations":
                (sp.Matrix([sp.diff(L, item) for item in v]).subs(substitution)-P).applyfunc(sp.factor),
            "five_shear_Hamiltonian": sp.factor(((P.T*v)[0]-L).subs(substitution)-H)}


@cache
def matter():
    N, U = sp.symbols("N U", positive=True)
    velocity, p = sp.symbols("advected_chi p_chi", real=True)
    L = U*velocity**2/(2*N)
    solution = N*p/U
    H = N*p**2/(2*U)
    return {"N": N, "U": U, "L": L, "p": p, "H": H,
            "free_matter_momentum_equation": sp.factor(sp.diff(L, velocity).subs(velocity, solution)-p),
            "free_matter_Hamiltonian": sp.factor((p*velocity-L).subs(velocity, solution)-H)}


@cache
def checks():
    out = {name: shear()[name] for name in ("five_shear_momentum_equations", "five_shear_Hamiltonian")}
    out.update({name: matter()[name] for name in ("free_matter_momentum_equation", "free_matter_Hamiltonian")})
    return out
