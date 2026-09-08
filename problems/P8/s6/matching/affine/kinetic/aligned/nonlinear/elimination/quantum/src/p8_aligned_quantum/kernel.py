"""Full constant-coefficient vector determinant and physical mode control."""
from functools import cache

import sympy as sp
from p8_affine_retuned import geometry

a, b, m2 = sp.symbols("a b m0_squared", positive=True)
omega, k = sp.symbols("omega k", real=True)


@cache
def determinant():
    momentum = sp.Matrix([omega, k, 0, 0])
    K = (momentum.T*momentum)[0]*sp.eye(4)-momentum*momentum.T+m2*sp.diag(a, b, b, b)
    expected = (omega**2+k**2+b*m2)**2*a*m2*(omega**2+b*k**2/a+b*m2)
    return {"K": K, "determinant": sp.factor(K.det()), "expected": expected,
            "full_four_component_determinant": sp.factor(K.det()-expected),
            "clock_three_polarization_determinant": sp.factor(K.det().subs({a: 1, b: 1})
                                                               -m2*(omega**2+k**2+m2)**3)}


@cache
def longitudinal():
    t, v, sigma = sp.symbols("temporal sigma_dot sigma", real=True)
    L = a*m2*t**2/2-b*m2*k**2*sigma**2/2+k**2*(v-t)**2/2
    solution = k**2*v/(a*m2+k**2)
    reduced = sp.factor(L.subs(t, solution))
    kinetic = a*m2*k**2/(a*m2+k**2)
    frequency = b*m2+b*k**2/a
    return {"L": L, "t": t, "v": v, "sigma": sigma, "solution": solution,
            "reduced": reduced, "kinetic": kinetic, "frequency": frequency,
            "full_temporal_Euler": sp.factor(sp.diff(L, t).subs(t, solution)),
            "physical_longitudinal_reduction": sp.factor(reduced-kinetic*(v**2-frequency*sigma**2)/2)}


@cache
def masses():
    data = geometry.update()
    return {"a": sp.factor(1/data["gamma_t"]), "b": sp.factor(1/data["gamma_s"])}


@cache
def checks():
    d, l = determinant(), longitudinal()
    return {name: (d if name in d else l)[name] for name in
            ("full_four_component_determinant", "clock_three_polarization_determinant",
             "full_temporal_Euler", "physical_longitudinal_reduction")}
