"""First-variation jet controls for counterterms vanishing on the clock."""
from functools import cache

import sympy as sp


@cache
def jet_check(maximum_order):
    if type(maximum_order) is not int or maximum_order not in (1, 2, 3, 4):
        raise ValueError("Require a controlled metric-jet order from one to four")
    u = sp.Symbol("clock_time", real=True)
    D = sp.Symbol("spatial_dimension", real=True)
    alpha, beta, H = [sp.Function(name)(u) for name in ("alpha", "beta", "H")]
    jets = sp.symbols("N0:"+str(2*maximum_order+2), real=True)
    N = jets[0]
    A = [sp.Function("A"+str(j))(u) for j in range(maximum_order+1)]
    B = [sp.Function("B"+str(j))(u) for j in range(maximum_order+1)]
    Ga = sum(value*jets[j]/N**2 for j, value in enumerate(A))
    Gb = sum(value*jets[j]/N**3 for j, value in enumerate(B))+jets[1]**2/N**4
    density = -N*(N-1)*(alpha*Ga+beta*Gb)

    def weighted_derivative(expression):
        return sp.diff(expression, u)+sum(jets[j+1]*sp.diff(expression, jet) for j, jet in enumerate(jets[:-1]))+D*H*expression

    clock = {jet: sp.Integer(1 if j == 0 else 0) for j, jet in enumerate(jets)}
    energy = sp.Integer(0)
    for j in range(maximum_order+1):
        term = sp.diff(density, jets[j])
        for _ in range(j):
            term = -weighted_derivative(term)
        energy -= term.subs(clock)
    expected = alpha*Ga.subs(clock)+beta*Gb.subs(clock)
    a = sp.Symbol("physical_scale_factor", positive=True)
    # Arbitrary scale dependence of the geometric coefficient does not
    # contribute at first order when its scalar mass deviation is zero.
    spatial_density = -N*a**D*(N-1)*(alpha*(Ga+a**2)+beta*(Gb+a**-3))
    pressure = sp.diff(spatial_density, a).subs(clock)
    return {"full_metric_jet_lapse_variation_"+str(maximum_order): sp.simplify(energy-expected),
            "zero_first_spatial_variation_of_linear_mass_terms_"+str(maximum_order): sp.simplify(pressure)}


@cache
def checks():
    out = {}
    for order in range(1, 5):
        out.update(jet_check(order))
    return out
