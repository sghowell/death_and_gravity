"""Exact symmetry identities and conditional functional-branch control.

No algebraic potential Hessian is relabelled as an inverse of a differential
operator. Boundary/initial data, constrained spaces and the operator norm
must be supplied separately for an application of the conditional bound.
"""

from fractions import Fraction as Q
from functools import cache

import sympy as sp


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Q)):
        raise TypeError("Only exact int/Fraction bounds are accepted")
    return Q(value)


def residual_bound(inverse_bound, contraction, residual):
    """K*r/(1-eta), CONDITIONAL on a constrained functional inverse and ball.

    For E(h)=Lh+N(h)=r with N(0)=0, ||L^-1||<=K and ||L^-1(N(h)-N(k))||
<=eta||h-k|| with eta<1 in the specified ball. Existence additionally
requires K||r||<=(1-eta)*radius and the corresponding complete space.
"""
    inverse_bound, contraction, residual = map(rational, (inverse_bound, contraction, residual))
    if inverse_bound <= 0 or not 0 <= contraction < 1 or residual < 0:
        raise ValueError("Require K>0, 0<=eta<1, residual>=0")
    return inverse_bound*residual/(1-contraction)


def ball_gate(inverse_bound, contraction, residual, radius):
    radius = rational(radius)
    if radius <= 0:
        raise ValueError("Positive ball radius required")
    return residual_bound(inverse_bound, contraction, residual) <= radius


def closed_heavy_graph_loop_bound(heavy_degrees):
    """Each nonempty connected heavy component with no heavy external legs
    has E-V+1>=1 if every vertex has positive even heavy degree.
    This is the diagrammatic parity counterpart of exact classical Delta0.
    """
    degrees = tuple(heavy_degrees)
    if not degrees or any(type(degree) is not int or degree < 2 or degree % 2 for degree in degrees):
        raise ValueError("A nonempty list of positive even heavy degrees is required")
    return sum(degrees)//2-len(degrees)+1


@cache
def checks():
    t = sp.Symbol("t", real=True)
    h = sp.Function("h")(t)
    kinetic, mass, quartic, light = (sp.Function(name)(t) for name in
                                    ("kinetic", "mass_squared", "quartic", "light_action"))
    hp = sp.diff(h, t)
    lagrangian = light+kinetic*hp**2/2-mass*h**2/2+quartic*h**2*hp**2
    euler = sp.diff(lagrangian, h)-sp.diff(sp.diff(lagrangian, hp), t)
    zero = {h: 0, sp.diff(h, t): 0, sp.diff(h, t, 2): 0}
    ell, relative, source = sp.symbols("light relative source", real=True)
    # A finite-jet chain-rule audit, not a formalization of the functional
    # implicit-function or well-posedness theorem in the written proof.
    full = ell**3+source*ell+(1+ell**2)*relative**2+relative**4
    wrong_map = ell**2
    substituted = full.subs(relative, wrong_map)
    return {
        "rolling_even_action_zero_branch_off_shell": sp.expand(euler.subs(zero, simultaneous=True)),
        "complete_linear_operator_includes_kinetic_derivative": sp.expand(
            euler.subs(quartic, 0).doit()+kinetic*sp.diff(h, t, 2)
            +sp.diff(kinetic, t)*hp+mass*h),
        "stationary_light_action_derivative": sp.diff(full.subs(relative, 0), ell)
                                                -sp.diff(full, ell).subs(relative, 0),
        "source_preserved_on_exact_branch": sp.diff(full.subs(relative, 0), source)-ell,
        "nonstationary_substitution_chain_rule": sp.expand(sp.diff(substituted, ell)
            -sp.diff(full, ell).subs(relative, wrong_map)
            -sp.diff(full, relative).subs(relative, wrong_map)*sp.diff(wrong_map, ell)),
    }


def controls():
    t, lam, h = sp.symbols("t lambda h", real=True)
    symmetry_breaking = h*(h*h-lam)
    resonant = sp.sin(t)
    loop_test = sp.log(1+lam*lam)
    return {
        "positive_mass_Dirichlet_resonant_equation": sp.diff(resonant, t, 2)+resonant,
        "positive_mass_Dirichlet_left": resonant.subs(t, 0),
        "positive_mass_Dirichlet_right": resonant.subs(t, sp.pi),
        "positive_mass_Dirichlet_nonzero_solution": resonant.subs(t, sp.pi/2),
        "nonzero_heavy_initial_state": sp.cos(t).subs(t, 0),
        "branch_switch_loses_linear_inverse": sp.diff(symmetry_breaking, h).subs({h: 0, lam: 0}),
        "broken_symmetry_linear_source": sp.diff(h, h),
        "even_heavy_quadratic_coefficient_can_produce_light_loop_dependence": sp.diff(loop_test, lam).subs(lam, 1),
        "conditional_residual_fixture": residual_bound(Q(2), Q(1, 4), Q(3, 100)),
        "conditional_ball_fixture": ball_gate(Q(2), Q(1, 4), Q(3, 100), Q(1, 10)),
    }
