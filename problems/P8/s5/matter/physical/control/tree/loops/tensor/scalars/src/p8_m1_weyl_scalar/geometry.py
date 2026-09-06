"""Physical-metric scalar Weyl geometry, without background field equations."""

import itertools
from functools import cache

import sympy as sp


@cache
def scalar_contraction_checks():
    """Contract the full linearized Riemann tensor for arbitrary scalar jets.

    By spatial rotation the dependence may be along z. The perturbation is
    h00=2A, h0z=-B_z, hij=-2*zeta*delta_ij-2*E_ij, on conformal Minkowski
    signature +---. No leading scalar equation or constraint is imposed.
    """
    signs = (1, -1, -1, -1)
    jets = sp.symbols("A_tt A_tz A_zz Z_tt Z_tz Z_zz B_ttz B_tzz B_zzz E_ttzz E_tzzz E_zzzz")
    pairs = ((0, 0), (0, 3), (3, 3))
    maps = []
    for start in range(0, 12, 3):
        mapping = dict(zip(pairs, jets[start:start+3]))
        mapping[3, 0] = mapping[0, 3]
        maps.append(mapping)
    a, zeta, shift, shear = maps

    def hessian(i, j, c, d):
        if i == j == 0:
            return 2*a.get((c, d), 0)
        if (i, j) in ((0, 3), (3, 0)):
            return -shift.get((c, d), 0)
        if i == j and i > 0:
            return -2*zeta.get((c, d), 0)-(2*shear.get((c, d), 0) if i == 3 else 0)
        return 0

    riemann = {(a, b, c, d): sp.expand((hessian(a, d, c, b)+hessian(b, c, d, a)
                - hessian(a, c, d, b)-hessian(b, d, c, a))/2)
               for a, b, c, d in itertools.product(range(4), repeat=4)}
    ricci = {(b, d): sum(signs[a]*riemann[a, b, a, d] for a in range(4))
             for b, d in itertools.product(range(4), repeat=2)}
    scalar = sum(signs[a]*ricci[a, a] for a in range(4))
    riemann2 = sum(sp.prod(signs[i] for i in indices)*value**2 for indices, value in riemann.items())
    ricci2 = sum(signs[b]*signs[d]*value**2 for (b, d), value in ricci.items())
    weyl2 = sp.factor(riemann2-2*ricci2+scalar**2/3)
    lensing_zz = jets[2]-jets[5]+jets[7]-jets[9]
    return {"full_four_index_scalar_Weyl": sp.expand(weyl2-sp.Rational(4, 3)*lensing_zz**2),
            "tracefree_Hessian_contraction": sp.Rational(1, 4)*(sp.Rational(2, 3)) * 8-sp.Rational(4, 3)}


def gauge_and_ADM_checks():
    a, zeta, shift_dot, shear_ddot, ht, tdot = sp.symbols("A zeta Bdot Eddot HcT Tdot")
    lensing = a-zeta+shift_dot-shear_ddot
    # B-E' changes by T, while A changes by -T'-Hc*T and zeta by -Hc*T.
    transformed = lensing.subs({a: a-tdot-ht, zeta: zeta-ht,
                               shift_dot: shift_dot+tdot}, simultaneous=True)
    t = sp.Symbol("t", real=True)
    scale, b = sp.Function("a")(t), sp.Function("b")(t)
    hubble = sp.diff(scale, t)/scale
    return {"Bardeen_time_gauge_invariance": sp.expand(transformed-lensing),
            "conformal_shift_to_cosmic": sp.simplify(scale*sp.diff(b/scale, t)-(sp.diff(b, t)-hubble*b))}


def controls():
    n, zeta, bdot, h, b, delta = sp.symbols("n zeta bdot H b delta")
    physical = (1-delta)*n-zeta+bdot-h*b
    return {"conformal_scalar_metric_has_no_Weyl": sp.expand((n-zeta).subs(zeta, n)),
            "static_lapse_unit_spatial_Hessian": sp.Rational(4, 3),
            "omitting_auxiliary_metric_map_changes_lensing": sp.expand((n-zeta+bdot-h*b)-physical),
            "omitting_shift_derivative_changes_lensing": -bdot}
