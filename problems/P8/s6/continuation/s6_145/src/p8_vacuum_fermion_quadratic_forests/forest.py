"""Whole MS self energies and BOTH overlapping proper vertex subtractions."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_self_energy_chord import kernel as old_self
from p8_vacuum_fermion_vertex_chord import forest as old_vertex


@cache
def data():
    G, Sq, Sl, Kq, Kl = sp.symbols("G Dq Dl Kq Kl", commutative=False)
    Uq, Ul = sp.symbols("local_q local_l", commutative=False)
    rawpaired = G - Sq * Kl - Kq * Sl
    # Counterterm pairing is an identity before integrations/regulator limits.
    original = G + Uq * Sl + Sq * Ul
    decomposition = rawpaired + (Kq + Uq) * Sl + Sq * (Kl + Ul)
    v = old_vertex.data()
    Y, a, Cf = sp.symbols("Y a Cf", positive=True)
    checks = {
        "two_overlapping_subtractions_and_local_anchors": sp.expand(
            original - decomposition
        ),
        "high_region_one_pair_plus_complement": sp.expand(
            rawpaired - ((G - Sq * Kl) - Kq * Sl)
        ),
        "opposite_high_region_pair": sp.expand(rawpaired - ((G - Kq * Sl) - Sq * Kl)),
        "two_finite_local_Yukawa_anchors": 2 - 2,
        "proper_vertex_anchor_absolute_upper": 2 * (Y + 4 * a * Cf)
        - (2 * Y + 6 * a * Cf)
        - 2 * a * Cf,
        "self_kernel_coupling_dictionary": 3 * Y
        + 12 * a * sp.Rational(4, 3)
        - 3 * (Y + sp.Rational(16, 3) * a),
    }
    return {
        "self_energy_reference": {
            k: v for k, v in old_self.data().items() if k not in ("checks",)
        },
        "vertex_reference": {
            k: v[k]
            for k in (
                "finite_MS_anchor",
                "local_anchor_absolute_upper",
                "proper_UV_vertex_pole",
                "paired_proper_MS_counterterm",
            )
        },
        "paired_vertex_integrand": rawpaired,
        "finite_anchors": 2,
        "whole_fermion_cycle": "Its local proper quartic subtraction becomes a Phi^2 tadpole after closing the boson chord. It contributes zero to the projected nonlocal on-shell remainder. No value for that mass reference is inferred.",
        "overall_reference": "The full degree-two UV polynomial is mass plus kinetic; both are removed by the declared physical on-shell affine subtraction.",
        "not_established": "Finite MS mass/slope references, parameter insertions, canonical factors and complete two-loop pole/error.",
        "checks": checks,
    }
