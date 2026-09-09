"""Exact center bridge to the already quantized scalar density phase."""

from functools import cache

import sympy as sp
from p8_proca_global_hadamard import canonical as prior
from p8_proca_global_support import model, phase

from . import center


@cache
def data():
    c = center.data()
    k = next(iter(c["center_packet_map"].free_symbols))
    old = model.old
    actual = model.substitution()
    E = prior.data()["packet_to_density"].subs(actual, simultaneous=True)
    E = E.subs({old.u: 0, prior.k: k}, simultaneous=True).applyfunc(sp.factor)
    ell = sp.Rational(1, 10)
    # Fixed spatial Q=(v,p_new,chi,Pchi_new), old density
    # Z=(v,chi,pi_v,pi_chi). Both shifts come from 3 P0 v chi.
    D = sp.Matrix([[1, 0, 0, 0], [0, 3 * ell, 1, 0], [0, 1, 0, 0], [3 * ell, 0, 0, 1]])
    C = (D * E).applyfunc(sp.factor)
    Omega = phase.data()["constant_symplectic_form"]
    OmegaQ = sp.diag(sp.Matrix([[0, 1], [-1, 0]]), sp.Matrix([[0, 1], [-1, 0]]))
    v, p, chi, pc = c["fields"][:4]
    lapse = (
        model.margin.scalar()["regular_lapse"]
        .subs(
            {
                old.v: v,
                old.matter: chi,
                old.shift: -(p - 3 * ell * chi) / (2 * k * k),
                old.pm: pc - 3 * ell * v,
                old.pi: 0,
                old.sigma: 0,
                old.q: k * k,
            },
            simultaneous=True,
        )
        .subs(actual, simultaneous=True)
        .subs(old.u, 0)
    )
    no_matter_shift = D.copy()
    no_matter_shift[3, 0] = 0
    no_metric_shift = D.copy()
    no_metric_shift[1, 1] = 0
    return {
        "old_density_to_fixed_spatial_phase": D,
        "actual_parent_packet_to_fixed_spatial_phase": C,
        "fixed_spatial_symplectic_form": OmegaQ,
        "missing_matter_boundary_shift_matrix": (no_matter_shift * E - C).applyfunc(
            sp.factor
        ),
        "missing_metric_boundary_shift_matrix": (no_metric_shift * E - C).applyfunc(
            sp.factor
        ),
        "checks": {
            "both_canonical_boundary_shifts_preserve_density_symplectic_form": (
                D * Omega * D.T - OmegaQ
            ).applyfunc(sp.factor),
            "full_parent_packet_map_is_independent_center_map": (
                C - c["center_packet_map"]
            ).applyfunc(sp.factor),
            "fixed_spatial_packet_map_has_inverse_k_symplectic_scale": (
                C * Omega * C.T - OmegaQ / k
            ).applyfunc(sp.factor),
            "actual_parent_linear_lapse_is_full_invariant_linear_lapse": sp.factor(
                lapse - c["n1"]
            ),
        },
    }
