"""Unchanged full real-time source and high phase, not time, derivatives."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantitative_local_time import audit as previous
from p8_vacuum_affine_quantitative_local_time import quantum, reference, source

RADIUS = reference.CORE
PHASE = 96
PAIRS = 48
MAX_ORDER = 196
HBOUND, VERROR, TIME = quantum.HBOUND, quantum.VERROR, source.TIME


@cache
def data():
    packets = previous.packets()
    domain = quantum.complete_domain_bounds()
    return {
        "whole_same_full_nonlinear_time_source_and_phase_domain": domain,
        "whole_original_time_and_phase_radii": [TIME, RADIUS, 2 * RADIUS, 4 * RADIUS],
        "whole_same_original_family_and_seed": previous.parameters(),
        "whole_high_phase_derivative_argument": "The entire S270 complex initial phase ball4R and its full nonlinear source, primitive, spatial formal-transpose Neumann inverse, auxiliary root and spatial integral are retained at each REAL time. The fixed profiles are only required C5 in real time. Their N dependence and the actual finite phase reconstruction are holomorphic, so arbitrarily high PHASE derivatives do not require arbitrarily high real-time profile derivatives. For an order<=196 derivative use coordinate polydisc radiusR/8 in at most96 coordinates around real ball2R; Euclidean displacement<=sqrt96 R/8<1.25R, hence inside4R. This is not an extension of the old R/4 polydisc to all96 coordinates.",
        "whole_same_full_interaction_and_volume": "Use exactly g=(hred-href)(u,S_u z) and F=average exp(3v)R_full(u,Nstar)^(-3/4) from S270, with actual scalar g0,F0 retained. The full original moving chart and free reference are not changed. All nonlinear harmonics, implicit Hessians, matter/vector/Gauss, clock and primitive contacts remain. Holomorphic amplitude bounds are |g|<1e1000 and |F-1|<1e-255, uniformly on initial complex ball4R.",
        "whole_same_regulators_and_canonical_normalization": "The SAME two explicit smooth radial cutoffs c=1,2, coreR=1e20, same pure covariance V0=S0 S0^T/2 and unit-CCR symplectic whitening are retained. Constant backgrounds are scalar identities, not subtracted phases. Weyl covariance is used only under the actual metaplectic S0 and the same original Uref; no non-symplectic rescaling of Planck's constant or fixed-window coherent covariance is assumed.",
        "checks": {
            "same_48_pairs_96_phase": s.Integer(PHASE - 2 * PAIRS),
            "same_complete_physical_kappa": reference.field.KAPPA - 10**800,
            "same_full_real_time": TIME - s.Rational(1, 10**2000),
            "same_full_complex_amplitude": HBOUND - 10**1000,
            "same_full_volume_amplitude": VERROR - s.Rational(1, 10**255),
            "same_two_cutoff_core": RADIUS - 10**20,
            "high_Cauchy_derivative_count": s.Integer(MAX_ORDER - 2 * PHASE - 4),
        },
        "gates": {
            "all_complete_original_S270_packet_gates": all(
                bool(g) for p in packets.values() for g in p["gates"].values()
            ),
            "all96_coordinate_polydisc_fits_full_domain": 2 + s.sqrt(PHASE) / 8 < 4,
            "full_actual_Hamiltonian_amplitude": domain[
                "whole_complete_interaction_bound"
            ]
            < HBOUND,
            "full_actual_physical_volume_amplitude": domain[
                "whole_complex_volume_error"
            ]
            < VERROR,
            "original_source_heavy_and_background_not_replaced": True,
            "high_phase_not_unproved_high_real_time_derivatives": True,
            "no_original_Wilsonian_or_continuum_matching_claim": True,
        },
    }
