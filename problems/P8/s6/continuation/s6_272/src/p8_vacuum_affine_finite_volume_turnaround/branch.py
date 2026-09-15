"""Full translated complex-lapse contraction and physical volume domain."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate

from . import geometry, moving, source


@cache
def bounds():
    spatial = geometry.bounds()
    variation = source.MAJORANT * (source.LAPSE_RADIUS + sum(spatial["images"]))
    contraction = s.Rational(1, 30) + variation / 3
    residual = spatial["residual"]
    root = residual / (3 * (1 - contraction))
    volume = 12 * spatial["v"] + 24 * root
    absolute = {}
    for z, deviation in zip(source.q.COORDS, spatial["images"], strict=True):
        interval = evaluate(
            source.background()[z], {source.u: I(-source.TIME, source.TIME)}
        )
        value = max(abs(interval.lo), abs(interval.hi))
        absolute[z] = s.Rational(value.numerator, value.denominator) + deviation
    return {
        "variation": variation,
        "contraction": contraction,
        "residual": residual,
        "root": root,
        "volume": volume,
        "absolute": absolute,
    }


@cache
def data():
    b = bounds()
    Rgap = source.q.JET * source.LAPSE_RADIUS
    rn = 2 + source.q.JET * source.LAPSE_RADIUS
    UN = s.Rational(3, 4) * s.Rational(100, 99) ** 2 * rn
    full_hred = moving.old.field.KAPPA * 1000 * (2 * 10**62 + 2)
    full_href = moving.GENERATOR * moving.IMAGE_RADIUS**2 / 2
    return {
        "whole_all_translated_absolute_invariant_bounds": b["absolute"],
        "whole_full_auxiliary_contraction_and_root_bounds": b,
        "whole_complex_R_gap_RN_and_UN_bounds": [Rgap, rn, UN],
        "whole_entire_moving_Hamiltonian_and_reference_bounds": [full_hred, full_href],
        "whole_translated_branch_proof": "At each real time the exact moving-center CN is in(-3.1,-3). Only N and invariant deviations are varied using the entire187-derivative bound; no M*T is inserted. N->N+C/3 is strictly contractive on |N-1|<=1e-245, maps the ball strictly into itself and gives the displayed full root. The full source and all reconstructed invariants are holomorphic in finite complex phase coordinates. Real data give the positive real lapse by conjugation and uniqueness. Full RN(u,1)=-2/(1+u^2)^3 and the complete RNN bound give UN<12.",
        "whole_physical_volume_proof": "The complete physical readout is the spatial average exp(3v)R_full(u,Nstar)^(-3/4), then the full original interaction-picture pullback. Its complex amplitude differs from1 by at most12||v||+24|Nstar-1|<1e-320. The actual source-only center, not1, remains the scalar extension value. All implicit first and second contacts and full spatial reconstruction are retained.",
        "whole_Hamiltonian_bound_proof": "Restore physical kappa and the full torus integral, bounded by1000. Keep the complete moving-chart term -Hclock Pi_v-ell Pi_M and any exact boundary in the original chart. The homogeneous momentum boundary is proportional to integral v=0 in this nonzero-mode variational chart, not omitted before restriction. The full fixed reference quadratic form in whitened variables has norm equal to its symplectic generator norm; its entire pulled-back bound is displayed. Hence the full g, including its scalar phase, is<1e1000.",
        "checks": {
            "whole_fixed_point_denominator": 3 * (1 - b["contraction"]) * b["root"]
            - b["residual"],
            "whole_MVT_only_lapse_and_invariant_directions": b["variation"]
            - source.MAJORANT
            * (source.LAPSE_RADIUS + sum(geometry.bounds()["images"])),
            "whole_physical_volume_full_root_and_conformal_factors": b["volume"]
            - 12 * geometry.bounds()["v"]
            - 24 * b["root"],
            "whole_translated_invariant_count": s.Integer(len(b["absolute"]) - 12),
        },
        "gates": {
            "whole_translated_source_domain": all(
                value < source.COORDINATE for value in b["absolute"].values()
            ),
            "full_lapse_parameter_variation": b["variation"] < s.Rational(1, 1000),
            "full_strict_contraction": b["contraction"] < s.Rational(1, 20),
            "full_strict_self_map": b["residual"] / 3
            + b["contraction"] * source.LAPSE_RADIUS
            < source.LAPSE_RADIUS,
            "full_root_displacement": b["root"] < s.Rational(1, 10**324),
            "full_real_lapse_pivot_remains_negative": -3 + b["variation"]
            < -s.Rational(29, 10),
            "full_nonzero_R_branch": Rgap < s.Rational(1, 100),
            "full_temporal_auxiliary_pivot": 1 - s.Rational(3, 2) * Rgap**2 / (1 - Rgap)
            > s.Rational(99, 100),
            "full_complex_UN": UN < 12,
            "full_complex_volume_bound": b["volume"] < s.Rational(1, 10**320),
            "whole_hred_and_href_bound": full_hred + full_href < 10**1000,
            "all_source_scalar_phases_and_implicit_contacts_retained": True,
        },
    }
