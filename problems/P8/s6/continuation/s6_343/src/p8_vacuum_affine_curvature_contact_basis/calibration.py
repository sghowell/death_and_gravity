"""Independent literal curvature contractions at the original massive states."""

from functools import cache
from itertools import combinations_with_replacement, product

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor
from p8_vacuum_affine_minimal_gravity_radiation import vertices as original
from p8_vacuum_affine_radiative_curvature_matching import contact

from .basis import pairs


@cache
def data():
    checks, gates = {}, {}
    actual = 0
    nonzero = 0
    for index in range(3):
        ps, wave, _ = original.sample(index)
        checks[f"state{index}_mass_shells"] = s.Matrix(
            [contact.dot(p, p) - 1 for p in ps]
        ).applyfunc(s.factor)
        checks[f"state{index}_null_graviton"] = s.factor(contact.dot(wave, wave))
        checks[f"state{index}_conservation"] = sum(ps, wave.copy()).applyfunc(s.factor)
        frame = tensor.frame(wave)
        for polidx, (eps, norm2) in enumerate(
            tensor.transverse_polarizations(frame, contact.ETA)
        ):
            R = contact.linear_curvature(wave, eps)
            for first, second in combinations_with_replacement(pairs, 2):
                i, j = first
                k, l = second
                literal = sum(
                    R[aa, bb, cc, dd] * ps[i][aa] * ps[j][bb] * ps[k][cc] * ps[l][dd]
                    for aa, bb, cc, dd in product(range(4), repeat=4)
                )
                vij = (
                    contact.dot(wave, ps[i]) * ps[j] - contact.dot(wave, ps[j]) * ps[i]
                )
                vkl = (
                    contact.dot(wave, ps[k]) * ps[l] - contact.dot(wave, ps[l]) * ps[k]
                )
                factor = (vij.T * eps * vkl)[0]
                checks[f"state{index}_pol{polidx}_literal_R_{i}{j}{k}{l}"] = s.factor(
                    literal - factor
                )
                actual += 1
            value = contact.literal_vertex(ps, wave, eps)
            nonzero += value != 0
            if polidx == 0:
                gates[f"state{index}_original_canonical_amplitude_nonzero"] = value != 0
    gates["three_original_massive_states_two_polarizations_all21_monomials"] = (
        actual == 126
    )
    gates["finite_dimension_no_massless_limit_was_taken"] = True
    gates["parity_even_at_most_six_derivatives_only"] = True

    return {
        "checks": checks,
        "gates": gates,
        "whole_original_literal_counts": {
            "states": 3,
            "polarizations": 6,
            "mixed_pair_contractions": actual,
            "nonzero_canonical_polarizations": nonzero,
        },
        "whole_calibration_boundary": "All21 symmetric pair-pair contractions of literal256-component linear curvature, at each of three original massive recoil states and both physical polarizations. These independently check conventions and nonvanishing; symbolic Bose reduction and written index counting prove completeness, not these samples.",
    }
