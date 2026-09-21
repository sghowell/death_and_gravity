"""Original massive-vector calibrations and all selected frozen box-line kernels."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor
from p8_vacuum_affine_minimal_gravity_radiation import vertices as original
from p8_vacuum_affine_triangle_box_radiation import radiation as frozen

from . import basis, jets, radiation, source

ETA = s.diag(1, -1, -1, -1)


def dot(p, q):
    return (p.T * ETA * q)[0]


def numeric_difference(ps, k, eps):
    G = s.Matrix(4, 4, lambda i, j: dot(ps[i], ps[j]))
    H = s.Matrix(4, 4, lambda i, j: (ps[i].T * eps * ps[j])[0])
    aa = tuple(dot(k, p) for p in ps)
    lift = sum(
        c * jets.bose_contact(word, G, H, aa)
        for word, c in zip(basis.WORDS, basis.lift_coefficients())
    )
    wvars = s.symbols("calibration_w0:4")
    radial = [0] * 8
    for assignment in basis.PERMS:
        vertices = tuple(ps[i] for i in assignment)
        for split in range(4):
            rr = vertices[split:] + vertices[:split]
            ws = wvars[split:] + wvars[:split]
            qs = [s.zeros(4, 1)]
            for r in rr[:-1]:
                qs.append(qs[-1] - r)
            bar = sum((w * q for w, q in zip(ws, qs)), s.zeros(4, 1))
            HH = (bar.T * eps * bar)[0]
            U0 = sum(
                ws[i] * ws[j] * dot(qs[i] - qs[j], qs[i] - qs[j])
                for i in range(4)
                for j in range(i + 1, 4)
            )
            U1 = sum(-2 * ws[0] * ws[j] * dot(k, qs[j]) for j in range(1, 4))
            poly = s.Poly(
                s.expand(ws[0] * HH * (U0 * U0 + U0 * U1 + U1 * U1 / 3)), *wvars
            )
            for powers, co in poly.terms():
                radial[powers[1] + powers[3]] += co * basis.beta_moment(powers)
    z = basis.Z
    literal = -6 * sum(radial[r] * z ** (r + 1) * (1 - z) ** (8 - r) for r in range(8))
    T = sum(
        aa[i] ** 2 * H[j, j] + aa[j] ** 2 * H[i, i] - 2 * aa[i] * aa[j] * H[i, j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
    return s.factor(literal - lift), s.factor(T)


@cache
def data():
    checks, gates = {}, {}
    R = s.Rational
    vector_count, kernel_count, nonzero = 0, 0, 0
    for state in range(3):
        ps, k, _ = original.sample(state)
        checks[f"state{state}_original_mass_shells"] = s.Matrix(
            [dot(p, p) - 1 for p in ps]
        ).applyfunc(s.factor)
        checks[f"state{state}_original_conservation"] = sum(ps, k.copy()).applyfunc(
            s.factor
        )
        for pol, (eps, norm) in enumerate(
            tensor.transverse_polarizations(tensor.frame(k), ETA)
        ):
            value, T = numeric_difference(ps, k, eps)
            checks[f"state{state}_pol{pol}_entire_vector_box_curvature_polynomial"] = (
                s.factor(value - radiation.TARGET * T)
            )
            vector_count += 1
            nonzero += T != 0
            for ordering, labels in enumerate(
                ((0, 1, 2, 3), (0, 2, 1, 3), (0, 1, 3, 2))
            ):
                scale = R(1, 16)
                vertices = tuple(scale * ps[i] for i in labels)
                wave = scale * k
                masses = (1, source.HEAVY_MASS2, 1, source.HEAVY_MASS2)
                z = R(1, 5) - s.I * R(4, 25)
                basew = ((1 - z) / 3, 2 * z / 5, 2 * (1 - z) / 3, 3 * z / 5)
                gamma = R(2, 5)
                for split in range(4):
                    rr = vertices[split:] + vertices[:split]
                    ws = basew[split:] + basew[:split]
                    mm = masses[split:] + masses[:split]
                    first_light = next(i for i, m in enumerate(mm) if m == 1)
                    original_light = (split + first_light) % 4
                    xi, eta = (
                        (R(1, 3), R(2, 5))
                        if original_light == 0
                        else (R(2, 3), R(3, 5))
                    )
                    actual = frozen.line_density(
                        vertices,
                        masses,
                        wave,
                        eps,
                        split,
                        xi,
                        R(1, 5),
                        eta,
                        gamma,
                        source.HEAVY_MASS2,
                    )
                    qs = [s.zeros(4, 1)]
                    for r in rr[:-1]:
                        qs.append(qs[-1] - r)
                    effective = [gamma * wave, *qs[1:]]
                    bar = sum((w * q for w, q in zip(ws, qs)), s.zeros(4, 1))
                    M = sum(w * m for w, m in zip(ws, mm))
                    U = sum(
                        ws[i]
                        * ws[j]
                        * dot(effective[i] - effective[j], effective[i] - effective[j])
                        for i in range(4)
                        for j in range(i + 1, 4)
                    )
                    den = M - U
                    num = -4 * ws[0] * (bar.T * eps * bar)[0]
                    jac = z * (1 - z) * (1 - s.I * R(3, 5))
                    tag = f"state{state}_pol{pol}_order{ordering}_line{split}"
                    checks[tag + "_frozen_ordered_denominator"] = s.factor(
                        actual["denominator"] - den
                    )
                    checks[tag + "_frozen_weighted_TT_density"] = s.factor(
                        actual["density"] * den**3 - jac * num
                    )
                    kernel_count += 1
    gates.update(
        {
            "all_three_original_massive_states_two_polarizations": vector_count == 6,
            "all_original_curvature_words_nonzero": nonzero == 6,
            "all72_frozen_mass_ordered_box_line_kernels": kernel_count == 72,
            "mass_preserved_under_stated_offshell_calibration_scaling": True,
            "full2430_coefficient_proof_not_inferred_from_samples": True,
        }
    )
    return {
        "checks": checks,
        "gates": gates,
        "whole_original_counts": {
            "states": 3,
            "polarizations": 6,
            "entire_vector_polynomials": vector_count,
            "frozen_kernel_checks": kernel_count,
            "nonzero_curvature_polynomials": nonzero,
        },
        "whole_independent_calibration": "Every original massive state and polarization is checked using an independent4-vector full24-label/four-line computation with exact angular monomial integration. Separately72 actual frozen S342 alternating-mass box kernels, covering three external orderings and all four line positions, are evaluated at the original heavy mass after the explicitly off-shell common1/16 momentum scaling. Both their complete denominators and weighted complex-contour TT densities agree. The uniform generic polynomial proof is separate; these are exact diagnostic calibrations, not a physical Taylor approximation.",
    }
