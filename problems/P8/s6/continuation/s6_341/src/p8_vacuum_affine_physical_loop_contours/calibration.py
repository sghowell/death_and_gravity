"""Exact finite controls at the original source, not uniform-proof substitutes."""

import itertools
from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree
from p8_vacuum_affine_minimal_gravity_radiation import vertices as actual
from p8_vacuum_affine_minimal_gravity_radiation import vertices as original

from .insertions import avg, dot, eta


@cache
def data():
    checks = {}
    R = s.Rational
    loop = s.Matrix(s.symbols("loop0:4", real=True))
    negative_cases = 0
    gap_cases = 0
    for vv, aa, bb in (
        (R(45, 8), 0, 0),
        (16, 2, 4),
        (R(25, 4), 1, -3),
        (-12, 2, 4),
        (0, 0, -3),
    ):
        for xx in (R(0), R(1, 32), R(1, 5), R(1, 2), R(4, 5), R(31, 32), R(1)):
            for tt in (R(0), R(1, 32 * 10**6), R(1, 10**6), R(1, 2), R(1)):
                for hh in (R(1, 4), R(1)):
                    XX = xx + s.I * hh * xx * (1 - xx) * (1 - 2 * xx) if vv > 4 else xx
                    ZZ = tt - s.I * hh * tt * (1 - tt)
                    LL = 1 - vv * XX * (1 - XX)
                    DD = (
                        (1 - ZZ) ** 2 * LL
                        + 10**6 * ZZ
                        + (1 - aa) * ZZ * (1 - ZZ)
                        - bb * ZZ**2
                    )
                    assert s.im(DD) <= -(10**6 - 20) * hh * tt * (1 - tt)
                    negative_cases += 1
                    if hh == 1:
                        assert (
                            s.expand(DD * s.conjugate(DD))
                            >= ((R(1, 32) + 10**6 * tt) / 64) ** 2
                        )
                        gap_cases += 1

    # Original recoil domain: shifting any one hard leg by k restores four-point
    # conservation; interpolation moves adjacent external momenta by gamma k.
    actual_virtualities = 0
    actual_pairs = 0
    for index in range(3):
        ps, k, _ = actual.sample(index)
        for leg in range(4):
            shifted = [p + k if j == leg else p for j, p in enumerate(ps)]
            checks[f"original_recoil{index}_leg{leg}_hard_conservation"] = sum(
                shifted, s.zeros(4, 1)
            )
            for p in shifted:
                a = actual.dot(p, p)
                assert 0 <= a <= 2
                actual_virtualities += 1
            for i, j in itertools.combinations(range(4), 2):
                vv = actual.dot(shifted[i] + shifted[j], shifted[i] + shifted[j])
                assert -12 <= vv <= 0 or R(45, 8) <= vv <= 16
                actual_pairs += 1
    cases = 0
    ttcases = 0
    endpointgates = 0
    for index in range(3):
        ps, k, _ = original.sample(index)
        for N, vertices, masses in (
            (3, [ps[0] + ps[1], ps[2], ps[3]], [1, 1, germs.MASS2]),
            (4, list(ps), [1, germs.MASS2, 1, germs.MASS2]),
            (4, list(ps), [germs.MASS2, 1, germs.MASS2, 1]),
        ):
            for split in range(N):
                r = vertices[split:] + vertices[:split]
                m = masses[split:] + masses[:split]
                qs = [s.zeros(4, 1)]
                for pp in r[:-1]:
                    qs.append(qs[-1] - pp)
                assert qs[-1] - r[-1] == k
                assert all(q.dot(q) <= 144 for q in qs)
                for gam in (R(0), R(1, 2), R(1)):
                    eff = [q + gam * k if j == 0 else q for j, q in enumerate(qs)]
                    rs = [eff[j] - eff[(j + 1) % N] for j in range(N)]
                    assert all(
                        rs[j]
                        == r[j]
                        + (
                            gam * k
                            if j == 0
                            else (1 - gam) * k
                            if j == N - 1
                            else s.zeros(4, 1)
                        )
                        for j in range(N)
                    )
                    assert sum(rs, s.zeros(4, 1)) == s.zeros(4, 1)
                    lights = [j for j, vv in enumerate(m) if vv == 1]
                    heavies = [j for j, vv in enumerate(m) if vv == germs.MASS2]
                    lv = dot(
                        eff[lights[0]] - eff[lights[1]], eff[lights[0]] - eff[lights[1]]
                    )
                    assert -12 <= lv <= 0 or R(45, 8) <= lv <= 16
                    if N == 4:
                        hv = dot(
                            eff[heavies[0]] - eff[heavies[1]],
                            eff[heavies[0]] - eff[heavies[1]],
                        )
                        assert -12 <= hv <= 16
                        assert all(0 <= dot(p, p) <= 2 for p in rs)
                        endpointgates += 6
                    else:
                        assert all(
                            0 <= dot(eff[a] - eff[b], eff[a] - eff[b]) <= 2
                            for a in lights
                            for b in heavies
                        )
                        endpointgates += 3
                    cases += 1
                # Independent literal source vertex after an actual complex
                # simplex shift, both physical transverse polarizations.
                X = R(1, 5) + s.I * R(12, 125)
                Z = R(1, 7) - s.I * R(6, 49)
                if N == 3:
                    weights = [None] * N
                    weights[lights[0]] = (1 - Z) * X
                    weights[lights[1]] = (1 - Z) * (1 - X)
                    weights[heavies[0]] = Z
                else:
                    weights = [None] * N
                    weights[lights[0]] = (1 - Z) * X
                    weights[lights[1]] = (1 - Z) * (1 - X)
                    weights[heavies[0]] = Z * R(2, 5)
                    weights[heavies[1]] = Z * R(3, 5)
                bar = sum((u * q for u, q in zip(weights, qs)), s.zeros(4, 1))
                yy = weights[0] * R(2, 3)
                RR = qs[0] - bar
                literal = tree.tensor(
                    loop + RR - yy * k, loop + RR + (1 - yy) * k, m[0], eta
                )
                for polidx, (ee, norm) in enumerate(
                    tensor.transverse_polarizations(tensor.frame(k), eta)
                ):
                    check = (
                        avg(
                            sum(
                                ee[a, b] * literal[a, b]
                                for a in range(4)
                                for b in range(4)
                            ),
                            loop,
                        )
                        - 2 * (RR.T * ee * RR)[0]
                    )
                    checks[
                        f"original{index}_N{N}_massstart{masses[0] == 1}_split{split}_pol{polidx}_literal_vertex"
                    ] = s.factor(check)
                    ttcases += 1

    return {
        "checks": {
            key: s.factor(value)
            if not isinstance(value, s.MatrixBase)
            else value.applyfunc(s.factor)
            for key, value in checks.items()
        },
        "gates": {
            "whole350_homotopy_controls": negative_cases == 350,
            "whole175_gap_controls": gap_cases == 175,
            "whole48_original_singleton_virtualities": actual_virtualities == 48,
            "whole72_original_pair_invariants": actual_pairs == 72,
            "whole99_cyclic_interpolated_routes": cases == 99,
            "whole513_actual_routing_domain_gates": endpointgates == 513,
            "whole66_actual_literal_TT_insertions": ttcases == 66,
            "finite_controls_do_not_replace_written_uniform_proof": True,
        },
        "whole_exact_calibration_counts": {
            "homotopy": negative_cases,
            "gap": gap_cases,
            "external_virtualities": actual_virtualities,
            "external_pairs": actual_pairs,
            "cyclic_routings": cases,
            "routing_domain_gates": endpointgates,
            "actual_TT_insertions": ttcases,
        },
        "whole_independent_branch_diagnostic_boundary": "A separate private60digit diagnostic compares the full primitive with the exact S235 cuts at both original mass and n=10^6, and direct2D rational contour quadrature with the primitive at n=10^6. Those numerical residuals are not exact certificates, rigorous errors or uniform estimates; the analytic contour proof supplies the bounds.",
    }
