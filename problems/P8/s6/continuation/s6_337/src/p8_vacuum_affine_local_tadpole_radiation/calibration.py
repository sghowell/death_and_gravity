"""Exact original recoil calibrations of the known local-loop radiation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor
from p8_vacuum_affine_minimal_gravity_radiation import vertices as old

from . import bounds, radiation, source, vertices

ETA = s.diag(1, -1, -1, -1)


@cache
def data():
    checks, gates, rows = {}, {}, []
    coefficients = radiation.local_coefficients()
    for index in range(3):
        ps, k, born_ps = old.sample(index)
        vectors = (*ps, k)
        gram = s.Matrix(
            5, 5, lambda i, j, vectors=vectors: (vectors[i].T * ETA * vectors[j])[0]
        )
        born_gram = s.Matrix(
            4, 4, lambda i, j, born_ps=born_ps: (born_ps[i].T * ETA * born_ps[j])[0]
        )
        four = -s.Add(
            *(
                co * vertices.quartic_vertex(name, born_gram)
                for name, co in coefficients.items()
            )
        )
        ss = 4 * ps[0][0] ** 2
        tt = (born_ps[0] + born_ps[2]).T * ETA * (born_ps[0] + born_ps[2])
        z = s.factor(1 + 2 * tt[0] / (ss - 4))
        A0 = bounds.born.original(ss, z, 0)
        checks[f"state{index}_original_mass_shells"] = s.Matrix(
            [gram[i, i] - 1 for i in range(4)]
        )
        checks[f"state{index}_emitted_null"] = gram[4, 4]
        checks[f"state{index}_momentum_conservation"] = sum(ps, k.copy()).applyfunc(
            s.factor
        )
        checks[f"state{index}_whole_known_fourpoint"] = s.factor(
            four + radiation.four_polynomial(ss, tt[0])
        )
        norm = s.S.Zero
        for polidx, (eps, norm2) in enumerate(
            tensor.transverse_polarizations(tensor.frame(k), ETA)
        ):
            pol = s.Matrix(
                5,
                5,
                lambda i, j, vectors=vectors, eps=eps: (
                    vectors[i].T * eps * vectors[j]
                )[0],
            )
            full = radiation.full_tt_without_common_loop_and_metric_factor(gram, pol)
            soft = s.Add(*((p.T * eps * p)[0] / (p.T * ETA * k)[0] for p in born_ps))
            remainder = s.factor(full - four * soft)
            checks[f"state{index}_polarization{polidx}_transverse"] = eps * k
            checks[f"state{index}_polarization{polidx}_trace"] = s.trace(ETA * eps)
            gates[f"state{index}_polarization{polidx}_nonzero_remainder"] = bool(
                remainder != 0
            )
            norm += remainder**2 / norm2
        norm = s.factor(norm)
        gates[f"state{index}_exact_original_Born_lower_bound"] = bool(
            A0 > s.Rational(1, 10**600)
        )
        gates[f"state{index}_whole_remainder_budget"] = bool(
            norm < (2 * bounds.UNSCALED_BUDGET * bounds.coefficient_budget()) ** 2
        )
        gates[f"state{index}_original_normalized_remainder_bound"] = bool(
            norm / (144**2 * A0**2) < bounds.REMAINDER_RATIO**2
        )
        gates[f"state{index}_original_fourpoint_relative_bound"] = bool(
            abs(four) / (144 * A0) < bounds.HARD_RATIO
        )
        rows.append(
            {
                "state": index,
                "omega": k[0],
                "s": ss,
                "z": z,
                "both_polarizations_nonzero": True,
                "original_normalized_remainder_upper": bounds.REMAINDER_RATIO
                / s.sqrt(source.KAPPA),
                "original_local_fourpoint_relative_upper": bounds.HARD_RATIO,
            }
        )
    gates["three_original_source_states_both_polarizations"] = len(rows) == 3
    gates["sample_calibration_not_substitute_for_uniform_written_bound"] = True
    return {
        "checks": checks,
        "gates": gates,
        "whole_original_recoil_calibrations": rows,
        "whole_calibration_scope": "Exact original n,g,C,kappa,lambda and original full positive Born. Both actual physical polarizations, nonzero local-loop remainder and the known S239 fourpoint term. The omitted common16pi^2 and sqrt(kappa) are restored in the stated bounds;16pi^2>144 is only a conservative comparison.",
    }
