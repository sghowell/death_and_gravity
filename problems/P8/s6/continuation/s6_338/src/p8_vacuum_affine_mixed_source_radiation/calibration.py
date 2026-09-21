"""Independent original-source component and both-polarization calibrations."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor
from p8_vacuum_affine_dimensional_real_remainder import born
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree
from p8_vacuum_affine_minimal_gravity_radiation import vertices as old

from . import radiation


@cache
def data():
    checks = {}
    eta = s.diag(1, -1, -1, -1)
    cases = []
    nn = germs.MASS2
    gg = germs.G
    kk = germs.KAPPA
    Q = 8 * (4 - nn) / kk  # common16pi^2 left out for exact rational comparisons
    for index in range(3):
        ps, k, born_ps = old.sample(index)
        vectors = (*ps, k)
        gram = s.Matrix(
            5, 5, lambda i, j, vectors=vectors: (vectors[i].T * eta * vectors[j])[0]
        )
        ak = tuple(gram[i, 4] for i in range(4))
        unit = tree.whole_tensor(
            ps,
            k,
            mass=1,
            heavy=nn,
            cubic=1,
            contact=-3 / (nn - s.Rational(4, 3)),
            kappa=1,
        )
        ss = 4 * ps[0][0] ** 2
        tt = (born_ps[0] + born_ps[2]).T * eta * (born_ps[0] + born_ps[2])
        z = s.factor(1 + 2 * tt[0] / (ss - 4))
        channels = (ss, tt[0], 4 - ss - tt[0])
        center = gg**2 * (
            sum(1 / (nn - v) for v in channels) - 3 / (nn - s.Rational(4, 3))
        )
        Am = born.forward.matter_born(ss, z, nn, gg)
        A0 = born.original(ss, z, 0)
        assert 0 < center < Am < A0
        norm = s.S.Zero
        for polidx, (eps, norm2) in enumerate(
            tensor.transverse_polarizations(tensor.frame(k), eta)
        ):
            pol = s.Matrix(
                5,
                5,
                lambda i, j, vectors=vectors, eps=eps: (
                    vectors[i].T * eps * vectors[j]
                )[0],
            )
            actual = s.S.Zero
            for j in (1, 2, 3):
                L = (0, j)
                R = tuple(i for i in range(4) if i not in L)
                value, _, _, _ = radiation.one_channel(L, R, gram, pol, ak, nn)
                actual += value
            quartic = sum(pol[i, i] / ak[i] for i in range(4))
            actual -= (6 + 6 * (4 - nn) / (nn - s.Rational(4, 3))) * quartic
            literal = sum(unit[i, j] * eps[i, j] for i in range(4) for j in range(4))
            checks[f"state{index}_pol{polidx}_literal42_plus_OS4_equals_complete26"] = (
                s.factor(actual - 2 * (4 - nn) * literal)
            )
            soft = sum((p.T * eps * p)[0] / (p.T * eta * k)[0] for p in born_ps)
            remainder = s.factor(Q * (gg**2 * literal - center * soft))
            assert remainder != 0
            norm += remainder**2 / norm2
        norm = s.factor(norm)
        assert not norm.has(s.Float)
        assert norm / (144**2 * A0**2) < s.Rational(1, 10**598) ** 2
        assert abs(Q * center) / (144 * A0) < s.Rational(1, 10**603)
        cases.append(
            {
                "index": index,
                "s": ss,
                "z": z,
                "centered_to_matter_ratio": s.factor(center / Am),
                "remainder_norm_squared_over_conservative_loop_Born": s.factor(
                    norm / (144**2 * A0**2)
                ),
            }
        )
    return {
        "checks": checks,
        "gates": {
            "three_actual_original_recoil_configurations": len(cases) == 3,
            "both_physical_TT_polarizations_per_configuration": len(checks) == 6,
            "source_remainder_nonzero_not_only_a_leading_soft_extension": True,
            "exact_no_float_norms_and_original_bound_margins": True,
            "independent_frozen_complete26_component_engine_used": True,
            "finite_calibration_not_substitute_for_uniform_proof": True,
        },
        "whole_original_recoil_calibrations": cases,
        "whole_normalization": "The exact rational comparison omits the common1/sqrt(kappa) and uses16pi^2>144. The same original n,g,kappa and complete positive A0 are used. Both physical polarizations are sewn with their actual tensor norms. The42 derivative-source graphs plus fixed OS4 contact are compared to the independent frozen full26 matter component engine.",
        "whole_uniform_boundary": "Three finite states verify independent formulas and numerical margins, not a uniform theorem. The analytic bound uses convexity, the original canonical current-change and matter-remainder proofs, and actual original source inequalities.",
    }
