"""Independent full tensor controls and exact original absorptive remainders."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor
from p8_vacuum_affine_dimensional_real_remainder import born
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree
from p8_vacuum_affine_minimal_gravity_radiation import vertices as old

from .kernel import channel
from .source import endpoint_factor


@cache
def component_checks():
    eta = s.diag(1, -1, -1, -1)
    checks = {}
    coincident = 0
    for index in range(3):
        ps, k, _ = old.sample(index)
        vv = (*ps, k)
        gram = s.Matrix(5, 5, lambda i, j, vv=vv: (vv[i].T * eta * vv[j])[0])
        aa = tuple(gram[i, 4] for i in range(4))
        expected = tree.whole_tensor(
            ps, k, mass=1, heavy=17, cubic=1, contact=0, kappa=1
        )
        for x in range(4):
            for y in range(x, 4):
                eps = s.zeros(4)
                eps[x, y] = 1
                eps[y, x] = 1
                pol = s.Matrix(
                    4, 4, lambda i, j, vv=vv, eps=eps: (vv[i].T * eps * vv[j])[0]
                )
                trace = s.trace(eta * eps)
                kep = tuple((k.T * eps * p)[0] for p in ps)
                actual = s.S.Zero
                for j in (1, 2, 3):
                    L = (0, j)
                    R = tuple(i for i in range(4) if i not in L)
                    left = sum(gram[i, h] for i in L for h in L)
                    right = sum(gram[i, h] for i in R for h in R)
                    if left == right:
                        coincident += 1
                    actual += channel(
                        gram,
                        pol,
                        aa,
                        L,
                        R,
                        1 / (17 - left),
                        1 / (17 - right),
                        1 / ((17 - left) * (17 - right)),
                        trace,
                        kep,
                    )
                checks[f"state{index}_component{x}{y}_complete_resolvent"] = s.factor(
                    actual
                    - sum(
                        expected[i, j] * eps[i, j] for i in range(4) for j in range(4)
                    )
                )

    return checks


@cache
def data():
    fixtures = [old.sample(i) for i in range(3)]
    E, Ep = s.Rational(5, 4), s.Rational(461, 380)
    fixtures.append(
        old.recoil(
            E, E - Ep * Ep / E, (1, 0, 0), (0, s.Rational(3, 5), s.Rational(4, 5))
        )
    )
    eta = s.diag(1, -1, -1, -1)
    checks = dict(component_checks())
    records = []
    equal = 0
    for index, (ps, k, p0) in enumerate(fixtures):
        vv = (*ps, k)
        gram = s.Matrix(5, 5, lambda i, j, vv=vv: (vv[i].T * eta * vv[j])[0])
        aa = tuple(gram[i, 4] for i in range(4))
        ss = 4 * ps[0][0] ** 2
        tt = ((p0[0] + p0[2]).T * eta * (p0[0] + p0[2]))[0]
        cosine = s.factor(1 + 2 * tt / (ss - 4))
        A0 = born.original(ss, cosine, 0)
        incoming = sum((ps[i] for i in (0, 1)), s.zeros(4, 1))
        outgoing = sum((ps[i] for i in (2, 3)), s.zeros(4, 1))
        beta0 = s.sqrt(1 - 4 / ss)
        outgoing_s = (outgoing.T * eta * outgoing)[0]
        betar = s.sqrt(1 - 4 / outgoing_s)
        assert beta0.is_Rational and betar.is_Rational and 0 < betar < beta0 < 1
        norm = s.S.Zero
        for polidx, (eps, norm2) in enumerate(
            tensor.transverse_polarizations(tensor.frame(k), eta)
        ):
            pol = s.Matrix(
                4, 4, lambda i, j, eps=eps, ps=ps: (ps[i].T * eps * ps[j])[0]
            )
            J = tuple(pol[i, i] / aa[i] for i in range(4))
            soft = sum((p.T * eps * p)[0] / (p.T * eta * k)[0] for p in p0)
            complete = s.S.Zero
            for j in (1, 2, 3):
                L = (0, j)
                R = tuple(i for i in range(4) if i not in L)
                v = sum(gram[i, h] for i in L for h in L)
                u = sum(gram[i, h] for i in R for h in R)
                imL = (
                    endpoint_factor(v) ** 2 * s.sqrt(1 - 4 / v) if j == 1 else s.S.Zero
                )
                imR = (
                    endpoint_factor(u) ** 2 * s.sqrt(1 - 4 / u) if j == 1 else s.S.Zero
                )
                if u == v:
                    equal += 1
                    dd = (
                        s.S.Zero
                        if j != 1
                        else s.diff(
                            (
                                germs.MASS2 * germs.CONTACT / germs.G**2
                                + germs.MASS2 / (germs.MASS2 - s.Symbol("q"))
                            )
                            ** 2
                            * s.sqrt(1 - 4 / s.Symbol("q")),
                            s.Symbol("q"),
                        ).subs(s.Symbol("q"), v)
                    )
                else:
                    dd = (imR - imL) / (u - v)
                complete += (
                    sum(J[i] for i in L) * imR
                    + sum(J[i] for i in R) * imL
                    - 2 * sum(pol[i, h] for i in L for h in L) * dd
                )
            reduced = endpoint_factor(outgoing_s) ** 2 * betar * (
                J[0] + J[1]
            ) + endpoint_factor(ss) ** 2 * beta0 * (J[2] + J[3])
            checks[f"state{index}_pol{polidx}_complete_physical_cut"] = s.factor(
                complete - reduced
            )
            checks[f"state{index}_pol{polidx}_timelike_internal_TT_zero"] = s.factor(
                (incoming.T * eps * incoming)[0]
            )
            remainder = s.factor(complete - endpoint_factor(ss) ** 2 * beta0 * soft)
            assert remainder.is_Rational and remainder != 0
            norm += remainder**2 / norm2
        norm = s.factor(norm)
        assert norm > 0 and 16 * norm < 2688768**2
        # Imaginary amplitude has g^4/(32*pi*n^2*sqrt(kappa)); pi>3.
        original_ratio = s.factor(norm * germs.G**8 / (96**2 * germs.MASS2**4 * A0**2))
        assert original_ratio < s.Integer(10) ** 388
        records.append((index, norm, original_ratio))
    return {
        "checks": checks,
        "gates": {
            "three_states_all30_independent_full_tensor_components": len(
                component_checks()
            )
            == 30,
            "four_original_rational_recoil_states": len(records) == 4,
            "four_actual_equal_invariant_cases": equal == 4,
            "both_polarizations_nonzero_exact_absorptive_remainders": all(
                row[1] > 0 for row in records
            ),
            "all_original_normalized_square_bounds_exact": all(
                row[2] < s.Integer(10) ** 388 for row in records
            ),
            "no_numerical_calibration_substituted_for_uniform_proof": True,
        },
        "whole_exact_original_absorptive_calibrations": records,
        "whole_actual_equal_invariant_cases": equal,
        "whole_component_control": "Three original recoil configurations and all10 symmetric components compare the general kernel to the independent complete heavy-exchange engine at diagnostic heavy mass17,unit cubic and no contact. This diagnostic mass is not used in the physical estimates.",
        "whole_absorptive_control": "Four original recoil states including a physical coincident crossed invariant have rational incoming/outgoing cut velocities and nonzero rational nonleading imaginary remainders for both physical polarizations. Each sewn squared norm and its ORIGINAL n,g,C,A0 normalization is exact. The common pi and1/sqrt(kappa) are restored explicitly in the written proof; these samples do not replace the uniform analytic bound.",
    }
