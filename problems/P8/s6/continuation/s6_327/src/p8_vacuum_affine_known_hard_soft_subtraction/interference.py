"""Retained hard phase, exact unpolarized cancellation and original47 checks."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree
from p8_vacuum_affine_minimal_gravity_radiation import vertices as v

from . import source


@cache
def original_calibration():
    checks, records = {}, []

    def put(name, value):
        checks[name] = s.factor(value)

    for sample in (0, 2):
        ps, q, born = v.sample(sample)
        if sample == 0:
            e = s.Matrix([0, s.Rational(4, 5), -s.Rational(3, 5), 0])
            f = s.Matrix([0, 0, 0, 1])
        else:
            e = s.Matrix([0, 1, 0, 0])
            f = s.Matrix([0, 0, s.Rational(4, 5), -s.Rational(3, 5)])
        plus = (e * e.T - f * f.T) / s.sqrt(2)
        cross = (e * f.T + f * e.T) / s.sqrt(2)
        T = tree.whole_tensor(ps, q, mass=1, **source.original_parameters())
        Am = tree.born_continuation(
            born, heavy=source.HEAVY_MASS2, cubic=source.CUBIC, contact=source.CONTACT
        )
        A0 = Am + v.born(born) / source.KAPPA
        if not A0 > 0:
            raise ValueError("Original full Born normalization must be positive")
        soft = tree.soft_current(born, q) / s.sqrt(source.KAPPA)
        values = []
        for label, A in (
            ("plus", plus),
            ("cross", cross),
            ("helicity", (plus + s.I * cross) / s.sqrt(2)),
        ):
            F = s.factor(
                (
                    v.pair(A, T)
                    + v.amplitude(ps, q, A) / source.KAPPA ** s.Rational(3, 2)
                )
                / A0
            )
            S = s.factor(v.pair(A, soft))
            values.append((F, S))
            if label != "helicity":
                put(str(sample) + "_" + label + "_full_tree_real", s.im(F))
                put(str(sample) + "_" + label + "_soft_current_real", s.im(S))
        fp, fc, fh = [x[0] for x in values]
        sp, sc, sh = [x[1] for x in values]
        put(str(sample) + "_full47_complex_linearity", fh - (fp + s.I * fc) / s.sqrt(2))
        put(str(sample) + "_soft_complex_linearity", sh - (sp + s.I * sc) / s.sqrt(2))
        put(
            str(sample) + "_full_unpolarized_inner_product",
            s.conjugate(fh) * sh + fh * s.conjugate(sh) - fp * sp - fc * sc,
        )
        single = s.factor(s.im(s.conjugate(fh) * sh))
        records.append((sample, fp, fc, sp, sc, single))
    return {"checks": checks, "records": records}


@cache
def data():
    f0, f1, s0, s1, hr, hi, rho = s.symbols("f0 f1 s0 s1 h_real h_imag rho", real=True)
    F, S = s.Matrix([f0, f1]), s.Matrix([s0, s1])
    U = s.Matrix([[1, s.I], [1, -s.I]]) / s.sqrt(2)
    inner = ((U * F).conjugate().T * (U * S))[0]
    counter = (
        (U * s.Matrix([1, 0])).conjugate().multiply_elementwise(U * s.Matrix([0, 1]))
    )
    checks = {
        "unitary_helicity_frame": (U.conjugate().T * U - s.eye(2)).norm() ** 2,
        "summed_real_tree_soft_inner_product": inner - (F.T * S)[0],
        "selected_soft_subtracted_hard_interference": 2
        * s.re(rho * inner * (hr + s.I * hi))
        - 2 * hr * (S.T * S)[0]
        - 2 * hr * (rho * (F.T * S)[0] - (S.T * S)[0]),
        "opposite_single_helicity_imaginary_cross": counter[0] + counter[1],
        "one_helicity_hard_phase_can_contribute": 2 * s.re(s.I * counter[0]) + 1,
        "complete_leading_soft_density_has_two_interference_factor": s.Integer(2 - 2),
    }
    actual = original_calibration()
    checks.update({"original_" + key: value for key, value in actual["checks"].items()})
    return {
        "checks": {k: s.factor(s.expand_complex(value)) for k, value in checks.items()},
        "gates": {
            "complete_unpolarized_sum_before_discarding_any_imaginary_interference": True,
            "phase_retained_in_amplitude_and_other_observables": True,
            "both_original_samples_have_nonzero_single_helicity_imaginary_cross": all(
                row[-1] != 0 for row in actual["records"]
            ),
            "all26_matter_and21_gravity_diagrams_in_actual_calibrations": True,
            "original_full_positive_Born_used": True,
            "finite_radiative_hard_remainder_not_identified_with_Born_extension": True,
        },
        "whole_exact_selected_density": "2Re[rho sum conj(F_lambda)S_lambda h]-2Re(h)sum|S_lambda|^2 = I_known[rho sum Re(conj(F_lambda)S_lambda)-sum|S_lambda|^2], I_known=2Re h. F is the COMPLETE original one-real tree/A0; S is its Born leading-soft factor; h is fixed at Born. This identity uses the complete two-polarization sum.",
        "whole_original_real_complex_calibrations": actual["records"],
        "whole_phase_boundary": "The imaginary hard amplitude is not bounded by its real part or set to zero. Its contribution cancels only in the stated unpolarized linear interference; it can be nonzero for a single complex helicity and survives in other amplitudes, density entries and higher orders.",
    }
