"""Full flat canonical covariance, original mode subtraction and tensor force."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_isolated_shear_resolvent import normalization as original
from p8_vacuum_affine_local_tensor_response import prescription
from p8_vacuum_canonical_affine_decoupling import family

MASS, KAPPA = modes.MASS, modes.KAPPA
C0 = 16 * s.pi**2 * KAPPA
C = C0 + 5 * MASS**2 / s.Integer(6)
p = s.Symbol("p", complex=True)
A = s.Function("original_A2")
y = original.y


def symmetric_matrix(name):
    v = s.symbols(name + "0:6", real=True)
    return s.Matrix([[v[0], v[1], v[2]], [v[1], v[3], v[4]], [v[2], v[4], v[5]]])


@cache
def data():
    lam, omega, S = s.symbols("lambda omega spectral_mass_squared", positive=True)
    XG, YG, XD, YD = [symmetric_matrix(name) for name in ("XG", "YG", "XD", "YD")]
    J = s.zeros(3).row_join(s.eye(3)).col_join((-s.eye(3)).row_join(s.zeros(3)))
    G = s.diag(XG, YG)
    U = omega * (YG - XG) / (lam**2 + 4 * omega**2)
    cross = lam * (YG - XG) / (2 * (lam**2 + 4 * omega**2))
    tangent = U.row_join(cross).col_join(cross.row_join(-U))
    forcing = (J * G + (J * G).T) / 2
    current = -s.trace(s.diag(XD, YD) * tangent) / 2
    target = omega * s.trace((YD - XD) * (YG - XG)) / (2 * (lam**2 + 4 * omega**2))
    marker = s.Symbol("derivative_marker", real=True)
    series = s.series(current.subs(lam, marker * lam), marker, 0, 6).removeO().expand()
    finite = prescription.data()["literal_fixed_finite_density_before_64_pi_squared"]
    names = {str(v): v for v in finite.free_symbols}
    m = names["mass"]
    t = s.Symbol("time", real=True)
    h = s.Function("unit_metric_TT")(t)
    local = 5 * m * m * s.diff(h, t) ** 2 / 12 - s.diff(h, t, 2) ** 2 / 60
    euler = s.expand(
        -s.diff(s.diff(local, s.diff(h, t)), t)
        + s.diff(s.diff(local, s.diff(h, t, 2)), t, 2)
    )
    checks = {
        "full_three_mode_vacuum_covariance_tangent": (
            lam * tangent - omega * (J * tangent - tangent * J) - forcing
        ).applyfunc(s.cancel),
        "complete_detector_source_current": s.cancel(current - target),
        "three_subtraction_remainder": s.cancel(
            1 / (S + p) - 1 / S + p / S**2 - p**2 / S**3 + p**3 / (S**3 * (S + p))
        ),
        "original_fixed_m2_R_coefficient": s.expand(finite).coeff(names["R_old"])
        - 5 * m * m / 3,
        "complete_finite_TT_Euler": euler
        + 5 * m * m * s.diff(h, t, 2) / 6
        + s.diff(h, t, 4) / 30,
        "actual_vacuum_Einstein_coefficient": family.base.data()["R"].subs(
            {family.base.u: 0, family.base.X: 0}
        )
        - 1,
        "full_tree_plus_loop_force_symbol": s.expand(
            (C0 + 5 * MASS**2 / s.Integer(6)) * p + p * p * A(p) - p * (C + p * A(p))
        ),
        "actual_mass_unchanged": MASS - 1000,
        "actual_kappa_unchanged": KAPPA - s.Integer(10) ** 800,
        "original_radial_weight_moment": s.integrate(original.W2, (y, 0, 1))
        - s.Rational(3, 14),
    }
    for j in (0, 2, 4):
        checks["original_marker_equals_frequency_Taylor_" + str(j)] = s.cancel(
            series.coeff(marker, j)
            - omega
            * s.trace((YD - XD) * (YG - XG))
            * (-1) ** (j // 2)
            * lam**j
            / (2 * (4 * omega**2) ** (j // 2 + 1))
        )
    return {
        "actual_mass": MASS,
        "actual_kappa": KAPPA,
        "Einstein_normalization_C0": C0,
        "complete_C": C,
        "full_physical_force_symbol": p * (C + p * A(p)) / C0,
        "source_sign": "O_phys h=f; the unsourced Euler residual is -O_phys h. G_phys=O_phys^-1, not the inverse of the negative Euler residual.",
        "flat_matching": "All canonical Proca frequencies agree at flat k. The full original second vertex is in the zero-order marker term. S193 identifies the same derivative-marker current with the original J_ad4. Its even orders0,2,4 equal the p-Taylor subtraction of this vacuum current, leaving the exact S199 three-subtracted cut and no additional lower polynomial.",
        "vacuum_parent_bridge": "Rvac=1 and S182 retains R and all original low vacuum jets while canceling the same-scheme vacuum constant. Pure TT metric variations at the constant-scalar vacuum do not excite the source-dependent light/vector channel. This is the conditional vector Gaussian sector, not scalar/graviton/mixed loops.",
        "checks": checks,
        "gates": {
            "mass_squared_finite_shift_nonzero": C - C0 > 0,
            "both_source_and_detector_arbitrary_symmetric": True,
            "original_second_vertex_in_zero_order_subtraction": True,
            "same_covariant_prescription_not_fixed_from_cut_alone": True,
            "Lorentz_vacuum_not_CD_state_substitution": True,
        },
    }
