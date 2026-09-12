"""Full constrained physical readouts on the joint complex domain."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference, stress

from . import domain


def transverse(ell, e, f, p, a, m):
    return s.Matrix([*(p * e), *(s.I * ell.cross(e) * f / a), 0, *(m * f * e)])


def longitudinal(n, rho, f, p, a, m, omega):
    return s.Matrix(
        [*(m * p * n / omega), 0, 0, 0, -s.I * rho * p / (a * omega), *(omega * f * n)]
    )


@cache
def data():
    a, m, k, omega = s.symbols("a m k omega", positive=True)
    f, p = s.symbols("f p")
    e1 = s.Matrix([1, 0, 0])
    e3 = s.Matrix([0, 0, 1])
    ell = k * e3
    previous = stress.data()["phase_stripped_mode_readouts"]
    T = transverse(ell, e1, f, p, a, m)
    L = longitudinal(e3, k, f, p, a, m, omega)
    checks = {
        "full_transverse_constraint_bridge": T - previous["one_transverse"],
        "full_longitudinal_constraint_bridge": L - previous["longitudinal"],
        "unchanged_joint_pair_upper_count": 2 * reference.FIELD_NORM**2 - 32000000,
        "same_real_volume_cancellation": a**3
        * a ** s.Rational(-3, 2)
        * a ** s.Rational(-3, 2)
        - 1,
    }
    return {
        "readouts": "In the ten-field order(E,B,mA0,mAsp), T=(p e,i(l cross e)f/a,0,m f e), L=((m/omega)p n,0,-i rho p/(a omega),omega f n), with rho=sqrt(l.l). The common a^(-3/2) is outside and cancels the stress-smearing volume exactly.",
        "genuine_complex_bounds": "The analytic frame vectors have Euclidean norm<2. The joint domain gives |f|<=nu^-1/2, |p|<64sqrt(nu), |m/omega|<2, ||l/a||<2nu and |rho/(a omega)|<3. Thus each of all ten components is below1000sqrt(nu) and the complete field norm is below4000sqrt(nu).",
        "complete_pair": "Every fixed full stress matrix has operator norm<=1/2, hence the entire sixteen-component pair tensor has Frobenius norm<=2||u_k||||u_l||<1e9nu. All nine polarization pairs use this same bound; no COM or longitudinal-only replacement is made.",
        "phase_sum": "Both W frequencies have real part>.5nu on the joint domain. Their inverse sum g=1/(W_k+W_l) is holomorphic with |g|<1/nu.",
        "Schwarz": "The conjugate-on-real-parameters readouts have analytic Schwarz continuations on the same domain and identical bounds. Holomorphy is asserted for these pre-current coefficients, not for taking an imaginary part as a complex operation.",
        "checks": checks,
        "gates": {
            "complete_transverse_electric_components": 2 * 64 < 1000,
            "complete_transverse_magnetic_components": 2 * 2 < 1000,
            "complete_longitudinal_electric_components": 2 * 64 * 2 < 1000,
            "complete_temporal_constraint_component": 3 * 64 < 1000,
            "complete_longitudinal_mass_components": 2 * 2 < 1000,
            "ten_field_norm_below_four_thousand": 10 * 1000**2
            < reference.FIELD_NORM**2,
            "full_pair_norm_below_existing_constant": 2 * reference.FIELD_NORM**2
            < reference.PAIR_REF,
            "same_joint_radius": domain.DELTA == s.Rational(1, 10**6),
        },
    }
