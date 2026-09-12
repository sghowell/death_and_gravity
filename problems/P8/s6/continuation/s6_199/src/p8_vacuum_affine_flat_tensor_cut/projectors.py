"""Full spin2/spin0 cut with the actual Wick and phase-space normalization."""

from functools import cache

import sympy as s

from . import polarizations as pol

S = s.Symbol("s", positive=True)
MASS = pol.MASS


def eigenvalues():
    p = pol.spectral_polynomials()
    return p["spin2_invariant"], p["spin0_invariant"]


def above_threshold_density(spin):
    a, b = eigenvalues()
    if (
        isinstance(spin, bool)
        or not isinstance(spin, (int, s.Integer))
        or spin not in (0, 2)
    ):
        raise ValueError("Only the conserved spin0 and spin2 sectors")
    return s.sqrt(1 - 4 * MASS**2 / S) * (a if spin == 2 else b) / (32 * s.pi**2)


def density(spin):
    value = above_threshold_density(spin)
    return s.Piecewise((s.Integer(0), S <= 4 * MASS**2), (value, True))


@cache
def data():
    I = s.eye(3)
    P0 = s.Matrix(9, 9, lambda u, v: I[u // 3, u % 3] * I[v // 3, v % 3] / 3)
    P2 = (
        s.Matrix(
            9,
            9,
            lambda u, v: (
                (
                    I[u // 3, v // 3] * I[u % 3, v % 3]
                    + I[u // 3, v % 3] * I[u % 3, v // 3]
                )
                / 2
            ),
        )
        - P0
    )
    x = s.Symbol("x", real=True)
    a, b = eigenvalues()
    beta = s.Symbol("beta", positive=True)
    phase = beta / (8 * s.pi)
    checks = {
        "spin0_projector": P0 * P0 - P0,
        "spin2_projector": P2 * P2 - P2,
        "orthogonal_spin_projectors": P0 * P2,
        "spin0_rank": s.trace(P0) - 1,
        "spin2_rank": s.trace(P2) - 5,
        "full_spin2_polynomial": s.expand(
            a - (13 * S * S + 56 * MASS * MASS * S + 48 * MASS**4) / 120
        ),
        "full_spin0_polynomial": s.expand(
            b - (S * S - 4 * MASS * MASS * S + 12 * MASS**4) / 12
        ),
        "phase_Wick_current_and_spectral_factors": 2 * phase / (4 * 2 * s.pi)
        - beta / (32 * s.pi**2),
        "spin0_strict_positive_square": s.expand(
            1 - 4 * x + 12 * x * x - 12 * (x - s.Rational(1, 6)) ** 2 - s.Rational(2, 3)
        ),
        "spin2_upper_on_threshold_interval": s.expand(
            30 - (13 + 56 * x + 48 * x * x) - (s.Rational(1, 4) - x) * (48 * x + 68)
        ),
        "spin0_upper_on_threshold_interval": s.expand(
            1 - (1 - 4 * x + 12 * x * x) - 4 * x * (1 - 3 * x)
        ),
    }
    return {
        "conserved_projectors": "theta_ab=eta_ab-P_a P_b/s; P0=theta_ab theta_cd/3, P2=(theta_ac theta_bd+theta_ad theta_bc)/2-P0. At timelike P the transverse space is three dimensional. The finite matrices check the COM representatives; tensor covariance gives arbitrary timelike momentum.",
        "Lorentz_phase_space": "The two-body integral with measure d^3k/((2pi)^3 2E) per leg and (2pi)^4 delta4(P-k-l) is beta/(8pi), beta=sqrt(1-4m^2/s).",
        "normalization": "The connected stress Wightman cut is2*phase_space*(a P2+b P0). For the actual current i theta<[T,T]>/4, its positive dispersion density is Wightman/(8pi), so rho2=beta*(13s^2+56m^2s+48m^4)/(3840pi^2) and rho0=beta*(s^2-4m^2s+12m^4)/(384pi^2).",
        "support": "Both densities vanish below and at4m^2, and are strictly positive above it. They are a flat external-metric Gaussian benchmark, not full interacting parent scattering cuts.",
        "retarded_sign": "With Fourier convention exp(iP.x), the current kernel i theta<[T,T]>/4 has denominators1/(sigma-z), z=(omega+i0)^2-|p|^2. Thus Im chi(s+i0)=+pi rho(s) for positive energy.",
        "densities": {"spin2": density(2), "spin0": density(0)},
        "checks": checks,
        "gates": {
            "two_physical_transverse_tensor_sectors": P2.shape == (9, 9),
            "no_Maxwell_only_replacement": a.subs(MASS, 0) == 13 * S * S / 120,
            "nonzero_longitudinal_scalar_UV": b.subs(MASS, 0) == S * S / 12,
        },
    }
