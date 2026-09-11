"""Pinned CD geometry, physical spinor normalization and modewise Moller limits."""

import hashlib
import json
from functools import cache

import sympy as s
from p8_affine import dictionary
from p8_exceptional_vacuum import analytic
from p8_vacuum_flat_dirac_production import calibration

WITNESS_SHA = "caf8c8e688a7565b9d00f921c099a28da00f97522ed26ad182a8227eb80cd4dd"


@cache
def data():
    path = dictionary.WITNESS
    if (
        path.name != "witness-CD_matter.json"
        or hashlib.sha256(path.read_bytes()).hexdigest() != WITNESS_SHA
    ):
        raise ValueError("The actual CD witness changed")
    background = json.loads(path.read_text())["background"]
    t, p, M, Md, q, H = s.symbols("t p M Mdot q H", real=True)
    a = (1 + t * t) ** 2
    hubble = s.diff(a, t) / a
    sigma1 = s.Matrix([[0, 1], [1, 0]])
    sigma3 = s.diag(1, -1)
    h = q * sigma1 + M * sigma3
    source = calibration.data()
    m = source["quadratic_subsystem_mean_mass"]
    d = source["profile_amplitude_rational_upper"]
    tau = source["same_profile_time_scale"]
    checks = {
        "pinned_scale_factor": a - s.sympify(background["a"], locals={"t": t}),
        "pinned_Hubble": s.simplify(
            hubble - s.sympify(background["H"], locals={"t": t})
        ),
        "actual_CD_tail_power": s.limit(t * hubble, t, s.oo) - 4,
        "actual_bounce_acceleration": s.diff(hubble, t).subs(t, 0) - 4,
        "physical_comoving_momentum": s.simplify(s.diff(p / a, t) + hubble * p / a),
        "volume_spinor_connection_cancelled": s.simplify(
            s.diff(a ** s.Rational(-3, 2), t) / a ** s.Rational(-3, 2) + 3 * hubble / 2
        ),
        "physical_mode_gap": s.trace(h * h) / 2 - q * q - M * M,
        "physical_squared_mode_offdiagonal": (h * h)[0, 1],
        "initial_rotation_connection": s.simplify(
            -(M * (-H * q) - q * Md) / (2 * (q * q + M * M))
            - q * (Md + H * M) / (2 * (q * q + M * M))
        ),
        "same_mass": m - 10**200,
        "same_amplitude_upper": d - 3 * 10**197,
        "same_transition_time": tau - s.Rational(1, 10**100),
        "same_clock_scale": analytic.KAPPA - 10**800,
        "inert_plus_active_color_flavor_copies": 36 + 6 - 42,
    }
    return {
        "witness_sha256": WITNESS_SHA,
        "witness_background": background,
        "scale_factor": a,
        "Hubble": hubble,
        "physical_scalar_curvature": s.factor(
            -6 * (s.diff(hubble, t) + 2 * hubble * hubble)
        ),
        "momentum": p / a,
        "physical_helicity_Hamiltonian": h,
        "initial_exact_connection": q * (Md + H * M) / (2 * (q * q + M * M)),
        "mean_mass": m,
        "amplitude_upper": d,
        "mass_transition_time": tau,
        "active_copies": 6,
        "inert_copies": 36,
        "all_copies": 42,
        "units": "Use the existing reference time unit mPhi=1 and original CD bounce scale one. The canonical clock is Phi=sqrt(kappa)t. tau=R/sqrt(kappa)=1e-100 is the mass transition time, NOT the geometric bounce time.",
        "ancestry": "The analytic affine-domain family inherits the retuned exceptional-vacuum tree; p8_affine.dictionary.target reads this pinned CD_matter witness. The retuning proportional to (X-1)^2 does not change its background jets.",
        "Cauchy_dictionary": "For signature(+---), psi_c=a^(3/2)psi identifies the physical Cauchy norm with flat L2. Its helicity Hamiltonian is q sigma1+M sigma3, q=p/a. No metric-dependent source variation is replaced by variation of this rescaling.",
        "Moller_bound": "For each fixed p and T>=max(1,tau), ||W_infinity-W(T)||<=p/(3T^3)+Delta*tau^8/(56T^7). W and adjoints converge in each finite mode; dominated convergence then gives strong one-particle unitary limits, NOT momentum-uniform norm convergence or a global Fock implementer.",
        "checks": checks,
    }
