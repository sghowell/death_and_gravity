"""Independent canonical-volume, York and lapse-Schur reconstruction."""

from functools import cache

import sympy as sp

from . import model


@cache
def data():
    d = model.generic()
    co = d["symbols"]
    epsilon, nu = sp.symbols("quadratic_order nonzero_mode_lapse", real=True)
    v, chi, Pv, Pc = model.FIELDS
    A = model.A
    p, l = model.p, model.l
    r, s = Pv / A**3, Pc / A**3
    q = model.k**2 / A**2

    def series(name):
        values = co[name]
        return (
            values[0] + epsilon * nu * values[1] + epsilon**2 * nu * nu * values[2] / 2
        )

    vol = sp.exp(3 * epsilon * v)
    trace = (p + epsilon * r / 3) / vol
    matter = (l + epsilon * s) / vol
    shear = (r - 9 * p * v - 3 * l * chi) ** 2 / 24
    local = vol * (
        series("a") * trace**2
        + series("b") * trace
        + series("c")
        + series("d") * matter**2
        + co["e"][0] * epsilon**2 * shear
        + co["g"][0] * epsilon**2 * q * chi**2
    )
    quadratic = sp.expand(sp.diff(local, epsilon, 2).subs(epsilon, 0) / 2)
    # Integrated curvature: volume*R3=e^zeta(-4 Delta zeta-2|grad zeta|^2)/A^2.
    # The constant-lapse second coefficient is 2q zeta^2.
    # A first-order lapse multiplies the unintegrated linear 4q zeta.
    quadratic += 2 * co["f"][0] * q * v * v + 4 * co["f"][1] * q * v * nu
    expected = (
        d["uneliminated_quadratic_density"]
        + nu * d["linear_inhomogeneous_lapse_force"]
        + nu * nu * d["independent_mean_lapse_pivot"] / 2
    )
    nstar = -d["linear_inhomogeneous_lapse_force"] / d["independent_mean_lapse_pivot"]
    kvec = sp.Matrix(
        sp.symbols("nonzero_wave_x nonzero_wave_y nonzero_wave_z", real=True)
    )
    projector = kvec * kvec.T / (kvec.dot(kvec)) - sp.eye(3) / 3
    sig = projector * (3 * l * chi - r + 9 * p * v) / 4
    york = (
        sig * kvec - kvec * (l * chi / 2 - r / 6 + sp.Rational(3, 2) * p * v)
    ).applyfunc(sp.factor)
    return {
        "natural_canonical_trace_density": trace,
        "natural_canonical_matter_charge_density": matter,
        "linear_York_tracefree_density": sig,
        "full_uneliminated_scalar_quadratic_density": quadratic,
        "nonzero_mode_lapse_solution": nstar,
        "checks": {
            "literal_canonical_volume_and_curvature_quadratic_expansion": sp.factor(
                quadratic - expected
            ),
            "linear_York_constraint_in_all_cartesian_directions": york,
            "linear_York_momentum_is_tracefree": sp.factor(sp.trace(sig)),
            "literal_linear_York_shear_square": sp.factor(sp.trace(sig * sig) - shear),
            "actual_nonzero_mode_lapse_stationarity": sp.factor(
                sp.diff(expected, nu).subs(nu, nstar)
            ),
            "actual_nonzero_mode_lapse_Schur_reconstruction": sp.factor(
                expected.subs(nu, nstar) - d["raw_density_quadratic_Hamiltonian"] / A**3
            ),
        },
    }


@cache
def canonical():
    bg = model.parent.coefficients()["background"]
    a = (1 + model.u**2) ** 2
    ell, H = bg["ell"], bg["H"]
    B = sp.eye(4)
    B[2, 0] = -18 * a**3 * H
    B[2, 1] = 3 * a**3 * ell
    B[3, 0] = 3 * a**3 * ell
    Omega = model.phase.data()["constant_symplectic_form"]
    o = model.global_model.old
    old = model.phase.data()["regular_density_generator"].subs(
        model.global_model.substitution(), simultaneous=True
    )
    old = old.subs(o.q, model.k**2 / a**2)
    raw = Omega * sp.hessian(
        model.on_clock()["raw_density_quadratic_Hamiltonian"], model.FIELDS
    )
    return {
        "old_to_natural_phase": B,
        "checks": {
            "both_natural_canonical_momentum_shifts_are_symplectic": (
                B * Omega * B.T - Omega
            ).applyfunc(sp.factor),
            "actual_natural_phase_generator_retains_time_dependent_boundary": (
                raw * B - B * old - sp.diff(B, model.u)
            ).applyfunc(sp.factor),
            "actual_natural_scalar_Hamiltonian_is_existing_global_scalar_Hamiltonian": model.bridge()[
                "difference"
            ],
        },
    }
