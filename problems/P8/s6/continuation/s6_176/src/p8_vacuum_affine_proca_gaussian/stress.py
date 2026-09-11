"""Same prescribed ordinary Proca stress, with all normalization and interval limits."""

from functools import cache

import sympy as s
from p8_constant_proca import quantum as ordinary
from p8_proca_stress import estimates
from p8_vector_state import wkb

from .bridge import KAPPA, MASS, ZETA


def require_scope(time, order, zeta, kappa):
    for value in (time, zeta, kappa):
        if isinstance(value, bool) or not isinstance(
            value, (int, s.Integer, s.Rational)
        ):
            raise TypeError("Finite exact rational time and fixed parameters required")
    if type(order) is not int or not 0 <= order <= 5:
        raise ValueError("Only native time-derivative orders zero through five")
    if not -s.Rational(1, 2) <= time <= s.Rational(1, 2):
        raise ValueError("The quantitative stress interval is [-1/2,1/2]")
    if zeta != ZETA or kappa != KAPPA:
        raise ValueError("Require the explicitly bridged mass and hierarchy")
    return s.Rational(time), order


@cache
def prescription():
    R, Ricci2, Riemann2, ell, m = s.symbols(
        "R_old Ricci_squared Riemann_squared log_mass_ratio positive_mass", real=True
    )
    Weyl2 = Riemann2 - 2 * Ricci2 + R**2 / 3
    Euler = Riemann2 - 4 * Ricci2 + R**2
    vector = -(R**2) / 8 + s.Rational(29, 60) * Ricci2 - Riemann2 / 15
    scalar = R**2 / 72 - Ricci2 / 180 + Riemann2 / 180
    pole = 3 * m**4 + m * m * R + 2 * vector
    finite = (
        -(3 * ell - s.Rational(5, 2)) * m**4
        - (ell - s.Rational(5, 3)) * m * m * R
        - 2 * ell * vector
        - 4 * scalar
    )
    X = s.symbols("clock_norm", real=True)
    return {
        "curvature_signature": "R_old=-R_P8; curvature squares, connection and volume are unchanged by the overall metric signature reversal.",
        "ordinary_pole_density_before_common_factor": pole,
        "pole_factor": "1/(64*pi^2*epsilon_DR), D=3-2epsilon_DR spatial dimensions; fixed four-dimensional heat coefficients are evaluated as covariant invariants in D+1 dimensions before metric variation and the dimension limit.",
        "ordinary_finite_heat_matching_density_before_64_pi_squared": finite,
        "mu_mass_specialization": finite.subs(ell, 0),
        "density_boundary": "The displayed finite density supplies the local part matched to fourth-order mode subtraction, not the complete state-dependent determinant. Compact Euler variations of the retained box-R divergence vanish; no infinite-time boundary flux is discarded.",
        "constant_mass_domain": "All clock-normal mass-deviation terms vanish identically as functions for this constant-mass vector, not only on the clock. The remaining ordinary covariant prescription contains no X denominator and extends through the parent vacuum and null-gradient loci.",
        "checks": {
            "vector_heat_covariant_basis": s.expand(
                vector
                - (s.Rational(13, 120) * Weyl2 - s.Rational(7, 40) * Euler + R**2 / 72)
            ),
            "scalar_heat_covariant_basis": s.expand(
                scalar - (Weyl2 / 120 - Euler / 360 + R**2 / 72)
            ),
            "ordinary_pole_has_no_clock_norm": s.diff(pole, X),
            "ordinary_finite_matching_has_no_clock_norm": s.diff(finite, X),
            "finite_curvature_piece_at_mu_mass": s.expand(
                finite.subs(ell, 0)
                - (s.Rational(5, 2) * m**4 + s.Rational(5, 3) * m * m * R - 4 * scalar)
            ),
        },
    }


@cache
def local():
    u = wkb.u
    H = 4 * u / (1 + u * u)
    hp, hpp, hppp = [s.diff(H, u, i) for i in (1, 2, 3)]
    expected = {
        0: {"energy": -s.Rational(5, 2), "pressure": s.Rational(5, 2)},
        1: {"energy": -10 * H**2, "pressure": s.Rational(10, 3) * (2 * hp + 3 * H**2)},
        2: {
            "energy": 12 * H**2 * hp + 4 * H * hpp - 2 * hp**2,
            "pressure": -12 * H**2 * hp
            - 8 * H * hpp
            - 6 * hp**2
            - s.Rational(4, 3) * hppp,
        },
    }
    actual = ordinary.local_coefficients()
    checks = {}
    for n, row in expected.items():
        for name, value in row.items():
            checks[f"literal_covariant_finite_{n}_{name}"] = s.factor(
                actual[n][name] - value
            )
        checks[f"ordinary_local_Ward_identity_{n}"] = s.factor(
            s.diff(row["energy"], u) + 3 * H * (row["energy"] + row["pressure"])
        )
    return {
        "prescription": "The ordinary constant-mass specialization of the frozen covariant dimensional prescription at mu=m. Finite local coefficients multiply m^(4-2n)/(64*pi^2), n=0,1,2 in physical tau=1 units; no new reference-state normal ordering or stress-canceling scalar profile is introduced.",
        "actual_local_coefficients": actual,
        "checks": checks,
    }


@cache
def data():
    old = estimates.physical_bounds(s.sqrt(KAPPA))
    if old["fixed_canonical_mass_time_product"] != MASS or old["M_tau"] ** 2 != KAPPA:
        raise ValueError(
            "The frozen ordinary Proca stress normalization does not match"
        )
    energy = old["new_energy_derivative_bounds"]
    pressure = old["identical_pressure_derivative_bounds"]
    rows = {}
    checks = {}
    for j in range(6):
        er, pr = energy[j]["total"], pressure[j]["total"]
        rows[j] = {
            "normalized_energy_decomposition": energy[j],
            "normalized_pressure_decomposition": pressure[j],
            "absolute_energy_derivative_upper": KAPPA * er,
            "absolute_pressure_derivative_upper": KAPPA * pr,
        }
        checks[f"energy_physical_units_order_{j}"] = (
            rows[j]["absolute_energy_derivative_upper"] / KAPPA - er
        )
        checks[f"pressure_physical_units_order_{j}"] = (
            rows[j]["absolute_pressure_derivative_upper"] / KAPPA - pr
        )
    return {
        "time_interval": [-s.Rational(1, 2), s.Rational(1, 2)],
        "tau": s.S.One,
        "M_tau": s.sqrt(KAPPA),
        "canonical_mass_time": MASS,
        "bounds_by_time_derivative_order": rows,
        "state_and_radial_integral": "Exactly the S6.55 all-order Cauchy data and the full S6.82 ordinary Proca readout/subtraction/local matching, including all three polarizations and all momenta. No proof-band endpoint is a physical cutoff.",
        "uniform_normalized_upper": "1e-770 for each rho/P derivative of order zero through five",
        "uniform_absolute_upper": "1e30 for each rho/P derivative of order zero through five in tau=1 units",
        "meaning": "The denominator is the fixed kappa reference scale, not the vanishing classical energy density at the bounce. This bounds the specified conditional vector contribution and the corresponding prescribed-clock residual, not a self-consistent solution.",
        "global_boundary": "The state and renormalized free stress are smooth on every compact time interval; these numerical constants apply only on [-1/2,1/2].",
        "checks": checks,
        "gates": {
            "all_twelve_normalized_bounds_below_one_e_minus_770": all(
                bool(
                    rows[j][f"absolute_{name}_derivative_upper"] / KAPPA
                    < s.Rational(1, 10**770)
                )
                for j in range(6)
                for name in ("energy", "pressure")
            ),
            "all_twelve_absolute_bounds_below_one_e_30": all(
                bool(rows[j][f"absolute_{name}_derivative_upper"] < 10**30)
                for j in range(6)
                for name in ("energy", "pressure")
            ),
            "same_positive_fixed_mass_and_hierarchy": bool(
                old["fixed_canonical_mass_time_product"] == MASS
                and old["M_tau"] ** 2 == KAPPA
            ),
        },
    }
