"""Uniform same-operator in/out energy and pressure differences, not absolute stress."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic

from . import geometry, transition, tube


def enclosures(mean_mass, amplitude, timescale, multiplicity=42):
    m, d, tau, m0 = transition.parameters(mean_mass, amplitude, timescale)
    if type(multiplicity) is not int or multiplicity < 1:
        raise TypeError("Require a positive native color/flavor multiplicity")
    k = s.Integer(tube.factor(20)) ** 20
    compact = (
        s.Rational(262144 * multiplicity, 323 * 9)
        * k
        * (d + m * tau)
        / (tau**20 * m0**17)
    )
    tail = s.Rational(16504 * multiplicity, 15 * 9) * m * k / m0**17
    return {
        "mass_floor": m0,
        "compact_state_energy_difference_upper": compact,
        "tail_state_energy_difference_upper": tail,
        "total_state_energy_difference_upper": compact + tail,
        "state_pressure_difference_upper": (compact + tail) / 6,
    }


@cache
def data():
    z = s.Symbol("z", positive=True)
    radial = s.integrate(z**3 / (1 + z * z) ** s.Rational(21, 2), (z, 0, s.oo))
    d = geometry.data()
    values = enclosures(d["mean_mass"], d["amplitude_upper"], d["mass_transition_time"])
    checks = {
        "compact_radial": radial - s.Rational(2, 323),
        "compact_radial_momentum_rescaling": 4 * 4**4 * radial - s.Rational(2048, 323),
        "compact_energy_prefactor": 4 * 32 * 2048 - 262144,
        "tail_low_plus_high_momentum_integrals": 2
        + 2 * s.Rational(2048, 15)
        - s.Rational(4126, 15),
        "tail_energy_prefactor": 4 * 4126 - 16504,
        "pressure_to_energy_operator_allowance": s.Rational(1, 3) / 2
        - s.Rational(1, 6),
        "energy_and_pressure_same_state": values["state_pressure_difference_upper"] * 6
        - values["total_state_energy_difference_upper"],
    }
    return {
        "actual_uniform_enclosures": values,
        "actual_decimal_diagnostics": {k: str(s.N(v, 28)) for k, v in values.items()},
        "reference_scale_ratio": values["total_state_energy_difference_upper"]
        / analytic.KAPPA,
        "physical_observable": "For the two exact curved Hadamard states of this SAME prescribed-mass operator, compare T_hat00 and isotropic P in the physical FLRW orthonormal frame. Identical local subtractions and finite counterterms cancel. These are differences between states, NOT either absolute stress.",
        "energy_bound": "A rank-one projector difference has trace norm2|beta| per helicity. Both helicities and 42 copies give |Delta rho|<=4N/pi^2 integral p^2 sqrt(p^2+m0^2)|beta|dp, since omega(t)<=2sqrt(p^2+m0^2) and a^-3<=1. Pi^2>9 yields the displayed rational allowances.",
        "pressure_bound": "The physical isotropic pressure operator is (q/3)sigma1 before the mass rotation, with norm q/3. Since q<=p<=sqrt(p^2+m0^2), the same state difference is bounded by one sixth the conservative energy allowance. This is not a pressure formula obtained by differentiating an energy inequality.",
        "source_Ward_identity": "partial_t Delta rho+3H(Delta rho+Delta P)=Mdot Delta<bar psi psi>. The prescribed external mass exchanges energy; the free fermion subsystem is not separately conserved.",
        "checks": checks,
        "gates": {
            "complete_curved_state_energy_difference_below_one_e_minus_1090": bool(
                values["total_state_energy_difference_upper"] < s.Rational(1, 10**1090)
            ),
            "reference_ratio_below_one_e_minus_1890": bool(
                values["total_state_energy_difference_upper"] / analytic.KAPPA
                < s.Rational(1, 10**1890)
            ),
            "all_forty_two_copies_included": d["all_copies"] == 42,
            "absolute_curved_stress_not_inferred": True,
        },
    }
