"""Exact positive-series majorants for finite-time hard-channel tree control.

No frozen propagator or cancellation is used in these bounds. They apply
only to the momentum fibers, frequency band, time window and tree orders stated
in notes/finite-time.md. In particular they are not all-orders UV control.
"""

from dataclasses import dataclass
from fractions import Fraction
from functools import cache

import sympy as sp
from p8_control import oscillator
from p8_s5 import compact

ORDER = 4
DELTA = Fraction(1, 100)
LOWER_K = 10**8
UPPER_K = 10**9
DERIVATIVE_K = 10**10
MODE_P = 10**11
Q_MAX = 10**20
F_MAX = 1000
SEED = 10**30


@dataclass(frozen=True)
class Series:
    coefficients: tuple

    @classmethod
    def at(cls, degree, value):
        return cls(tuple(Fraction(value) if i == degree else Fraction(0) for i in range(ORDER+1)))

    def __add__(self, other):
        if not isinstance(other, Series):
            other = Series.at(0, other)
        return Series(tuple(a+b for a, b in zip(self.coefficients, other.coefficients)))

    __radd__ = __add__

    def __mul__(self, other):
        if not isinstance(other, Series):
            other = Series.at(0, other)
        return Series(tuple(sum(self.coefficients[i]*other.coefficients[n-i] for i in range(n+1))
                            for n in range(ORDER+1)))

    __rmul__ = __mul__

    def __pow__(self, power):
        if not isinstance(power, int) or power < 0:
            raise ValueError("Use a nonnegative integer power")
        result = Series.at(0, 1)
        for _ in range(power):
            result *= self
        return result

    def from_degree(self, degree):
        return Series(tuple(c if i >= degree else Fraction(0) for i, c in enumerate(self.coefficients)))

    def analytic(self, exponent):
        """Majorant of (1+self)^exponent when self has zero constant term."""
        if self.coefficients[0] != 0:
            raise ValueError("Analytic remainder must vanish at degree zero")
        return sum(abs(Fraction(str(sp.binomial(exponent, n))))*self**n for n in range(ORDER+1))


def domain_checks():
    delta = sp.Rational(DELTA.numerator, DELTA.denominator)
    switch = sp.Rational(9, 50)
    drift = delta/(1-delta)
    assert (switch+drift)**2 < sp.Rational(1, 17)
    assert (switch-drift)**2 > sp.Rational(1, 65)
    # e^v <= 1/(1-v) for 0<=v<1 supplies rational scale/energy enclosures.
    a_exponent = 2*delta/(1-delta)
    energy_exponent = 20*delta/(1-delta)
    assert 1/(1-a_exponent) < 2
    assert 1/(1-energy_exponent) < 2
    assert LOWER_K**2//16 >= 10**14
    assert (8*UPPER_K)**2 < Q_MAX
    assert 4*UPPER_K < DERIVATIVE_K
    assert 5*UPPER_K < MODE_P/4
    assert Fraction(8*MODE_P, LOWER_K**2) < 1
    f_bounds = {}
    for chart in ("unitary", "gamma"):
        f = oscillator.derive(chart)["f"].subs(oscillator.q, 1/oscillator.r)
        report = oscillator.coefficient_bound(f, chart)
        assert sp.Rational(report["absolute_bound"]) < F_MAX
        f_bounds[chart] = report
    # K in [1/32,130] in fixed local units; a^-3/2/sqrt(2K)<12.
    scalar_p = 12*(260*(MODE_P+2*F_MAX)+260*(8*Q_MAX+4))
    gamma_p = 96*Q_MAX
    tensor_p = 2*(MODE_P+6)
    assert max(12, scalar_p, gamma_p, tensor_p) < SEED
    return {"chart_switch_abs_x": str(switch), "x_drift": str(drift),
            "a_log_drift": str(a_exponent), "energy_log_growth": str(energy_exponent),
            "energy_growth_upper": str(1/(1-energy_exponent)),
            "normalization_f_bounds": f_bounds,
            "canonical_seed_bounds": {"scalar_coordinate": "12", "scalar_momentum": str(scalar_p),
                                      "gamma_unitary_momentum": str(gamma_p), "tensor_momentum": str(tensor_p),
                                      "common_loose_seed": str(SEED)}}


@cache
def build():
    checks = domain_checks()
    # Four possible labelled legs per field/component; the degree norm sums
    # absolute coefficients of all masks. This also bounds any three-leg case.
    field = Series.at(1, 4*SEED)
    perturbation = 3*field  # 2*zeta*I+gamma, entrywise
    metric = 1+perturbation
    inverse = 1+sum(3**(n-1)*perturbation**n for n in range(1, ORDER+1))
    # Six determinant permutations, with the exact background value removed.
    determinant_remainder = (6*metric**3).from_degree(1)
    volume = determinant_remainder.analytic(sp.Rational(1, 2))
    inverse_volume = determinant_remainder.analytic(sp.Rational(-1, 2))
    connection = Fraction(9, 2)*DERIVATIVE_K*inverse*perturbation
    curvature = 9*inverse*(6*DERIVATIVE_K*connection+18*connection**2)
    pi = Series.at(0, 4)+(4+Fraction(1, 6)+1+4)*field
    shift_orders = []
    for degree in (1, 2, 3):
        source = 3*DERIVATIVE_K*pi+9*connection*pi
        # Each entry of W <=2*source/k_min^2 <=2*source and |LW|<=4*K*|W|.
        correction = Series.at(degree, 8*DERIVATIVE_K*source.coefficients[degree])
        pi += correction
        shift_orders.append(str(correction.coefficients[degree]))
    mixed = 3*pi*metric*inverse_volume
    sigma = (3*mixed).from_degree(1)  # the exact background cancels
    shear2 = (12*mixed**2).from_degree(2)  # exact background/linear cancellations
    density = Series.at(0, 0)
    entries = []
    for order in compact.coefficient_report()["hamiltonian"]:
        for powers, proof in order.items():
            p, r, s = map(int, powers.split(','))
            # ell(t)/ell0 in [1/2,2], weight p+2r+2s-2 lies in [-2,6].
            weight = p+2*r+2*s-2
            assert -2 <= weight <= 6
            bound = 64*Fraction(proof["absolute_bound"])
            density += bound*sigma**p*curvature**r*shear2**s
            entries.append({"powers": powers, "compact_bound": proof["absolute_bound"],
                            "fixed_unit_bound": str(bound), "unit_weight": weight})
    # The explicit Adot and gamma generating terms are at most quadratic.
    # a^3<=8 converts the reduced physical-volume density to fixed coordinates.
    h = 8*(volume*density+72*pi*metric)
    # A loose 6! factor covers normalized repeated-mode occupation factors
    # in the connected tree Wick contractions (four external plus two ends).
    B3, B4 = 720*h.coefficients[3], 720*h.coefficients[4]
    phase_measure = 3*(4*UPPER_K+1)**3  # exceeds three-species d^3k/(2*pi)^3 volume
    # At fixed total momentum, one-particle species measure is 3 and the
    # two-particle relative-momentum/species measure is <=3*phase_measure.
    # This deliberately loose common Schur majorant bounds both. Vacuum
    # production and a physical periodic-universe state count are NOT used.
    multiplicity = 4*phase_measure**3
    C3 = multiplicity*B3
    C4 = multiplicity*(B4+18*B3**2)
    scale = 1
    while scale < 1000*C3 or scale*scale < 1000*C4:
        scale *= 10
    return {"domain_checks": checks,
            "window": {"half_duration_over_ell0": str(DELTA), "external_k_min": str(LOWER_K),
                       "external_k_max": str(UPPER_K), "proper_subset_k_min": str(LOWER_K),
                       "Fourier_measure": "d^3k/(2*pi)^3; fixed total momentum fibers on R^3",
                       "momentum_species_measure_upper": str(phase_measure),
                       "Schur_measure_majorant": str(multiplicity)},
            "stationary_coefficient_majorants": entries,
            "momentum_correction_majorants": shift_orders,
            "series_coefficients": {"curvature": list(map(str, curvature.coefficients)),
                                    "sigma": list(map(str, sigma.coefficients)),
                                    "shear2": list(map(str, shear2.coefficients)),
                                    "Hamiltonian": list(map(str, h.coefficients))},
            "cubic_kernel_bound": str(B3), "quartic_Hamiltonian_kernel_bound": str(B4),
            "cubic_transition_block_bound_numerator": str(C3),
            "quartic_connected_tree_block_bound_numerator": str(C4),
            "sufficient_M_tau": str(scale), "sufficient_M_tau_power10": len(str(scale))-1,
            "target_block_norm": "1/1000",
            "not_a_necessary_scale_or_all_orders_cutoff": True}
