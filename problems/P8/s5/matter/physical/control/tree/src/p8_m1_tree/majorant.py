"""Exact positive-recursion bounds for coupled M1 cubic/quartic trees.

Only local hard-band one-to-two/two-to-one and hard-masked two-to-two blocks
are addressed.  Bounds use exact free canonical columns, not frozen scalar
frequencies or a scalar-only propagator.  No infinite-time vacuum, full Fock
space norm, loop estimate, or optimal physical cutoff follows from this code.
"""

from dataclasses import dataclass
from fractions import Fraction
from functools import cache

import sympy as sp
from p8_m1 import series as invariant_series
from p8_m1_control import bounds as free_bounds
from p8_m1_control import model, oscillator

ORDER = 4
DELTA = Fraction(1, 100)
LOWER_K = 10**11
UPPER_K = 10**12
DERIVATIVE_K = 10**13
MODE_P = 10**14
Q_MAX = 10**27
SPECIES = 4  # Two coupled scalar columns and two normalized TT columns.
EXCHANGE_FACTOR = 3*SPECIES*2
WICK_FACTOR = 720
TARGET_DENOMINATOR = 1000


@dataclass(frozen=True)
class Series:
    """Nonnegative degree norms, with exact convolution truncated after H4.

    A degree norm sums absolute coefficients of all labelled masks.  Products
    with repeated nilpotent labels vanish, so ordinary positive convolution
    is an upper bound even though it also counts those forbidden products.
    """

    coefficients: tuple

    def __post_init__(self):
        values = tuple(Fraction(value) for value in self.coefficients)
        if len(values) != ORDER+1 or any(value < 0 for value in values):
            raise ValueError("Use five nonnegative exact degree coefficients")
        object.__setattr__(self, "coefficients", values)

    @classmethod
    def at(cls, degree, value):
        if not isinstance(degree, int) or not 0 <= degree <= ORDER:
            raise ValueError("Degree must be between zero and four")
        return cls(tuple(Fraction(value) if i == degree else Fraction(0)
                         for i in range(ORDER+1)))

    def __add__(self, other):
        if not isinstance(other, Series):
            other = Series.at(0, other)
        return Series(tuple(a+b for a, b in zip(self.coefficients, other.coefficients)))

    __radd__ = __add__

    def __mul__(self, other):
        if not isinstance(other, Series):
            other = Series.at(0, other)
        return Series(tuple(sum(self.coefficients[i]*other.coefficients[n-i]
                                for i in range(n+1)) for n in range(ORDER+1)))

    __rmul__ = __mul__

    def __pow__(self, power):
        if not isinstance(power, int) or power < 0:
            raise ValueError("Use a nonnegative integer power")
        result = Series.at(0, 1)
        for _ in range(power):
            result *= self
        return result

    def from_degree(self, degree):
        if not 0 <= degree <= ORDER+1:
            raise ValueError("Invalid degree truncation")
        return Series(tuple(value if i >= degree else 0
                            for i, value in enumerate(self.coefficients)))

    def analytic(self, exponent):
        """Absolute finite binomial expansion of (1+self)^exponent."""
        if self.coefficients[0] != 0:
            raise ValueError("An analytic remainder must have zero constant term")
        return sum(abs(Fraction(str(sp.binomial(exponent, n))))*self**n
                   for n in range(ORDER+1))


def _power_of_ten_at_least(value):
    result = 1
    while result < value:
        result *= 10
    return result


@cache
def canonical_boundary_bounds():
    """Prove ||Sbar||<=q*CS using exact finite-band Laurent coefficients."""
    rows = {}
    for chart in ("unitary", "gamma"):
        data = oscillator.derive(chart)
        coefficients = {}
        for suffix in ("11", "22", "12_factor"):
            value = model.canonical(model.z*data["momentum_boundary"+suffix])
            model.no_high_frequency_pole(value)
            coefficients[suffix] = free_bounds.coefficient_bound(value, chart)
        values = {key: Fraction(row["absolute_bound"]) for key, row in coefficients.items()}
        bound = max(values["11"], values["22"])+10*values["12_factor"]
        rows[chart] = {"z_times_boundary_coefficients": coefficients,
                       "boundary_operator_bound_over_q": str(bound)}
    return rows


@cache
def domain_checks(lower_k=LOWER_K, upper_k=UPPER_K, derivative_k=DERIVATIVE_K,
                  mode_p=MODE_P, q_max=Q_MAX):
    """Exact sufficient inequalities; a smaller unsafe band fails closed.

    Initial external magnitudes lie in [L,U]; all proper subset magnitudes
    are at least L.  Internal tree magnitudes are at most 2U.  Momenta here
    are dimensionless in the fixed centre unit ell0, with a(t0)=1.
    """
    lower_k, upper_k, derivative_k, mode_p, q_max = map(
        Fraction, (lower_k, upper_k, derivative_k, mode_p, q_max))
    if lower_k <= 0 or upper_k < lower_k:
        raise ValueError("Require a positive ordered hard band")
    threshold = Fraction(free_bounds.build_bounds()["q_threshold"])
    drift = DELTA/(1-DELTA)
    energy_ratio = Fraction(11, 7)
    # Each canonical free column starts with E_j=R_jj/2, R=sqrt(W(left)).
    # W<=3*k_fixed^2/2 and |k_internal|<=4U give lambda_max(R)<5U.
    energy_upper = energy_ratio*5*upper_k/2
    margins = {
        "gamma_window": Fraction(1, 4)-Fraction(9, 50)-drift,
        "unitary_window": Fraction(9, 50)-drift-Fraction(1, 9),
        "scale_factor_exponential": 1-4*drift,
        "scale_factor_upper_below_two": 2-1/(1-4*drift),
        "ell_ratio_lower_above_half": 1-DELTA-Fraction(1, 2),
        "ell_ratio_upper_below_two": 2-(1+DELTA),
        "centre_high_q": lower_k**2-2*threshold,
        "whole_window_high_q": lower_k**2/16-threshold,
        "whole_window_q_upper": q_max-(8*upper_k)**2,
        "derivative_band": derivative_k-8*upper_k,
        "York_denominator_above_one": lower_k**2/4-1,
        "initial_frequency_upper": 25-24,
        "mode_energy_below_half_MODE_P": mode_p/2-energy_upper,
        "mode_coordinate_squared_below_one": 1-16*energy_upper/lower_k**2,
        "mode_momentum_squared_below_MODE_P_squared": mode_p**2-2*energy_upper,
        "full_window_duration_below_one": 1-2*DELTA,
    }
    if any(value <= 0 for value in margins.values()):
        failed = [name for name, value in margins.items() if value <= 0]
        raise ValueError("Unsafe tree domain: "+", ".join(failed))
    boundary = canonical_boundary_bounds()
    CS = max(Fraction(row["boundary_operator_bound_over_q"]) for row in boundary.values())
    # alpha_bar is in [1/100,512].  s=ell/ell0 in [1/2,2] gives
    # ||T_fixed||<46, ||T_fixed^-1||<20, and a^-3/2<3.
    chart_Q = Fraction(60)
    chart_p = 138*(mode_p+2*q_max*CS)
    # q_fixed>=L^2/4>1; retaining only >=1 is a conservative gamma seed.
    zeta = max(chart_Q, chart_p/2)
    pv = max(chart_p, 2*q_max*chart_Q)
    # l_fixed<=1/5.  Both mixed-boundary shifts must be included.
    old_p = pv+Fraction(3, 5)*chart_Q
    matter_density_perturbation = chart_p+Fraction(3, 5)*zeta
    tensor_Q, tensor_p = Fraction(6), 2*(mode_p+12)
    terminal = {"chart_coordinate_vector_norm": chart_Q,
                "chart_momentum_vector_norm": chart_p, "zeta": zeta,
                "pre_mixed_metric_momentum": old_p,
                "matter_field": chart_Q, "matter_density_perturbation": matter_density_perturbation,
                "tensor_coordinate_entry": tensor_Q, "tensor_momentum_entry": tensor_p}
    seed = _power_of_ten_at_least(max(terminal.values()))
    return {"margins": {key: str(value) for key, value in margins.items()},
            "S5_7_q_threshold": str(threshold), "chart_switch_abs_x": "9/50",
            "ell_and_scale_factor_ratio_enclosure": ["1/2", "2"],
            "physical_momentum_enclosure": [str(lower_k/2), str(4*upper_k)],
            "proper_subset_squared_lower": str(lower_k**2/4),
            "local_and_fixed_squared_momentum_upper": str(q_max),
            "derivative_bound": str(derivative_k), "free_energy_ratio_upper": "11/7",
            "free_column_energy_upper": str(energy_upper),
            "free_column_coordinate_norm_upper": "1", "free_column_momentum_norm_upper": str(mode_p),
            "canonical_boundary_bounds": boundary, "common_boundary_bound_over_q": str(CS),
            "terminal_phase_seed_bounds": {key: str(value) for key, value in terminal.items()},
            "common_terminal_phase_seed": str(seed)}


def geometry_majorants(field, derivative_bound):
    """Full coordinate geometry with explicit matrix/determinant counts."""
    perturbation = 3*field
    metric = 1+perturbation
    inverse = 1+sum(3**(n-1)*perturbation**n for n in range(1, ORDER+1))
    determinant_remainder = (6*metric**3).from_degree(1)
    volume = determinant_remainder.analytic(sp.Rational(1, 2))
    inverse_volume = determinant_remainder.analytic(sp.Rational(-1, 2))
    connection = Fraction(9, 2)*derivative_bound*inverse*perturbation
    curvature = 9*inverse*(6*derivative_bound*connection+18*connection**2)
    return {"metric": metric, "inverse": inverse, "volume": volume,
            "inverse_volume": inverse_volume, "connection": connection, "curvature": curvature}


def physical_series(seed, derivative_bound, inverse_squared_momentum_bound, *,
                    include_matter_source=True, york_order=3):
    """Inductive degree majorant for all M1 phase invariants through H4.

    The seed includes both mixed-boundary shifts.  The two optional omissions
    exist only for negative controls; build() uses the complete recursion.
    """
    if not 0 <= york_order <= 3:
        raise ValueError("Use zero through three York orders")
    field = Series.at(1, 4*Fraction(seed))
    K = Fraction(derivative_bound)
    inverse_k2 = Fraction(inverse_squared_momentum_bound)
    geo = geometry_majorants(field, K)
    pi = 8+18*field
    matter_density = Fraction(1, 5)+field
    corrections = []
    for degree in range(1, york_order+1):
        source = 3*K*pi+9*geo["connection"]*pi
        if include_matter_source:
            source += Fraction(3, 2)*K*matter_density*geo["inverse"]*field
        coefficient = 8*K*inverse_k2*source.coefficients[degree]
        pi += Series.at(degree, coefficient)
        corrections.append(str(coefficient))
    mixed = 3*pi*geo["metric"]*geo["inverse_volume"]
    sigma = (3*mixed).from_degree(1)
    shear2 = (12*mixed**2).from_degree(2)
    # Only exact zero background/linear terms are removed, as certified by S5.6.
    eta = (matter_density*geo["inverse_volume"]).from_degree(1)
    matter_gradient = (9*K*K*geo["inverse"]*field**2).from_degree(2)
    return {**geo, "momentum": pi, "matter_density": matter_density,
            "sigma": sigma, "rho": geo["curvature"], "eta": eta,
            "shear2": shear2, "matter_gradient": matter_gradient,
            "York_correction_majorants": corrections}


def scale_checks():
    """Homogeneous M,tau and fixed-patch power counting, before any bound."""
    mass, tau, ell0 = sp.symbols("M tau ell0", positive=True)
    h3, h4 = sp.symbols("H3 H4")
    coupling = mass*ell0
    checks = {
        "covariant_action_prefactor": sp.expand(tau**4*mass**2/tau**2-(mass*tau)**2),
        "canonical_matter_action_prefactor": sp.expand(tau**4*(mass/tau)**2-(mass*tau)**2),
        "cubic_canonical_power": sp.cancel(coupling**2*h3/coupling**3-h3/coupling),
        "quartic_canonical_power": sp.cancel(coupling**2*h4/coupling**4-h4/coupling**2),
        "physical_time_cancels_Hamiltonian_unit": sp.cancel(ell0*(h3/(ell0*coupling))-h3/coupling),
    }
    e = sp.Symbol("field_order")
    generic = sum(sp.Symbol(f"h{n}")*e**n for n in range(5))
    scaled = sp.expand(coupling**2*generic.subs(e, e/coupling))
    checks.update({f"homogeneous_order_{n}": sp.cancel(scaled.coeff(e, n)
                  -coupling**(2-n)*sp.Symbol(f"h{n}")) for n in range(5)})
    return checks


def negative_controls():
    """Actual nonzero matter subalgebra terms that vacuum-only reuse omits."""
    e, zeta, pm, matter_l, gradient = sp.symbols("e zeta P l gradient", real=True)
    density = matter_l+e*pm+3*matter_l*e*zeta
    inverse_volume = sp.series((1+2*e*zeta)**sp.Rational(-3, 2), e, 0, 3).removeO()
    eta = sp.expand(density*inverse_volume-matter_l)
    return {
        "matter_York_quadratic_source": -sp.I*pm*gradient/2,
        "matter_eta_nonlinear_second_order": sp.expand(eta.coeff(e, 2)),
        "matter_gradient_metric_cubic": -2*zeta*gradient**2,
        "metric_mixed_boundary_shift": 3*matter_l*sp.Symbol("s", real=True),
        "matter_mixed_boundary_shift": 3*matter_l*zeta,
    }


@cache
def build():
    domain = domain_checks()
    data = physical_series(Fraction(domain["common_terminal_phase_seed"]), DERIVATIVE_K,
                           1/Fraction(domain["proper_subset_squared_lower"]))
    density, entries = Series.at(0, 0), []
    invariants = (data["sigma"], data["rho"], data["eta"], data["shear2"], data["matter_gradient"])
    for order in invariant_series.coefficient_report()["hamiltonian"]:
        for key, proof in order.items():
            powers = tuple(map(int, key.split(",")))
            p, r, eta, shear, gradient = powers
            weight = p+2*r+eta+2*shear+2*gradient-2
            if not -2 <= weight <= 6:
                raise ValueError("Unproved fixed-unit coefficient scaling weight")
            bound = 64*Fraction(proof["uniform_absolute_upper"])
            term = Series.at(0, bound)
            for invariant, power in zip(invariants, powers):
                term *= invariant**power
            density += term
            entries.append({"powers": key, "compact_bound": proof["uniform_absolute_upper"],
                            "fixed_unit_bound": str(bound), "unit_weight": weight})
    # The -2H Pi:g boundary is nonlinear; H_fixed<=8 and 9 contractions give
    # 144.  Adot, -l*pi_chi, gamma and normalization generators have degree<=2.
    # A^3<=8 converts physical-volume density to fixed spatial coordinates.
    h = 8*(data["volume"]*density+144*data["momentum"]*data["metric"])
    B3, B4 = WICK_FACTOR*h.coefficients[3], WICK_FACTOR*h.coefficients[4]
    phase_measure = SPECIES*(4*UPPER_K+1)**3
    # Fixed total momentum: one-particle species measure=4; two-particle
    # relative-momentum/species measure<=4*phase_measure.  No delta is bounded.
    schur = 8*phase_measure**3
    if schur <= SPECIES+SPECIES*phase_measure:
        raise ValueError("Insufficient fixed-fiber Schur measure majorant")
    C3, C4 = schur*B3, schur*(B4+EXCHANGE_FACTOR*B3**2)
    scale = 1
    while scale < TARGET_DENOMINATOR*C3 or scale*scale < TARGET_DENOMINATOR*C4:
        scale *= 10
    return {"domain_checks": domain,
            "window": {"half_duration_over_ell0": str(DELTA), "external_k_min": str(LOWER_K),
                       "external_k_max": str(UPPER_K), "proper_subset_k_min": str(LOWER_K),
                       "Fourier_measure": "d^3k/(2*pi)^3; fixed total momentum fibers on R^3",
                       "propagating_columns": SPECIES, "momentum_species_measure_upper": str(phase_measure),
                       "Schur_measure_majorant": str(schur)},
            "stationary_coefficient_majorants": entries,
            "momentum_correction_majorants": data["York_correction_majorants"],
            "series_coefficients": {key: list(map(str, data[key].coefficients)) for key in
                                    ("rho", "sigma", "eta", "shear2", "matter_gradient", "momentum")},
            "Hamiltonian_series_coefficients": list(map(str, h.coefficients)),
            "cubic_kernel_bound": str(B3), "quartic_Hamiltonian_kernel_bound": str(B4),
            "exchange_factor": EXCHANGE_FACTOR, "Wick_factor": WICK_FACTOR,
            "cubic_transition_block_bound_numerator": str(C3),
            "quartic_connected_tree_block_bound_numerator": str(C4),
            "scale_restoration": {"action_prefactor": "(M*ell0)^2",
                                  "dimensionless_H3_factor": "1/(M*ell0)",
                                  "dimensionless_H4_factor": "1/(M*ell0)^2",
                                  "physical_Hamiltonian_extra_factor": "1/ell0",
                                  "ell0_lower_bound": "tau"},
            "sufficient_M_tau": str(scale), "sufficient_M_tau_power10": len(str(scale))-1,
            "cubic_block_bound_at_sufficient_scale": str(C3/scale),
            "quartic_tree_block_bound_at_sufficient_scale": str(C4/(scale*scale)),
            "target_block_norm": f"1/{TARGET_DENOMINATOR}",
            "not_a_necessary_scale_or_all_orders_cutoff": True}
