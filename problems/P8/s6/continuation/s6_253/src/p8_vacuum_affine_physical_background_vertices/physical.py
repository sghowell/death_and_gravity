"""Complete homogeneous physical probe jets in the fixed affine-clock chart.

Jets store the value, first derivative and SECOND DERIVATIVE, not Taylor
coefficients. The actual profile bindings are source-pinned in parent.py.
"""

from dataclasses import dataclass
from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase as previous_phase
from p8_vacuum_affine_scalar_tame_propagator import charts
from p8_vector_state import wkb

from . import coupled, parent

u = parent.u
a, zeta, kappa, heavy_mass2 = s.symbols(
    "positive_reference_scale positive_zeta positive_kappa positive_heavy_mass_squared",
    positive=True,
)
eta = s.symbols("physical_lapse_probe_0:3", real=True)
vp = s.symbols("physical_log_scale_probe_0:3", real=True)
psi = s.symbols("physical_matter_probe_0:3", real=True)


def dtime(expression):
    return (
        s.diff(expression, u)
        + parent.H * a * s.diff(expression, a)
        + sum(
            chain[j + 1] * s.diff(expression, chain[j])
            for chain in (eta, vp, psi)
            for j in range(2)
        )
    )


@dataclass(frozen=True)
class Jet:
    values: tuple

    @classmethod
    def constant(cls, value):
        return cls((s.sympify(value), s.S.Zero, s.S.Zero))

    @classmethod
    def wrap(cls, value):
        return value if isinstance(value, cls) else cls.constant(value)

    def __add__(self, other):
        other = self.wrap(other)
        return Jet(tuple(x + y for x, y in zip(self.values, other.values)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(tuple(-x for x in self.values))

    def __sub__(self, other):
        return self + -self.wrap(other)

    def __rsub__(self, other):
        return self.wrap(other) + -self

    def __mul__(self, other):
        y = self.wrap(other).values
        x = self.values
        return Jet(
            (
                x[0] * y[0],
                x[1] * y[0] + x[0] * y[1],
                x[2] * y[0] + 2 * x[1] * y[1] + x[0] * y[2],
            )
        )

    __rmul__ = __mul__

    def __pow__(self, power):
        power = s.sympify(power)
        if power == 0:
            return Jet.constant(1)
        if power == 1:
            return self
        x, xd, xdd = self.values
        return Jet(
            (
                x**power,
                power * x ** (power - 1) * xd,
                power * x ** (power - 1) * xdd
                + power * (power - 1) * x ** (power - 2) * xd**2,
            )
        )

    def __truediv__(self, other):
        return self * self.wrap(other) ** -1

    def __rtruediv__(self, other):
        return self.wrap(other) * self**-1

    def time_derivative(self):
        return Jet(tuple(dtime(x) for x in self.values))

    def factored(self):
        return Jet(tuple(s.factor(x) for x in self.values))


@cache
def coefficient_map():
    """All homogeneous lapse, scale and matter probe jets of the FULL block."""
    rows = parent.lapse_jets()
    functions = rows["rows"]
    lapse = Jet((s.S.One, eta[0], s.S.Zero))
    physical_scale = Jet((a, a * vp[0], a * vp[0] ** 2))
    rjets = rows["R_lapse_jets_zero_through_four"]
    R = Jet((s.S.One, rjets[1] * eta[0], rjets[2] * eta[0] ** 2))
    RN = Jet((rjets[1], rjets[2] * eta[0], rjets[3] * eta[0] ** 2))
    ahat = (physical_scale * R ** s.Rational(1, 4)).factored()
    Hhat = (ahat.time_derivative() / ahat).factored()
    matter_rate = Jet((parent.ell, psi[1], s.S.Zero))

    def f(name, derivative=0):
        row = functions[name]
        return Jet(
            (
                row[derivative],
                row[derivative + 1] * eta[0],
                row[derivative + 2] * eta[0] ** 2,
            )
        )

    D, Z = f("M") / lapse, f("U") / lapse
    DN = f("M", 1) / lapse - f("M") / lapse**2
    DNN = f("M", 2) / lapse - 2 * f("M", 1) / lapse**2 + 2 * f("M") / lapse**3
    ZN = f("U", 1) / lapse - f("U") / lapse**2
    ZNN = f("U", 2) / lapse - 2 * f("U", 1) / lapse**2 + 2 * f("U") / lapse**3
    B, BN, BNN = f("B"), f("B", 1), f("B", 2)
    F, FN, FNN = f("Fhat"), f("Fhat", 1), f("Fhat", 2)
    D, Z, DN, DNN, ZN, ZNN = [x.factored() for x in (D, Z, DN, DNN, ZN, ZNN)]
    L0 = -3 * D * Hhat**2 + 3 * B * Hhat + lapse * F + Z * matter_rate**2 / 2
    theta = (-Hhat * DN + BN / 2).factored()
    w = matter_rate * ZN
    Cnn = (
        -3 * Hhat**2 * DNN
        + 3 * Hhat * BNN
        + 2 * FN
        + lapse * FNN
        + ZNN * matter_rate**2 / 2
    ) / 2
    Qn = -3 * Hhat**2 * DN + 3 * Hhat * BN + F + lapse * FN + ZN * matter_rate**2 / 2
    cv = (-18 * D * Hhat + 9 * B).factored()
    c = (Z * matter_rate).factored()
    Vvv = s.Rational(9, 2) * L0 - (cv.time_derivative() + 3 * Hhat * cv) / 2
    Vvs = -3 * (c.time_derivative() + 3 * Hhat * c)
    qhat = coupled.P**2 / ahat**2
    params = {
        coupled.D: D,
        coupled.Z: Z,
        coupled.J: Cnn + 3 * theta**2 / D - w**2 / (2 * Z),
        coupled.th: theta,
        coupled.w: w,
        coupled.c: c,
        coupled.Lnv: 3 * Qn + 4 * (f("C3") + lapse * f("C3", 1)) * qhat,
        coupled.Vvv: Vvv,
        coupled.Vvs: Vvs,
        coupled.C: lapse * f("C3"),
        coupled.Y: lapse * f("Cchi"),
        coupled.K: zeta * f("Cchi") / (lapse * ahat**2),
        coupled.Yv: lapse * f("Cchi") / ahat**2,
        coupled.r: R - 1,
        coupled.rN: RN,
        coupled.dH: Hhat - parent.H,
        coupled.q: qhat,
        a: ahat,
    }
    return {
        "parameters": {key: value.factored().values for key, value in params.items()},
        "lapse": lapse.values,
        "physical_scale": physical_scale.values,
        "hat_scale": ahat.values,
        "physical_to_hat_Hubble": Hhat.values,
        "matter_rate": matter_rate.values,
        "R": R.values,
        "M": f("M").values,
        "C3": f("C3").values,
        "scalar_time_boundary_cv": cv.values,
        "scalar_time_boundary_matter_charge": c.values,
        "aligned_temporal_vector_background": (3 * (R - 1) * (Hhat - parent.H))
        .factored()
        .values,
        "parent_lapse_rows": functions,
        "convention": "N=1+epsilon eta; physical a=a_ref exp(epsilon v_phys); chi=chi_ref+epsilon psi; a_hat=R^(1/4) a_physical. Every probe and its needed derivatives vanishes near the fixed initial state slice. P is held comoving. Fixed profiles, H_clock(u), kappa, zeta and all masses are NOT varied.",
        "scope": "Complete homogeneous probe jets of this gauge-fixed quadratic background functional, in its stated fixed affine-clock fluctuation chart. Not arbitrary inhomogeneous covariant vertices, a nonlinear quantum field-map measure or a Ward-restored physical stress.",
    }


def coefficient_directional_matrices(hessian, parameters):
    """Entire matrix directional derivatives, prior to physical substitution."""
    first = {p: s.Symbol(str(p) + "_physical_first", real=True) for p in parameters}
    second = {p: s.Symbol(str(p) + "_physical_second", real=True) for p in parameters}
    d1 = sum((hessian.diff(p) * first[p] for p in parameters), s.zeros(*hessian.shape))
    d2 = sum((hessian.diff(p) * second[p] for p in parameters), s.zeros(*hessian.shape))
    for i, p in enumerate(parameters):
        d2 += hessian.diff(p, 2) * first[p] ** 2
        for r in parameters[i + 1 :]:
            d2 += 2 * hessian.diff(p, r) * first[p] * first[r]
    return d1, d2, first, second


@cache
def vertices():
    C = s.sqrt(kappa) * s.diag(1, 1, s.sqrt(zeta), a**3, a**3, a**3 / s.sqrt(zeta))
    H = (
        kappa
        * a**3
        * C.inv().T
        * s.hessian(coupled.reduced_hamiltonian(), coupled.PHASE)
        * C.inv()
    )
    parameters = (
        coupled.D,
        coupled.Z,
        coupled.J,
        coupled.K,
        coupled.th,
        coupled.w,
        coupled.c,
        coupled.Lnv,
        coupled.Vvv,
        coupled.Vvs,
        coupled.C,
        coupled.Y,
        coupled.Yv,
        coupled.r,
        coupled.rN,
        coupled.dH,
        coupled.q,
        a,
    )
    d1, d2, first, second = coefficient_directional_matrices(H, parameters)
    return {
        "whole_canonical_Hessian": H,
        "whole_first_directional_Hessian": d1,
        "whole_second_directional_Hessian": d2,
        "coefficient_first_symbols": first,
        "coefficient_second_symbols": second,
        "parameters": parameters,
        "physical_substitution_rule": "For each parameter p, substitute p, p_physical_first, p_physical_second simultaneously by the three derivative entries of coefficient_map.parameters[p]. This applies to EVERY entry of the displayed 6x6 Hessians, including q=P^2/a_hat^2, the canonical density, all off-reference Ward contacts and the complete source square. Cross directional probes follow by polarization of the entire second derivative.",
    }


@cache
def other_modes():
    mapping = coefficient_map()
    N = Jet(mapping["lapse"])
    ahat, ap = Jet(mapping["hat_scale"]), Jet(mapping["physical_scale"])
    M, C3 = Jet(mapping["M"]), Jet(mapping["C3"])
    blocks = {
        "each_of_two_TT": (2 * ahat * N * C3 * coupled.P**2, N / (ahat**3 * M)),
        "each_of_two_transverse_Proca": (N * (coupled.P**2 / ap + ap / zeta), N / ap),
        "one_heavy_scalar": (N * (heavy_mass2 * ap**3 + coupled.P**2 * ap), N / ap**3),
    }
    matrices = {
        name: tuple(
            s.diag(left.values[j], right.values[j]).applyfunc(s.factor)
            for j in range(3)
        )
        for name, (left, right) in blocks.items()
    }
    q = coupled.P**2 / a**2
    half_omega = s.Matrix([[0, 1], [-1, 0]])
    normalizers = {
        "longitudinal_Proca": a / (1 + zeta * q),
        "transverse_Proca": a,
        "heavy_scalar": a**3,
    }
    source_hessians = {
        "longitudinal_Proca": s.diag(a / zeta, (1 + zeta * q) / a),
        "transverse_Proca": matrices["each_of_two_transverse_Proca"][0],
        "heavy_scalar": matrices["one_heavy_scalar"][0],
    }
    transforms, checks = {}, {}
    for name, g in normalizers.items():
        rate = s.factor(dtime(g) / (2 * g))
        transform = s.Matrix([[1 / s.sqrt(g), 0], [-s.sqrt(g) * rate, s.sqrt(g)]])
        frequency2 = q + (heavy_mass2 if name == "heavy_scalar" else 1 / zeta)
        potential = s.factor(dtime(rate) + rate**2)
        oscillator_hessian = s.diag(frequency2 - potential, 1)
        transformed_generator = transform.inv() * (
            half_omega * source_hessians[name] * transform - transform.applyfunc(dtime)
        )
        checks[name + "_entire_fixed_state_canonical_map"] = (
            transform * half_omega * transform.T - half_omega
        )
        checks[name + "_entire_fixed_state_time_connection"] = (
            transformed_generator - half_omega * oscillator_hessian
        )
        if name != "heavy_scalar":
            old_kind = (
                "longitudinal" if name.startswith("longitudinal") else "transverse"
            )
            old_potential = wkb.frequency(old_kind)["U"].subs(wkb.z, q / (q + 1 / zeta))
            checks[name + "_unchanged_source_pinned_all_order_state_operator"] = (
                potential - old_potential
            )
        else:
            checks["heavy_unchanged_minimal_KG_operator"] = (
                potential - 3 * s.diff(parent.H, u) / 2 - 9 * parent.H**2 / 4
            )
        transforms[name] = {
            "oscillator_to_physical_canonical_map": transform,
            "normalizer_squared": g,
            "half_log_normalizer_rate": rate,
            "fixed_reference_oscillator_Hessian": oscillator_hessian,
        }
    checks["both_TT_reference_operators_unchanged"] = matrices["each_of_two_TT"][
        0
    ] - s.diag(a * coupled.P**2, a**-3)
    return {
        "whole_other_five_mode_Hessian_and_first_second_derivatives": matrices,
        "fixed_existing_state_canonical_transports": transforms,
        "same_eight_physical_modes": "Two coupled light scalars, three Proca polarizations, two TT gravitons, one heavy H scalar. The homogeneous probed Gaussian splits into one coupled three-mode block and five single modes. At the clock the three-mode block is the S251 two-scalar state times the existing S55/S176 longitudinal state; all other states and the S240 H state are unchanged.",
        "fixed_state_rule": "The mode maps act on oscillator (v,v_dot) data with their entire shear and time connection. A covariance already expressed in (v,v_dot-rate*v) must NOT receive the rate shear twice. Initial data are fixed and every probe vanishes near the common initial slice. This is a product free reference, not an interacting state or a new sampling prescription.",
        "heavy_source_rule": "The exact current H source starts at clock field degree8. Every background/fluctuation derivative used here has total degree at most4, so its source terms vanish, but the entire source-free massive H determinant and its metric vertices remain.",
        "checks": {
            name: value.applyfunc(s.factor)
            if isinstance(value, s.MatrixBase)
            else s.factor(value)
            for name, value in checks.items()
        },
        "gates": {
            "all_eight_physical_modes_included": 3 + 2 + 2 + 1 == 8,
            "no_fixed_state_or_mass_reset": True,
            "full_H_source_free_determinant_retained": True,
            "no_uniform_perturbed_high_frequency_hyperbolicity_inferred": True,
        },
    }


@cache
def data():
    mapping, packet = coefficient_map(), vertices()
    rows = mapping["parameters"]
    Rjets = parent.lapse_jets()["R_lapse_jets_zero_through_four"]
    r1, r2 = Rjets[1:3]
    logscale1 = vp[0] + r1 * eta[0] / 4
    logscale2 = (r2 - r1**2) * eta[0] ** 2 / 4
    reference = parent.reference_data()["whole_reference_coefficients"]
    E, ell = reference["E"], reference["ell"]
    q0 = coupled.P**2 / a**2
    reference_targets = {
        coupled.D: 1,
        coupled.Z: 1,
        coupled.J: reference["Jc"],
        coupled.K: zeta / a**2,
        coupled.th: reference["Theta"],
        coupled.w: -ell * E,
        coupled.c: ell,
        coupled.Lnv: 2 * E * q0 + 3 * reference["T"],
        coupled.Vvv: s.Rational(9, 2) * reference["A"],
        coupled.Vvs: 0,
        coupled.C: s.Rational(1, 2),
        coupled.Y: 1,
        coupled.Yv: a**-2,
        coupled.r: 0,
        coupled.rN: r1,
        coupled.dH: 0,
        coupled.q: q0,
        a: a,
    }
    checks = {
        "physical_hat_scale_first_chain": mapping["hat_scale"][1] - a * logscale1,
        "physical_hat_scale_second_chain": mapping["hat_scale"][2]
        - a * (logscale1**2 + logscale2),
        "physical_hat_Hubble_first_chain": mapping["physical_to_hat_Hubble"][1]
        - dtime(logscale1),
        "physical_hat_Hubble_second_chain": mapping["physical_to_hat_Hubble"][2]
        - dtime(logscale2),
        "physical_fixed_comoving_q_first_chain": rows[coupled.q][1]
        + 2 * q0 * logscale1,
        "physical_fixed_comoving_q_second_chain": rows[coupled.q][2]
        - q0 * (4 * logscale1**2 - 2 * logscale2),
        "physical_aligned_vector_first_embedding": mapping[
            "aligned_temporal_vector_background"
        ][1],
        "physical_aligned_vector_second_embedding": mapping[
            "aligned_temporal_vector_background"
        ][2]
        - 6 * r1 * eta[0] * dtime(logscale1),
    }
    checks.update(
        {
            "entire_physical_reference_parameter_" + str(p): rows[p][0] - target
            for p, target in reference_targets.items()
        }
    )
    for name in (
        "whole_canonical_Hessian",
        "whole_first_directional_Hessian",
        "whole_second_directional_Hessian",
    ):
        checks[name + "_symmetric"] = packet[name] - packet[name].T
    symbolic_ref = {
        coupled.D: 1,
        coupled.Z: 1,
        coupled.J: charts.J,
        coupled.K: zeta / a**2,
        coupled.th: charts.th,
        coupled.w: charts.w,
        coupled.c: charts.l,
        coupled.Lnv: 2 * charts.E * coupled.P**2 / a**2 + 3 * charts.T,
        coupled.Vvv: s.Rational(9, 2) * charts.A,
        coupled.Vvs: 0,
        coupled.C: s.Rational(1, 2),
        coupled.Y: 1,
        coupled.Yv: a**-2,
        coupled.r: 0,
        coupled.dH: 0,
        coupled.q: coupled.P**2 / a**2,
    }
    whole_ref = packet["whole_canonical_Hessian"].subs(symbolic_ref, simultaneous=True)
    index = [0, 1, 3, 4]
    scalar_ref = (-previous_phase.OMEGA * previous_phase.canonical_generator()).subs(
        {
            previous_phase.a: a,
            previous_phase.kappa: kappa,
            charts.q: coupled.P**2 / a**2,
        },
        simultaneous=True,
    )
    checks["entire_six_phase_reference_retains_whole_S251_scalar"] = (
        whole_ref.extract(index, index) - scalar_ref
    )
    checks["entire_six_phase_reference_retains_longitudinal_Proca"] = whole_ref.extract(
        [2, 5], [2, 5]
    ) - s.diag(a / zeta, (1 + zeta * coupled.P**2 / a**2) / a)
    checks["reference_product_has_no_spurious_light_vector_bilinear"] = (
        whole_ref.extract(index, [2, 5])
    )
    first_mixed = (
        s.diff(packet["whole_canonical_Hessian"][3, 5], coupled.r).subs(
            symbolic_ref, simultaneous=True
        )
        * r1
        * eta[0]
    )
    Cphase = s.sqrt(kappa) * s.diag(1, 1, s.sqrt(zeta), a**3, a**3, a**3 / s.sqrt(zeta))
    normal_hessian = coupled.held_vector_data()["whole_clock_normal_vector_Hessian"]
    normal_hessian = kappa * a**3 * Cphase.inv().T * normal_hessian * Cphase.inv()
    normal_substitution = {p: value[0] for p, value in rows.items()}
    normal_substitution[s.Symbol("Z_N", real=True)] = (
        parent.lapse_jets()["rows"]["U"][1] - 1
    )
    return {
        "whole_physical_homogeneous_probe_map": mapping,
        "whole_canonical_Hamiltonian_vertices": packet,
        "first_mixed_lapse_source_vertex_example": first_mixed,
        "whole_clock_normal_vector_canonical_Hessian": normal_hessian,
        "normal_vector_reference_coefficient_substitution": normal_substitution,
        "held_second_contact_to_subtract": normal_hessian
        * mapping["aligned_temporal_vector_background"][2],
        "complete_held_W0_second_vertex_rule": "Take the ENTIRE aligned second 6x6 vertex after the simultaneous physical substitution, and subtract the canonical Hessian of coupled.held_vector_data.whole_clock_normal_vector_Hamiltonian_vertex times the displayed actual S_second. First vertices coincide because S_first=0. No normal onepoint or measure contact is set to zero.",
        "canonical_boundary_rule": "The raw-to-normalized scalar action differs by d_t[kappa*a_hat^3*(cv*v^2/2+3*c*v*sigma)]. Restore this boundary and transport the same canonical state in any alternate chart. The a_hat^3 phase density and its full time connection are already in the displayed canonical Hessian. Probes vanish near the initial state slice; an endpoint phase cannot be used to reset the state.",
        "checks": {
            name: value.applyfunc(s.factor)
            if isinstance(value, s.MatrixBase)
            else s.factor(value)
            for name, value in checks.items()
        },
        "gates": {
            "all_eighteen_coefficient_and_density_jets_present": len(rows) == 18,
            "complete_first_and_second_six_phase_vertices_present": packet[
                "whole_first_directional_Hessian"
            ].shape
            == (6, 6)
            and packet["whole_second_directional_Hessian"].shape == (6, 6),
            "physical_lapse_induces_nonzero_light_vector_first_vertex": first_mixed
            != 0,
            "off_reference_matter_Ward_first_contact_nonzero": rows[coupled.Vvs][1]
            != 0,
            "physical_lapse_time_derivatives_retained": rows[coupled.Vvv][1].has(
                eta[2]
            ),
            "physical_second_embedding_contact_nonzero": mapping[
                "aligned_temporal_vector_background"
            ][2]
            != 0,
        },
    }


def omega(dimension):
    if type(dimension) is not int or dimension not in (2, 4, 6, 8, 16):
        raise ValueError("Require an exact supported even physical phase dimension")
    half = dimension // 2
    return (
        s.zeros(half)
        .row_join(s.eye(half))
        .col_join((-s.eye(half)).row_join(s.zeros(half)))
    )


def exact_matrix(value, dimension=None, symmetric=False):
    matrix = s.Matrix(value)
    omega(matrix.rows)
    if (
        matrix.rows != matrix.cols
        or dimension is not None
        and matrix.shape != (dimension, dimension)
    ):
        raise ValueError(
            "Require a complete square phase matrix of the stated dimension"
        )
    if any(
        x.has(s.Float)
        or x.free_symbols
        or x.is_real is not True
        or x.is_finite is not True
        for x in matrix
    ):
        raise ValueError("Require fixed provably finite real exact entries")
    if symmetric and matrix != matrix.T:
        raise ValueError("Require an entire symmetric Weyl coefficient")
    return matrix


def pure_covariance(value):
    matrix = exact_matrix(value, symmetric=True)
    if any(
        s.factor(matrix[:j, :j].det()).is_positive is not True
        for j in range(1, matrix.rows + 1)
    ):
        raise ValueError("Require a strictly positive covariance")
    O = omega(matrix.rows)
    if (matrix * O * matrix - O / 4).applyfunc(s.cancel) != s.zeros(matrix.rows):
        raise ValueError("Require the full pure-state canonical uncertainty equality")
    return matrix


def symplectic(value, dimension):
    matrix = exact_matrix(value, dimension)
    O = omega(dimension)
    if (matrix * O * matrix.T - O).applyfunc(s.cancel) != s.zeros(dimension):
        raise ValueError("Require the complete symplectic state transport")
    return matrix


def assemble_product(covariances):
    if not isinstance(covariances, (tuple, list)) or len(covariances) != 7:
        raise ValueError(
            "Require the complete two-scalar block and six single-mode blocks"
        )
    dimensions = (4, 2, 2, 2, 2, 2, 2)
    blocks = [
        pure_covariance(exact_matrix(value, d))
        for value, d in zip(covariances, dimensions)
    ]
    # Input: scalar(v,s,pv,ps), L, TT1, TT2, PT1, PT2, H. Output all Q then all P.
    order = (0, 1, 4, 6, 8, 10, 12, 14, 2, 3, 5, 7, 9, 11, 13, 15)
    return pure_covariance(s.diag(*blocks).extract(order, order))


def kernels(A, B, left, right, initial, contact=None):
    V = pure_covariance(initial)
    dim = V.rows
    A, B = exact_matrix(A, dim, True), exact_matrix(B, dim, True)
    L, R = symplectic(left, dim), symplectic(right, dim)
    C = s.zeros(dim) if contact is None else exact_matrix(contact, dim, True)
    O = omega(dim)
    W = L * (V + s.I * O / 2) * R.T
    VA, VB = L * V * L.T, R * V * R.T
    PA, PB = L.T * A * L, R.T * B * R
    connected = s.expand(s.trace(A * W * B * W.T) / 2)
    chi = s.cancel(s.trace((PA * O * PB - PB * O * PA) * V) / 2)
    return {
        "mean_A": s.cancel(s.trace(A * VA) / 2),
        "mean_B": s.cancel(s.trace(B * VB) / 2),
        "whole_unequal_time_W": W,
        "ordered_connected_Wick": connected,
        "symmetric_noise": s.re(connected),
        "retarded_observable_before_step": chi,
        "effective_retarded_before_step": -chi,
        "effective_second_contact": s.cancel(-s.trace(C * VA) / 2),
    }


def symplectic_fixture(case, dim):
    if type(case) is not int or not 0 <= case < 12:
        raise ValueError("Unknown finite fixture")
    omega(dim)
    n = dim // 2
    linear = s.eye(n)
    linear[0, 1 if n > 1 else 0] += s.Rational(case + 1, 17)
    lower, upper = s.eye(dim), s.eye(dim)
    lower[n:, :n] = s.diag(*(s.Rational(case + j + 1, 19 + j) for j in range(n)))
    upper[:n, n:] = s.diag(*(s.Rational(case + 2 * j + 1, 29 + j) for j in range(n)))
    return lower * s.diag(linear, linear.inv().T) * upper


@cache
def response_data():
    checks = {}
    for dim in (6, 16):
        O = omega(dim)
        prepared = symplectic_fixture(0, dim)
        V = prepared * prepared.T / 2
        L, R = symplectic_fixture(1, dim), symplectic_fixture(2, dim)
        A = s.diag(*(s.Rational(j + 1, j + 3) for j in range(dim)))
        B = s.diag(*(s.Rational(j + 3, j + 7) for j in range(dim)))
        A[0, dim - 1] = A[dim - 1, 0] = s.Rational(1, 13)
        B[1, dim // 2] = B[dim // 2, 1] = -s.Rational(1, 11)
        contact = A + B
        packet = kernels(A, B, L, R, V, contact)
        reverse = kernels(B, A, R, L, V)
        variation = L * R.inv() * O * B * R
        dV = variation * V * L.T + L * V * variation.T
        checks[str(dim) + "_whole_Wick_commutator"] = (
            2 * s.im(packet["ordered_connected_Wick"])
            - packet["retarded_observable_before_step"]
        )
        checks[str(dim) + "_whole_Duhamel_kick_response"] = (
            s.trace(A * dV) / 2 - packet["retarded_observable_before_step"]
        )
        checks[str(dim) + "_whole_reverse_contraction"] = reverse[
            "ordered_connected_Wick"
        ] - s.conjugate(packet["ordered_connected_Wick"])
        checks[str(dim) + "_whole_fixed_state_CCR"] = (
            packet["whole_unequal_time_W"]
            - s.conjugate(packet["whole_unequal_time_W"])
            - s.I * L * O * R.T
        )
        checks[str(dim) + "_entire_second_vertex_contact"] = (
            packet["effective_second_contact"] + s.trace(contact * L * V * L.T) / 2
        )
    product = assemble_product([s.eye(4) / 2] + [s.eye(2) / 2 for _ in range(6)])
    checks["whole_eight_mode_product_normalization"] = product - s.eye(16) / 2
    return {
        "entire_finite_CTP_and_Wick_rule": "Use the complete S252 fixed-state Gaussian trace with the full16-dimensional relative symplectic evolution and its continued metaplectic logarithm. At the reference the covariance is the seven-block product, reordered into eight Q and eight P. The homogeneous perturbed evolution has a coupled6-dimensional block plus five2-dimensional blocks; their common-regulator metaplectic traces multiply with continuous phases. Pure centered covariance is not a zero-point-force subtraction.",
        "complete_physical_probe_rule": "Insert every entry of physical.data's whole first and second Hessians, including the fixed profile jets, both tensor modes, all Proca modes and H. For a held-W probe use the complete embedding contact. The finite Gaussian current is minus Tr(H_first V)/2, the seagull minus Tr(H_second V)/2, and the retarded connected kernel minus theta times the displayed susceptibility. The full Wick contraction uses W transpose, never its adjoint.",
        "continuum_boundary": "The unchanged reference is the tensor product of the existing positive states with the common time orientation. Its finite-jet Wick products and mixed-mode noise use the same S252 contraction-graph argument. This does not select a physical covariant counterfunctional, extend every arbitrary off-reference background, prove a uniform nonlinear high-frequency neighborhood, or evaluate a renormalized mean or omitted-loop norm.",
        "formal_grading": "The two previously fixed stress profiles are coefficients of formal loop grade1, held fixed under every probe. The complete physical jets contain their variations through the coefficient functions, not derivatives of live state means. Expand the whole current inverse lapse pivot and determinant with their matching countervertices at one common regulator. Source-square mixed vertices belong to the same full light/vector loop functional; a separate source-free Proca determinant plus a light determinant omits them.",
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "whole_six_and_sixteen_dimensional_interfaces_supported": omega(6).shape
            == (6, 6)
            and omega(16).shape == (16, 16),
            "whole_product_same_reference_preparations_required": True,
            "continued_CTP_phase_and_common_regulator_required": True,
            "physical_loop_mean_and_omitted_loop_norm_not_established": True,
        },
    }
