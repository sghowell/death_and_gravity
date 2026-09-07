"""Fixed-pulse own-f endpoint separation without assuming an outer limit.

The actual central remainder need not converge at its shrinking corner.
Uniform finite-window comparison is sufficient to separate two exact phase
sequences.  The source is already loaded from zero data by loading.py.
"""

from functools import cache

import sympy as sp
from p8_own_scattering import connection

from . import loading

RHO = connection.RHO


def _exact_real(value, name, *, nonnegative=False):
    if isinstance(value, bool) or value is sp.true or value is sp.false:
        raise TypeError(f"{name} excludes booleans")
    value = sp.sympify(value)
    if not isinstance(value, sp.Expr):
        raise TypeError(f"{name} must be scalar")
    if value.has(sp.Float) or value.is_finite is not True or value.is_real is not True:
        raise ValueError(f"{name} must be exact, finite and real")
    if nonnegative and value.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    return value


def phase_value(delta):
    delta = _exact_real(delta, "delta")
    if delta.is_positive is not True:
        raise ValueError("delta must be positive; zero is not a literal action")
    a = sp.Rational(loading.A)
    return RHO * sp.asinh(sp.sqrt(8) * a / sp.sqrt(delta))


def phase_sequences(n, offset=0):
    """Exact phase-separated deltas; n>=2 guarantees the certified range.

    Any nonnegative fixed offset is admitted. For the original-Q readout,
    choose offset=arg(A*(w1+i*w2))/2 modulo pi in [0,pi). This depends on
    the one fixed loaded limit, not on delta. No numerical phase fitting is
    required for the existence and separation theorem.
    """
    if isinstance(n, bool) or not isinstance(n, (int, sp.Integer)):
        raise TypeError("sequence index must be an actual integer")
    n = int(n)
    if n < 2:
        raise ValueError("n>=2 is the stated common delta-domain guard")
    offset = _exact_real(offset, "offset", nonnegative=True)
    theta_plus = offset + sp.pi * n
    theta_minus = theta_plus + sp.pi / 2
    numerator = 8 * sp.Rational(loading.A)**2
    return {
        "n": n, "offset": offset,
        "theta_plus": theta_plus, "theta_minus": theta_minus,
        "delta_plus": numerator / sp.sinh(theta_plus / RHO)**2,
        "delta_minus": numerator / sp.sinh(theta_minus / RHO)**2,
    }


def reference_cauchy(theta):
    """Reference Cauchy map V0 D(theta) T_t D(theta) V0^*, not T_t alone."""
    theta = _exact_real(theta, "theta")
    v0 = connection.plane_wave_frame(0)
    rotation = sp.diag(sp.exp(sp.I * theta), sp.exp(-sp.I * theta))
    return v0 * rotation * connection.time_transfer() * rotation * v0.conjugate().T


@cache
def generic_reference():
    ar, ai, beta, theta = sp.symbols("A_real A_imag beta theta", real=True)
    w1, w2 = sp.symbols("w1 w2", real=True)
    aa = ar + sp.I * ai
    a_theta = ar * sp.cos(2 * theta) + ai * sp.sin(2 * theta)
    b_theta = ai * sp.cos(2 * theta) - ar * sp.sin(2 * theta)
    oscillatory = sp.Matrix([[a_theta, -b_theta], [b_theta, a_theta]])
    constant = sp.Matrix([[0, -beta], [-beta, 0]])
    return {
        "ar": ar, "ai": ai, "beta": beta, "theta": theta,
        "A": aa, "w1": w1, "w2": w2, "w": sp.Matrix([w1, w2]),
        "oscillatory": oscillatory, "constant": constant,
        "reference": oscillatory + constant,
        "alignment_complex": aa * (w1 + sp.I * w2),
        "alignment_squared_norm": (ar**2 + ai**2) * (w1**2 + w2**2),
    }


def endpoint_limit(k0, k_u0, *, side=1):
    """Limiting actual map at u=+/-a, not an action at delta=0.

    Inputs are the punctured coefficient limits there.  Q_u=tau*Q_T.
    The first row implies Q(+a)=sqrt(r_delta/k_delta)*Z_out,1.
    """
    if type(side) is not int or side not in (-1, 1):
        raise ValueError("side must be the actual integer -1 or +1")
    k0 = _exact_real(k0, "k0")
    k_u0 = _exact_real(k_u0, "k_u0")
    if k0.is_positive is not True:
        raise ValueError("the actual kinetic limit must be positive")
    a = sp.Rational(loading.A)
    whitening = sp.Matrix([[sp.sqrt(k0), 0], [k_u0 / (2 * sp.sqrt(k0)), sp.sqrt(k0)]])
    radial = sp.Matrix([[1 / sp.sqrt(a), 0], [-side / (2 * RHO * sp.sqrt(a)), sp.sqrt(a) / RHO]])
    return radial * whitening


def calibration():
    f = loading.rational
    a = loading.A
    central = loading.REMAINDER_CAP * (a**2 + 2 * a * f(loading.L / 20))
    # L/20=1/2000 is a strict upper bound for sqrt(delta/8).
    tails = 3 * loading.DELTA_MAX / (32 * a**2)
    budget = (central + tails) * f(sp.Rational(4, 5))
    error = f(sp.Rational(4, 399))
    separation = 2 - 2 * error
    scalar_gap = separation * f(sp.Rational(9, 1280)) * f(sp.Rational(2, 5))
    sequence_delta_upper = 8 * a**2 / f(sp.Rational(315, 16))**2
    return {
        "a": a, "delta_max": loading.DELTA_MAX,
        "remainder_cap": loading.REMAINDER_CAP,
        "central_integral_upper": central,
        "tail_integral_upper": tails,
        "perturbation_budget_upper": budget,
        "budget_gate": f(sp.Rational(1, 200)),
        "transfer_error_upper": error,
        "normalized_pair_separation_factor": separation,
        "scalar_Q_gap_lower_over_eta": scalar_gap,
        "scalar_Q_gap_gate_over_eta": f(sp.Rational(1, 200)),
        "sequence_delta_upper": sequence_delta_upper,
        "sequence_index_min": 2,
        "pi_lower_from_finite_integral": f(sp.Rational(135904, 45045)),
        "reference_A_abs_squared": connection.coefficients()["A_abs_squared"],
        "endpoint_vector_scope": "actual (Q,Q_u), with exact canonical/asinh maps retained",
        "normalized_separation": "liminf ||Z(delta_n+)-Z(delta_n-)|| > (790/399)||w|| for the paired phases",
        "scalar_readout": "limsup_delta->0 Q(+a)-liminf_delta->0 Q(+a) > eta/200",
        "subsequential_reference_limits_asserted": False,
        "central_remainder_convergence_assumed": False,
        "matter_source_or_EFT_verdict": False,
    }


@cache
def identities():
    data = generic_reference()
    ar, ai, beta, theta = (data[name] for name in ("ar", "ai", "beta", "theta"))
    aa = data["A"]
    v0 = connection.plane_wave_frame(0)
    rotation = sp.diag(sp.exp(sp.I * theta), sp.exp(-sp.I * theta))
    tt = sp.Matrix([[sp.conjugate(aa), sp.I * beta], [-sp.I * beta, aa]])
    direct = v0 * rotation * tt * rotation * v0.conjugate().T
    output = {}

    def add_matrix(name, matrix):
        for row in range(matrix.rows):
            for col in range(matrix.cols):
                output[f"{name}_{row}{col}"] = sp.simplify(sp.expand_complex(matrix[row, col]))

    add_matrix("direct_real_reference", direct - data["reference"])
    shifted = data["reference"].subs(theta, theta + sp.pi / 2)
    add_matrix("quarter_phase_changes_only_A", data["reference"] - shifted - 2 * data["oscillatory"])
    add_matrix("oscillatory_norm", data["oscillatory"].T * data["oscillatory"] - (ar**2 + ai**2) * sp.eye(2))
    z = data["alignment_complex"]
    output["alignment_modulus"] = sp.expand(z * sp.conjugate(z) - data["alignment_squared_norm"])
    output["alignment_rotated_imaginary_cleared"] = sp.simplify(sp.im(z * sp.conjugate(z)))
    k0 = sp.Symbol("k0", positive=True)
    ku0 = sp.Symbol("k_u0", real=True)
    for side in (-1, 1):
        matrix = endpoint_limit(k0, ku0, side=side)
        output[f"endpoint_limit_determinant_{side}"] = sp.simplify(matrix.det() - k0 / RHO)
        output[f"Q_readout_{side}"] = sp.simplify(matrix[0, 0] * sp.sqrt(sp.Rational(loading.A) / k0) - 1)
        output[f"Q_readout_velocity_{side}"] = matrix[0, 1]
    x = sp.Symbol("x", real=True)
    partial = sum((-1)**n * x**(2 * n) for n in range(8))
    output["pi_lower_positive_remainder"] = sp.expand(1 - (1 + x**2) * partial - x**16)
    output["pi_lower_integral"] = 4 * sp.integrate(partial, (x, 0, 1)) - sp.Rational(135904, 45045)
    theta_positive = sp.Symbol("theta_positive", positive=True)
    delta_sequence = 8 * sp.Rational(loading.A)**2 / sp.sinh(theta_positive / RHO)**2
    output["sequence_positive_sqrt"] = sp.simplify(
        sp.sqrt(8) * sp.Rational(loading.A) / sp.sqrt(delta_sequence) - sp.sinh(theta_positive / RHO)
    )
    inverse_argument = theta_positive / RHO
    # For a positive real argument, cosh>0 and this squared identity fix
    # sqrt(1+sinh^2)=cosh without imposing a complex branch rewrite.
    output["sequence_positive_cosh_square"] = sp.simplify(
        sp.cosh(inverse_argument)**2 - sp.sinh(inverse_argument)**2 - 1
    )
    inverse_log = sp.log(sp.sinh(inverse_argument) + sp.cosh(inverse_argument))
    output["sequence_phase_identity"] = sp.simplify(RHO * inverse_log.rewrite(sp.exp) - theta_positive)
    return output


def checks():
    data = calibration()
    f = loading.rational
    return {
        "half_width_inside_branch": 0 < loading.A < loading.L,
        "same_canonical_remainder_cap": data["remainder_cap"] == 44,
        "central_integral_value": data["central_integral_upper"] == f(sp.Rational(33, 25000)),
        "tail_integral_value": data["tail_integral_upper"] == f(sp.Rational(3, 800)),
        "budget_value": data["perturbation_budget_upper"] == f(sp.Rational(507, 125000)),
        "budget_below_gate": data["perturbation_budget_upper"] < data["budget_gate"],
        "error_from_geometric_exponential": 4 * (1 / (1 - data["budget_gate"] / 2) - 1) == data["transfer_error_upper"],
        "pair_separation_value": data["normalized_pair_separation_factor"] == f(sp.Rational(790, 399)),
        "pair_separation_positive": data["normalized_pair_separation_factor"] > 0,
        "scalar_gap_value": data["scalar_Q_gap_lower_over_eta"] == f(sp.Rational(237, 42560)),
        "scalar_gap_exceeds_gate": data["scalar_Q_gap_lower_over_eta"] > data["scalar_Q_gap_gate_over_eta"],
        "reciprocal_kinetic_sqrt_margin": loading.source_bounds()["kinetic_upper"] < f(sp.Rational(25, 4)),
        "pi_lower_above_three": data["pi_lower_from_finite_integral"] > 3,
        "rho_below_four_thirds": f(sp.Rational(7, 4)) < f(sp.Rational(16, 9)),
        "sinh_lower_at_n_two": f(sp.Rational(9, 2)) + f(sp.Rational(9, 2))**3 / 6 == f(sp.Rational(315, 16)),
        "real_inverse_cosh_positive": sp.cosh(sp.Symbol("positive_argument", positive=True)).is_positive is True,
        "sequence_delta_in_domain": 0 < data["sequence_delta_upper"] < loading.DELTA_MAX,
    }
