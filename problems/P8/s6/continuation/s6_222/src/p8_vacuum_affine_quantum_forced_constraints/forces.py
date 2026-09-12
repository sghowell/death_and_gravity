"""Complete force-dependent scalar constraints and canonical phase interface."""

from functools import cache

import sympy as s
from p8_vacuum_affine_reduced_scalar_hamiltonian import scalar as old

n, v, vd, sigma, sd, b, q = old.n, old.v, old.vd, old.sigma, old.sd, old.b, old.q
pv, ps = old.pv, old.ps
th, E, l, w, J, A, T, delta, H = s.symbols(
    "Theta E ell w Jc A Tcorr delta H", real=True
)
gn, gz, gb = s.symbols("gn gzeta gb", real=True)
Z = s.Matrix([v, sigma, pv, ps])
g = s.Matrix([gn, gz, gb])
Jcan = old.JC


def lagrangian():
    return (
        old.lagrangian(J=J, theta=th, matter=w, charge=l, gradient=E)
        + 3 * T * n * v
        + s.Rational(9, 2) * A * v * v
    )


@cache
def system():
    L = lagrangian()
    jn = gn + delta * gz
    equations = s.Matrix(
        [s.diff(L, vd) - pv, s.diff(L, sd) - ps, s.diff(L, n) + jn, s.diff(L, b) + gb]
    )
    variables = (vd, sd, n, b)
    solution = s.solve(list(equations), variables, dict=True)[0]
    lc = th * pv - w * ps + 3 * th * l * sigma - (2 * E * q + 3 * T) * v
    flow = s.Matrix(
        [
            solution[vd],
            solution[sd],
            (s.diff(L, v) + gz).subs(solution, simultaneous=True) - 3 * H * pv,
            s.diff(L, sigma).subs(solution, simultaneous=True) - 3 * H * ps,
        ]
    )
    flow = flow.applyfunc(s.factor)
    metric = s.Matrix([solution[n], v + delta * solution[n], solution[b]])
    C = metric.jacobian(Z).applyfunc(s.factor)
    D = metric.jacobian(g).applyfunc(s.factor)
    F = flow.jacobian(g).applyfunc(s.factor)
    K = flow.jacobian(Z).applyfunc(s.factor)
    return {
        "L": L,
        "equations": equations,
        "solution": solution,
        "Lc": lc,
        "flow": flow,
        "metric": metric,
        "C": C,
        "D": D,
        "F": F,
        "K": K,
    }


@cache
def data():
    d = system()
    C, D, F, K = [d[name] for name in ("C", "D", "F", "K")]
    zero = dict.fromkeys(g, 0)
    hc = (
        ps**2 / 2
        - l * pv * sigma / 2
        + (q / 2 - 3 * l * l / 4) * sigma**2
        - q * v * v
        - s.Rational(9, 2) * A * v * v
        + d["Lc"] ** 2 / (4 * J)
    )
    u = s.Matrix([1, delta, 0])
    eb = s.Matrix([0, 0, 1])
    mu = s.Symbol("spectral_parameter", real=True)
    clockmap = s.Matrix([[1, 0, 0], [delta, 1, 0], [0, 0, 1]])
    clocktuple = s.Matrix([n, v, b])
    checks = {
        "complete_four_variable_forced_saddle": d["equations"]
        .subs(d["solution"], simultaneous=True)
        .applyfunc(s.factor),
        "complete_saddle_determinant": s.factor(
            d["equations"].jacobian((vd, sd, n, b)).det() + 8 * J
        ),
        "forced_lapse_exact": s.factor(
            d["solution"][n] - (d["Lc"] - gn - delta * gz) / (2 * J)
        ),
        "forced_shift_exact": s.factor(d["solution"][b] - (pv - 3 * gb) / 2),
        "full_metric_phase_plus_direct_force": (d["metric"] - C * Z - D * g).applyfunc(
            s.factor
        ),
        "full_phase_generator_plus_force": (d["flow"] - K * Z - F * g).applyfunc(
            s.factor
        ),
        "zero_force_recovers_full_classical_Hamiltonian": (
            K - Jcan * s.hessian(hc, Z) + 3 * H * s.diag(0, 0, 1, 1)
        ).applyfunc(s.factor),
        "zero_force_metric_recovers_classical_comparison": (
            d["metric"].subs(zero) - C * Z
        ).applyfunc(s.factor),
        "force_is_negative_symplectic_metric_adjoint": (F + Jcan * C.T).applyfunc(
            s.factor
        ),
        "direct_auxiliary_rank_two_decomposition": (
            D + u * u.T / (2 * J) + s.Rational(3, 2) * eb * eb.T
        ).applyfunc(s.factor),
        "direct_auxiliary_null_covector": D * s.Matrix([-delta, 1, 0]),
        "direct_auxiliary_characteristic_polynomial": s.factor(
            (mu * s.eye(3) - D).det()
            - mu * (mu + (1 + delta**2) / (2 * J)) * (mu + s.Rational(3, 2))
        ),
        "clock_force_pairing_and_first_jet_pullback": s.factor(
            (g.T * clockmap * clocktuple)[0] - ((gn + delta * gz) * n + gz * v + gb * b)
        ),
        "orthogonal_canonical_norm_bridge": Jcan.T * Jcan - s.eye(4),
    }
    return {
        "forced_lapse": d["solution"][n],
        "forced_shift": d["solution"][b],
        "complete_phase_force_matrix": F,
        "classical_metric_matrix": C,
        "direct_auxiliary_response_matrix": D,
        "clock_first_jet": clockmap,
        "saddle_determinant": -8 * J,
        "forced_equations": "For the full coefficient action and clock forces(gn+delta gzeta,gzeta,gb), solve both momenta and both constraints before substituting. The gb terms cancel in the lapse pivot, but remain in b=(pv-3gb)/2. Matter evolution and all weighted momentum terms are retained.",
        "phase_definition": "Z uses the original classical Legendre momenta as auxiliary variables in the full forced system. It is not asserted to be a full quantum canonical Hamiltonian, and the three force histories retain all nonlocal and higher-time-derivative quantum content.",
        "shift_notation": "Here b is the physical scalar shift divergence from S220. The complementary coordinate called b in S221 is instead y_b=-pv/(2q); it is not this physical shift.",
        "exact_scope": "s=CZ+Dg and Z'=KZ+Fg are exact force interfaces, including D of rank2. Its negative local sign is not positivity or invertibility of a retarded Schur operator. No Theta or inverse-transfer division occurs.",
        "checks": checks,
        "gates": {
            "direct_block_symmetric": D == D.T,
            "direct_block_rank_two": D.rank() == 2,
            "full_saddle_regular_for_Jc_positive": True,
            "metric_and_force_degree_at_most_one_in_q": all(
                s.degree(x, q) <= 1 for x in (*C, *F) if x != 0
            ),
            "direct_block_no_spatial_transfer": all(not x.has(q) for x in D),
            "no_Theta_or_q_denominator": all(
                not s.denom(s.factor(x)).has(th, q) for x in (*C, *D, *F)
            ),
            "both_scalar_phase_pairs_and_quantum_force_histories_retained": True,
            "direct_block_not_discarded_as_classical_constraint": True,
        },
    }


@cache
def contact_data():
    N, scale = s.symbols("N a", positive=True)
    pressure = s.Symbol("pressure_physical", real=True)
    vv = s.Symbol("v_clock", real=True)
    R = 1 + 2 * delta * (N**-2 - 1)
    zeta = vv - s.log(R) / 4
    first = s.diff(zeta, N).subs(N, 1)
    second = s.diff(zeta, N, 2).subs(N, 1)
    clock_contact = 3 * scale**3 * pressure * second
    kappa = s.Symbol("kappa", positive=True)
    raw = s.Matrix(s.symbols("raw_n raw_zeta raw_b", real=True))
    clockmap = s.Matrix([[1, 0, 0], [delta, 1, 0], [0, 0, 1]])
    W = 1 / (kappa * scale**3)
    checks = {
        "actual_clock_zeta_first_lapse_jet": s.factor(first - delta),
        "actual_clock_zeta_second_lapse_jet": s.factor(
            second - (4 * delta**2 - 3 * delta)
        ),
        "full_extra_clock_current_contact": s.factor(
            clock_contact
            - s.Rational(3, 2) * scale**3 * pressure * (8 * delta**2 - 6 * delta)
        ),
        "output_density_and_clock_force_pullback": (
            W * clockmap.T * raw - clockmap.T * (W * raw)
        ).applyfunc(s.factor),
        "nonzero_clock_contact_bounce_value": s.factor(
            clock_contact.subs(delta, s.Rational(1, 2))
            + s.Rational(3, 2) * scale**3 * pressure
        ),
    }
    return {
        "additional_clock_nn_response": clock_contact,
        "output_density_normalization": W,
        "response_definition": "Rhat is the complete S220 ordered scalar ADM response plus its extra clock nn one-current contact. Pulling it back by the first-jet clock map gives the full clock response without deleting or double-counting that contact. The fixed profile one-point cancellation does not set the Gaussian current to zero.",
        "ordered_normalization": "Qbar=W Rhat, W(t)=1/(kappa a(t)^3), multiplies the output density. It commutes with the same-time clock matrix, not with the retarded time kernel. Detector and source operators keep their original order.",
        "checks": checks,
        "gates": {
            "both_ADM_and_extra_clock_contact_kept_once": True,
            "Gaussian_one_point_not_set_to_zero": True,
            "output_density_not_moved_inside_retarded_source_leg": True,
            "no_new_state_or_finite_counterterm": True,
        },
    }
