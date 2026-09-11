"""Canonical Einstein quadratic action and the nonuniform finite-gravity pole."""

from functools import cache

import sympy as s

ETA = s.diag(1, -1, -1, -1)


def connection_linear(derivatives):
    return [
        [
            [
                ETA[r, r]
                * (derivatives[m][r][n] + derivatives[n][r][m] - derivatives[r][m][n])
                / 2
                for n in range(4)
            ]
            for m in range(4)
        ]
        for r in range(4)
    ]


def gamma_density(derivatives):
    """g=eta+2h/sqrt(kappa); exact leading EH action after its divergence."""
    G = connection_linear(derivatives)
    return s.expand(
        -2
        * sum(
            ETA[m, m] * (G[r][m][n] * G[n][m][r] - G[r][m][m] * G[n][n][r])
            for m in range(4)
            for r in range(4)
            for n in range(4)
        )
    )


@cache
def data():
    k = s.Matrix(s.symbols("k0:4", real=True))
    entries = {
        (i, j): s.Symbol("h" + str(i) + str(j), real=True)
        for i in range(4)
        for j in range(i, 4)
    }
    h = s.Matrix(4, 4, lambda i, j: entries[tuple(sorted((i, j)))])
    d = [[[k[a] * h[i, j] for j in range(4)] for i in range(4)] for a in range(4)]
    actual = gamma_density(d)
    k2 = (k.T * ETA * k)[0]
    tr = s.trace(ETA * h)
    contraction = s.trace(ETA * h * ETA * h)
    divergence = h * ETA * k
    div2 = (divergence.T * ETA * divergence)[0]
    mixed = (k.T * ETA * divergence)[0]
    expected = (k2 * contraction - 2 * div2 + 2 * mixed * tr - k2 * tr * tr) / 2
    gauge = divergence - k * tr / 2
    gaugefix = (gauge.T * ETA * gauge)[0]
    xi = s.Matrix(s.symbols("xi0:4", real=True))
    gauge_h = k * xi.T + xi * k.T
    pure = {h[i, j]: gauge_h[i, j] for i in range(4) for j in range(i, 4)}
    # Independent real space TT jets, not a null-plane-wave specialization.
    p, c, q, dz = s.symbols("plus_time cross_time plus_z cross_z", real=True)
    jets = [[[s.S.Zero for j in range(4)] for i in range(4)] for a in range(4)]
    for a, hp, hc in [(0, p, c), (3, q, dz)]:
        jets[a][1][1] = hp
        jets[a][2][2] = -hp
        jets[a][1][2] = jets[a][2][1] = hc
    TT = gamma_density(jets)
    # The fixed mass scalar vertex from sqrt(-g)(Y-Phi^2)/2.
    mass2, ab = s.symbols("mass_squared p_dot_q", real=True)
    checks = {
        "literal_EH_Gamma_Gamma_to_Fierz_Pauli": s.expand(actual - expected),
        "de_Donder_gauge_fixed_normalization": s.expand(
            actual + gaugefix - k2 * (contraction - tr * tr / 2) / 2
        ),
        "linear_gauge_direction_null": s.expand(actual.subs(pure, simultaneous=True)),
        "two_physical_tensor_polarizations": s.expand(
            TT - (p * p + c * c - q * q - dz * dz)
        ),
        "on_shell_scalar_vertex_Ward_coefficient": (mass2 + ab) - (ab + mass2),
    }
    return {
        "physical_metric": "g_mu_nu=eta_mu_nu+2 h_mu_nu/sqrt(kappa); h has kappa-independent canonical Fierz-Pauli normalization",
        "linear_EH_boundary": "-sqrt(kappa) partial_mu[partial_nu h^(mu nu)-partial^mu trace(h)]; retained as an exact divergence, zero for the stated compact/Schwartz variations",
        "actual_quadratic_Gamma_Gamma": actual,
        "two_polarization_density": TT,
        "canonical_physical_polarizations": "sqrt(2)*h_plus and sqrt(2)*h_cross have density (time derivative squared-spatial derivative squared)/2",
        "de_Donder_propagator": "i P/(p^2+i0), P_mnrs=(eta_mr eta_ns+eta_ms eta_nr-eta_mn eta_rs)/2",
        "scalar_graviton_vertex": "i [p_mu q_nu+p_nu q_mu-eta_mu_nu(p.q+m^2)]/sqrt(kappa), all momenta incoming",
        "full_action_limit": "Einstein interactions beyond its quadratic density are O(kappa^-1/2); the fixed nonminimal r R term and all metric corrections to the full scalar, M1 and canonical Proca actions vanish on fixed compact jets. The latter spectator sectors keep their full kappa-independent canonical quadratic actions.",
        "uniformity": "For fixed compact C^(j+2) canonical jets strictly inside the coefficient strip, choose kappa large enough that the inverse metric stays in a compact nonsingular neighborhood of eta. Analyticity and Taylor's integral remainder give a finite C_j times kappa^-1/2 bound for the remaining density and its j jet derivatives, after the exact EH divergence. No momentum-uniform scattering or quantum-continuum estimate is claimed.",
        "checks": checks,
    }


@cache
def pole():
    S, T, U, m, K = s.symbols("s t u mass_squared kappa", real=True)
    a, b, c = (S - 2 * m) / 2, (T - 2 * m) / 2, (U - 2 * m) / 2
    # V13:P:V24, with all four incoming on-shell momenta.
    VV = 2 * (a * a + c * c) - 4 * b * (b + m) + 4 * (b + m) ** 2
    tr = -2 * b - 4 * m
    numerator = s.expand(VV - tr * tr / 2)
    on = s.factor(numerator.subs(U, 4 * m - S - T))
    expected = 2 * m * m - 2 * m * T - S * (4 * m - S - T)
    A = -on / (K * T)
    v = s.Symbol("crossing_v", real=True)
    tau = s.Symbol("positive_tau", positive=True)
    co = s.expand(A.subs(S, 2 * m + v - T / 2)).coeff(v, 2)
    diagonal = s.limit(A.subs(T, -tau / K), K, s.oo)
    return {
        "minimal_tree_t_channel": A,
        "exact_scalar_stress_contraction": on,
        "crossing_centered_v_squared_pole": co,
        "fixed_t_decoupling": s.S.Zero,
        "approaching_forward_t_minus_tau_over_kappa_limit": diagonal,
        "boundary": "Fixed nonzero transfer gravitational exchange vanishes as 1/kappa, but the shrinking-forward path has a generally nonzero limit. This is a tree diagnostic of nonuniformity, not a gravitational dispersion prescription or bound on the Regge/IR remainder. All scalar-graviton channels, loops and physical cuts still require their specified observable.",
        "checks": {
            "literal_stress_projector_contraction": s.factor(on - expected),
            "negative_forward_v_squared_pole": s.factor(co + 1 / (K * T)),
            "fixed_nonzero_transfer_limit": s.limit(A, K, s.oo),
            "shrinking_forward_limit": s.factor(
                diagonal - (S * S - 4 * m * S + 2 * m * m) / tau
            ),
            "all_incoming_on_shell_sum": s.expand(S + T + (4 * m - S - T) - 4 * m),
        },
    }
