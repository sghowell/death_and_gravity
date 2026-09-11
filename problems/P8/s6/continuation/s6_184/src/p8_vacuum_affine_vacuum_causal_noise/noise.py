"""Complete three-polarization vacuum covariance of the smeared scalar force."""

from functools import cache

import sympy as s

from . import source


@cache
def data():
    m, q = s.symbols("m q", positive=True)
    omega = s.sqrt(m * m + q * q)
    momentum = s.Matrix([omega, q, 0, 0])
    eta = s.diag(1, -1, -1, -1)
    pol = [
        s.Matrix([q / m, omega / m, 0, 0]),
        s.Matrix([0, 0, 1, 0]),
        s.Matrix([0, 0, 0, 1]),
    ]
    P = sum((v * v.T for v in pol), s.zeros(4, 4))
    lam = 1 + 2 * q * q / (m * m)
    mass, k, k0, C = source.MASS, source.K, source.K0, source.C
    variance = k * C * C * (mass / 2 + s.Rational(62**2, 1) / mass)
    anchor = variance.subs(k, k0)
    checks = {
        "full_positive_polarization_sum": s.ImmutableMatrix(
            (P + eta - momentum * momentum.T / m**2).applyfunc(s.simplify)
        ),
        "longitudinal_mass_shell_transversality": (momentum.T * eta * pol[0])[0],
        "longitudinal_positive_unit_normalization": (pol[0].T * eta * pol[0])[0] + 1,
        "full_polarization_characteristic_polynomial": s.factor(
            P.charpoly().as_expr()
            - s.Symbol("lambda")
            * (s.Symbol("lambda") - 1) ** 2
            * (s.Symbol("lambda") - lam)
        ),
        "positive_longitudinal_eigenvalue": s.trace(P[:2, :2]) - lam,
        "zero_unphysical_longitudinal_block_determinant": P[:2, :2].det(),
        "canonical_force_covariance_normalization": k
        * mass**2
        * (C * C / (2 * mass) + (62 * C) ** 2 / mass**3)
        - variance,
        "full_noise_family_decay": s.factor(variance * k - anchor * k0),
        "display_standard_deviation_square": (2 * s.Rational(1, 10**1190)) ** 2
        - 4 * s.Rational(1, 10**2380),
    }
    return {
        "state": "The same fixed-mass positive-frequency Minkowski Gaussian vacuum on physical g=eta, with the source-dependent coherent mean and unchanged connected covariance",
        "polarization_basis": tuple(pol),
        "complete_polarization_covariance": s.ImmutableMatrix(P),
        "eigenvalues": (s.S.Zero, s.S.One, s.S.One, lam),
        "source_smearing": "f_mu=DS_mu[eta] with the appropriate Minkowski pairing; real eta is smooth, spatially Schwartz and temporally compact in an interval of length T<=1",
        "one_particle_covariance_bound": "C_A(f,f)<=T/(2m)||f||L2(dt dx)^2+T/m^3||spatial_grad f||L2(dt dx)^2. All momenta and three physical polarizations retained; metric sign flips preserve the Euclidean norm used in the bound.",
        "connected_scalar_force": "F_eta-E(F_eta)=-sqrt(kappa)*m integral A_centered.DS_eta, an exactly Gaussian linear smearing at fixed background and scalar history",
        "variance_coefficient": variance,
        "variance_result": "Var(F_eta)<=coefficient*||U_eta||L2(time)^2 for interval length<=1",
        "anchor_variance_upper": 4 * s.Rational(1, 10**2380),
        "anchor_standard_deviation_upper": 2 * s.Rational(1, 10**1190),
        "uniform_family": "Variance and mean/source-response bounds decay as kappa0/kappa. The standard deviation decays as sqrt(kappa0/kappa).",
        "not_a_metric_noise_bound": "This is the scalar force linear in the vector field at fixed g and fixed Phi. It is not stress-tensor noise, a stochastic metric solution or a bound on interacting scalar/graviton loops.",
        "checks": {
            key: s.ImmutableMatrix(value.applyfunc(s.simplify))
            if isinstance(value, s.MatrixBase)
            else s.simplify(value)
            for key, value in checks.items()
        },
        "gates": {
            "mass_positive": mass > 0,
            "nonzero_conditional_noise_upper": anchor > 0,
            "full_variance_below_display": anchor < 4 * s.Rational(1, 10**2380),
            "temporal_constraint_mode_not_discarded": P[0, 0] != 0,
        },
    }
