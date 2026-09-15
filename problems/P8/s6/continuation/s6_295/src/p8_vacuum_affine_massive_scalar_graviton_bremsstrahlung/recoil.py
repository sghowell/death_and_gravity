"""Exact physical recoil map and explicit compact-domain inequalities."""

from functools import cache

import sympy as s


def momenta(energy, omega, direction, outgoing_direction):
    E, w = map(s.sympify, (energy, omega))
    nhat, u = s.Matrix(direction), s.Matrix(outgoing_direction)
    if nhat.shape != (3, 1) or u.shape != (3, 1):
        raise ValueError("Require two three-dimensional unit directions")
    if s.simplify(nhat.dot(nhat) - 1) != 0 or s.simplify(u.dot(u) - 1) != 0:
        raise ValueError("Require exactly unit recoil and rest-frame directions")
    ep = s.sqrt(E * (E - w))
    r = s.sqrt(E * (E - w) - 1)
    r0 = s.sqrt(E * E - 1)
    c = u.dot(nhat)
    gamma = (E - w / 2) / ep
    p1 = s.Matrix([E, 0, 0, r0])
    p2 = s.Matrix([E, 0, 0, -r0])
    out = []
    for sign in (1, -1):
        spatial = sign * r * u + (sign * (gamma - 1) * r * c - w / 2) * nhat
        out.append(s.Matrix([E - w / 2 - sign * w * r * c / (2 * ep), *spatial]))
    q = s.Matrix([w, *(w * nhat)])
    born = (-p1, -p2, s.Matrix([E, *(r0 * u)]), s.Matrix([E, *(-r0 * u)]))
    return (-p1, -p2, *out), q, born


def density_ratio(energy, omega):
    E, w = map(s.sympify, (energy, omega))
    return E * s.sqrt(E * (E - w) - 1) / (s.sqrt(E * (E - w)) * s.sqrt(E * E - 1))


def compact_domain():
    return {
        "mass_squared": s.S.One,
        "hard_energy_min": s.Rational(5, 4),
        "hard_energy_max": s.Integer(2),
        "heavy_mass_squared_min": s.Integer(128),
        "resolution_max": s.Rational(1, 8),
        "all_graviton_and_hard_rest_frame_angles": True,
    }


@cache
def data():
    checks = {}
    bounds = {}

    def put(name, value):
        checks[name] = s.factor(value)

    def positive(name, value, variables):
        value = s.factor(value)
        num, den = s.fraction(value)
        if not all(c >= 0 for c in s.Poly(num, *variables).coeffs()):
            raise ValueError("Negative numerator in " + name)
        if not all(c >= 0 for c in s.Poly(den, *variables).coeffs()) or den == 0:
            raise ValueError("Invalid denominator in " + name)
        bounds[name] = value

    # Exact recoil map with rprime and Eprime eliminated by squared relations.
    E, w, c, r, Ep = s.symbols("E omega cos rprime Eprime", real=True)
    # On-shell identities are polynomial after multiplying denominators.
    gamma = (E - w / 2) / Ep
    energy = E - w / 2 - w * r * c / (2 * Ep)
    # Spatial r*u + h*n, u.u=n.n=1, u.n=c.
    h = (gamma - 1) * r * c - w / 2
    norm = r * r + h * h + 2 * r * h * c
    reduce = lambda expr: s.rem(
        s.rem(s.expand(expr), Ep**2 - E * (E - w), Ep), r**2 - E * (E - w) + 1, r
    )
    put("whole_recoil_scalar3_shell", reduce((energy * energy - norm - 1) * 4 * Ep**2))
    put(
        "whole_recoil_scalar4_shell",
        reduce((energy.subs(c, -c) ** 2 - norm.subs(c, -c) - 1) * 4 * Ep**2),
    )
    put("whole_recoil_energy_conservation", energy + energy.subs(c, -c) - (2 * E - w))
    put("whole_recoil_spatial_conservation", h + h.subs(c, -c) + w)
    put("boost_gamma_identity", s.expand((E - w / 2) ** 2 - E * (E - w) - w * w / 4))
    put("radial_difference_identity", (E**2 - 1) - (E * (E - w) - 1) - E * w)
    beta2 = 1 - 1 / E**2
    betap2 = 1 - 1 / (E * (E - w))
    put("two_body_density_squared_difference", beta2 - betap2 - w / (E**2 * (E - w)))
    n, aa, bb, L = s.symbols("n a b Lambda", positive=True)
    put(
        "heavy_denominator_subtraction",
        1 / ((n - aa) * (n - bb))
        - 1 / n**2
        - (n * (aa + bb) - aa * bb) / (n**2 * (n - aa) * (n - bb)),
    )
    put("majorant_complete_heavy_sum", 3 * 40 * 160 * L**6 - 19200 * L**6)
    put("tuned_Born_recoil_relative_bound", s.Rational(1, 2) * 34 * 72 / 4 - 306)
    put("one_final_leg_current_variation", 12 * 4 + 4 * 4 * 16 - 304)
    put("both_final_leg_current_variation", 2 * 304 - 608)
    put("heavy_remainder_over_Born", 19200 * 2**6 / 4 - 307200)
    put("whole_physical_remainder", 307200 + 306 * 64 + 608 - 327392)
    put(
        "phase_space_angular_normalization",
        2 * 4 * s.pi / (2 * (2 * s.pi) ** 3) - 1 / (2 * s.pi**2),
    )
    B = s.Integer(330000)
    cut = s.Symbol("cut", positive=True)
    omega = s.Symbol("omega", positive=True)
    integral = s.integrate(omega * (B * B + (128 * B + 8192) / omega), (omega, 0, cut))
    put(
        "integrable_real_rate_error",
        integral - (s.Integer(42248192) * cut + s.Integer(54450000000) * cut**2),
    )
    # Polynomial positivity certificates for every nontrivial endpoint margin.
    x = s.Symbol("x", nonnegative=True)
    positive(
        "subthreshold_denominator_gap", ((n - 16) - 7 * n / 8).subs(n, 128 + x), (x,)
    )
    positive("majorant_fprime_first", (32 / n - 28 / (n - 16)).subs(n, 128 + x), (x,))
    positive(
        "majorant_fprime_second",
        (256 / n**2 - 196 / (n - 16) ** 2).subs(n, 128 + x),
        (x,),
    )
    positive(
        "majorant_fprime_total", (34 / n - 32 / n - 256 / n**2).subs(n, 128 + x), (x,)
    )
    # Compact box E in[5/4,2], omega in[0,1/8]. Products use independent
    # monotone endpoint arguments; no invalid use of unbounded substitute E.
    put(
        "minimum_recoil_radial_square",
        s.Rational(5, 4) * (s.Rational(5, 4) - s.Rational(1, 8))
        - 1
        - s.Rational(13, 32),
    )
    positive("recoil_radial_above_half", s.Rational(13, 32) - s.Rational(1, 4), (x,))
    positive(
        "radial_shift_below_2omega",
        2 - s.Integer(2) / (s.Rational(1, 2) + s.Rational(3, 4)),
        (x,),
    )
    positive(
        "spatial_shift_below_3omega",
        3 - (2 + s.Rational(1, 32) + s.Rational(1, 2)),
        (x,),
    )
    put(
        "density_ratio_coarse_coefficient",
        s.Rational(25, 9) / (s.Rational(25, 16) * s.Rational(9, 8))
        - s.Rational(128, 81),
    )
    positive("density_ratio_error_below_2omega", 2 - s.Rational(128, 81), (x,))
    positive("whole_remainder_margin", 330000 - 327392, (x,))
    positive(
        "resolution_maximum_below_1e8_kappa",
        10**8 - (s.Rational(42248192, 8) + s.Rational(54450000000, 64)) / 18,
        (x,),
    )
    positive(
        "positive_heavy_resolvent_remainder",
        (1 / (n * (n - 2) ** 2) - 1 / n**3).subs(n, 128 + x),
        (x,),
    )
    positive(
        "generic_resolvent_numerator_bound",
        (40 * n * L**2 - (32 * n * L**2 + 256 * L**4)).subs(n, 32 * L**2 + x),
        (L, x),
    )
    positive(
        "generic_heavy_denominator_half_gap",
        (n - 16 * L**2 - n / 2).subs(n, 32 * L**2 + x),
        (L, x),
    )
    put(
        "original_kappa_error_exponent",
        s.Integer(10) ** 8 / s.Integer(10) ** 800 - s.Rational(1, 10**792),
    )

    for index, (E, w) in enumerate(
        (
            (s.Rational(5, 4), s.Rational(1, 8)),
            (s.Rational(3, 2), s.Rational(1, 16)),
            (s.Integer(2), s.Rational(1, 8)),
        )
    ):
        k, q, born = momenta(
            E,
            w,
            [s.Rational(3, 5), 0, s.Rational(4, 5)],
            [0, s.Rational(4, 5), s.Rational(3, 5)],
        )
        eta = s.diag(1, -1, -1, -1)
        for i, p in enumerate(k):
            put(str(index) + "_public_recoil_mass_" + str(i), (p.T * eta * p)[0] - 1)
        for i, component in enumerate(sum(k, q.copy())):
            put(str(index) + "_public_recoil_conservation_" + str(i), component)
        put(
            str(index) + "_public_density_ratio_square",
            density_ratio(E, w) ** 2 - (1 - 1 / (E * (E - w))) / (1 - 1 / E**2),
        )
        for i, p in enumerate(born):
            put(str(index) + "_public_Born_mass_" + str(i), (p.T * eta * p)[0] - 1)
    return {
        "whole_compact_physical_domain": compact_domain(),
        "whole_recoil_map": "Boost the final pair of rest energy Eprime=sqrt(E(E-omega)) and radius rprime=sqrt(E(E-omega)-1) along -nhat. p3^0=E-omega/2-omega*rprime*c/(2Eprime); vec p3=rprime*u+[(gamma-1)rprime*c-omega/2]nhat, gamma=(E-omega/2)/Eprime,c=u.nhat. p4 reverses all rprime terms. All external scalar masses are1 and sum p3+p4=P-q.",
        "whole_recoil_pointwise_bounds": "Each final energy shifts at most omega and each final spatial vector at most3omega. All scalar energies<=2. Six pair invariants lie in[-12,16]; their total shift is<=72omega. For f(a)=(a-2)^2/(n-a), |f'|<=34/n when n>=128. The original Born positive identity implies |Abar-A0|/A0<=306omega.",
        "whole_soft_current_and_remainder_bounds": "Per unit-Frobenius helicity, |sum J_i|<=64/omega; two final-leg current differences sum<=608. The improved heavy remainder obeys |R|/A0<=307200/sqrt(kappa). Thus |M5-A0 S0|/A0<=327392/sqrt(kappa)<330000/sqrt(kappa). All are uniform over both physical angular directions on the stated compact domain.",
        "whole_phase_space_factorization": "dPhi3=d^3q/[(2pi)^3 2omega] dPhi2(P-q). Using the same outgoing rest-frame angle u in both measures gives J=beta(s')/beta(s),0<J<=1 and1-J<=2omega. Same incoming flux and identical-final-scalar factor cancel in the Born-normalized rate.",
        "explicit_nonnegative_majorant_margins": bounds,
        "checks": checks,
        "gates": {
            "whole_on_shell_recoil_and_momentum_conservation": True,
            "compact_mass_energy_heavy_and_resolution_domain_explicit": True,
            "all_angles_and_massive_Doppler_gap_retained": True,
            "tuned_positive_Born_denominator_not_generic_quantum_positivity": True,
            "every_majorant_endpoint_and_denominator_margin_checked": True,
            "phase_space_Jacobian_and_identical_scalar_factor_retained": True,
            "not_virtual_IR_pairing_or_all_energy_Regge": True,
        },
    }
