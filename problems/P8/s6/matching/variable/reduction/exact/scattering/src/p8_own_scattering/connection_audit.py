"""Independent Acb/Arb corroboration of the coefficient-16 connection.

No formulas or evaluated expressions are imported from ``connection.py``.
The generic Gamma connection coefficients and ordinary (not regularized)
Gauss functions are evaluated independently at 128 and 256 bits.  All inputs,
enclosure thresholds, and reported outer intervals are exact rationals.

The primary analytic arguments are DLMF 15.10.21, 15.8.2, and 5.5.3:
https://dlmf.nist.gov/15.10.E21
https://dlmf.nist.gov/15.8.E2
https://dlmf.nist.gov/5.5.E3
Finite ball overlaps do NOT prove these analytic identities.  They corroborate
their branches, normalizations, and phase dictionary.  In contrast, the ball
inequality |B| > 1/40 is itself a rigorous numerical inequality at this fixed
rho.  No finite-window ODE estimate or physical scattering claim is added.
"""

from fractions import Fraction

from flint import acb, arb, ctx, fmpq

PRECISIONS = (128, 256)
DELTA_FIXTURES = (Fraction(1), Fraction(1, 10**6))
X_FIXTURES = (Fraction(-1), Fraction(-1, 4), Fraction(0), Fraction(1, 4), Fraction(1))
ENCLOSURE_DENOMINATOR = 10**9


def _precision(value):
    if type(value) is not int:
        raise TypeError("precision must be an integer, not bool or a rounded value")
    if value not in PRECISIONS:
        raise ValueError("this audit specifies only 128 and 256 bits")
    return value


def _q(numerator, denominator=1):
    return acb(fmpq(numerator, denominator))


def _matrix_product(left, right):
    return tuple(tuple(sum((left[i][k] * right[k][j] for k in range(2)), acb(0))
                       for j in range(2)) for i in range(2))


def _matrix_determinant(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def _matrix_inverse(matrix):
    determinant = _matrix_determinant(matrix)
    if determinant.contains(0):
        raise ArithmeticError("independent parity basis is not enclosed as invertible")
    return ((matrix[1][1] / determinant, -matrix[0][1] / determinant),
            (-matrix[1][0] / determinant, matrix[0][0] / determinant))


def _coefficients():
    rho = arb(7).sqrt() / 2
    irho = acb(0, rho)
    a, b, c = _q(-1, 2), _q(3, 2), 1 - irho
    # The two generic coefficients in DLMF 15.10.21, before reflection.
    aa = c.gamma() * (c - a - b).gamma() / ((c - a).gamma() * (c - b).gamma())
    bb = c.gamma() * (a + b - c).gamma() / (a.gamma() * b.gamma())
    closed_b = acb(0, 1 / (arb.pi() * rho).sinh())

    pa, pb = _q(-1, 4) + irho / 2, _q(-1, 4) - irho / 2
    # Coefficient of x**(1/2+i*rho): the (-z)**(-pb) term in 15.8.2.
    # The Gamma(c) factors convert the source's regularized F to ordinary F.
    ee = (_q(1, 2).gamma() * (pa - pb).gamma() * acb(8)**(-pb)
          / (pa.gamma() * (_q(1, 2) - pb).gamma()))
    oa, ob, oc = pa + _q(1, 2), pb + _q(1, 2), _q(3, 2)
    oo = (oc.gamma() * (oa - ob).gamma() * acb(8)**(-ob)
          / (oa.gamma() * (oc - ob).gamma()))
    return {
        "rho": rho, "irho": irho, "a": a, "b": b, "c": c,
        "pa": pa, "pb": pb, "oa": oa, "ob": ob,
        "A": aa, "B_gamma": bb, "B_closed": closed_b,
        "E_plus": ee, "O_plus": oo,
    }


def _finite_point(x_fraction, data, alpha, beta):
    """Cauchy values evaluated on real hypergeometric arguments off the cut."""
    x = _q(x_fraction.numerator, x_fraction.denominator)
    pa, pb, oa, ob = (data[key] for key in ("pa", "pb", "oa", "ob"))
    z = -8 * x**2
    even = z.hypgeom_2f1(pa, pb, _q(1, 2))
    even_x = -16 * x * (2 * pa * pb) * z.hypgeom_2f1(pa + 1, pb + 1, _q(3, 2))
    odd_factor = z.hypgeom_2f1(oa, ob, _q(3, 2))
    odd = x * odd_factor
    odd_x = odd_factor - 16 * x**2 * (oa * ob / _q(3, 2)) * z.hypgeom_2f1(
        oa + 1, ob + 1, _q(5, 2)
    )

    cosh_t = (1 + 8 * x**2).sqrt()
    tanh_t = arb(8).sqrt() * x / cosh_t
    t = (arb(8).sqrt() * x).asinh()
    w = (1 - tanh_t) / 2
    y = 1 - w
    a, b, c, irho = (data[key] for key in ("a", "b", "c", "irho"))
    plus, minus = (irho * t).exp(), (-irho * t).exp()
    ff = w.hypgeom_2f1(a, b, c)
    psi = plus * ff
    psi_t = plus * (irho * ff - 2 * w * (1 - w) * a * b / c
                    * w.hypgeom_2f1(a + 1, b + 1, c + 1))
    jost_y = cosh_t.sqrt() * psi
    jost_y_x = arb(8).sqrt() / cosh_t.sqrt() * (psi_t + tanh_t * psi / 2)

    # Euler-transform the second left factor BEFORE extracting its phase.
    # Both small-y Gauss factors then have a=-1/2,b=3/2 and c=1+/-i*rho.
    left_plus = y.hypgeom_2f1(a, b, 1 + irho)
    left_minus = y.hypgeom_2f1(a, b, 1 - irho)
    y_t = 2 * y * (1 - y)
    left_plus_t = (irho * left_plus + y_t * a * b / (1 + irho)
                   * y.hypgeom_2f1(a + 1, b + 1, 2 + irho))
    left_minus_t = (-irho * left_minus + y_t * a * b / (1 - irho)
                    * y.hypgeom_2f1(a + 1, b + 1, 2 - irho))
    connected = data["A"] * plus * left_plus + data["B_gamma"] * minus * left_minus
    connected_t = (data["A"] * plus * left_plus_t
                   + data["B_gamma"] * minus * left_minus_t)
    wronskian_x = even * odd_x - even_x * odd
    wronskian_t = psi * psi_t.conjugate() - psi_t * psi.conjugate()
    return {
        "even": even, "even_x": even_x, "odd": odd, "odd_x": odd_x,
        "jost_y": jost_y, "jost_y_x": jost_y_x,
        "Wronskian_x": wronskian_x, "Wronskian_t": wronskian_t,
        "residuals": {
            "parity_Wronskian": wronskian_x - 1,
            "Jost_parity_value": jost_y - alpha * even - beta * odd,
            "Jost_parity_x_derivative": jost_y_x - alpha * even_x - beta * odd_x,
            "Jost_time_Wronskian": wronskian_t / (-2 * irho) - 1,
            "left_right_Jost_value": psi - connected,
            "left_right_Jost_t_derivative": psi_t - connected_t,
        },
    }


def _compute():
    """Compute at the caller's active precision; no global cached balls."""
    data = _coefficients()
    rho, irho, aa, bb, ee, oo = (data[key] for key in
                                ("rho", "irho", "A", "B_gamma", "E_plus", "O_plus"))
    a_squared = aa.real**2 + aa.imag**2
    b_squared = bb.real**2 + bb.imag**2
    absolute_b = abs(bb)
    parity_wronskian = 4 * rho * (ee * oo.conjugate()).imag
    time = ((aa.conjugate(), -bb.conjugate()), (-bb, aa))
    residuals = {
        "Gamma_B_reflection": bb - data["B_closed"],
        "Gamma_A_modulus": acb(a_squared - (arb.pi() * rho).coth()**2),
        "unit_flux": acb(a_squared - b_squared - 1),
        "scattering_flux_sum": acb(1 / a_squared + b_squared / a_squared - 1),
        "asymptotic_parity_Wronskian": acb(parity_wronskian - 1),
        "time_transfer_determinant": _matrix_determinant(time) - 1,
    }

    parity_basis = ((ee, oo), (ee.conjugate(), oo.conjugate()))
    parity_inverse = _matrix_inverse(parity_basis)
    parity_reflection = ((acb(1), acb(0)), (acb(0), acb(-1)))
    radial = {}
    for delta_fraction in DELTA_FIXTURES:
        delta = _q(delta_fraction.numerator, delta_fraction.denominator)
        p_plus = _q(1, 2) + irho
        p_minus = _q(1, 2) - irho
        scale_plus, scale_minus = delta**(-p_plus / 2), delta**(-p_minus / 2)
        scaled = ((scale_plus * ee, scale_plus * oo),
                  (scale_minus * ee.conjugate(), scale_minus * oo.conjugate()))
        parity = _matrix_product(_matrix_product(scaled, parity_reflection),
                                 _matrix_inverse(scaled))
        phase = (2 * irho * (4 * arb(2).sqrt() / delta.sqrt()).log()).exp()
        jost = ((-bb.conjugate(), aa.conjugate() * phase), (aa / phase, -bb))
        label = str(delta_fraction)
        radial[label] = {"parity": parity, "Jost": jost, "phase": phase}
        for row in range(2):
            for col in range(2):
                residuals[f"radial_{label}_{row}{col}"] = parity[row][col] - jost[row][col]
        residuals[f"radial_{label}_determinant"] = _matrix_determinant(parity) + 1

    # y_R=sqrt(cosh(t))*psi_R has right x-power coefficient
    # 8**(1/4)*(4*sqrt(2))**(i*rho), not 1, before solving in the parity basis.
    right_coefficient = arb(8).sqrt().sqrt() * (irho * (4 * arb(2).sqrt()).log()).exp()
    alpha = parity_inverse[0][0] * right_coefficient
    beta = parity_inverse[1][0] * right_coefficient
    points = {}
    for x_fraction in X_FIXTURES:
        label = str(x_fraction)
        fixture = _finite_point(x_fraction, data, alpha, beta)
        points[label] = fixture
        residuals.update({f"x_{label}_{key}": value
                          for key, value in fixture["residuals"].items()})

    # Deliberate mistakes are separated by rigorous balls, not float deltas.
    negative_controls = {
        "wrong_B_sign_disjoint": not bb.overlaps(-data["B_closed"]),
        "missing_radial_phase_disjoint": not radial["1"]["parity"][0][1].overlaps(aa.conjugate()),
        "wrong_Wronskian_orientation_disjoint": not acb(-parity_wronskian).contains(1),
        "reflected_flux_not_B_squared": b_squared - b_squared / a_squared > 0,
        "wrong_unit_Jost_x_amplitude_disjoint": not points["0"]["jost_y"].overlaps(parity_inverse[0][0]),
    }
    data.update({
        "A_abs_squared": a_squared, "B_abs_squared": b_squared,
        "B_abs": absolute_b, "asymptotic_parity_Wronskian": parity_wronskian,
        "time_transfer": time, "radial": radial, "points": points,
        "residuals": residuals, "negative_controls": negative_controls,
    })
    return data


def balls(precision=128):
    """Return fresh rigorous balls; restores the caller's global precision.

    For new arithmetic on the returned balls use a matching ``ctx.workprec``
    context.  ``report`` evaluates all comparisons inside the proper context.
    This API accepts only the two explicitly audited integer precisions.
    """
    precision = _precision(precision)
    with ctx.workprec(precision):
        return _compute()


def _outer_real_interval(value):
    if not value.is_finite():
        raise ArithmeticError("nonfinite ball cannot be serialized as an enclosure")
    lower = Fraction(str(value.lower().fmpq()))
    upper = Fraction(str(value.upper().fmpq()))
    scaled_lower, scaled_upper = lower * ENCLOSURE_DENOMINATOR, upper * ENCLOSURE_DENOMINATOR
    low_integer = scaled_lower.numerator // scaled_lower.denominator
    high_integer = -((-scaled_upper.numerator) // scaled_upper.denominator)
    return [str(Fraction(low_integer, ENCLOSURE_DENOMINATOR)),
            str(Fraction(high_integer, ENCLOSURE_DENOMINATOR))]


def _summarize(data, precision):
    cap_denominator = 1 << (precision // 2)
    cap = arb(fmpq(1, cap_denominator))
    residual_checks = {
        key: {"contains_zero": value.is_finite() and value.contains(0),
              "absolute_value_below_cap": value.is_finite() and value.abs_upper() < cap}
        for key, value in data["residuals"].items()
    }
    checks = {
        "all_residuals_contain_zero": all(item["contains_zero"] for item in residual_checks.values()),
        "all_residuals_below_rational_cap": all(item["absolute_value_below_cap"]
                                                for item in residual_checks.values()),
        "rho_strict_rational_interval": arb(fmpq(5, 4)) < data["rho"] < arb(fmpq(4, 3)),
        "B_strict_lower_1_over_40": data["B_abs"] > arb(fmpq(1, 40)),
        "B_strict_upper_1_over_30": data["B_abs"] < arb(fmpq(1, 30)),
        **data["negative_controls"],
    }
    enclosures = {key: {"real": _outer_real_interval(data[key].real),
                        "imag": _outer_real_interval(data[key].imag)}
                  for key in ("A", "B_gamma", "B_closed", "E_plus", "O_plus")}
    enclosures["B_abs"] = _outer_real_interval(data["B_abs"])
    return {
        "precision_bits": precision,
        "residual_absolute_cap": str(Fraction(1, cap_denominator)),
        "residual_count": len(residual_checks),
        "residual_checks": residual_checks,
        "outer_rational_enclosures": enclosures,
        "checks": checks,
    }


def report():
    """Deterministic, JSON-compatible metadata from two fresh ball replays.

    No rounded decimal display is read back as an exact value.  Outer rational
    intervals are formed from exact binary endpoints by directed integer
    rounding.  Overlap, width, and strict inequality assertions use live balls.
    """
    runs, computed = [], {}
    for precision in PRECISIONS:
        with ctx.workprec(precision):
            computed[precision] = _compute()
            runs.append(_summarize(computed[precision], precision))
    keys = ("A", "B_gamma", "B_closed", "E_plus", "O_plus")
    cross_precision = {key: computed[128][key].overlaps(computed[256][key]) for key in keys}
    checks = {f"precision_{run['precision_bits']}_passes": all(run["checks"].values()) for run in runs}
    checks["cross_precision_coefficients_overlap"] = all(cross_precision.values())
    return {
        "engine": "python-flint Acb/Arb rigorous balls",
        "role": "numerical corroboration, not proof of the analytic connection",
        "precisions_bits": list(PRECISIONS),
        "delta_fixtures": [str(value) for value in DELTA_FIXTURES],
        "x_fixtures": [str(value) for value in X_FIXTURES],
        "outer_enclosure_grid_denominator": ENCLOSURE_DENOMINATOR,
        "runs": runs,
        "cross_precision_coefficient_overlap": cross_precision,
        "checks": checks,
        "boundaries": [
            "Gamma and hypergeometric identity proofs remain analytic and source-backed",
            "finite point checks are not a uniform ODE or asymptotic remainder estimate",
            "no physical S-matrix, matter forcing, vacuum, quantum production, or UV claim",
        ],
    }


def checks():
    """Top-level pass flags; see report for individual residual enclosures."""
    return report()["checks"]
