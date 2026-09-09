"""Exact nonminimal stress, conformal bridge and flat spectral moments."""

from functools import cache

import sympy as sp


@cache
def data():
    xi = sp.Symbol("xi", real=True)
    z = sp.Symbol("field", real=True)
    grad = sp.Matrix(sp.symbols("field_d0:4", real=True))
    h = sp.symbols("field_h0:10", real=True)
    hess = sp.Matrix(
        [
            [h[0], h[1], h[2], h[3]],
            [h[1], h[4], h[5], h[6]],
            [h[2], h[5], h[7], h[8]],
            [h[3], h[6], h[8], h[9]],
        ]
    )
    eta = sp.diag(1, -1, -1, -1)
    norm = (grad.T * eta * grad)[0]
    box = sum(eta[i, i] * hess[i, i] for i in range(4))
    square_hess = 2 * grad * grad.T + 2 * z * hess
    square_box = sum(eta[i, i] * square_hess[i, i] for i in range(4))
    flat = grad * grad.T - eta * norm / 2 + xi * (eta * square_box - square_hess)
    tr = sum(eta[i, i] * flat[i, i] for i in range(4))
    a = sp.Symbol("a", positive=True)
    ap, app = sp.symbols("a_prime a_second", real=True)
    A = ap / a
    phi = z / a
    dphi = grad / a
    dphi[0] = (grad[0] - A * z) / a
    ddphi = hess / a
    ddphi[0, 0] = (hess[0, 0] - 2 * A * grad[0] + (2 * A * A - app / a) * z) / a
    for i in range(1, 4):
        ddphi[0, i] = ddphi[i, 0] = (hess[0, i] - A * grad[i]) / a

    def Gamma(k, i, j):
        if k == 0 and i == j:
            return A
        if k > 0 and ((i == 0 and j == k) or (j == 0 and i == k)):
            return A
        return 0

    cov_hess = sp.Matrix(
        4, 4, lambda i, j: ddphi[i, j] - sum(Gamma(k, i, j) * dphi[k] for k in range(4))
    )
    physical_box = sum(eta[i, i] * cov_hess[i, i] for i in range(4)) / a**2
    R = 6 * app / a**3
    G = sp.diag(-3 * A * A, *[2 * app / a - A * A] * 3)
    sq_cov = 2 * dphi * dphi.T + 2 * phi * cov_hess
    sq_box = sum(eta[i, i] * sq_cov[i, i] for i in range(4)) / a**2
    curved = dphi * dphi.T - eta * (dphi.T * eta * dphi)[0] / 2
    curved += xi * (a * a * eta * sq_box - sq_cov - G * phi * phi)
    curved = curved.applyfunc(sp.factor)

    time, x, y, w, tau, amp = sp.symbols("time x y z tau amplitude", real=True)
    coords = (time, x, y, w)
    example = amp * (3 + (time / tau) ** 2 + (x * x + y * y + w * w) / (3 * tau * tau))
    sub = {z: example}
    sub.update({grad[i]: sp.diff(example, coords[i]) for i in range(4)})
    sub.update(
        {
            hess[i, j]: sp.diff(example, coords[i], coords[j])
            for i in range(4)
            for j in range(i, 4)
        }
    )
    exrho = sp.factor(
        flat[0, 0].subs(sub, simultaneous=True).subs(xi, sp.Rational(1, 6))
    )
    exeed = sp.factor(
        (flat[0, 0] - tr / 2).subs(sub, simultaneous=True).subs(xi, sp.Rational(1, 6))
    )
    exline = sp.factor(exrho.subs({x: 0, y: 0, w: 0}))
    f, fp = sp.symbols("sampler sampler_prime", real=True)
    # A total derivative relates the weighted stress to positive squares.
    original_rho = (
        grad[0] ** 2 / 2
        + (sp.Rational(1, 2) - 2 * xi) * sum(grad[i] ** 2 for i in range(1, 4))
        - 2 * xi * z * hess[0, 0]
    )
    positive = f * f * grad[0] ** 2 / 2 + (sp.Rational(1, 2) - 2 * xi) * f * f * sum(
        grad[i] ** 2 for i in range(1, 4)
    )
    positive += 2 * xi * (f * grad[0] + fp * z) ** 2 - 2 * xi * fp * fp * z * z
    boundary_derivative = (
        -2 * xi * (2 * f * fp * z * grad[0] + f * f * (grad[0] ** 2 + z * hess[0, 0]))
    )
    theta, k = sp.symbols("theta k", positive=True)
    spectral = sp.integrate(
        (1 - 2 * xi) * k**3 + 2 * xi * k * (theta - k) ** 2, (k, 0, theta)
    )
    constant = sp.factor(spectral / theta**4 / (4 * sp.pi**2))
    checks = {
        "flat_offshell_trace": sp.factor(tr - (6 * xi - 1) * norm - 6 * xi * z * box),
        "flat_on_shell_energy_before_sampler": sp.factor(
            flat[0, 0] - original_rho - 2 * xi * z * box
        ),
        "positive_sampler_squares_exact_boundary": sp.factor(
            f * f * original_rho - positive - boundary_derivative
        ),
        "all_component_conformal_stress_weight": (
            curved.subs(xi, sp.Rational(1, 6)) - flat.subs(xi, sp.Rational(1, 6)) / a**2
        ).applyfunc(sp.factor),
        "conformal_field_equation": sp.factor(physical_box + R * phi / 6 - box / a**3),
        "conformal_offshell_trace_is_equation_term": sp.factor(
            tr.subs(xi, sp.Rational(1, 6)) - z * box
        ),
        "coherent_local_example_solves_massless_equation": sp.factor(
            box.subs(sub, simultaneous=True)
        ),
        "coherent_conformal_energy_is_EED": sp.factor(exrho - exeed),
        "coherent_line_exact_negative_formula": sp.factor(
            exline - amp**2 / tau**2 * (sp.Rational(4, 3) * time**2 / tau**2 - 2)
        ),
        "spectral_coefficient_exact": sp.factor(
            constant - (3 - 4 * xi) / (48 * sp.pi**2)
        ),
        "conformal_spectral_coefficient": sp.factor(
            constant.subs(xi, sp.Rational(1, 6)) - 7 / (144 * sp.pi**2)
        ),
    }
    return {
        "xi": xi,
        "flat_stress": flat,
        "flat_trace": tr,
        "conformal_stress": curved,
        "coherent_local_field": example,
        "coherent_line_energy": exline,
        "flat_quantum_coefficient_without_hbar": constant,
        "checks": checks,
    }


if __name__ == "__main__":
    d = data()
    for name, v in d["checks"].items():
        assert all(
            x == 0 for x in (list(v) if isinstance(v, sp.MatrixBase) else [v])
        ), (name, v)
    print("EXACT_NONMINIMAL_STRESS_CONFORMAL_AND_COHERENT_IDENTITIES_PASS", flush=True)
    print(d["coherent_line_energy"], flush=True)
    print(d["flat_quantum_coefficient_without_hbar"], flush=True)
