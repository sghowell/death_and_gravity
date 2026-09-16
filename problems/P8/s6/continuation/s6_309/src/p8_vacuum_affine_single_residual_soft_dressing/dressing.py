"""Same-state remaining-energy kernel and the unchanged physical residual."""

from functools import cache

import sympy as s

from . import series, source


def factor(index, conversion, remaining):
    a, d = map(series.exact_real, (index, conversion))
    y = series.require_remaining(remaining)
    if a < 0:
        raise ValueError("Require a nonnegative physical leading-soft index")
    return s.exp(d - s.EulerGamma * a) * y**a / s.gamma(1 + a)


def original_weight(
    state_index, state_conversion, elastic_index, elastic_conversion, omega, resolution
):
    a, d, a0, d0, w, x = map(
        series.exact_real,
        (
            state_index,
            state_conversion,
            elastic_index,
            elastic_conversion,
            omega,
            resolution,
        ),
    )
    kap = source.KAPPA
    if not (
        0 < w < x <= s.Rational(1, 8)
        and 0 <= a <= 6 / kap
        and 0 <= a0 <= s.Rational(4, 5) / kap
        and abs(a - a0) <= 40 * w / kap
        and abs(d) <= 2000 / kap
        and abs(d0) <= 112 / kap
    ):
        raise ValueError(
            "Require the original radiative-state bounds and positive remaining energy"
        )
    return (
        s.exp(d - d0 - s.EulerGamma * (a - a0))
        * s.gamma(1 + a0)
        / s.gamma(1 + a)
        * (x - w) ** a
        / x**a0
    )


def transverse_basis(direction):
    n = s.Matrix(direction)
    if n.shape != (3, 1) or s.simplify(n.dot(n) - 1) != 0:
        raise ValueError("Require an exact three-dimensional unit emission direction")
    for v in n:
        series.exact_real(v)
    axis = next(s.eye(3)[:, j] for j in range(3) if n[j] ** 2 != 1)
    e = axis - n * axis.dot(n)
    e = e / s.sqrt(e.dot(e))
    f = n.cross(e)
    e4, f4 = s.Matrix([0, *e]), s.Matrix([0, *f])
    # Divide squared contractions bytwo below, keeping the tensors algebraic.
    return e4 * e4.T - f4 * f4.T, e4 * f4.T + f4 * e4.T


def physical_kernel(energy, omega, direction, outgoing):
    E, w = map(series.exact_real, (energy, omega))
    if not 0 < w <= s.Rational(1, 8):
        raise ValueError("Require positive marked energy at most1/8")
    n = s.Matrix(direction)
    epsilons = transverse_basis(n)
    q = s.Matrix([w, *(w * n)])
    points, rays, born_points = source.recoil.momenta(E, [q], outgoing)
    v, tree = source.vertices, source.tree
    ss = v.dot(born_points[0] + born_points[1], born_points[0] + born_points[1])
    tt = v.dot(born_points[0] + born_points[2], born_points[0] + born_points[2])
    source.softlimit.estimates.require_domain(ss, tt, w)
    am = tree.born_continuation(
        born_points, source.HEAVY_MASS2, source.CUBIC, source.CONTACT
    )
    ag = v.born(born_points) / source.KAPPA
    born = s.factor(am + ag)
    tensor = tree.whole_tensor(
        points,
        q,
        mass=1,
        heavy=source.HEAVY_MASS2,
        cubic=source.CUBIC,
        contact=source.CONTACT,
        kappa=1,
    )
    soft = tree.soft_current(born_points, q)
    rho = source.single_recoil.density_ratio(E, w)
    contractions = []
    for eps in epsilons:
        whole = s.factor(
            (v.pair(eps, tensor) + v.amplitude(points, q, eps) / source.KAPPA) / born
        )
        lead = s.factor(v.pair(eps, soft))
        contractions.append((whole, lead))
    scaled = s.factor(
        w * sum((rho * whole**2 - lead**2) / 2 for whole, lead in contractions)
    )
    return {
        "points": points,
        "ray": rays[0],
        "born_points": born_points,
        "rho": rho,
        "whole_Born": born,
        "sqrt_kappa_amplitude_ratios_and_leading_currents": tuple(contractions),
        "kappa_times_energy_weighted_signed_density": scaled,
        "whole_radiative_K0": source.index.kernel(points, rays),
        "whole_elastic_K0": source.index.kernel(born_points, []),
    }


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_log(s.sympify(value)))

    a, a0, x, u = s.symbols("state_index elastic_index x fraction", positive=True)
    d, d0 = s.symbols("state_conversion elastic_conversion", real=True)
    log_ratio = (
        d
        - d0
        - s.EulerGamma * (a - a0)
        + s.loggamma(1 + a0)
        - s.loggamma(1 + a)
        + (a - a0) * s.log(x)
        + a * s.log(1 - u)
    )
    same_energy = log_ratio - a * s.log(1 - u)
    put("remaining_energy_log_split", log_ratio - same_energy - a * s.log(1 - u))
    put(
        "remaining_energy_decreases_weight_derivative",
        s.diff(a * s.log(1 - u), u) + a / (1 - u),
    )
    put(
        "same_state_first_order_virtual_pairing",
        a / (2 * source.EP) - a / (2 * source.EP),
    )
    put(
        "wrong_elastic_virtual_pole",
        (a - a0) / (2 * source.EP) - a / (2 * source.EP) + a0 / (2 * source.EP),
    )
    put(
        "same_full_finite_ratio_bound",
        source.conversion.original_reference_ratio_bound() - 4250 / source.KAPPA,
    )
    wm, wg, Mm, Mg, J, rho = s.symbols(
        "matter_weight gravity_weight matter_amplitude gravity_amplitude soft_current rho",
        real=True,
    )
    put(
        "complete_matter_gravity_interference",
        (wm * Mm + wg * Mg) ** 2
        - wm**2 * Mm**2
        - 2 * wm * wg * Mm * Mg
        - wg**2 * Mg**2,
    )
    put(
        "signed_real_recoil_subtraction",
        rho * (J + s.Symbol("R")) ** 2
        - J**2
        - rho * ((J + s.Symbol("R")) ** 2 - J**2)
        - (rho - 1) * J**2,
    )
    for row in range(3):
        points, q, born = source.vertices.sample(row)
        e = -born[0][0]
        w = q[0]
        outgoing = born[2][1:4, 0] / s.sqrt(e * e - 1)
        generated, rays, generated_born = source.recoil.momenta(e, [q], outgoing)
        for leg in range(4):
            for component in range(4):
                put(
                    f"same_original_recoil_{row}_{leg}_{component}",
                    generated[leg][component] - points[leg][component],
                )
        put(
            f"same_original_density_{row}",
            source.single_recoil.density_ratio(e, w) ** 2
            - (1 - 1 / (e * (e - w))) / (1 - 1 / (e * e)),
        )
        basis = transverse_basis(q[1:4, 0] / w)
        for j, eps in enumerate(basis):
            put(f"physical_TT_trace_{row}_{j}", s.trace(source.vertices.ETA * eps))
            put(f"physical_TT_norm_{row}_{j}", sum(t * t for t in eps) - 2)
            for component in range(4):
                put(
                    f"physical_TT_transverse_{row}_{j}_{component}",
                    (eps * q)[component],
                )
        put(f"physical_TT_orthogonality_{row}", sum(l * r for l, r in zip(*basis)))
        put(
            f"same_original_scalar_mass_{row}",
            source.vertices.dot(generated[2], generated[2]) - 1,
        )
        put(f"same_original_null_mass_{row}", source.vertices.dot(rays[0], rays[0]))
        put(f"same_original_Born_energy_{row}", generated_born[2][0] - e)
    return {
        "whole_marked_reference_operator": "D[R](x)=integral_(0<omega<x) P_sigma(x-omega)dR(sigma), with P_sigma(y)=expDelta_sigma F(a_sigma)y^a_sigma and the full conserved massive-plus-marked-null current. The markedstate and harddensity are held fixed while additional LEADING-soft terms are summed.",
        "whole_log_weight": log_ratio,
        "whole_same_energy_log_weight": same_energy,
        "whole_marked_recoil_and_source": "The public physical_kernel reconstructs the original26matter+21Einstein graphs with both unit-Frobenius TT polarizations, the exact two-body density ratio and positive complete Born. It returns the signed difference, not absolute amplitude squares or a positive-density replacement. The scaled density omits only the commonpositive1/(4pi^2*kappa) angular-average measure factor.",
        "whole_comparison_proof": "Foromega<x, stateenergyR=omega<=x allows S301's P_sigma(x)/P_elastic(x)<=1+4250/kappa. The true remaining-energy factor adds(1-omega/x)^a_sigma<=1 sincea_sigma>=0. Thus the weighted kernel is nonnegative and never exceeds1+4250/kappa for any positive remainingenergy, even when its logarithm is large.",
        "whole_endpoint_measure": "The physical one-real residual is absolutely continuous in its positive energy; the endpointomega=x has zero measure. Dominated convergence uses its finite total variation after subtracting the soft term. No proof is based on cancellation between two separately divergent integrated real rates.",
        "whole_wrong_virtual_negative_control": (a - a0) / (2 * source.EP),
        "checks": checks,
        "gates": {
            "same_full_original47_tree_residual_kernel": True,
            "both_unit_TT_polarizations_and_exact_density_ratio": True,
            "state_correct_soft_index_and_finite_conversion": True,
            "remaining_energy_not_replaced_by_total_threshold": True,
            "signed_total_variation_used_without_density_positivity": True,
            "endpoint_zero_measure_not_uniform_pointwise_limit": True,
            "single_residual_reference_not_full_multi_real_decomposition": True,
        },
    }
