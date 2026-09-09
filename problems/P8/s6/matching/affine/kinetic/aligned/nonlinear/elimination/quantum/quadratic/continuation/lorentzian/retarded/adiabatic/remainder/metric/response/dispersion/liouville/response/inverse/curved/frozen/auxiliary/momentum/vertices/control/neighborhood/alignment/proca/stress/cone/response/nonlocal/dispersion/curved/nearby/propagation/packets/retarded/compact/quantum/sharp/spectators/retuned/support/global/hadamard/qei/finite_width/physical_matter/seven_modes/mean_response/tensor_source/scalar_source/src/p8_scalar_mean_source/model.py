"""Candidate off-background scalar quadratic action; native bridge under development."""

from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model as parent
from p8_proca_global_support import model as global_model
from p8_proca_global_support import phase
from p8_proca_physical_matter import coefficients as actual

u = parent.u
A = sp.Symbol("independent_positive_hat_scale", positive=True)
k = sp.Symbol("nonzero_comoving_scalar_momentum", positive=True)
p, l = sp.symbols(
    "independent_mean_trace independent_mean_matter_charge_density", real=True
)
v, chi, Pv, Pc = sp.symbols(
    "natural_log_curvature natural_matter natural_curvature_density_momentum natural_matter_density_momentum",
    real=True,
)
FIELDS = (v, chi, Pv, Pc)


def unit(i, n=1):
    return tuple(n if j == i else 0 for j in range(9))


@cache
def coefficients():
    rows = actual.data()["rows"]
    bg = parent.coefficients()["background"]
    p0 = -2 * bg["H"]
    l0 = bg["ell"]
    aa = rows[unit(0, 2)]
    dd = rows[unit(1, 2)]
    bb = tuple(sp.factor(rows[unit(0)][j] - 2 * aa[j] * p0) for j in range(5))
    cc = tuple(
        sp.factor(rows[(0,) * 9][j] - aa[j] * p0 * p0 - bb[j] * p0 - dd[j] * l0 * l0)
        for j in range(5)
    )
    return {
        "a": aa,
        "b": bb,
        "c": cc,
        "d": dd,
        "e": rows[unit(3)],
        "g": rows[unit(7)],
        "f": rows[unit(8)],
    }


@cache
def generic():
    co = {name: sp.symbols("scalar_" + name + "0:4", real=True) for name in "abcdefg"}
    a, b, c, d, e, f, g = (co[name] for name in "abcdefg")
    r, s = Pv / A**3, Pc / A**3
    q = k * k / A**2
    Q = (
        a[0] * (sp.Rational(9, 2) * p * p * v * v - 2 * p * r * v + r * r / 9)
        + sp.Rational(9, 2) * c[0] * v * v
    )
    Q += d[0] * (sp.Rational(9, 2) * l * l * v * v - 6 * l * v * s + s * s)
    Q += (
        e[0] * (r - 9 * p * v - 3 * l * chi) ** 2 / 24
        + g[0] * q * chi * chi
        + 2 * f[0] * q * v * v
    )
    force = (
        (-3 * a[1] * p * p + 3 * c[1] - 3 * d[1] * l * l + 4 * f[1] * q) * v
        + (2 * a[1] * p + b[1]) * r / 3
        + 2 * d[1] * l * s
    )
    pivot = a[2] * p * p + b[2] * p + c[2] + d[2] * l * l
    shift = {co[name][j]: co[name][j + 1] for name in co for j in range(3)}
    deriv = lambda val: sum(sp.diff(val, key) * value for key, value in shift.items())
    H = A**3 * (Q - force * force / (2 * pivot))
    nsource = A**3 * (
        deriv(Q)
        - force * deriv(force) / pivot
        + force * force * deriv(pivot) / (2 * pivot * pivot)
    )
    return {
        "symbols": co,
        "raw_density_quadratic_Hamiltonian": H,
        "uneliminated_quadratic_density": Q,
        "linear_inhomogeneous_lapse_force": force,
        "independent_mean_lapse_pivot": pivot,
        "actual_mean_lapse_source": nsource,
        "mean_hat_scale_direct_source": sp.diff(H, p) / 3,
        "mean_trace_direct_source": -A * sp.diff(H, A) / 3 + l * sp.diff(H, l),
        "mean_matter_field_direct_source": sp.diff(H, l),
    }


@cache
def substitution():
    bg = parent.coefficients()["background"]
    co = coefficients()
    g = generic()
    return {
        **{g["symbols"][name][j]: co[name][j] for name in co for j in range(4)},
        p: -2 * bg["H"],
        l: bg["ell"],
        A: (1 + u * u) ** 2,
    }


@cache
def on_clock():
    return {
        key: sp.factor(value.subs(substitution(), simultaneous=True))
        for key, value in generic().items()
        if key != "symbols"
    }


@cache
def bridge():
    bg = parent.coefficients()["background"]
    a = (1 + u * u) ** 2
    pv, ps = phase.pv, phase.ps
    ell, H = bg["ell"], bg["H"]
    change = {
        Pv: pv + 3 * a**3 * ell * chi - 18 * a**3 * H * v,
        Pc: ps + 3 * a**3 * ell * v,
    }
    original = global_model.old
    target = phase.data()["regular_density_Hamiltonian"].subs(
        global_model.substitution(), simultaneous=True
    )
    target = target.subs(
        {original.q: k * k / a**2, original.v: v, original.matter: chi},
        simultaneous=True,
    )
    generating = -9 * a**3 * (sp.diff(H, u) + 3 * H * H) * v * v
    raw = on_clock()["raw_density_quadratic_Hamiltonian"]
    difference = sp.factor(raw.subs(change, simultaneous=True) + generating - target)
    return {
        "difference": difference,
        "generating_time_derivative": generating,
        "raw_to_existing_density_phase": change,
        "target": target,
        "raw_clock_Hamiltonian": raw,
    }


@cache
def checks():
    co = coefficients()
    rows = actual.data()["rows"]
    bg = parent.coefficients()["background"]
    p0 = -2 * bg["H"]
    l0 = bg["ell"]
    checks = {
        "actual_scalar_background_potential_jets_reconstruct_literal_retuned_rows": sp.Matrix(
            [
                sp.factor(
                    co["a"][j] * p0 * p0
                    + co["b"][j] * p0
                    + co["c"][j]
                    + co["d"][j] * l0 * l0
                    - rows[(0,) * 9][j]
                )
                for j in range(5)
            ]
        ),
        "actual_matter_charge_coefficient_jets_reconstruct_literal_rows": sp.Matrix(
            [sp.factor(2 * co["d"][j] * l0 - rows[unit(1)][j]) for j in range(5)]
        ),
        "actual_mean_lapse_pivot_matches_original_retuned_positive_J": sp.factor(
            on_clock()["independent_mean_lapse_pivot"] + 2 * actual.data()["Jnew"]
        ),
    }
    g = generic()
    H = g["raw_density_quadratic_Hamiltonian"]
    offset = sp.Symbol("independent_background_lapse_offset", real=True)
    translated = {
        values[j]: sum(
            values[h] * offset ** (h - j) / sp.factorial(h - j) for h in range(j, 4)
        )
        for values in g["symbols"].values()
        for j in range(3)
    }
    direct = sp.diff(H.subs(translated, simultaneous=True), offset).subs(offset, 0)
    checks["literal_background_lapse_differentiation_retains_third_pivot_jet"] = (
        sp.factor(direct - g["actual_mean_lapse_source"])
    )
    PB, Q = sp.symbols(
        "independent_canonical_homogeneous_trace independent_canonical_homogeneous_matter_charge",
        real=True,
    )
    canonical = H.subs({p: PB / (3 * A**3), l: Q / A**3}, simultaneous=True)
    back = {PB: 3 * A**3 * p, Q: A**3 * l}
    scale = sp.diff(canonical, PB).subs(back, simultaneous=True)
    trace = (-A * sp.diff(canonical, A) / (3 * A**3)).subs(
        back, simultaneous=True
    ) - 3 * p * scale
    field = sp.diff(canonical, Q).subs(back, simultaneous=True)
    checks.update(
        {
            "actual_scalar_mean_scale_source_from_canonical_homogeneous_momentum": sp.factor(
                scale - g["mean_hat_scale_direct_source"] / A**3
            ),
            "actual_scalar_mean_trace_source_keeps_volume_and_charge_chain": sp.factor(
                trace - g["mean_trace_direct_source"] / A**3
            ),
            "actual_scalar_mean_matter_field_source_from_fixed_canonical_charge": sp.factor(
                field - g["mean_matter_field_direct_source"] / A**3
            ),
        }
    )
    return checks
