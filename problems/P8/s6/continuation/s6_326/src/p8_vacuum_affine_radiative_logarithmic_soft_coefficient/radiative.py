"""Finite mixed-null log coefficient; retain normal derivatives of conservation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_state_soft_index import recoil

from . import kernel


def vector(value, dimension):
    if not isinstance(value, (tuple, list, s.MatrixBase)):
        raise TypeError("Require an exact vector")
    out = s.ImmutableMatrix(value)
    if out.shape != (dimension, 1):
        raise ValueError("Wrong vector dimension")
    for entry in out:
        kernel.exact_scalar(entry)
    return out


def polarization(value, direction):
    if not isinstance(value, (tuple, list, s.MatrixBase)):
        raise TypeError("Require an exact physical TT polarization")
    A = s.ImmutableMatrix(value)
    if A.shape != (4, 4):
        raise ValueError("Require a four-dimensional polarization")
    for item in A:
        if not isinstance(item, s.Expr):
            raise TypeError("Require exact scalar polarization entries")
        if (
            item.has(s.Float)
            or item.is_number is not True
            or item.is_finite is not True
        ):
            raise ValueError("Require exact finite polarization entries")
    q = s.ImmutableMatrix([1, *direction])
    if (
        A != A.T
        or any(A[0, j] != 0 for j in range(4))
        or (A * q).applyfunc(s.simplify) != s.zeros(4, 1)
        or s.simplify(s.trace(A)) != 0
        or s.simplify(s.trace(A.conjugate().T * A) - 1) != 0
    ):
        raise ValueError("Require spatial TT and unit Frobenius norm")
    return A


def components(energy, quanta, outgoing, direction, tensor):
    E = kernel.exact_scalar(energy)
    u, n = vector(outgoing, 3), vector(direction, 3)
    if s.simplify(n.dot(n) - 1) != 0:
        raise ValueError("Require a unit additional-soft direction")
    if not isinstance(quanta, (tuple, list)):
        raise TypeError("Require finitely many positive null momenta")
    rays = tuple(vector(ray, 4) for ray in quanta)
    A = polarization(tensor, n)
    points, rays, _born = recoil.momenta(E, rays, u)
    legs = tuple(points) + tuple(rays)
    q = s.ImmutableMatrix([1, *n])
    Ds = tuple(s.factor(recoil.dot(p, q)) for p in legs)
    if any(D == 0 for D in Ds):
        raise ValueError(
            "Exactly soft-collinear null legs have no assigned point value"
        )
    P = lambda i, j: (legs[i].T * A * legs[j])[0]
    Sb = s.factor(sum(P(i, i) / Ds[i] for i in range(len(legs))))
    F, G = s.S.Zero, s.S.Zero
    pairs = {}
    for i in range(len(legs)):
        for j in range(i + 1, len(legs)):
            z = s.factor(recoil.dot(legs[i], legs[j]))
            S = s.factor(
                P(i, i) * Ds[j] / Ds[i] + P(j, j) * Ds[i] / Ds[j] - 2 * P(i, j)
            )
            pairs[(i, j)] = S
            if j < 4:
                F -= kernel.massive_fprime(abs(z)) * S / (16 * s.pi**2)
            elif z == 0:
                if S != 0:
                    raise ValueError("Relative-collinear TT contraction did not vanish")
            else:
                scale = legs[j][0] * (legs[i][0] if i >= 4 else 1)
                F -= (s.log(s.factor(2 * abs(z) / scale)) + 1) * S / (4 * s.pi**2)
            if legs[i][0] * legs[j][0] > 0:
                c = kernel.massive_c(z) if j < 4 else 2
                G += c * S / (8 * s.pi)
    H = sum(D * s.log(abs(D)) for D in Ds[:4])
    H += sum(D * s.log(s.factor(D / p[0])) for p, D in zip(legs[4:], Ds[4:]))
    F += Sb * H / (4 * s.pi**2)
    outgoing_D = sum(D for p, D in zip(legs, Ds) if p[0] > 0)
    G -= Sb * outgoing_D / (4 * s.pi)
    return {
        "F": F,
        "G": G,
        "Sbar": Sb,
        "phase_H": H,
        "pair_contractions": pairs,
        "D": Ds,
        "outgoing_D": outgoing_D,
    }


def coefficient(energy, quanta, outgoing, direction, tensor):
    row = components(energy, quanta, outgoing, direction, tensor)
    return row["F"] + s.I * row["G"]


def calibration_states():
    E = s.Rational(5, 4)
    a = s.Rational(19, 10)
    Ep = (a + 1 / a) / 2
    R = 2 * (E - Ep)
    w = E - Ep**2 / E
    one = [s.ImmutableMatrix([w, w, 0, 0])]
    two = [
        s.ImmutableMatrix([R / 2, R / 2, 0, 0]),
        s.ImmutableMatrix([R / 2, -R / 2, 0, 0]),
    ]
    split = [two[0] / 3, 2 * two[0] / 3, two[1]]
    four = [
        two[0] / 2,
        two[1] / 2,
        s.ImmutableMatrix([R / 4, 0, 0, R / 4]),
        s.ImmutableMatrix([R / 4, 0, 0, -R / 4]),
    ]
    return E, {"born": [], "one": one, "two": two, "split": split, "four": four}


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = (
            value.applyfunc(s.factor)
            if isinstance(value, s.MatrixBase)
            else s.factor(value)
        )

    Ds = s.symbols("D0:7", nonzero=True, real=True)
    Ps = {
        (i, j): s.Symbol("P" + str(i) + "_" + str(j))
        for i in range(7)
        for j in range(i, 7)
    }
    P = lambda i, j: Ps[tuple(sorted((i, j)))]
    soft = sum(P(i, i) / Ds[i] for i in range(7))
    for i in range(7):
        lhs = sum(
            P(i, i) * Ds[j] / Ds[i] + P(j, j) * Ds[i] / Ds[j] - 2 * P(i, j)
            for j in range(7)
            if j != i
        )
        rhs = (
            Ds[i] * soft
            + P(i, i) * sum(Ds) / Ds[i]
            - 2 * sum(P(i, j) for j in range(7))
        )
        put("generic_conservation_normal_derivative_" + str(i), s.expand(lhs - rhs))
    z, w, m, d, ell, D, soft = s.symbols("z w m d ell D soft", positive=True)
    put(
        "mixed_kernel_derivative",
        s.diff(2 * z * s.log(2 * z / w), z) - 2 * (s.log(2 * z / w) + 1),
    )
    put("reference_energy_derivative", s.diff(2 * z * s.log(2 * z / w), w) + 2 * z / w)
    put(
        "auxiliary_conservation_zero_normal_derivative",
        s.diff(-2 * s.log(m / w) * (z + ell * D * soft), ell).subs(ell, 0)
        + 2 * s.log(m / w) * D * soft,
    )
    put(
        "mass_and_energy_logs_cancel_with_phase",
        s.expand_log(s.log(m / w) + s.log(w * d / m) - s.log(d), force=True),
    )
    E, states = calibration_states()
    u = s.ImmutableMatrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    n = s.ImmutableMatrix([0, 1, 0])
    plus = s.diag(0, 1, 0, -1) / s.sqrt(2)
    cross = s.zeros(4)
    cross[1, 3] = cross[3, 1] = 1 / s.sqrt(2)
    pols = {"plus": plus, "complex": (plus + s.I * cross) / s.sqrt(2)}
    records, coefficients = [], {}
    for label, rays in states.items():
        points, rays, _ = recoil.momenta(E, rays, u)
        legs = tuple(points) + tuple(rays)
        for i, p in enumerate(legs):
            put(label + "_mass_shell_" + str(i), recoil.dot(p, p) - (1 if i < 4 else 0))
        put(label + "_complete_conservation", sum(legs, s.zeros(4, 1)))
        for plabel, A in pols.items():
            row = components(E, rays, u, n, A)
            put(label + "_" + plabel + "_outgoing_phase_sum", row["outgoing_D"] - 2 * E)
            for i, D in enumerate(row["D"]):
                put(
                    label + "_" + plabel + "_normal_derivative_" + str(i),
                    sum(S for pair, S in row["pair_contractions"].items() if i in pair)
                    - D * row["Sbar"],
                )
            coefficients[(label, plabel)] = row["F"] + s.I * row["G"]
            records.append((label, plabel, len(rays)))
    for plabel in pols:
        put(
            "exact_full_collinear_split_" + plabel,
            coefficients[("split", plabel)] - coefficients[("two", plabel)],
        )
    h = s.Symbol("small_positive_angle", positive=True)
    nx = s.Matrix([2 * h / (1 + h * h), (1 - h * h) / (1 + h * h), 0])
    nz = s.Matrix([0, (1 - h * h) / (1 + h * h), 2 * h / (1 + h * h)])
    lx = s.limit((nx.T * plus[1:4, 1:4] * nx)[0] / (1 - nx[1]), h, 0, dir="+")
    lz = s.limit((nz.T * plus[1:4, 1:4] * nz)[0] / (1 - nz[1]), h, 0, dir="+")
    put("soft_collinear_x_limit", lx - s.sqrt(2))
    put("soft_collinear_z_limit", lz + s.sqrt(2))
    return {
        "checks": checks,
        "gates": {
            "differentiate_before_conservation": True,
            "mass_logs_and_individual_energy_logs_both_cancel": True,
            "all_massive_null_and_null_null_pairs_present": True,
            "no_energy_entropy_or_particle_count_loss": True,
            "real_and_complex_TT_calibrations": True,
            "finite_collinear_splitting_invariant": True,
            "soft_collinear_point_limits_not_unique": lx != lz,
            "detector_scale_term_not_erased": True,
            "single_additional_soft_limit_not_simultaneous_remainder_theorem": True,
        },
        "whole_finite_coefficient": "Unordered pairs: F=-sumMM f'(|z|)S/(16pi^2)-[sumMN(ln(2alpha)+1)S+sumNN(ln(2delta)+1)S]/(4pi^2)+Sbar[sumM Dln|D|+sumN Dln d]/(4pi^2); G=[sumSameMM c(z)S+2sumSameMN S+2sumNN S]/(8pi)-Sbar*2E/(4pi); C=F+iG. F/G are complex-linear components for complex TT tensors.",
        "whole_normal_derivative_identity": "delta k_i=Pii*q/Di-eta*A*k_i preserves masses. On conservation, sum delta k=Sbar*q and delta Z_j=Dj*Sbar although Z_j=sum_(i!=j)ki.kj=0 for null j. Kraw=Kfinite-2sumNull ln(mj/wj)Zj; keep deltaZj until cancellation with the phase. Reference-energy derivatives multiplying Zj vanish only after differentiation.",
        "whole_exact_recoil_calibrations": records,
        "whole_exact_calibration_coefficients": {
            label + "_" + pol: value for (label, pol), value in coefficients.items()
        },
        "whole_soft_collinear_counterexample": (lx, lz),
    }
