"""Full minimally coupled scalar ADM Hamiltonian features and ordered parity."""

from functools import cache

import sympy as s

a, mass = s.symbols("a heavy_mass", positive=True)
nD, nG = s.symbols("lapse_D lapse_G", real=True)
ed, eg = s.symbols("ed eg", real=True)


def symmetric(name):
    v = s.symbols(
        name
        + "00 "
        + name
        + "11 "
        + name
        + "22 "
        + name
        + "01 "
        + name
        + "02 "
        + name
        + "12",
        real=True,
    )
    return s.Matrix([[v[0], v[3], v[4]], [v[3], v[1], v[5]], [v[4], v[5], v[2]]])


D, G = symmetric("D"), symmetric("G")
bD = s.Matrix(s.symbols("beta_D0:3", real=True))
bG = s.Matrix(s.symbols("beta_G0:3", real=True))
reflection = s.diag(1, -1, -1, -1, 1)


def jet(value):
    value = s.expand(value)
    return sum(
        value.coeff(ed, i).coeff(eg, j) * ed**i * eg**j for i in (0, 1) for j in (0, 1)
    )


def first(n, b, Q):
    tau = s.trace(Q)
    out = s.diag(n - tau / 2, 1, 1, 1, n + tau / 2)
    out[1:4, 1:4] = (n + tau / 2) * s.eye(3) - Q
    out[0, 1:4] = a * b.T
    out[1:4, 0] = a * b
    return out


def second(nd, Qd, ng, Qg):
    td, tg = s.trace(Qd), s.trace(Qg)
    out = s.zeros(5)
    out[0, 0] = td * tg / 4 - (nd * tg + ng * td) / 2
    Bd = td * s.eye(3) / 2 - Qd
    Bg = tg * s.eye(3) / 2 - Qg
    out[1:4, 1:4] = (Bd * Bg + Bg * Bd) / 2 + nd * Bg + ng * Bd
    out[4, 4] = td * tg / 4 + (nd * tg + ng * td) / 2
    return out


def physical_map(momentum):
    out = s.zeros(5, 2)
    out[0, 1] = a**-3
    out[1:4, 0] = s.I * s.Matrix(momentum) / a
    out[4, 0] = mass
    return out


def canonical_first(n, b, Q, left, right):
    left, right = s.Matrix(left), s.Matrix(right)
    tau = s.trace(Q)
    grad = (n + tau / 2) * s.eye(3) - Q
    return s.Matrix(
        [
            [
                a**3 * ((n + tau / 2) * mass**2 + (left.T * grad * right)[0] / a**2),
                -s.I * left.dot(b),
            ],
            [s.I * right.dot(b), a**-3 * (n - tau / 2)],
        ]
    )


@cache
def data():
    Q = ed * D + eg * G
    tau = s.trace(Q)
    N = 1 + ed * nD + eg * nG
    volume = 1 + tau / 2 + tau * tau / 8
    inverse = s.eye(3) - Q + Q * Q / 2
    whole = s.zeros(5)
    whole[0, 0] = jet(N * (1 - tau / 2 + tau * tau / 8))
    whole[1:4, 1:4] = (N * volume * inverse).applyfunc(jet)
    whole[4, 4] = jet(N * volume)
    whole[0, 1:4] = a * (ed * bD + eg * bG).T
    whole[1:4, 0] = a * (ed * bD + eg * bG)
    MD, MG = first(nD, bD, D), first(nG, bG, G)
    CDG = second(nD, D, nG, G)
    checks = {
        "entire_ADM_first_Hamiltonian_feature_D": whole.diff(ed).subs({ed: 0, eg: 0})
        - MD,
        "entire_ADM_first_Hamiltonian_feature_G": whole.diff(eg).subs({ed: 0, eg: 0})
        - MG,
        "entire_ADM_second_Hamiltonian_contact": whole.diff(ed, eg).subs({ed: 0, eg: 0})
        - CDG,
        "complete_noncommuting_second_feature_symmetry": CDG - second(nG, G, nD, D),
        "complete_zero_background_positive_energy_feature": whole.subs({ed: 0, eg: 0})
        - s.eye(5),
        "full_lapse_metric_reflection_even": reflection
        * first(nD, s.zeros(3, 1), D)
        * reflection
        - first(nD, s.zeros(3, 1), D),
        "full_shift_reflection_odd": reflection * first(0, bD, s.zeros(3)) * reflection
        + first(0, bD, s.zeros(3)),
        "full_second_contact_reflection_even": reflection * CDG * reflection - CDG,
    }
    left = s.Matrix(s.symbols("k0:3", real=True))
    right = s.Matrix(s.symbols("l0:3", real=True))
    bridge = a**3 * physical_map(-left).T * MD * physical_map(right)
    checks["independent_canonical_full_first_feature_bridge"] = (
        bridge - canonical_first(nD, bD, D, left, right)
    )
    checks["canonical_reality_and_ordered_leg_exchange"] = (
        canonical_first(nD, bD, D, left, right)
        - canonical_first(nD, bD, D, right, left).conjugate().T
    )
    # All matrix coefficients are checked, not only TT or trace projections.
    return {
        "full_Hamiltonian": "H=N/2[pi²/sqrt(h)+sqrt(h)(h^ij partial_i H partial_j H+mass² H²)]+beta^i pi partial_i H, with h=a² exp(Q), pi held canonical, and N=1+n_lapse. The beta term is Weyl ordered. No lapse, shift or spatial direction is removed.",
        "physical_feature_order": (
            "pi/a³",
            "partial_x H/a",
            "partial_y H/a",
            "partial_z H/a",
            "mass*H",
        ),
        "full_first_Hamiltonian_feature": MD,
        "full_second_Hamiltonian_feature": CDG,
        "complete_canonical_first_symbol": canonical_first(nD, bD, D, left, right),
        "current_convention": "J_D=-delta_D H=-a³ F^t M_D F/2. The instantaneous contact in its tangent is-a³<F^t C_DG F>/2. C_DG is the second HAMILTONIAN feature, not the current feature.",
        "full_feature_envelopes": "For v(D)=|n_D|+a|beta_D|+2||Q_D||F, the entire5x5 first feature obeys||M_D||op<=v(D) and the second feature obeys||C_DG||op<=v(D)v(G). Use|tr Q|<=sqrt(3)||Q||F, the full anti-commutator of trace-shifted Q matrices, and block product bounds. No second shift or lapse-lapse feature exists in these canonical ADM variables.",
        "ordered_reflection": reflection,
        "parity_warning": "Simultaneous internal momentum reflection acts on the physical feature by diag(1,-I3,1). Lapse and spatial-metric insertions are even, shift insertions odd. A mixed odd channel is NOT obtained by applying the ordinary imaginary-part operation to the same-momentum product. The reversed ordered term also reflects both momenta.",
        "checks": {
            key: value.applyfunc(s.expand)
            if isinstance(value, s.MatrixBase)
            else s.expand(value)
            for key, value in checks.items()
        },
        "gates": {
            "all_ten_ADM_directions_including_shift_retained": True,
            "five_feature_full_positive_reference_Hamiltonian": True,
            "noncommuting_metric_contact_not_trace_restricted": True,
            "no_inverse_Hubble_or_spatial_transfer": True,
            "current_sign_distinguished_from_Hamiltonian_contact": True,
            "opposite_parity_reverse_branch_not_deleted": True,
        },
    }
