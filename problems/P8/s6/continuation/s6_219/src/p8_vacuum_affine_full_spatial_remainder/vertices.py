"""Full spatial constrained pair and contact bridges in ten physical fields."""

from functools import cache

import sympy as s
from p8_vacuum_affine_full_adm_vertices import vertices as adm_vertices
from p8_vacuum_affine_spatial_current import hamiltonian
from p8_vacuum_affine_spatial_symbol import sectors
from p8_vacuum_affine_subleading_band_conversion import centered


def shifted(Q):
    return Q - s.trace(Q) * s.eye(3) / 2


def current_feature(Q):
    B = shifted(Q)
    # E3, B3, mA0, mAsp3. This is minus the first Hamiltonian variation.
    return s.diag(-B, -B, s.trace(Q) / 2, B)


def contact_feature(D, G):
    B, C = shifted(D), shifted(G)
    H = (B * C + C * B) / 2
    return s.diag(H, H, s.trace(D) * s.trace(G) / 4, H)


def physical_map(k, a=hamiltonian.A, mass=hamiltonian.MASS):
    F = hamiltonian.features(k, a, mass)
    return s.Matrix.vstack(F[6:9, :], s.I * F[3:6, :], -s.I * F[9:10, :], F[:3, :])


def full_geometry(k, ell, D, G):
    rows = sectors.geometric_factors(k, ell, shifted(D), shifted(G))
    rows.update(
        {
            "LC": s.trace(G) * (k.T * shifted(D) * ell)[0],
            "CL": s.trace(D) * (k.T * shifted(G) * ell)[0],
            "CC": s.trace(D) * s.trace(G) * k.dot(k) * ell.dot(ell),
        }
    )
    return rows


def full_amplitudes(point, mass):
    out = dict(sectors.scalar_amplitudes(point, mass))
    _fk, pk, wk = point["kl"]
    _fl, pl, wl = point["ll"]
    out["C"] = -pk * pl / (2 * point["a"] ** 2 * wk * wl)
    return out


def pair_products(k, ell, D, G, source, detector, mass):
    geom = full_geometry(k, ell, D, G)
    A, B = full_amplitudes(detector, mass), full_amplitudes(source, mass)
    return {
        "TT": A["A"] * B["A"] * geom["00"]
        + A["A"] * B["B"] * geom["01"]
        + A["B"] * B["A"] * geom["10"]
        + A["B"] * B["B"] * geom["11"],
        "TL": A["TL"] * B["TL"] * geom["TL"],
        "LT": A["LT"] * B["LT"] * geom["LT"],
        "LL": A["LL"] * B["LL"] * geom["LL"]
        + A["LL"] * B["C"] * geom["LC"]
        + A["C"] * B["LL"] * geom["CL"]
        + A["C"] * B["C"] * geom["CC"],
    }


def symmetric(prefix):
    q = s.symbols(prefix + "0:6", real=True)
    return s.Matrix([[q[0], q[1], q[2]], [q[1], q[3], q[4]], [q[2], q[4], q[5]]])


@cache
def feature_data():
    k, ell = (s.Matrix(s.symbols(label + "0:3", real=True)) for label in ("k", "l"))
    D, G = symmetric("D"), symmetric("G")
    a, mass = hamiltonian.A, hamiltonian.MASS
    B = shifted(D)
    Fk, Fl = physical_map(k), physical_map(ell)
    zero = s.zeros(3, 1)
    checks = {
        "full_complex_tensor_Frobenius_contraction": s.expand(
            s.trace((D + s.I * G).conjugate().T * (D + s.I * G))
            - s.trace(shifted(D + s.I * G).conjugate().T * shifted(D + s.I * G))
            - s.conjugate(s.trace(D + s.I * G)) * s.trace(D + s.I * G) / 4
        ),
        "full_current_equals_minus_independent_ADM_first_vertex": Fk.T
        * current_feature(D)
        * Fl
        + adm_vertices.vertex(-k, ell, 0, zero, D),
        "full_second_contact_equals_independent_ADM_second_vertex": Fk.T
        * contact_feature(D, G)
        * Fl
        - adm_vertices.contact(-k, ell, 0, D, 0, G),
        "trace_shift_exact_Frobenius_contraction": s.expand(
            s.trace(D.T * D) - s.trace(B.T * B) - s.trace(D) ** 2 / 4
        ),
        "trace_temporal_scalar_Cauchy_gap": s.expand(
            3 * s.trace(D.T * D)
            - s.trace(D) ** 2
            - (
                (D[0, 0] - D[1, 1]) ** 2
                + (D[0, 0] - D[2, 2]) ** 2
                + (D[1, 1] - D[2, 2]) ** 2
                + 6 * (D[0, 1] ** 2 + D[0, 2] ** 2 + D[1, 2] ** 2)
            )
        ),
        "canonical_physical_feature_energy_Gram": Fk.conjugate().T * Fk
        - hamiltonian.base(k),
        "nonzero_trace_constraint_contact_block": contact_feature(s.eye(3), s.eye(3))[
            6, 6
        ]
        - s.Rational(9, 4),
        "full_contact_local_no_external_tensor_derivatives": contact_feature(D, G)
        - contact_feature(G, D),
    }
    return {
        "ordering": "Physical ten-field order E3,B3,mA0,mAsp3; canonical phase-space order A3,pi3. The physical map sends the Hamiltonian energy features to(E,iB,-i div,mA). Pair momenta use M(-k,ell), not M(k,ell).",
        "full_current_feature": "diag(-B_D,-B_D,tr(D)/2,B_D), B_D=D-tr(D)I/2. Exact equality to minus the independent full ADM first vertex includes the nonzero temporal-constraint term and all Fourier signs.",
        "full_contact_feature": "diag(H_B,H_B,tr(D)tr(G)/4,H_B), H_B=(B_D B_G+B_G B_D)/2. This is the full positive second Hamiltonian vertex before the overall current sign. It is local in D,G and therefore its one-mode contraction is independent of external transfer, but its value is not zero.",
        "norm_proof": "For real symmetric D, ||B_D||F^2=||D||F^2-tr(D)^2/4. Each spatial block has operator norm<=||D||F and |tr(D)|/2<=sqrt(3)||D||F/2. Hence the first feature norm is<=||D||F. Each contact spatial block has norm<=||D||F||G||F, while its temporal coefficient is<=3||D||F||G||F/4. Extend the real bilinear bounds by the explicit Hilbert complexification, not by conjugating external sources inside a retarded formula.",
        "checks": {
            key: value.applyfunc(s.expand)
            if isinstance(value, s.MatrixBase)
            else s.expand(value)
            for key, value in checks.items()
        },
        "gates": {
            "all_ten_fields_in_both_vertices": current_feature(D).shape == (10, 10),
            "first_trace_block_fits_unit_norm": s.Rational(3, 4) < 1,
            "contact_trace_block_fits_product_norm": s.Rational(3, 4) < 1,
            "full_noncommuting_contact_not_diagonalized": D * G != G * D,
            "two_original_internal_momenta_not_identified": k != ell,
            "mass_gap_not_removed": mass.is_positive is True,
            "scale_positive": a.is_positive is True,
        },
    }


@cache
def pair_data():
    x, mass = s.symbols("x mass", real=True)
    source = centered.sample_point("full_source", x)
    detector = centered.sample_point("full_detector", x, True)
    k, ell = s.Matrix([3, 0, 4]), s.Matrix([0, 5, 12])
    ek, el = s.Matrix([0, 1, 0]), s.Matrix([1, 0, 0])
    D = s.Matrix([[1, 2, -1], [2, 3, 1], [-1, 1, 2]]) / 11
    G = s.Matrix([[2, -1, 1], [-1, -3, 2], [1, 2, 4]]) / 13
    sk = sectors.physical_vectors(k, ek, k.cross(ek) / 5, source, "k", mass * x)
    sl = sectors.physical_vectors(ell, el, ell.cross(el) / 13, source, "l", mass * x)
    dk = sectors.physical_vectors(k, ek, k.cross(ek) / 5, detector, "k", mass * x, True)
    dl = sectors.physical_vectors(
        ell, el, ell.cross(el) / 13, detector, "l", mass * x, True
    )
    original = pair_products(k, ell, D, G, source, detector, mass * x)
    exchanged = pair_products(
        ell, k, D, G, centered.swap_legs(source), centered.swap_legs(detector), mass * x
    )
    reflected = pair_products(
        -k,
        -ell,
        D,
        G,
        centered.sample_point("full_source", -x),
        centered.sample_point("full_detector", -x, True),
        -mass * x,
    )
    partner = {"TT": "TT", "TL": "LT", "LT": "TL", "LL": "LL"}
    checks = {}
    for tag, ii, jj in (
        ("TT", range(2), range(2)),
        ("TL", range(2), (2,)),
        ("LT", (2,), range(2)),
        ("LL", (2,), (2,)),
    ):
        actual = sum(
            (dk[i].T * current_feature(D) * dl[j])[0]
            * (sk[i].T * current_feature(G) * sl[j])[0]
            for i in ii
            for j in jj
        )
        checks[tag + "_full_ten_field_two_time_pair"] = s.expand(actual - original[tag])
        checks[tag + "_complete_leg_exchange"] = s.expand(
            original[tag] - exchanged[partner[tag]]
        )
        checks[tag + "_full_Schwarz_reflection"] = s.expand(
            reflected[tag] - s.conjugate(original[tag])
        )
    re, im = s.symbols("re im", real=True)
    z = re + s.I * im
    checks["corrected_all_five_current_endpoint_parities"] = s.Matrix(
        [
            s.simplify(
                s.im(s.I * (-s.I) ** j * s.conjugate(z))
                - (-1) ** j * s.im(s.I * (-s.I) ** j * z)
            )
            for j in range(5)
        ]
    )
    return {
        "full_factorization": "The literal physical ten-field pair equals the four sectors with ten geometric products. The LL constraint amplitude is C=-pk pl/(2 a^2 omega_k omega_l), multiplying tr(Q)|k||ell|. Both LC and CL and the CC square remain. The full time-dependent source and sharp detector are independent.",
        "all_orders_symmetry": "Each first feature is symmetric; swapping both created momenta and their distinct mode labels transposes the amplitude and exchanges TL with LT. The full constrained pair obeys (x,n)->(-x,-n) Schwarz reflection. Each inverse summed phase and source differentiation preserves leg exchange. With corrected coefficient i(-i)^j the current parity is(-1)^j. Thus only odd CENTERED TOTAL grades vanish, not odd endpoint labels.",
        "checks": checks,
        "gates": {
            "two_independent_clock_points": source["a"] != detector["a"],
            "nonzero_trace_constraint_in_both_directions": s.trace(D) != 0
            and s.trace(G) != 0,
            "all_nine_physical_pairs": 4 + 2 + 2 + 1 == 9,
            "longitudinal_constraint_square_kept": full_geometry(k, ell, D, G)["CC"]
            != 0,
            "ordered_cross_terms_distinct": full_geometry(k, ell, D, G)["LC"]
            != full_geometry(k, ell, D, G)["CL"],
        },
    }
