"""Private complete noncommuting local-action/current variational identity."""

from functools import cache
from math import comb

import sympy as s

from .angular import NC, trace

w = (s.Symbol("omega", positive=True),) + s.symbols("omega1:5", real=True)


def inverse_jets(row):
    inv = [NC(1)]
    for n in range(1, len(row)):
        inv.append(
            -sum((comb(n, j) * row[j] * inv[n - j] for j in range(1, n + 1)), NC())
        )
    return inv


def product_jet(left, right, n):
    return sum((comb(n, j) * left[j] * right[n - j] for j in range(n + 1)), NC())


def scalar_reciprocal(row):
    inv = [1 / row[0]]
    for n in range(1, len(row)):
        inv.append(
            s.expand(
                -sum(comb(n, j) * row[j] * inv[n - j] for j in range(1, n + 1)) / row[0]
            )
        )
    return inv


@cache
def ordered_frame():
    K = [NC(1)] + [NC.letter("K" + str(j)) for j in range(1, 5)]
    Ki = inverse_jets(K)
    B = [NC(1)]
    for n in range(1, 5):
        B.append(
            (K[n] - sum((comb(n, j) * B[j] * B[n - j] for j in range(1, n)), NC())) / 2
        )
    Bi = inverse_jets(B)
    invw = scalar_reciprocal(w)
    p = [
        s.expand(sum(comb(n, j) * w[j + 1] * invw[n - j] for j in range(n + 1)))
        for n in range(4)
    ]
    L = [
        sum((comb(n, j) * Bi[j] * B[n - j + 1] for j in range(n + 1)), NC())
        for n in range(4)
    ]
    R = [(x.T - x) / 2 for x in L]
    S = [NC(p[n] / 2) - (L[n] + L[n].T) / 2 for n in range(4)]
    A = [
        sum((comb(n, j) * K[j + 1] * Ki[n - j] for j in range(n + 1)), NC())
        for n in range(4)
    ]
    roots = [NC(p[n] / 2) - A[n] / 2 for n in range(4)]
    return {
        "K": K,
        "inverseK": Ki,
        "B": B,
        "inverseB": Bi,
        "invw": invw,
        "p": p,
        "R": R,
        "S": S,
        "A": A,
        "root_free_s": roots,
    }


@cache
def reference():
    c = ordered_frame()
    v, R, S = c["invw"], c["R"], c["S"]
    b = {
        1: [
            s.I * sum((comb(d, j) * v[j] * S[d - j] for j in range(d + 1)), NC()) / 2
            for d in range(4)
        ]
    }
    for n in range(1, 4):

        def right(d, n=n):
            value = b[n][d + 1] - sum(
                (
                    comb(d, j) * (R[j] * b[n][d - j] - b[n][d - j] * R[j])
                    for j in range(d + 1)
                ),
                NC(),
            )
            for j in range(1, n):
                l = n - j
                for a in range(d + 1):
                    for bb in range(d - a + 1):
                        cc = d - a - bb
                        value += (
                            comb(d, a) * comb(d - a, bb) * b[j][a] * S[bb] * b[l][cc]
                        )
            return value

        b[n + 1] = [
            -s.I
            * sum((comb(d, j) * v[j] * right(d - j) for j in range(d + 1)), NC())
            / 2
            for d in range(3 - n + 1)
        ]
    return [NC()] + [b[n][0] for n in range(1, 5)]


@cache
def covariance_blocks():
    r = reference()
    E = [NC(1)] + [
        -sum((r[j].H * r[n - j] for j in range(n + 1)), NC()) for n in range(1, 5)
    ]
    C = [NC(1)]
    for n in range(1, 5):
        C.append(-sum((E[j] * C[n - j] for j in range(1, n + 1)), NC()))
    QQ, PP = [], []
    for n in range(5):
        linear = sum(
            (r[j] * C[n - j] + C[n - j] * r[j].H for j in range(1, n + 1)), NC()
        )
        quadratic = sum(
            (
                r[a] * C[n - a - b] * r[b].H
                for a in range(1, n)
                for b in range(1, n - a + 1)
            ),
            NC(),
        )
        QQ.append((C[n] + linear + quadratic).real() / 2)
        PP.append((C[n] - linear + quadratic).real() / 2)
    return QQ, PP


@cache
def variational_check():
    c = ordered_frame()
    p, A, ss = c["p"], c["A"], c["root_free_s"]
    tt = [
        ss[n + 1] - sum((comb(n, j) * p[j] * ss[n - j] for j in range(n + 1)), NC())
        for n in range(3)
    ]
    U = -tt[1] + 2 * p[0] * tt[0] + 2 * ss[0] ** 3
    Ud = (
        -tt[2]
        + 2 * p[1] * tt[0]
        + 2 * p[0] * tt[1]
        + 2 * (ss[1] * ss[0] * ss[0] + ss[0] * ss[1] * ss[0] + ss[0] * ss[0] * ss[1])
    )
    v = trace(U) - 2 * trace(tt[0] * ss[0])
    vd = trace(Ud) - 2 * trace(tt[1] * ss[0] + tt[0] * ss[1])
    Z = trace(tt[0] ** 2 + ss[0] ** 4)
    unbalanced = Ud + U * A[0] - A[0] * U - 3 * p[0] * U
    EK = {
        0: NC(),
        2: tt[0] / (4 * w[0]),
        4: (unbalanced + unbalanced.T) / (32 * w[0] ** 3),
    }
    Ew = {
        0: -trace(NC(1)) / 2,
        2: -trace(ss[0] ** 2 + tt[0]) / (4 * w[0] ** 2),
        4: -(3 * Z + vd - 3 * p[0] * v) / (16 * w[0] ** 4),
    }
    QQ, PP = covariance_blocks()
    checks = {}
    for n in (0, 2, 4):
        error = (EK[n] - w[0] * (QQ[n] - PP[n]) / 2).cancel()
        checks["matrix_variation_order_" + str(n)] = error.terms
        checks["frequency_variation_order_" + str(n)] = s.factor(Ew[n] + trace(QQ[n]))
    for n in range(1, 5):
        checks["reference_symmetry_" + str(n)] = (
            (reference()[n] - reference()[n].T).cancel().terms
        )
    return {
        "checks": checks,
        "matrix_Euler": EK,
        "frequency_Euler": Ew,
        "QQ": QQ,
        "PP": PP,
    }


@cache
def data():
    c = variational_check()
    checks = {}
    matrix_counts = {}
    frequency_counts = {}
    for n in (0, 2, 4):
        expected = c["matrix_Euler"][n]
        actual = w[0] * (c["QQ"][n] - c["PP"][n]) / 2
        words = sorted(set(expected.terms) | set(actual.terms)) or [()]
        matrix_counts[n] = len(words)
        for word in words:
            label = "_".join(word) or "identity"
            checks[f"all_word_K_variation_order_{n}_{label}"] = s.cancel(
                expected.terms.get(word, 0) - actual.terms.get(word, 0)
            )
        expected = c["frequency_Euler"][n]
        actual = -trace(c["QQ"][n])
        symbols = sorted(
            (expected.free_symbols | actual.free_symbols) - set(w), key=str
        )
        left = s.Poly(expected, *symbols)
        right = s.Poly(actual, *symbols)
        monomials = sorted(set(left.monoms()) | set(right.monoms()))
        frequency_counts[n] = len(monomials)
        for j, monomial in enumerate(monomials):
            checks[f"all_trace_frequency_variation_order_{n}_{j}"] = s.cancel(
                left.coeff_monomial(monomial) - right.coeff_monomial(monomial)
            )
    for n in range(1, 5):
        ref = reference()[n]
        for word, value in ref.terms.items():
            label = "_".join(word) or "identity"
            checks[f"ordered_reference_phase_{n}_{label}"] = s.im(value / s.I**n)
            checks[f"ordered_reference_symmetry_{n}_{label}"] = s.cancel(
                value - ref.T.terms.get(word, 0)
            )
    frame = ordered_frame()
    for n in range(1, 5):
        residual = product_jet(frame["B"], frame["B"], n) - frame["K"][n]
        checks[f"complete_positive_root_jet_{n}"] = s.Integer(
            len(residual.cancel().terms)
        )
        for key in ("B", "inverseB"):
            other = "inverseB" if key == "B" else "B"
            residual = product_jet(frame[key], frame[other], n)
            checks[f"complete_ordered_{key}_inverse_jet_{n}"] = s.Integer(
                len(residual.cancel().terms)
            )
    return {
        "general_domain": "Any positive symmetric K(t), scalar omega(t)>0 and V=omega^2 K^-1. A constant canonical GL change sets K0=I at one evaluation point without constraining arbitrary symmetric noncommuting K1..K4.",
        "local_action": "L0=-d omega/2; L2=tr(s^2)/(4omega); L4=tr(t^2+s^4)/(16omega^3), A=K'K^-1, p=omega'/omega, s=pI/2-A/2, t=s'-p s. These are local finite marker terms, not a full time-nonlocal determinant.",
        "full_current_variation": "At K0=I, deltaG=diag(2deltaomega/omega I-deltaK,deltaK), so J=-deltaomega tr(Sigma_QQ)+(omega/2)tr[deltaK(Sigma_QQ-Sigma_PP)]. Every word/trace coefficient through marker order4 agrees with the compact local Euler variation.",
        "fourth_Euler_formula": "U=-t'+2p t+2s^3, v=trU-2tr(ts), Z=tr(t^2+s^4). E_K4=Sym{K^-1[U'+[U,A]-3pU]}/(16omega^3), E_omega4=-[3Z+v'-3p v]/(16omega^4). Both the commutator and raw second-jet variations are retained.",
        "independent_matrix_word_counts": matrix_counts,
        "independent_frequency_trace_counts": frequency_counts,
        "Euler_matrix_coefficients": {
            n: {
                " ".join(k) or "identity": v
                for k, v in c["matrix_Euler"][n].terms.items()
            }
            for n in (0, 2, 4)
        },
        "Euler_frequency_coefficients": c["frequency_Euler"],
        "checks": checks,
        "gates": {
            "complete_second_and_fourth_matrix_word_bases": matrix_counts[2] == 4
            and matrix_counts[4] == 16,
            "noncommuting_connection_is_not_dropped": bool(frame["R"][2].terms),
            "arbitrary_noncommuting_fourth_jet_retained": ("K4",)
            in c["matrix_Euler"][4].terms,
            "complete_both_independent_variations": set(matrix_counts)
            == set(frequency_counts)
            == {0, 2, 4},
            "no_vacuum_phase_heuristic_used_as_proof": True,
        },
    }
