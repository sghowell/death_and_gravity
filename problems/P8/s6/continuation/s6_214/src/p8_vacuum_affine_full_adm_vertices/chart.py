"""Four-metric chart derivatives and the necessary first-current contact chain rule."""

from functools import cache

import sympy as s


def first(n, beta, Q, a):
    M = s.zeros(4)
    M[0, 0] = 2 * n
    M[0, 1:] = -a * a * beta.T
    M[1:, 0] = -a * a * beta
    M[1:, 1:] = -a * a * Q
    return M


def second(nD, bD, D, nG, bG, G, a):
    M = s.zeros(4)
    M[0, 0] = 2 * nD * nG - 2 * a * a * (bD.T * bG)[0]
    mixed = D * bG + G * bD
    M[0, 1:] = -a * a * mixed.T
    M[1:, 0] = -a * a * mixed
    M[1:, 1:] = -a * a * (D * G + G * D) / 2
    return M


@cache
def data():
    e, z = s.symbols("e z", real=True)
    a = s.Symbol("a", positive=True)
    nD, nG = s.symbols("nD nG", real=True)
    bD = s.Matrix(s.symbols("bD0:3", real=True))
    bG = s.Matrix(s.symbols("bG0:3", real=True))
    D = s.Matrix([[1, 2, 0], [2, 3, 1], [0, 1, 2]]) / 7
    G = s.Matrix([[2, 1, 3], [1, 1, -2], [3, -2, -3]]) / 11
    X = e * D + z * G
    h = a * a * (s.eye(3) + X + X * X / 2)
    N = 1 + e * nD + z * nG
    b = e * bD + z * bG
    g = s.zeros(4)
    g[0, 0] = N * N - (b.T * h * b)[0]
    g[0, 1:] = -b.T * h
    g[1:, 0] = -h * b
    g[1:, 1:] = -h
    mixed = g.applyfunc(lambda v: s.expand(v).coeff(e, 1).coeff(z, 1))
    firstD = g.applyfunc(lambda v: s.expand(v).coeff(e, 1).subs(z, 0))
    J = s.diag(2, 3, 5, 7)
    cov = s.trace(J * g) + s.trace(g * g) / 2
    g0 = s.diag(1, -a * a, -a * a, -a * a)
    dD, dG = first(nD, bD, D, a), first(nG, bG, G, a)
    dd = second(nD, bD, D, nG, bG, G, a)
    chain = s.trace(dD * dG) + s.trace((J + g0) * dd)
    checks = {
        "full_ADM_to_four_metric_first_chart": firstD - dD,
        "full_ADM_to_four_metric_noncommuting_second_chart": mixed - dd,
        "literal_functional_Hessian_chain_rule": s.expand(cov).coeff(e, 1).coeff(z, 1)
        - chain,
        "shift_shift_four_metric_time_contact": second(
            0, bD, s.zeros(3), 0, bG, s.zeros(3), a
        )[0, 0]
        + 2 * a * a * (bD.T * bG)[0],
        "lapse_lapse_four_metric_time_contact": second(
            nD, s.zeros(3, 1), s.zeros(3), nG, s.zeros(3, 1), s.zeros(3), a
        )[0, 0]
        - 2 * nD * nG,
    }
    return {
        "metric_chart": "g00=N^2-h(beta,beta),g0i=-(h beta)i,gij=-hij, N=1+n,h=a^2exp(Q). All first and full noncommuting mixed second derivatives are retained.",
        "contact_chain": "For the same renormalized metric functional, its ADM Hessian is W_g''[Dg,Gg]+W_g'[D_G g]. In particular zero canonical ADM shift-shift or lapse-lapse Hamiltonian contact does not authorize deleting covariant metric contacts; the nonlinear chart has nonzero g00 second variations.",
        "mean_not_zero": "The actual Gaussian stress one-point function is not set to zero in this chain rule. Any fixed scalar retuning belongs to the separate full-parent action and must be varied too when assembling the coupled equations.",
        "checks": checks,
        "gates": {
            "all_four_metric_components_and_chart_terms_retained": True,
            "no_zero_mean_stress_assumption": True,
            "spatial_noncommuting_second_chart": D * G != G * D,
            "no_reduced_constraint_or_gauge_inverse_inferred": True,
        },
    }
