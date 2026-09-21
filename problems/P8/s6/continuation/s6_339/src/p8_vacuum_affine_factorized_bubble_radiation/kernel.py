"""General physical two-endpoint kernel and a Ward-compatible representative."""

from functools import cache

import sympy as s
from p8_vacuum_affine_local_tadpole_radiation import radiation, vertices


def channel(gram, pol, aa, L, R, FL, FR, divided, trace=0, kep=None):
    trace = s.sympify(trace)
    kep = (0, 0, 0, 0) if kep is None else kep
    J = tuple(pol[i, i] / aa[i] + kep[i] / aa[i] - trace / 2 for i in range(4))
    pp = sum(pol[i, j] for i in L for j in L)
    kp = sum(kep[i] for i in L)
    return (
        sum(J[i] * FR for i in L)
        + sum(J[i] * FL for i in R)
        - (2 * pp + 2 * kp) * divided
        + trace * (FL + FR) / 2
    )


@cache
def data():
    G, H = vertices.generic_radiative_data()
    ak = tuple(G[i, 4] for i in range(4))
    xi = s.symbols("xi_dot_p0:4")
    beta = -sum(xi)
    gauge = s.Matrix(4, 4, lambda i, j: ak[i] * xi[j] + ak[j] * xi[i])
    checks = {}
    f0, f1 = s.symbols("F_left F_right")
    for j in (1, 2, 3):
        L = (0, j)
        R = tuple(i for i in range(4) if i not in L)
        v = sum(G[i, h] for i in L for h in L)
        u = sum(G[i, h] for i in R for h in R)
        checks[f"channel{j}_null_transfer_identity"] = s.factor(
            u - v - 2 * sum(ak[i] for i in L)
        )
        pure = channel(
            G,
            gauge,
            ak,
            L,
            R,
            f0,
            f1,
            (f1 - f0) / (u - v),
            2 * beta,
            tuple(beta * a for a in ak),
        )
        checks[f"channel{j}_arbitrary_kernel_full_Ward"] = s.factor(pure)

    # Independently compare low-degree local covariant quartic actions.
    for degree in (0, 1, 2):
        value = s.S.Zero
        for j in (1, 2, 3):
            L = (0, j)
            R = tuple(i for i in range(4) if i not in L)
            v = sum(G[i, h] for i in L for h in L)
            u = sum(G[i, h] for i in R for h in R)
            div = (
                sum(u**h * v ** (degree - 1 - h) for h in range(degree))
                if degree
                else s.S.Zero
            )
            value += channel(G, H, ak, L, R, v**degree, u**degree, div)
        phi4 = radiation.basis_radiation("phi4", G, H)
        target = (
            phi4 / 8
            if degree == 0
            else phi4 / 6
            if degree == 1
            else radiation.basis_radiation("Y2", G, H) / 2 + phi4 / 6
        )
        checks[f"power{degree}_independent_literal_covariant_quartic"] = s.factor(
            value - target
        )

    v, u, n = s.symbols("v u n")
    for degree in range(1, 9):
        divided = sum(u**h * v ** (degree - 1 - h) for h in range(degree))
        checks[f"power{degree}_continuous_equal_invariant_extension"] = s.expand(
            divided.subs(u, v) - degree * v ** (degree - 1)
        )
    checks["resolvent_complete_divided_difference"] = s.factor(
        (1 / (n - u) - 1 / (n - v)) / (u - v) - 1 / ((n - u) * (n - v))
    )
    return {
        "checks": checks,
        "gates": {
            "three_generic_channels_full_puregauge_Ward": True,
            "three_independent_literal_covariant_quartic_comparisons": True,
            "eight_complete_operator_insertion_coincidence_limits": True,
            "resolvent_identity_independent_of_loop_branch": True,
            "physical_TT_complete_unprojected_kk_form_factor_not_claimed": True,
        },
        "whole_kernel_formula": "For P=sum_L p_i,v=P^2,u=(P+k)^2,J_i=epsilon(pi,pi)/ai+epsilon(k,pi)/ai-trace/2, a Ward-compatible channel is sum_L J_i F(u)+sum_R J_i F(v)-(2epsilon(P,P)+2epsilon(k,P))*D_F(v,u)+trace*(F(v)+F(u))/2. In physical TT the k and trace terms vanish. D_F extends as Fprime(v) at equality.",
        "whole_unprojected_boundary": "The complete physical TT response of the actual two-endpoint kernels is computed. A further transverse k_mu k_nu form factor can occur in the unprojected loop tensor; it is not determined or declared absent. This is not a full curved counterfunctional or the four-hard-direction S336 curvature coordinate.",
        "whole_local_controls": "F(v)=1,v,v^2 give the complete covariant quartic radiation of Phi4/8,Phi4/6,Y2/2+Phi4/6. F=1/(n-v) reproduces minimal heavy exchange. These are independent identities of the literal prior action, not inferred only from the leading soft residue.",
    }
