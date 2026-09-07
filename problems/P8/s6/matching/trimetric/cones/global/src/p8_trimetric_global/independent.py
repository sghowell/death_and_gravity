"""Independent coefficientwise Fraction and exact source/dynamics replay.

No SymPy or new primary modules are imported. Only pinned immutable Poly
and Dual arithmetic primitives are reused, without ancestor physics verdicts.
"""

from fractions import Fraction as Q

from p8_composite_modes.independent import Poly
from p8_trimetric.independent import Dual


def identities():
    r, c, p, h, hi, dhi, g = (Poly.variable(name) for name in ("y", "c", "p", "x", "u", "v", "m4"))
    balance_cleared = 6*p*r**3*(r*h-c*hi)+6*p*hi*r**3*(c-1)
    expected = 6*p*r**3*(r*h-hi)
    friedmann_prime_cleared = 6*g*hi*dhi*c-6*p*r**3*(r*h-c*hi)
    lapse_minus_space_cleared = -2*g*c*dhi-2*p*r**3*c+2*p*r**3
    kg, kf, cg, cf, pg, pf, hp, nh = (Poly.variable(name) for name in
        ("A", "B", "c", "u", "alpha", "beta", "z", "rho"))
    k = kg+kf
    kp = -2*h*(kg*(1-cg)+kf*(1-cf))
    ray = kg*(hp+(1-cg)*h**2)-pg*(1-cg)+kf*(hp+(1-cf)*h**2)-pf*(1-cf)
    unull = pg*(cg-1)+pf*(cf-1)-nh*Q(1, 2)
    left, right, ym, y0, yp = (Poly.variable(name) for name in ("A", "B", "x", "z", "u"))
    secant_cleared = right*(y0-ym)-left*(yp-y0)
    chord = right*ym+left*yp-(left+right)*y0
    hs = [Poly.variable(name) for name in ("alpha", "beta", "m4", "p")]
    js = [Poly.variable(name) for name in ("rho", "beta0", "beta1", "beta2")]
    pe, pv = Poly.variable("A"), Poly.variable("B")
    def e2(values):
        return sum((values[i]*values[j] for i in range(4) for j in range(i+1, 4)), Poly())
    difference = [a-b for a, b in zip(hs, js, strict=True)]
    determinant_numerator = ((pe+pv)*(pe*e2(hs)+pv*e2(js))-e2([pe*a+pv*b for a, b in zip(hs, js, strict=True)]))*Q(1, 2)
    fp_numerator = pe*pv*Q(1, 4)*(sum(difference, Poly())**2-sum((x*x for x in difference), Poly()))
    return {"source_Bianchi_factorization_cleared": balance_cleared-expected,
            "source_balance_from_full_Einstein_residuals": friedmann_prime_cleared+3*hi*lapse_minus_space_cleared+expected,
            "combined_monotonic_equation_residual": 2*ray-2*unull-(2*k*hp-h*kp+nh),
            "unequal_proper_duration_secant_chord": secant_cleared+chord,
            "full_traceful_FP_stationary_determinant": determinant_numerator-fp_numerator}


def phase_fixtures():
    """Off-shell full residual/normalization identities, not fake solutions.

    K=4 in these fixtures, so normalized derivative is independently obtained
    by literal quotient differentiation of H/sqrt(K) in the Dual ring.
    Arbitrary residuals need not vanish; the summed identity must hold.
    """
    specs = [
        [(Q(4), Q(1), Q(1, 2), Q(2))],
        [(Q(16), Q(2), Q(3, 2), Q(-3))],
        [(Q(1), Q(1), Q(1, 2), Q(1)), (Q(3), Q(1), Q(3, 2), Q(2))],
        [(Q(1), Q(1), Q(2), Q(-2)), (Q(3), Q(1), Q(1, 3), Q(1))],
        [(Q(1), Q(1), Q(3, 4), Q(1)), (Q(3), Q(1), Q(5, 4), Q(-1))],
    ]
    h, hp, null = Q(-2, 3), Q(5, 7), Q(11, 13)
    result = []
    for spec in specs:
        kdual = Dual(0)
        ray, weighted = Q(0), -null/2
        for g, r, c, p in spec:
            rprime = r*(1-c)*h
            rdual = Dual(r, rprime)
            kdual += g/(rdual*rdual)
            ray += g/r**2*(hp+(1-c)*h**2)-p/r*(1-c)
            weighted += p/r*(c-1)
        if kdual.value != 4:
            raise ValueError("Fixture normalization expected exact sqrt(K)=2")
        root = Dual(Q(2), kdual.tangent/4)
        normalized = Dual(h, hp)/root
        combined = 2*kdual.value*hp-h*kdual.tangent+null
        if 2*ray-2*weighted != combined or normalized.tangent+null/16 != combined/16:
            raise ValueError("Independent physical dynamics/normalization fixture failed")
        result.append({"links": len(spec),
                       "inputs": {"links": [{key: str(value) for key, value in zip(("G", "R", "c", "p"), entry, strict=True)}
                                            for entry in spec],
                                  "H_u": str(h), "H_u_prime": str(hp), "n_h": str(null)},
                       "K": str(kdual.value), "K_prime": str(kdual.tangent),
                       "combined_residual": str(combined), "Z_prime": str(normalized.tangent)})
    return result


def actual_rolling_fixtures():
    """Actual local free-scalar solution family, including negative links.

    G_i=epsilon=1,p_g=p_f=q,B=-6q,b_i=-q. Require A>5/6 and
    q(A-1)>0; the analytic positive-root flow is the literal continued
    source/Einstein family, not a ghost or vacuum-health assertion.
    """
    result = []
    for a, q in ((Q(2), Q(1)), (Q(3, 2), Q(1)), (Q(9, 10), Q(-1)), (Q(11, 12), Q(-1))):
        c, n = 6*a-5, a/(6*a-5)
        rho = 12*q*(a-1)/a
        hr2 = 2*q*(a**3-1)/3
        a_r_rate_over_hr = -6*a*(a-1)/c
        hrprime = -6*q*a**3*(a-1)/c
        hu2 = hr2/a**2
        huprime = (hrprime/a-hr2*a_r_rate_over_hr/a**2)/n
        k, kp_over_h = 2/a**2, -4*(1-c)/a**2
        residuals = (-6*q+6*q/a+rho/2, -6*q+2*q/n+4*q/a-rho/2,
                     3*hr2-2*q*a**3+2*q, 2*hrprime+3*hr2-2*q*n*a**2+2*q,
                     3*(1+a_r_rate_over_hr/a)+(12*q/a**2)*a_r_rate_over_hr/(2*rho),
                     2*k*huprime-hu2*kp_over_h+2*rho)
        if any(residuals) or not (a > 0 and n > 0 and rho > 0 and hr2 > 0 and k > 0):
            raise ValueError("Actual independent positive-NEC rolling solution fixture failed")
        result.append({"A": str(a), "q": str(q), "N": str(n), "rho": str(rho),
                       "H_u_squared": str(hu2), "H_u_prime": str(huprime),
                       "K": str(k), "K_prime_over_H_u": str(kp_over_h),
                       "actual_combined_residual": "0"})
    return result


def calibrated_weight_fixtures():
    result = []
    for rg, rf in ((Q(1), Q(1)), (Q(1, 2), Q(3, 2)), (Q(3, 4), Q(2))):
        pg, pf, rg0, rf0 = Q(-2), Q(1), Q(1), Q(1)
        for p, r, r0 in ((pg, rg, rg0), (pf, rf, rf0)):
            h2 = 2*p*(r**3-r0**3)/3
            loss = p/r-p/r0
            if h2 < 0 or loss != -3*h2/(2*r*r0*(r*r+r*r0+r0*r0)) or loss > 0:
                raise ValueError("Same-action calibrated weight loss failed")
        total = pg/rg+pf/rf
        if total > -1:
            raise ValueError("Mixed-sign positive-FP denominator bound failed")
        result.append({"R_g": str(rg), "R_f": str(rf), "S": str(total), "S0": "-1"})
    return result


def checks():
    polynomials = identities()
    if any(not value.is_zero() for value in polynomials.values()):
        raise ValueError("Independent coefficientwise identity failed")
    return {"coefficientwise_identities": dict.fromkeys(polynomials, "0"),
            "off_shell_source_dynamics_normalization_fixtures": phase_fixtures(),
            "actual_positive_NEC_rolling_fixtures": actual_rolling_fixtures(),
            "same_action_mixed_weight_fixtures": calibrated_weight_fixtures()}
