"""Independent Fraction replay of the coefficient and source norm arguments.

No SymPy and no primary-module imports. Finite polynomial inverse checks
are controls on the all-degree proof; they are not its replacement.
"""

from fractions import Fraction as Q
from math import comb, factorial, prod


def fraction(value):
    if isinstance(value, bool) or not isinstance(value, (int, Q)):
        raise TypeError("Independent arithmetic admits integers and Fractions only")
    return Q(value)


def inverse_monomial(j, n):
    if (isinstance(j, bool) or isinstance(n, bool) or not isinstance(j, int)
            or not isinstance(n, int) or min(j, n) < 0):
        raise ValueError("Require nonnegative integer monomial powers")
    # Closed factorial/product formula, not the production recursive inverse.
    return {(j+shift, n-2*shift): Q((-1)**shift*factorial(n),
            factorial(n-2*shift)*prod(8*(m*m+3*m+12)
                                    for m in range(n, n-2*shift-1, -2)))
            for shift in range(n//2+1)}


def apply_leading(polynomial):
    output = {}
    for (j, n), coefficient in polynomial.items():
        output[j, n] = output.get((j, n), Q(0))+8*(n*n+3*n+12)*coefficient
        if n >= 2:
            output[j+1, n-2] = output.get((j+1, n-2), Q(0))+n*(n-1)*coefficient
    return {key: value for key, value in output.items() if value}


def analytic_constants():
    r, v = Q(1, 20), Q(1, 400)
    cp, cm, dp, dm = 2+v, 2-v, 1+v, 1-v
    power_error = cp*dp**12-2
    denominator = 10-power_error
    power_derivative = 12*cp*dp**11
    t = 6*(6+power_error)/(dm*denominator)
    tv = 6*(power_derivative/(dm*denominator)
            +(6+power_error)/(dm**2*denominator)
            +(6+power_error)*power_derivative/(dm*denominator**2))
    o = 48*dp**5/((1-v/2)*denominator)
    ov = 48/(1-v/2)*(5*dp**4/denominator+dp**5*power_derivative/denominator**2)
    nt, op = t+2*v*tv+v*t*t, o+2*v*ov
    weights = (cp*dp**12/denominator, 8/denominator)
    cf, sqw, qs = cp**2*dp**8/4, 4*dp**6/((1-v/2)*denominator), 4/dm**4
    cross = qs*sqw*(cf-1)
    dn = cp*dp**4-2
    udp = 8*v*cp*dp**3
    dpp = 8*cp*dp**3+48*v*cp*dp**2-16
    ne = 16*((1+v)*(8/(cm*dm**14)+1/dm**2)-5)
    tq = (dn-9*v)/(7*v)+(16*cp*dp**3-32)/70+(dpp+ne+9*dn)/84
    aa, bb, cc = 9*v/2, Q(23, 3)*dn+5*udp, (1+5*v)/7
    col_l, col_q = aa+v*cc, bb+tq
    le, lo = 2*(aa+5*v*v/7), 2*col_l*r
    return {"R": r, "v": v, "denominator_lower": denominator,
            "theta_over_u": t, "theta_v": tv, "omega_over_u": o, "omega_v": ov,
            "normalization": nt, "omega_prime": op,
            "A": qs*(weights[0]+cf*weights[1])+nt,
            "C": cross+2*op+2*v*o*t, "E": cross+2*v*o*t,
            "b_analytic": qs*(weights[1]+cf*weights[0])+nt+4*v*o*o,
            "D": dn, "u_Dprime": udp, "Dpp_minus16": dpp, "N_minus80": ne,
            "q_self": tq, "light_self": aa, "light_from_q_scaled": bb,
            "q_from_light": cc, "column_light": col_l, "column_q": col_q,
            "even_l_error": le, "even_q": (5*v/7+cc*le)/(1-tq),
            "odd_l_error": lo, "odd_q": cc*(r+lo)/(1-tq),
            "inverse_norm": Q(1, 84), "inverse_u_derivative_norm": Q(1, 70),
            "inverse_u2_second_derivative_norm": Q(1, 7),
            "inverse_spring": dn*(cp*dp**12+8)/(dm**6*(80-ne)),
            "K_g": dp**6/8, "K_f": 1/(cm*dm**6), "G_g": dp**2/8,
            "physical_ag": dp**3/((1-v/2)*(1-power_error/10)),
            "physical_bg": 2/(dm**3*(1-power_error/10)),
            "physical_bf": cp*dp**9/(3*(1-power_error/10))}


def physical_constants():
    delta, v, r = Q(1, 10**9), Q(1, 400), Q(1, 100)
    a0, cp, dp = Q(7, 2640), 2+delta, 1+r*r
    even = (Q(2, 5), Q(80), Q(9, 2500), Q(74403, 100000))
    odd = (Q(1, 250), Q(6, 5), Q(3, 625), Q(3201, 5000))
    def column(values):
        l, lp, rel, relp = values
        return 2*(l+lp+10*rel+(r*relp+rel/2)/Q(3, 10))
    return {"delta_max": delta, "v": v, "omega_error": 3*delta**3/5,
            "ag_squared_lower": 8/(cp*dp**12+8), "ag_squared_upper": 4*cp*dp**6/10,
            "w1_upper": cp*dp**12/10, "D_pointwise": cp*dp**4-2,
            "Dprime_pointwise": 8*cp*r*dp**3,
            "D_delta_coefficient": dp**4, "Dprime_delta_coefficient": 8*r*dp**3,
            "g_even_lower": Q(79, 100), "g_oddprime_lower": Q(79, 100),
            "g_evenprime_upper": Q(27, 100), "g_odd_upper": Q(11, 1000),
            "W_lower": Q(79, 100)**2-Q(27, 100)*Q(11, 1000),
            "W_upper": Q(501, 500)*Q(101, 100)+Q(27, 100)*Q(11, 1000),
            "kinetic_lower": Q(9, 20), "kinetic_upper": Q(11, 6),
            "q_even_K_delta_coefficient": -a0, "q_even_K_remainder_coefficient": Q(300),
            "B_loss_per_delta": 104+(20+1280*a0)*delta+(640*a0+384000)*delta**2+192000*delta**3,
            "center_linear_coefficient": Q(14, 33), "center_remainder_coefficient": 48000+105*a0,
            "endpoint_even_column": column(even), "endpoint_odd_column": column(odd),
            "endpoint_matrix_coefficient": Q(200), "old_full_transfer_norm": Q(42),
            "fixed_light_coefficient": Q(8600), "fixed_source_L1_coefficient": Q(126000000)}


def source_constants():
    # Independent derivative-jet propagation for S1..S4. A jet has ordinary
    # derivatives, so product convolution includes its binomial coefficient.
    def jet(norm):
        return tuple(Q(norm)*factorial(j)*40**j for j in range(4))
    def multiply(left, right):
        return tuple(sum(Q(comb(n, j))*left[j]*right[n-j] for j in range(n+1))
                     for n in range(min(len(left), len(right))))
    def add(left, right):
        return tuple(a+b for a, b in zip(left, right))
    one, field, inv = jet(1), jet(5), jet(Q(1, 100))
    n = multiply(one, field)
    first = multiply(one, field[1:])
    second = multiply(one[1:], field)
    m = add(tuple(2*item for item in first), second)
    aa, bb = multiply(inv, m), multiply(inv, n)
    coeff = {
        1: 5*(2*one[0]*field[1]+one[1]*field[0]+one[1]*aa[1]
              +one[0]*aa[2]+4*one[0]*aa[0]+m[0]),
        2: 5*(one[0]*field[0]+one[0]*aa[1]+one[1]*(aa[0]+bb[1])
              +one[0]*(aa[1]+bb[2])+4*one[0]*bb[0]+n[0]),
        3: 5*(one[0]*(aa[0]+bb[1])+one[1]*bb[0]+one[0]*bb[1]),
        4: 5*one[0]*bb[0],
    }
    polynomial = {0: Q(1)}
    bumps = []
    for _ in range(5):
        bumps.append(sum(abs(value)*factorial(n) for n, value in polynomial.items()))
        new = {}
        for n, value in polynomial.items():
            new[n+2] = new.get(n+2, Q(0))+value
            if n:
                new[n+1] = new.get(n+1, Q(0))-n*value
        polynomial = {key: value for key, value in new.items() if value}
    inverse = [Q(9)]
    for n in range(1, 5):
        inverse.append(9*sum(Q(comb(n, j))*2*bumps[j]*inverse[n-j] for j in range(1, n+1)))
    cutoff = {0: Q(1)}
    for n in range(1, 5):
        cutoff[n] = 200**n*sum(Q(comb(n, j))*bumps[j]*inverse[n-j] for j in range(n+1))
    return {"coefficients": coeff, "bump_jets": bumps, "reciprocal_jets": inverse,
            "zeta_jets": cutoff, "source_sup": sum(coeff[n]*cutoff[n] for n in range(1, 5))}


def checks():
    cases = 0
    for j in range(4):
        for n in range(13):
            image = apply_leading(inverse_monomial(j, n))
            if image != {(j, n): Q(1)}:
                raise ValueError("Independent leading inverse failed")
            cases += 1
    a, p, s = analytic_constants(), physical_constants(), source_constants()
    inequalities = {"light_column_below_half": a["column_light"] < Q(1, 2),
                    "q_column_below_half": a["column_q"] < Q(1, 2),
                    "even_l_bound": a["even_l_error"] < 10*a["v"],
                    "even_q_bound": a["even_q"] < 3*a["v"],
                    "odd_l_bound": a["odd_l_error"] < a["R"]/40,
                    "odd_q_bound": a["odd_q"] < a["R"]/5,
                    "W_lower": p["W_lower"] > Q(3, 5), "W_upper": p["W_upper"] < 2,
                    "center_remainder": p["center_remainder_coefficient"] < 48001,
                    "endpoint_sum": p["endpoint_even_column"]+p["endpoint_odd_column"] < 200,
                    "fixed_source_identity": p["fixed_source_L1_coefficient"] == 42*3000000,
                    "source_finite_calibration": s["source_sup"] == 321137485366608000}
    if not all(inequalities.values()):
        raise ValueError("An independent rational analytic/source bound failed")
    return {"exact_inverse_polynomials": cases, "strict_and_exact_controls": inequalities}
