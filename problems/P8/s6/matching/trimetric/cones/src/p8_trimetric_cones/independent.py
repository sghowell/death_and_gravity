"""Independent Fraction/polynomial checks; no SymPy or new core imports.

Only immutable exact arithmetic primitives are reused from pinned ancestors.
Their cosmological, vacuum and principal-cone verdicts are not imported.
"""

from fractions import Fraction as Q

from p8_composite_modes.independent import Poly
from p8_trimetric.independent import Dual


def identities():
    pg, pf, ge, gf, gu, j = (Poly.variable(name) for name in ("A", "B", "x", "z", "u", "p"))
    den = pg+pf
    numerator = pg*ge+pf*gf+2*j
    lag = -(pg*(ge-gu)**2+pf*(gf-gu)**2)*Q(1, 4)+j*gu
    stationary_numerator = -(pg*(ge*den-numerator)**2+pf*(gf*den-numerator)**2)*Q(1, 4)+j*numerator*den
    reduced_numerator = -pg*pf*den*(ge-gf)**2*Q(1, 4)+j*(pg*ge+pf*gf)*den+j**2*den
    b, eps, rho = (Poly.variable(name) for name in ("C", "v", "rho"))
    lapse = b+3*den+eps*rho*Q(1, 2)
    space = b+pg*(ge+2)+pf*(gf+2)-eps*j*Q(1, 2)
    return {
        "weighted_cone_from_two_full_u_equations": space-lapse-pg*(ge-1)-pf*(gf-1)+eps*(rho+j)*Q(1, 2),
        "auxiliary_stationary_action_cleared": stationary_numerator-reduced_numerator,
        "full_sourced_auxiliary_completed_square": lag*den+pg*pf*(ge-gf)**2*Q(1, 4)-j*(pg*ge+pf*gf)-j**2
                                                  +(den*gu-numerator)**2*Q(1, 4),
        "exact_speed_square_excess": (1+ge)**2-1-2*ge-ge**2,
    }


def clock_fixtures():
    """Reconstruct from physical measure/derivatives, not primary formulas."""
    result = []
    for a, n, g in ((Q(2), Q(2, 7), Q(1)), (Q(3, 2), Q(5, 4), Q(2))):
        # dT=N dt_r, a_h=A a_r; the r EH action has coefficient 2G.
        # Convert its measure and BOTH derivative types independently.
        measure = 1/(n*a**3)
        time_derivative = n
        spatial_derivative = a
        kinetic = 2*g*measure*time_derivative**2
        gradient = 2*g*measure*spatial_derivative**2
        result.append({"A": a, "N": n, "G": g, "physical_measure_ratio": measure,
                       "G_T_h": kinetic, "F_T_h": gradient,
                       "common_speed_squared": gradient/kinetic})
    if tuple(result[0][key] for key in ("G_T_h", "F_T_h", "common_speed_squared")) != (Q(1, 14), Q(7, 2), Q(49)):
        raise ValueError("First independent physical-clock pullback fixture failed")
    if tuple(result[1][key] for key in ("G_T_h", "F_T_h", "common_speed_squared")) != (Q(40, 27), Q(32, 15), Q(36, 25)):
        raise ValueError("Second independent physical-clock pullback fixture failed")
    return result


def literal_lapse_scale_jets():
    ne, ae, nv, av, nu0, au0, pg, pf, b, eps, velocity, potential = map(Q, (2, 3, 5, 7, 11, 13, 17, 19, -23, 2, 29, 31))
    rho = velocity**2/(2*nu0**2)+potential
    pressure = velocity**2/(2*nu0**2)-potential
    pe, pv = pg*ae/au0, pf*av/au0
    ce, cf = au0*ne/(nu0*ae), au0*nv/(nu0*av)
    expected = (-2*au0**3*(b+3*(pe+pv)+eps*rho/2),
                -6*nu0*au0**2*(b+pe*(ce+2)+pv*(cf+2)-eps*pressure/2))
    for index in range(2):
        nu, au = Dual(nu0, int(index == 0)), Dual(au0, int(index == 1))
        volume = nu*au*au*au
        action = -2*volume*(b+pg*(ne/nu+3*ae/au)+pf*(nv/nu+3*av/au))
        action += eps*au*au*au*(velocity**2/(2*nu)-nu*potential)
        if action.tangent != expected[index]:
            raise ValueError("Literal independent lapse/scale source variation failed")
    return 2


def rolling_fixtures():
    result = {}
    for a in (Q(7, 6), Q(5, 4), Q(2), Q(101, 100)):
        n = a/(6*a-5)
        rho = 12*(a-1)/a
        h2 = 2*(a**3-1)/3
        a_rate_over_h = -6*a*(a-1)/(6*a-5)
        h_rate = -6*a**3*(a-1)/(6*a-5)
        rho_derivative = 12/a**2
        residuals = (-6+6/a+rho/2, -6+2/n+4/a-rho/2,
                     3*h2-2*a**3+2, 2*h_rate+3*h2-2*n*a**2+2,
                     2*a**2*a_rate_over_h-2*h_rate,
                     3*(1+a_rate_over_h/a)+rho_derivative*a_rate_over_h/(2*rho),
                     a/n-1-a*2*rho/4)
        if any(residuals) or not (n > 0 and rho > 0 and h2 > 0 and a/n > 1):
            raise ValueError("A full rolling-background independent identity failed")
        result[str(a)] = {"N": str(n), "rho": str(rho), "H_r_squared": str(h2), "c_T": str(a/n)}
    return result


def checks():
    polynomials = identities()
    if any(not poly.is_zero() for poly in polynomials.values()):
        raise ValueError("An independent polynomial identity failed")
    jets = literal_lapse_scale_jets()
    fixtures = rolling_fixtures()
    # Independent exact mixed-sign controls with physically distinct scope.
    pg, pf, ce, cf, nh = Q(-2), Q(1), Q(1, 2), Q(1), Q(2)
    if pg*(ce-1)+pf*(cf-1) != nh/2:
        raise ValueError("Mixed-sign algebraic countercontrol failed")
    spring = 2*pg*pf/(pg+pf)
    if spring != 4:
        raise ValueError("Mixed-sign flat relative spring control failed")
    return {"coefficientwise_identities": dict.fromkeys(polynomials, "0"),
            "literal_source_lapse_scale_jet_directions": jets,
            "independent_physical_clock_pullbacks": [{key: str(value) for key, value in fixture.items()}
                                                     for fixture in clock_fixtures()],
            "full_rolling_Fraction_fixtures": fixtures,
            "mixed_sign_auxiliary_ce_cf": [str(ce), str(cf)],
            "mixed_sign_flat_spring": str(spring)}
