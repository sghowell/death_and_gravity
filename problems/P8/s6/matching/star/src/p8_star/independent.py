"""Separate Fraction polynomial/coframe-jet replay; no SymPy or core imports.

Only exact arithmetic primitives are reused from immutable S6.13 ancestry.
The no-crossing topology/Gronwall argument is a written proof, not sampling.
"""

from fractions import Fraction as Q
from itertools import combinations

from p8_composite_modes.independent import NAMES, Poly
from p8_trimetric.independent import Dual


def product(values):
    result = 1
    for value in values:
        result *= value
    return result


def derivative(poly, name):
    index = NAMES.index(name)
    terms = {}
    for powers, coefficient in poly.terms.items():
        if powers[index]:
            next_powers = list(powers)
            next_powers[index] -= 1
            terms[tuple(next_powers)] = coefficient*powers[index]
    return Poly(terms)


def coefficient(poly, name, degree):
    index = NAMES.index(name)
    terms = {}
    for powers, value in poly.terms.items():
        if powers[index] == degree:
            reduced = list(powers)
            reduced[index] = 0
            terms[tuple(reduced)] = value
    return Poly(terms)


def polynomial_identities():
    r, n, hi, hu = (Poly.variable(name) for name in ("y", "c", "x", "z"))
    betas = [Poly.variable(f"beta{i}") for i in range(5)]
    roots = (n, r, r, r)
    full = sum((beta*sum((product(group) for group in combinations(roots, degree)), Poly())
                for degree, beta in enumerate(betas)), Poly())
    a, b = coefficient(full, "c", 0), coefficient(full, "c", 1)
    j = betas[1]+2*betas[2]*r+betas[3]*r**2
    radial = derivative(a, "y")+n*derivative(b, "y")
    leaf_null = 2*a-2*(a+n*b-r*radial*Q(1, 3))
    center_null_cleared = 2*n*b-2*r*radial*Q(1, 3)
    divergence = 2*derivative(a, "y")*r*(n*hu-hi)+6*hi*(r-n)*j
    k1, k2, g0, h, hp, s1, s2, null = (Poly.variable(name) for name in
                                        ("A", "B", "C", "x", "z", "u", "v", "p"))
    cone1, cone2 = n, r
    k = g0+k1+k2
    kp = -2*h*(k1*(1-cone1)+k2*(1-cone2))
    leaf1 = k1*(hp+(1-cone1)*h**2)-s1*(1-cone1)
    leaf2 = k2*(hp+(1-cone2)*h**2)-s2*(1-cone2)
    central = g0*hp+null*Q(1, 2)+s1*(1-cone1)+s2*(1-cone2)
    return {
        "four_eigenvalue_expansion": full-a-n*b,
        "A_from_elementary_subsets": a-betas[0]-3*betas[1]*r-3*betas[2]*r**2-betas[3]*r**3,
        "B_from_elementary_subsets": b-betas[1]-3*betas[2]*r-3*betas[3]*r**2-betas[4]*r**3,
        "A_derivative": derivative(a, "y")-3*j,
        "B_reciprocity_derivative": 3*b-r*derivative(b, "y")-3*j,
        "unfactored_leaf_Bianchi": divergence-6*n*j*(r*hu-hi),
        "leaf_null_from_scale_variation": leaf_null-2*(r-n)*j,
        "central_null_from_scale_variation_cleared": center_null_cleared-2*(n-r)*j,
        "weighted_null_cancellation_cleared": center_null_cleared+leaf_null,
        "dynamic_star_residual_combination": leaf1+leaf2+central-k*hp+h*kp*Q(1, 2)-null*Q(1, 2),
    }


def literal_jet_fixture(beta, leaf, center):
    """Differentiate four independent positive coframe entries before isotropy."""
    beta = tuple(Q(value) for value in beta)
    ni, ai = map(Q, leaf)
    nu, au = map(Q, center)
    entries = (ni, ai, ai, ai, nu, au, au, au)
    jets = []
    for direction in range(8):
        fields = [Dual(value, int(index == direction)) for index, value in enumerate(entries)]
        roots = [fields[index+4]/fields[index] for index in range(4)]
        polynomial = sum((beta[degree]*sum((product(group) for group in combinations(roots, degree)), Dual(0))
                          for degree in range(5)), Dual(0))
        lag = -2*product(fields[:4])*polynomial
        jets.append(lag.tangent)
    leaf_pressures = [value/(ni*ai**2) for value in jets[1:4]]
    center_pressures = [value/(nu*au**2) for value in jets[5:8]]
    if len(set(leaf_pressures)) != 1 or len(set(center_pressures)) != 1:
        raise ValueError("Literal three-direction spatial pressures disagree")
    rho_i, pressure_i = -jets[0]/ai**3, leaf_pressures[0]
    rho_u, pressure_u = -jets[4]/au**3, center_pressures[0]
    r, n = au/ai, nu/ni
    if rho_i+pressure_i+n*r**3*(rho_u+pressure_u):
        raise ValueError("Independent literal two-metric null cancellation failed")
    return {"beta": beta, "R": r, "N": n, "jets": tuple(jets),
            "rho_i": rho_i, "pressure_i": pressure_i,
            "rho_u": rho_u, "pressure_u": pressure_u,
            "null_i": rho_i+pressure_i, "null_u": rho_u+pressure_u}


def stress_fixtures():
    return [literal_jet_fixture(beta, leaf, center) for beta, leaf, center in (
        ((1, -2, 3, -4, 5), (2, 3), (5, 7)),
        ((0, 1, Q(-1, 2), 0, Q(1, 2)), (Q(1, 3), 1), (1, 1)),
        ((2, 0, 0, -2, 3), (Q(3, 2), Q(5, 3)), (Q(2, 3), Q(7, 4))),
    )]


def dynamic_fixtures():
    fixtures = []
    for gu, gs, rs, cs, h, null in (
        (Q(0), (Q(2),), (Q(3),), (Q(7, 4),), Q(-2), Q(5)),
        (Q(3), (Q(2), Q(5)), (Q(3), Q(2)), (Q(7, 4), Q(1, 3)), Q(0), Q(7)),
        (Q(1), (Q(2), Q(3), Q(4)), (Q(1), Q(2), Q(3)),
         (Q(1, 2), Q(4, 3), Q(2)), Q(3, 5), Q(1, 7)),
        (Q(2), (), (), (), Q(-1), Q(3)),
    ):
        rp = tuple(r*(1-c)*h for r, c in zip(rs, cs, strict=True))
        kjet = Dual(gu)
        for g, r, rate in zip(gs, rs, rp, strict=True):
            rjet = Dual(r, rate)
            kjet += g/(rjet*rjet)
        k, kp = kjet.value, kjet.tangent
        # Solve the summed Raychaudhuri equation using the independent K jet.
        hp = (h*kp-null)/(2*k)
        if k <= 0 or k*hp-h*kp/2 != -null/2:
            raise ValueError("Independent finite-star dynamic identity failed")
        rates = tuple(rate/r for rate, r in zip(rp, rs, strict=True))
        cap = max((abs(rate) for rate in rates), default=Q(0))
        if abs(kp/(2*k)) > cap:
            raise ValueError("Positive-weight compact comparison failed")
        fixtures.append({"G_u": gu, "Gs": gs, "Rs": rs, "cs": cs, "H": h, "nh": null,
                         "K": k, "Kprime": kp, "Hprime": hp, "Rprimes": rp,
                         "log_rate_cap": cap, "comparison_coefficient": kp/(2*k)})
    return fixtures


def actual_branch_fixtures():
    fixtures = []
    beta = (0, 1, Q(-1, 2), 0, Q(1, 2))
    for scale in (Q(1, 2), Q(1), Q(2), Q(3)):
        t = scale**3
        ni, nu = 1/(3*t), Q(1)
        fields = literal_jet_fixture(beta, (ni, scale), (nu, scale))
        # a=T^(1/3), its first two jets, and the proper-time conversion.
        adot, addot = scale/(3*t), -2*scale/(9*t*t)
        nip = -1/(3*t*t)
        hu, hi = adot/scale, adot/(ni*scale)
        hup = addot/scale-hu**2
        hip = addot/(ni*scale)-adot*nip/(ni**2*scale)-adot**2/(ni*scale**2)
        rho, pressure = Q(1, 3)/t**2, Q(1, 3)/t**2
        residuals = (3*hi**2-fields["rho_i"],
                     -2*hip/ni-3*hi**2-fields["pressure_i"],
                     3*hu**2-rho-fields["rho_u"],
                     -2*hup-3*hu**2-pressure-fields["pressure_u"],
                     -1/t**2+3*hu/t, hup+(rho+pressure)/2)
        wrong_all_k = 4*hup+rho+pressure
        if any(residuals) or wrong_all_k != -Q(2, 3)/t**2:
            raise ValueError("Actual algebraic Einstein/scalar or wrong-K control failed")
        fixtures.append({"T": t, "a_u": scale, "n_i": ni, "H_i": hi, "H_u": hu,
                         "rho_m": rho, "null_m": rho+pressure,
                         "wrong_all_K_undivided_defect": wrong_all_k})
    return fixtures


def disconnected_fixtures():
    fixtures = []
    for t in (Q(-1), Q(0), Q(2)):
        scale = 1+t*t
        fields = literal_jet_fixture((0, 0, 0, 0, 0), (1, 1), (1, scale))
        hu, hup = 2*t/scale, 2*(1-t*t)/scale**2
        # All eight independent potential variations are zero; G_u=0,
        # actual Minkowski leaf H_i=0, and zero scalar stress/current.
        residuals = fields["jets"]+(3*Q(1)*0**2-fields["rho_i"],
                                     -2*Q(1)*0-3*Q(1)*0**2-fields["pressure_i"],
                                     3*Q(0)*hu**2-fields["rho_u"],
                                     -2*Q(0)*hup-3*Q(0)*hu**2-fields["pressure_u"])
        if any(residuals):
            raise ValueError("Disconnected literal Euler control failed")
        fixtures.append({"T": t, "a_u": scale, "H_u": hu, "H_u_prime": hup,
                         "added_Gu1_constraint_defect": 3*hu**2})
    if fixtures[1]["H_u_prime"] != 2 or fixtures[2]["added_Gu1_constraint_defect"] == 0:
        raise ValueError("The necessary disconnected-central-kinetic exception was lost")
    return fixtures


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    return str(value)


def checks():
    polynomials = polynomial_identities()
    if any(not value.is_zero() for value in polynomials.values()):
        raise ValueError("Independent coefficientwise star identity failed")
    stresses = stress_fixtures()
    return {"coefficientwise_identities": dict.fromkeys(polynomials, "0"),
            "literal_coframe_jet_directions": 8*len(stresses),
            "literal_stress_fixtures": serialize(stresses),
            "finite_dynamic_star_fixtures": serialize(dynamic_fixtures()),
            "actual_algebraic_branch_fixtures": serialize(actual_branch_fixtures()),
            "literal_disconnected_euler_fixtures": serialize(disconnected_fixtures())}
