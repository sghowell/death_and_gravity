"""Fraction polynomial/full-coframe replay; no SymPy or new primary imports."""

from fractions import Fraction as Q

from p8_composite_modes.independent import Poly
from p8_star.independent import derivative
from p8_trimetric.independent import Dual, determinant


def literal_fixture(beta, scales, lapses, lam=Q(1)):
    beta, scales, lapses = (tuple(Q(value) for value in values) for values in (beta, scales, lapses))
    matrices = [[[Q(n if row == 0 else a) if row == column else Q(0) for column in range(4)]
                 for row in range(4)] for a, n in zip(scales, lapses, strict=True)]
    gradients = []
    directions = 0
    for which in range(len(beta)):
        gradient = [[Q(0) for _ in range(4)] for _ in range(4)]
        for row in range(4):
            for col in range(4):
                summed = [[sum((b*Dual(matrix[i][j], int(index == which and (i, j) == (row, col)))
                                for index, (b, matrix) in enumerate(zip(beta, matrices, strict=True))), Dual(0))
                           for j in range(4)] for i in range(4)]
                gradient[row][col] = (-lam*determinant(summed)).tangent
                directions += 1
        gradients.append(gradient)
    rho, pressure = [], []
    for a, n, gradient in zip(scales, lapses, gradients, strict=True):
        if any(gradient[row][col] for row in range(4) for col in range(4) if row != col):
            raise ValueError("Aligned off-diagonal first variation failed")
        p = [gradient[i][i]/(n*a*a) for i in range(1, 4)]
        if len(set(p)) != 1:
            raise ValueError("Independent spatial directions disagree")
        rho.append(-gradient[0][0]/a**3)
        pressure.append(p[0])
    nulls = tuple(r+p for r, p in zip(rho, pressure, strict=True))
    weighted = sum(n*a**3*null for n, a, null in zip(lapses, scales, nulls, strict=True))
    if weighted:
        raise ValueError("Literal full-coframe weighted interaction null cancellation failed")
    sa, sn = sum(b*a for b, a in zip(beta, scales, strict=True)), sum(b*n for b, n in zip(beta, lapses, strict=True))
    return {"beta": beta, "a": scales, "n": lapses, "lambda": lam, "S_a": sa, "S_N": sn,
            "rho": tuple(rho), "pressure": tuple(pressure), "nulls": nulls,
            "directions": directions, "weighted_null": weighted}


def stress_fixtures():
    return [literal_fixture(beta, aa, nn) for beta, aa, nn in (
        ((1, 2, 3), (2, 3, 4), (1, 2, 3)),
        ((1, -1, 1), (1, 2, 3), (1, 2, 3)),
        ((1, 1, -1), (1, 2, 1), (1, 1, 2)),
        ((1, 0, 2), (2, 3, 4), (1, 2, 3)),
    )]


def polynomial_identities():
    aa = tuple(Poly.variable(name) for name in ("x", "z", "u"))
    nn = tuple(Poly.variable(name) for name in ("A", "B", "C"))
    vv = tuple(Poly.variable(name) for name in ("v", "accel", "root"))
    bb = tuple(Poly.variable(f"beta{i}") for i in range(3))
    lam = Poly.variable("p")
    sa = sum((b*a for b, a in zip(bb, aa, strict=True)), Poly())
    sn = sum((b*n for b, n in zip(bb, nn, strict=True)), Poly())
    sadot = sum((b*v for b, v in zip(bb, vv, strict=True)), Poly())
    lag = -lam*sn*sa**3
    weighted = Poly()
    out = {}
    for index, (a, n, v, b, aname, nname) in enumerate(zip(
            aa, nn, vv, bb, ("x", "z", "u"), ("A", "B", "C"), strict=True)):
        rho_numerator = -derivative(lag, nname)
        pressure_numerator = derivative(lag, aname)*Q(1, 3)
        weighted += n*rho_numerator+a*pressure_numerator
        rdot = sum((derivative(rho_numerator, name)*vel
                    for name, vel in zip(("x", "z", "u"), vv, strict=True)), Poly())
        # Literal continuity equation multiplied by n_i^2 a_i^4;
        # the density-volume derivative terms cancel before any sum division.
        cleared_bianchi = a*(n*rdot+3*v*pressure_numerator)
        out[f"density_numerator_{index}"] = rho_numerator-lam*b*sa**3
        out[f"pressure_numerator_{index}"] = pressure_numerator+lam*b*sn*sa**2
        out[f"unfactored_Bianchi_cleared_{index}"] = cleared_bianchi-3*lam*b*sa**2*a*(n*sadot-sn*v)
    out["literal_weighted_null_sum"] = weighted
    return out


def rate_fixtures():
    fixtures = []
    for gs, ys, cs, h, nulls in (
        ((1, 2, 3), (1, 2, 3), (1, Q(3, 2), Q(1, 3)), Q(-2), (1, 2, 3)),
        ((2, 1), (1, Q(3, 2)), (1, 7), Q(0), (3, 0)),
        ((3,), (1,), (1,), Q(4), (0,)),
    ):
        gs, ys, cs, nulls = (tuple(Q(v) for v in values) for values in (gs, ys, cs, nulls))
        yp = tuple(y*(c-1)*h for y, c in zip(ys, cs, strict=True))
        kjet = Dual(0)
        for g, y, rate in zip(gs, ys, yp, strict=True):
            yjet = Dual(y, rate)
            kjet += g*yjet*yjet
        k, kp = kjet.value, kjet.tangent
        # Pull the null sources using the actual measure n_i a_i^3 / n_r a_r^3.
        ni, ai = tuple(c*y for c, y in zip(cs, ys, strict=True)), ys
        source = sum(n*a**3*value for n, a, value in zip(ni, ai, nulls, strict=True))
        hp = (h*kp-source)/(2*k)
        if k <= 0 or source < 0 or k*hp-h*kp/2+source/2:
            raise ValueError("Independent proper-clock rate identity failed")
        fixtures.append({"Gs": gs, "ys": ys, "cs": cs, "H": h, "nulls": nulls,
                         "K": k, "Kprime": kp, "Hprime": hp, "yprimes": yp,
                         "weighted_null": source})
    return fixtures


def actual_stiff_fixtures():
    ys, gs, beta = (Q(1), Q(2), Q(3)), (Q(1), Q(2), Q(3)), (1, -1, 1)
    fixtures = []
    for base in (Q(1, 2), Q(1), Q(2)):
        t = base**3
        scales = tuple(y*base for y in ys)
        stress = literal_fixture(beta, scales, ys)
        cc = tuple(-rho/g for rho, g in zip(stress["rho"], gs, strict=True))
        hs = tuple(Q(1, 3)/(y*t) for y in ys)
        h_rates = tuple(-Q(1, 3)/(y*y*t*t) for y in ys)
        matter = tuple(g/(3*y*y*t*t) for g, y in zip(gs, ys, strict=True))
        equations = []
        for y, g, h, hp, cosmological, rho, ri, pi in zip(
                ys, gs, hs, h_rates, cc, matter, stress["rho"], stress["pressure"], strict=True):
            equations += [g*(3*h*h-cosmological)-rho-ri,
                          g*(-2*hp-3*h*h+cosmological)-rho-pi,
                          -1/(y*y*t*t)+3*h/(y*t)]
        if any(equations) or cc != (Q(-8), Q(1, 2), Q(-8, 81)):
            raise ValueError("Actual mixed-sign proportional scalar/EH solution failed")
        fixtures.append({"T": t, "S_a": stress["S_a"], "S_N": stress["S_N"],
                         "Hs": hs, "matter_rho": matter, "Lambda": cc})
    return fixtures


def auxiliary_full_matrix_jets():
    diagonal = (Q(2), Q(3), Q(4), Q(5))
    u = [[diagonal[i] if i == j else Q(0) for j in range(4)] for i in range(4)]
    # det(w) tr(w^-1 U) is sum over columns of det(w with that
    # column replaced by U), so no inverse is imported into this audit.
    directions = 0
    for row in range(4):
        for col in range(4):
            w = [[Dual(u[i][j], int((i, j) == (row, col))) for j in range(4)] for i in range(4)]
            trace_adjugate = Dual(0)
            for column in range(4):
                replaced = [[u[i][j] if j == column else w[i][j] for j in range(4)] for i in range(4)]
                trace_adjugate += determinant(replaced)
            lag = 3*determinant(w)-trace_adjugate
            if lag.tangent or lag.value != -determinant(u):
                raise ValueError("Independent source-preserving full auxiliary stationary map failed")
            directions += 1
    return directions


def actual_vacuum_fixtures():
    gs, ys = (Q(1), Q(2), Q(3)), (Q(1), Q(2), Q(3))
    stress = literal_fixture((1, -1, 1), ys, ys)
    fixtures = []
    for h0 in (Q(-2), Q(0), Q(3, 2)):
        hs = tuple(h0/y for y in ys)
        cc = tuple(3*h*h-r/g for h, r, g in zip(hs, stress["rho"], gs, strict=True))
        residuals = []
        for g, h, cosmological, ri, pi in zip(gs, hs, cc, stress["rho"], stress["pressure"], strict=True):
            residuals += [g*(3*h*h-cosmological)-ri, g*(-3*h*h+cosmological)-pi]
        if any(residuals):
            raise ValueError("Independent actual proportional vacuum Einstein fixture failed")
        fixtures.append({"H0": h0, "Hs": hs, "Lambda": cc, "scaled_H_prime": Q(0)})
    return fixtures


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    return str(value)


def checks():
    polynomials, stresses = polynomial_identities(), stress_fixtures()
    if any(not value.is_zero() for value in polynomials.values()):
        raise ValueError("Independent determinant polynomial replay failed")
    return {"coefficientwise_identities": dict.fromkeys(polynomials, "0"),
            "full_coframe_first_jet_directions": sum(row["directions"] for row in stresses),
            "literal_stress_fixtures": serialize(stresses),
            "physical_clock_fixtures": serialize(rate_fixtures()),
            "actual_scalar_Einstein_fixtures": serialize(actual_stiff_fixtures()),
            "actual_vacuum_Einstein_fixtures": serialize(actual_vacuum_fixtures()),
            "full_auxiliary_stationary_directions": auxiliary_full_matrix_jets()}
