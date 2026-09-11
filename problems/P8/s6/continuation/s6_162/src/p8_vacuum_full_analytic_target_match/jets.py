"""Exact sixth/eighth scalar-field germs of the actual rank-regular target."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import family
from p8_offshell_vacuum import jets as old_jets

PHI, Y, L3, L4, L5 = s.symbols("Phi Y L3 L4 L5", real=True)
N, KAPPA, LAM = s.symbols("n kappa lambda", positive=True)
VARIABLES = (PHI, Y, L3, L4, L5)
WEIGHTS = (1, 2, 4, 4, 6)


def literal(n, kappa, coupling):
    if type(n) is not int or n < 6 or n % 2:
        raise ValueError(
            "Require an even native order at least six for these local germs"
        )
    u, X = family.u, family.X
    T = X**n / (X**n + (1 - X) ** n)
    w = n * X**2 * s.exp(-n * X**2)
    B = T + (1 - T) * w
    R = 1 + B * (X - 1) / (1 + u * u) ** 3
    RX = s.diff(R, X)
    tree = family.data()["original_retuned_tree_scalar"]
    F0 = tree.subs({u: 0, X: 0})
    Fv = (
        (X - u * u) / 2
        + (coupling * kappa - n * F0) * X * X
        + (-s.Rational(n, 3) - F0) * u**4
    )
    F = tree + s.exp(-(u**4)) * (1 - T) * (1 - w) * (Fv - tree)
    return {
        "F": F,
        "R": R,
        "RX": RX,
        "A3": RX / X,
        "A4": -RX / X - s.Rational(7, 4) * RX**2 / R,
        "A5": RX**2 / (R * X),
    }


@cache
def data():
    u, X = family.u, family.X
    tree = family.data()["original_retuned_tree_scalar"]
    treejet = sum(
        s.diff(tree, u, a, X, b).subs({u: 0, X: 0})
        * u**a
        * X**b
        / (s.factorial(a) * s.factorial(b))
        for a in range(5)
        for b in range(3)
        if a + 2 * b <= 4
    )
    F0 = tree.subs({u: 0, X: 0})
    Fv = (X - u * u) / 2 + (LAM * KAPPA - N * F0) * X * X + (-N / 3 - F0) * u**4
    defect = u**4 + N * X * X - u**8 / 2 - N * u**4 * X * X - N * N * X**4
    lower = Fv + defect * (treejet - Fv)
    # T starts at X^n; for n>=6 it affects no displayed weighted degree.
    R = 1 + (-N * X * X + N * X**3 + N * N * X**4) * (1 - 3 * u * u + 6 * u**4)
    RX = s.diff(R, X)
    a3 = RX / X
    a4 = -a3 - s.Rational(7, 4) * RX * RX * (2 - R)
    a5 = RX * RX * (2 - R) / X
    t = s.Symbol("field_degree", real=True)
    sub = {u: t * PHI / s.sqrt(KAPPA), X: t * t * Y / KAPPA}
    scaled = s.expand(KAPPA * lower.subs(sub, simultaneous=True))
    scaled += s.expand(a3.subs(sub, simultaneous=True)) * t**4 * L3 / KAPPA
    scaled += s.expand(a4.subs(sub, simultaneous=True)) * t**4 * L4 / KAPPA
    scaled += s.expand(a5.subs(sub, simultaneous=True)) * t**6 * L5 / KAPPA**2
    pieces = {
        degree: s.factor(s.expand(scaled).coeff(t, degree)) for degree in (2, 4, 6, 8)
    }
    e6 = (
        -(2681 * PHI**2 + 853 * Y) * (PHI**4 + N * Y**2) / (200 * KAPPA**2)
        + 3 * N * (2 * PHI**2 + Y) * (L3 - L4) / KAPPA**2
    )
    e8 = (
        200 * PHI**8 * N
        - 29073 * PHI**8
        - 4878 * PHI**6 * Y
        - 600 * PHI**4 * Y**2 * KAPPA * LAM
        + 200 * PHI**4 * Y**2 * N**2
        - 29409 * PHI**4 * Y**2 * N
        - 1872 * PHI**4 * Y**2
        - 4878 * PHI**2 * Y**3 * N
        - 600 * Y**4 * KAPPA * LAM * N
        - 1872 * Y**4 * N
    ) / (600 * KAPPA**3)
    e8 += (
        N
        * (
            -12 * L3 * PHI**4
            - 9 * L3 * PHI**2 * Y
            + 4 * L3 * Y**2 * N
            + 12 * L4 * PHI**4
            + 9 * L4 * PHI**2 * Y
            - 11 * L4 * Y**2 * N
            + 4 * L5 * Y * N
        )
        / KAPPA**3
    )
    j = old_jets.data()
    # Substitute n=gamma*kappa explicitly rather than a structural ratio
    # match inside factored expressions.
    quartic = s.expand(
        pieces[4].subs(N, KAPPA * s.Symbol("quartic_gamma", real=True))
    ).subs(
        {
            PHI: j["phi"],
            Y: j["X"],
            L3: j["L3"],
            L4: j["L4"],
            LAM: s.Symbol("quartic_lambda", real=True),
        },
        simultaneous=True,
    )
    checks = {
        "literal_retuned_vacuum_constant": F0 + s.Rational(28, 25),
        "literal_retuned_weighted_tree_jet": s.expand(
            treejet
            + s.Rational(28, 25)
            + s.Rational(2781, 200) * u * u
            + s.Rational(753, 200) * X
            + s.Rational(9579, 200) * u**4
            + s.Rational(813, 100) * X * u * u
            + s.Rational(78, 25) * X * X
        ),
        "canonical_quadratic_mass_and_residue": s.expand(pieces[2] - (Y - PHI**2) / 2),
        "canonical_quartic_matches_literal_four_dimensional_target": s.expand(
            quartic - j["target_quartic"]
        ),
        "complete_sixth_density": s.expand(pieces[6] - e6),
        "complete_eighth_density": s.expand(pieces[8] - e8),
        "finite_kappa_sextic_not_assumed_zero": s.expand(
            s.diff(e6, L3) - 3 * N * (2 * PHI**2 + Y) / KAPPA**2
        ),
    }
    for degree, piece in pieces.items():
        checks[f"weighted_homogeneity_{degree}"] = s.expand(
            piece.subs(
                {v: t**w * v for v, w in zip(VARIABLES, WEIGHTS)}, simultaneous=True
            )
            - t**degree * piece
        )
    return {
        "retuned_tree_weighted_degree_four_jet": treejet,
        "canonical_density_field_degrees": pieces,
        "canonical_normalization": "u=Phi/sqrt(kappa), X=Y/kappa. The flat scalar density is kappa F+(A3 L3+A4 L4)/kappa+A5 L5/kappa^2, with L5=Z^2. Curvature terms vanish only for the stipulated fixed flat metric; their metric variations are not being set to zero.",
        "rational_step_jet_order": "T begins at weighted scalar-field degree 2n. Its removal is exact only for these finite germs, not for the full action or a loop integral.",
        "checks": checks,
    }
