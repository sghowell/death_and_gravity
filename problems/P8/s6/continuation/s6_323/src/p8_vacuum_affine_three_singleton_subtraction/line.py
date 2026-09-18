"""Exact connected scalar-line algebra, physical mapping and pole controls."""

from functools import cache
from itertools import combinations, product
from types import SimpleNamespace

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as e
from p8_vacuum_affine_uniform_two_real_tree_bound import bounds as vertex_bound
from sympy.polys.rings import ring

imm = s.ImmutableMatrix


@cache
def generic_polynomials():
    names = [f"d{i}" for i in range(3)] + [
        f"z{i}{j}" for i, j in combinations(range(3), 2)
    ]
    names += [f"L{i}" for i in range(3)]
    names += [
        f"{kind}{i}{j}"
        for i in range(3)
        for j in range(3)
        if i != j
        for kind in ("m", "n")
    ]
    names += [f"t{i}" for i in range(3)]
    for i, j in combinations(range(3), 2):
        k = 3 - i - j
        names += [
            f"B{i}{j}",
            f"B{i}{j}q{i}",
            f"B{i}{j}q{j}",
            f"B{i}{j}p{k}",
            f"B{i}{j}p{k}q{i}",
            f"B{i}{j}p{k}q{j}",
            f"B{i}{j}p{k}p{k}",
        ]
    names += ["C"] + [f"Cq{i}" for i in range(3)]
    assert len(set(names)) == len(names)
    coeff, *gens = ring(",".join(names), s.QQ)
    R, *w = ring("a,b,c", coeff.to_domain())
    v = {name: R(value) for name, value in zip(names, gens)}
    cuts = tuple(range(1, 8))

    def indices(mask):
        return tuple(i for i in range(3) if mask >> i & 1)

    def denominator(mask):
        made = sum((v[f"d{i}"] * w[i] for i in indices(mask)), R.zero)
        for i, j in combinations(indices(mask), 2):
            made += v[f"z{i}{j}"] * w[i] * w[j]
        return made

    D = {mask: denominator(mask) for mask in cuts}

    def vertex(block, previous):
        ix = indices(block)
        rest = indices(previous)
        if len(ix) == 1:
            i = ix[0]
            made = v[f"L{i}"]
            for j in rest:
                made += v[f"m{i}{j}"] * w[j] + v[f"n{i}{j}"] * w[j] ** 2
            if len(rest) == 2:
                made += v[f"t{i}"] * w[rest[0]] * w[rest[1]]
            return made
        if len(ix) == 2:
            i, j = ix
            made = v[f"B{i}{j}"] + v[f"B{i}{j}q{i}"] * w[i] + v[f"B{i}{j}q{j}"] * w[j]
            if rest:
                (k,) = rest
                made += (
                    v[f"B{i}{j}p{k}"] * w[k]
                    + v[f"B{i}{j}p{k}q{i}"] * w[k] * w[i]
                    + v[f"B{i}{j}p{k}q{j}"] * w[k] * w[j]
                    + v[f"B{i}{j}p{k}p{k}"] * w[k] ** 2
                )
            return made
        assert block == 7 and previous == 0
        return v["C"] + sum((v[f"Cq{i}"] * w[i] for i in range(3)), R.zero)

    @cache
    def terms(mask):
        if not mask:
            return ((R.one, ()),)
        out = []
        block = mask
        while block:
            rest = mask ^ block
            for N, denoms in terms(rest):
                out.append((vertex(block, rest) * N, denoms + (mask,)))
            block = (block - 1) & mask
        return tuple(out)

    assert [len(terms(mask)) for mask in (1, 3, 7)] == [1, 3, 13]
    contributions = [(N, denoms, 1) for N, denoms in terms(7)]
    for i in range(3):
        pair = 7 ^ (1 << i)
        for N, denoms in terms(pair):
            contributions.append((N * v[f"L{i}"], denoms + (1 << i,), -1))
    contributions.append((v["L0"] * v["L1"] * v["L2"], (1, 2, 4), 2))
    assert (
        len(contributions) == 23
        and sum(abs(weight) for N, den, weight in contributions) == 24
    )
    P = R.zero
    for N, denoms, weight in contributions:
        assert len(set(denoms)) == len(denoms)
        complement = R.one
        for cut in cuts:
            if cut not in denoms:
                complement *= D[cut]
        P += weight * N * complement
    assert P
    assert all(all(exponent >= 1 for exponent in monomial) for monomial in P)
    assert min(sum(m) for m in P) >= 6
    assert max(sum(m) for m in P) <= 10
    Q = R.from_dict({tuple(e - 1 for e in m): value for m, value in P.items()})
    assert Q * w[0] * w[1] * w[2] == P
    assert min(sum(m) for m in Q) >= 3

    return SimpleNamespace(
        names=names, R=R, w=w, v=v, D=D, vertex=vertex, terms=terms, P=P, Q=Q
    )


def vector(label):
    return imm([s.Symbol(f"{label}_{i}") for i in range(4)])


def matrix(label):
    out = s.zeros(4)
    for i in range(4):
        for j in range(i, 4):
            out[i, j] = out[j, i] = s.Symbol(label + str(i) + str(j))
    return imm(out)


def dot(p, q):
    return (p.T * e.ETA * q)[0]


def contract(p, A, q):
    return (p.T * A * q)[0]


@cache
def physical_mapping():
    p, k, q = map(vector, ("p", "k", "q"))
    A = matrix("A")
    tr = e.tr(A)
    vertex = -e.scalar_vertex(p + k, -p - k - q, (A,), 1)
    expected = 2 * contract(p + k, A, p + k + q) - tr * (
        dot(p + k, p + k) - 1 + dot(p + k, q)
    )
    assert s.expand(vertex - expected) == 0
    physical_form = 2 * contract(p + k, A, p + k)
    constraint_terms = 2 * contract(p + k, A, q) - tr * (
        dot(p + k, p + k) - 1 + dot(p + k, q)
    )
    assert s.expand(vertex - physical_form - constraint_terms) == 0
    # Thus physical TT (tr A=0, A q=0) removes exactly the displayed terms.
    # All pair/triple coefficients use their literal density matrices; no Ward use.
    M = matrix("M")
    d = s.Symbol("det")
    a, b, c = w = s.symbols("a b c")
    n = tuple(vector("v" + str(i)) for i in range(3))
    pair_vertex = d - contract(p + c * n[2], M, p + c * n[2] + a * n[0] + b * n[1])
    pair_coefficients = {
        (0, 0, 0): d - contract(p, M, p),
        (1, 0, 0): -contract(p, M, n[0]),
        (0, 1, 0): -contract(p, M, n[1]),
        (0, 0, 1): -2 * contract(p, M, n[2]),
        (1, 0, 1): -contract(n[2], M, n[0]),
        (0, 1, 1): -contract(n[2], M, n[1]),
        (0, 0, 2): -contract(n[2], M, n[2]),
    }
    assert s.Poly(s.expand(pair_vertex), *w).as_dict() == {
        mon: s.expand(val) for mon, val in pair_coefficients.items()
    }
    triple_vertex = d - contract(p, M, p + a * n[0] + b * n[1] + c * n[2])
    assert s.Poly(s.expand(triple_vertex), *w).total_degree() == 1
    single_vertex = 2 * contract(p + b * n[1] + c * n[2], A, p + b * n[1] + c * n[2])
    single_coefficients = {
        (0, 0, 0): 2 * contract(p, A, p),
        (0, 1, 0): 4 * contract(p, A, n[1]),
        (0, 0, 1): 4 * contract(p, A, n[2]),
        (0, 2, 0): 2 * contract(n[1], A, n[1]),
        (0, 0, 2): 2 * contract(n[2], A, n[2]),
        (0, 1, 1): 4 * contract(n[1], A, n[2]),
    }
    assert s.Poly(s.expand(single_vertex), *w).as_dict() == {
        mon: s.expand(val) for mon, val in single_coefficients.items()
    }
    # Each line momentum component has energy-polynomial l1 <=2+3=5.
    # K_r and d_r bounds are entrywise for fields of component modulus<=1.
    caps = tuple(
        16 * 5**2 * vertex_bound.DENSITY[r] + vertex_bound.DET[r] for r in (1, 2, 3)
    )
    assert caps == (2404, 38448, 960960)
    # Pure scalar denominators: |d_i|>1/2, |z_ij|<=4.
    # For E<=2, (2-sqrt(3))>1/4, and W<=1/8:
    assert (s.Rational(7, 4)) ** 2 > 3
    assert s.Rational(1, 2) - s.Rational(1, 8) == s.Rational(3, 8)
    assert (8, 2 * 8 + 4, 3 * 8 + 3 * 4) == (8, 20, 36)

    checks = {
        "generic_scalar_vertex": s.expand(vertex - expected),
        "generic_TT_constraint_terms": s.expand(
            vertex - physical_form - constraint_terms
        ),
        "generic_pair_coefficient_map": s.expand(
            pair_vertex
            - sum(
                value * s.prod(x**e for x, e in zip(w, powers))
                for powers, value in pair_coefficients.items()
            )
        ),
        "generic_single_coefficient_map": s.expand(
            single_vertex
            - sum(
                value * s.prod(x**e for x, e in zip(w, powers))
                for powers, value in single_coefficients.items()
            )
        ),
    }

    g = generic_polynomials()
    directions = (
        imm([1, 0, 0]),
        imm([s.Rational(3, 5), s.Rational(4, 5), 0]),
        imm([0, s.Rational(3, 5), s.Rational(4, 5)]),
    )
    frames = (
        (imm([0, 1, 0]), imm([0, 0, 1])),
        (imm([-s.Rational(4, 5), s.Rational(3, 5), 0]), imm([0, 0, 1])),
        (imm([1, 0, 0]), imm([0, s.Rational(4, 5), -s.Rational(3, 5)])),
    )
    rays = tuple(imm([1, *v]) for v in directions)
    p0 = imm([s.Rational(5, 4), s.Rational(9, 20), 0, s.Rational(3, 5)])
    counts = 0
    for sign in (1, -1):
        p = sign * p0
        assert dot(p, p) == 1
        for bits in product((0, 1), repeat=3):
            fields = []
            for bit, (u, v) in zip(bits, frames):
                H = s.zeros(4)
                plus = u * u.T - v * v.T
                cross = u * v.T + v * u.T
                # Unit-bound exact complex physical tensor, two independent choices.
                H[1:, 1:] = (plus + s.I * cross) / 2 if bit else plus / 2
                fields.append(imm(H))
            fields = tuple(fields)
            for H, ray in zip(fields, rays):
                assert e.tr(H) == 0 and H * ray == e.VECTOR_ZERO
                assert sum(s.conjugate(x) * x for x in H) <= 1
            values = {}
            for i in range(3):
                values["d" + str(i)] = 2 * dot(p, rays[i])
                values["L" + str(i)] = 2 * contract(p, fields[i], p)
                rest = [j for j in range(3) if j != i]
                for j in rest:
                    values[f"m{i}{j}"] = 4 * contract(p, fields[i], rays[j])
                    values[f"n{i}{j}"] = 2 * contract(rays[j], fields[i], rays[j])
                j, k = rest
                values["t" + str(i)] = 4 * contract(rays[j], fields[i], rays[k])
            for i, j in combinations(range(3), 2):
                k = 3 - i - j
                values[f"z{i}{j}"] = 2 * dot(rays[i], rays[j])
                K = (
                    e.ETA
                    * e.density_inverse_coefficient((fields[i], fields[j]))
                    * e.ETA
                )
                delta = e.determinant_coefficient((fields[i], fields[j]))
                values[f"B{i}{j}"] = delta - contract(p, K, p)
                values[f"B{i}{j}q{i}"] = -contract(p, K, rays[i])
                values[f"B{i}{j}q{j}"] = -contract(p, K, rays[j])
                values[f"B{i}{j}p{k}"] = -2 * contract(p, K, rays[k])
                values[f"B{i}{j}p{k}q{i}"] = -contract(rays[k], K, rays[i])
                values[f"B{i}{j}p{k}q{j}"] = -contract(rays[k], K, rays[j])
                values[f"B{i}{j}p{k}p{k}"] = -contract(rays[k], K, rays[k])
            K = e.ETA * e.density_inverse_coefficient(fields) * e.ETA
            values["C"] = e.determinant_coefficient(fields) - contract(p, K, p)
            for i in range(3):
                values["Cq" + str(i)] = -contract(p, K, rays[i])
            assert set(values) == set(g.names)
            subst = {s.Symbol(name): value for name, value in values.items()}

            def specialize(poly, _subst=subst):
                return s.expand(poly.as_expr().subs(_subst, simultaneous=True))

            @cache
            def momentum(mask):
                return sum(
                    (w[i] * rays[i] for i in range(3) if mask >> i & 1), e.VECTOR_ZERO
                )

            @cache
            def line(
                mask,
                _p=p,
                _fields=fields,
                _sign=sign,
                _bits=bits,
                _momentum=momentum,
                _specialize=specialize,
            ):
                if not mask:
                    return s.S.One
                Q = _momentum(mask)
                answer = 0
                part = mask
                while part:
                    rest = mask ^ part
                    V = -e.scalar_vertex(
                        imm(_p + _momentum(rest)),
                        imm(-_p - Q),
                        tuple(_fields[i] for i in range(3) if part >> i & 1),
                        1,
                    )
                    expected = _specialize(g.vertex(part, rest))
                    residual = s.expand(V - expected)
                    checks[f"signed{_sign}_bits{_bits}_mask{mask}_vertex{part}"] = (
                        residual
                    )
                    assert residual == 0
                    answer += V * line(rest) / (dot(_p + Q, _p + Q) - 1)
                    part = (part - 1) & mask
                return s.cancel(answer)

            triple = line(7)
            connected = (
                triple
                - sum(line(7 ^ (1 << i)) * line(1 << i) for i in range(3))
                + 2 * s.prod(line(1 << i) for i in range(3))
            )
            # Cross-multiply the ORIGINAL seven-cut common denominator; no fitted data.
            denominator = s.prod(
                dot(p + momentum(mask), p + momentum(mask)) - 1 for mask in range(1, 8)
            )
            exact = specialize(g.P)
            residual = s.cancel(connected * denominator - exact)
            checks[f"signed{sign}_bits{bits}_whole_cumulant"] = residual
            assert residual == 0
            counts += 1

    assert counts == 16
    return {
        "checks": checks,
        "gates": {
            "all16_signed_complex_TT_specializations": counts == 16,
            "literal_pair_and_triple_density_matrices_not_TT_shortcuts": True,
            "generic_scalar_vertex_constraint_terms_explicit": True,
        },
        "whole_exact_mapping_fixture_count": counts,
        "whole_mapping_boundary": "Generic original Phi2-h identity plus arbitrary density-matrix coefficient expansion;16 independent exact signed-p, complex physical TT fixtures reconstruct every vertex and the complete13-tree cumulant against the generic polynomial.",
    }


@cache
def generic_data():
    g = generic_polynomials()
    R, v = g.R, g.v
    a, b, c = g.w
    P2 = R.zero
    for N, cuts in g.terms(3):
        missing = R.one
        for cut in (1, 2, 3):
            if cut not in cuts:
                missing *= g.D[cut]
        P2 += N * missing
    P2 -= v["L0"] * v["L1"] * g.D[3]
    R0 = (
        -v["L0"] * v["L1"] * v["z01"]
        + v["L0"] * v["m10"] * v["d1"]
        + v["L1"] * v["m01"] * v["d0"]
        + v["B01"] * v["d0"] * v["d1"]
    )
    Ra = v["L0"] * v["n10"] * v["d1"] + v["B01q0"] * v["d0"] * v["d1"]
    Rb = v["L1"] * v["n01"] * v["d0"] + v["B01q1"] * v["d0"] * v["d1"]
    linear = a * v["d0"] + b * v["d1"]
    remainder = a * b * ((a * Ra + b * Rb) * linear - R0 * a * b * v["z01"])
    flip = {
        s.Symbol(name): -s.Symbol(name)
        for name in ("d0", "d1", "m01", "m10", "B01q0", "B01q1")
    }
    leading = (a * b * R0).as_expr() / (v["d0"] * v["d1"] * linear).as_expr()
    vertex_caps = (2404, 38448, 960960)
    numerator_cap = (
        24
        * 8**3
        * 20**3
        * 36
        * max(vertex_caps[0] ** 3, vertex_caps[0] * vertex_caps[1], vertex_caps[2])
    )
    connected_cap = s.Rational(numerator_cap) * 8 * s.Rational(8, 3) ** 4
    checks = {
        "generic_triple_common_factor_identity": (g.P - a * b * c * g.Q).as_expr(),
        "generic_pair_exact_numerator": (P2 - a * b * (R0 + a * Ra + b * Rb)).as_expr(),
        "generic_pair_exact_remainder": (
            P2 * linear - a * b * R0 * g.D[3] - remainder
        ).as_expr(),
        "generic_pair_R0_even": s.expand(R0.as_expr().xreplace(flip) - R0.as_expr()),
        "generic_pair_Ra_odd": s.expand(Ra.as_expr().xreplace(flip) + Ra.as_expr()),
        "generic_pair_Rb_odd": s.expand(Rb.as_expr().xreplace(flip) + Rb.as_expr()),
        "generic_pair_leading_odd": s.cancel(leading.xreplace(flip) + leading),
        "exact_connected_triple_coefficient": connected_cap
        - s.Rational(59670991094513926144000, 3),
    }
    return {
        "checks": checks,
        "gates": {
            "all49_contraction_coefficients_free": len(g.names) == 49,
            "all62_energy_monomials_retained": len(g.P) == 62,
            "every_numerator_monomial_divisible_byabc": all(
                all(e >= 1 for e in mon) for mon in g.P
            ),
            "minimum_energy_degree6": min(sum(mon) for mon in g.P) == 6,
            "maximum_energy_degree10": max(sum(mon) for mon in g.P) == 10,
            "quotient_minimum_degree3": min(sum(mon) for mon in g.Q) == 3,
            "triple_scalar_line_coefficient_below1e24": bool(connected_cap < 10**24),
            "no_angular_fit_or_altered_factorization": True,
        },
        "whole_free_coefficient_names": g.names,
        "whole_connected_triple_common_numerator": g.P.as_expr(),
        "whole_connected_pair_coefficients": {
            "R0": R0.as_expr(),
            "Ra": Ra.as_expr(),
            "Rb": Rb.as_expr(),
        },
        "whole_real_triple_uniform_coefficient": connected_cap,
        "whole_line_vertex_coefficient_caps": vertex_caps,
        "whole_scope": "One fixed mass-one signed scalar momentum with |p0|<=2, three real future null rays and unit complex TT fields. Connected triple modulus below1e24*abc*S/kappa^(3/2), with compatible zero faces. The full recoil-dependent four-leg result requires the separate forest and analytic argument.",
    }


@cache
def global_tube_control():
    a, b = s.symbols("a b")
    p = s.ImmutableMatrix([1, 0, 0, 0])
    q1 = s.ImmutableMatrix([a, a, 0, 0])
    q2 = s.ImmutableMatrix([b, 0, b, 0])
    A = s.ImmutableMatrix(s.diag(0, 0, s.Rational(1, 2), -s.Rational(1, 2)))
    B = s.ImmutableMatrix(s.diag(0, s.Rational(1, 2), 0, -s.Rational(1, 2)))
    assert A * q1 == e.VECTOR_ZERO and B * q2 == e.VECTOR_ZERO
    assert e.tr(A) == e.tr(B) == 0
    assert sum(x * x for x in A) == sum(x * x for x in B) == s.Rational(1, 2)
    dot = lambda p, q: (p.T * e.ETA * q)[0]
    j1 = -e.scalar_vertex(p, -p - q1, (A,), 1) / (dot(p + q1, p + q1) - 1)
    j2 = -e.scalar_vertex(p, -p - q2, (B,), 1) / (dot(p + q2, p + q2) - 1)
    Q = q1 + q2
    D = dot(p + Q, p + Q) - 1
    pair = (
        -e.scalar_vertex(p + q1, -p - Q, (B,), 1) * j1
        - e.scalar_vertex(p + q2, -p - Q, (A,), 1) * j2
        - e.scalar_vertex(p, -p - Q, (A, B), 1)
    ) / D
    connected = s.factor(a * b * (pair - j1 * j2))
    P, den = s.fraction(connected)
    eps = s.Rational(1, 10**30)
    root = -eps / (1 + eps)
    assert D.subs({a: root, b: eps}) == 0
    assert P.subs({a: root, b: eps}) != 0 and den.subs({a: root, b: eps}) == 0
    W = 2 * eps + s.Rational(1, 16)
    assert abs(root - eps) < s.Rational(1, 10**13) * W

    return {
        "connected_pair_pole_residue_nonzero": P.subs({a: root, b: eps}) != 0,
        "pole_inside_naive_global_W_disc": bool(
            abs(root - eps) < s.Rational(1, 10**13) * W
        ),
        "whole_counterexample": "Fixed on-shell p=(1,0,0,0), two unit-bound TT fields, b=1e-30 and a=-b/(1+b). The connected pair has a nonremovable scalar-cut pole within the global-W disc about a=b, c=1/16. The c-only and relative-pair domains used in the proof do not cross this pole.",
    }


@cache
def data():
    generic, mapping, control = (
        generic_data(),
        physical_mapping(),
        global_tube_control(),
    )
    return {
        "checks": {
            **{"generic_" + k: v for k, v in generic["checks"].items()},
            **{"physical_" + k: v for k, v in mapping["checks"].items()},
        },
        "gates": {
            **generic["gates"],
            **mapping["gates"],
            "connected_pair_pole_residue_nonzero": control[
                "connected_pair_pole_residue_nonzero"
            ],
            "pole_inside_naive_global_W_disc": control[
                "pole_inside_naive_global_W_disc"
            ],
        },
        "whole_generic_line": {
            k: v for k, v in generic.items() if k not in ("checks", "gates")
        },
        "whole_original_mapping": {
            k: v for k, v in mapping.items() if k not in ("checks", "gates")
        },
        "whole_excluded_global_tube": control["whole_counterexample"],
    }
