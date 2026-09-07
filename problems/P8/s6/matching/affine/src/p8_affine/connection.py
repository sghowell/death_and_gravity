"""Literal, unrestricted affine-connection elimination in a timelike frame.

Source convention: (-+++), X=-s**2 and Gamma^a_{bc}, with the derivative
index LAST.  No symmetry of the lower connection indices is imposed.
The functions below concern an algebraic auxiliary field, not a propagator.
"""

from functools import cache
from itertools import permutations, product

import sympy as sp

P, S = sp.symbols("p s", positive=True)
C, F3, F4 = sp.symbols("c F3 F4", real=True)
PP, PX, CP, CX = sp.symbols("p_phi p_X c_phi c_X", real=True)
SIGN = (-1, 1, 1, 1)
INDICES = tuple(product(range(4), repeat=3))
KAPPA = sp.symbols("kappa_0:64", real=True)
H_PAIRS = tuple((a, b) for a in range(4) for b in range(a, 4))
H_SYMBOLS = sp.symbols("H00 H01 H02 H03 H11 H12 H13 H22 H23 H33", real=True)
H = sp.zeros(4)
for _pair, _symbol in zip(H_PAIRS, H_SYMBOLS, strict=True):
    H[_pair[0], _pair[1]] = H[_pair[1], _pair[0]] = _symbol
H = sp.ImmutableMatrix(H)


def index(a, b, c):
    """Index of kappa^a_{bc}; accept only literal integer tensor indices."""
    for value in (a, b, c):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError("Tensor indices must be integers")
        if not 0 <= value < 4:
            raise ValueError("Tensor index outside 0,...,3")
    return 16*a + 4*b + c


def _epsilon(indices):
    if len(set(indices)) != 4:
        return 0
    inversions = sum(indices[i] > indices[j] for i in range(4) for j in range(i+1, 4))
    return (-1)**inversions


def _factor_matrix(matrix):
    def normalized(value):
        result = sp.factor(value)
        return sp.S.Zero if result.is_zero is True else result

    return sp.ImmutableMatrix(matrix.applyfunc(normalized))


@cache
def literal():
    """Build the full 64-component polynomial from the curvature contractions.

    ``hessian`` is d_k d_k L, not the matrix multiplying k.k without 1/2.
    ``source`` follows by integrating the two derivatives of kappa in
    R^a_{bcd} by parts BEFORE imposing the rest frame.  The double-dual
    Einstein coupling uses its literal two-epsilon contraction.
    """
    matrix = sp.zeros(64)
    source = sp.zeros(64, 1)

    def add(left, right, coefficient):
        matrix[left, right] += coefficient
        matrix[right, left] += coefficient

    # p*g^{bd} R^a_{bad}, using the full affine curvature KK product.
    for a, b, e in product(range(4), repeat=3):
        weight = P*SIGN[b]
        add(index(a, e, a), index(e, b, b), weight)
        add(index(a, e, b), index(e, b, a), -weight)
    for a, b in product(range(4), repeat=2):
        dp_a = S*(PP*int(a == 0)-2*PX*H[a, 0])
        dp_b = S*(PP*int(b == 0)-2*PX*H[b, 0])
        source[index(a, b, b)] -= SIGN[b]*dp_a
        source[index(a, b, a)] += SIGN[b]*dp_b

    # At v_a=(s,0,0,0), c*G_Gamma^{ab}v_a v_b is c*s²/2 sum_ij R^i_{jij}.
    for i, j, e in product(range(1, 4), range(1, 4), range(4)):
        weight = C*S**2/2
        add(index(i, e, i), index(e, j, j), weight)
        add(index(i, e, j), index(e, j, i), -weight)

    # Derivative of c*v_alpha*v_beta in the unrestricted epsilon contraction.
    # eps^{gamma alpha a b} eps_gamma^{beta c d} g_aa / 4.
    for gamma, alpha, a, b in permutations(range(4)):
        for beta, c, d in permutations(i for i in range(4) if i != gamma):
            weight = sp.Rational(SIGN[a]*SIGN[gamma]
                                 *_epsilon((gamma, alpha, a, b))
                                 *_epsilon((gamma, beta, c, d)), 4)

            def derivative(e, alpha=alpha, beta=beta):
                coefficient = S*(CP*int(e == 0)-2*CX*H[e, 0])
                return (coefficient*S**2*int(alpha == beta == 0)
                        + C*S*(H[e, alpha]*int(beta == 0)
                               + int(alpha == 0)*H[e, beta]))

            source[index(a, b, d)] -= weight*derivative(c)
            source[index(a, b, c)] += weight*derivative(d)

    trace_spatial = sum(H[i, i] for i in range(1, 4))
    # Literal L3=2s² tr(H_ij-s*kappa^0_{ji}); L4=s²[(tr F)^2-tr(F²)].
    for i in range(1, 4):
        source[index(0, i, i)] -= 2*F3*S**3 + 2*F4*S**3*trace_spatial
    for i, j in product(range(1, 4), repeat=2):
        source[index(0, i, j)] += 2*F4*S**3*H[i, j]
        add(index(0, i, i), index(0, j, j), F4*S**4)
        add(index(0, j, i), index(0, i, j), -F4*S**4)

    tr_h = sum(SIGN[a]*H[a, a] for a in range(4))
    h_squared = sum(SIGN[a]*SIGN[b]*H[a, b]**2 for a, b in product(range(4), repeat=2))
    vhv = S**2*H[0, 0]
    vh2v = S**2*sum(SIGN[a]*H[0, a]**2 for a in range(4))
    # The remaining Levi-Civita c*R_ab*v^a*v^b is reduced by its full IBP.
    bare = (2*F3*S**2*trace_spatial
            + F4*S**2*(trace_spatial**2-sum(H[i, j]**2
                                           for i, j in product(range(1, 4), repeat=2)))
            + C*(tr_h**2-h_squared)
            + CP*(-S**2*tr_h-vhv) + 2*CX*(vhv*tr_h-vh2v))
    vector = sp.Matrix(KAPPA)
    return {
        "symbols": {"p": P, "s": S, "X": -S**2, "c": C, "F3": F3, "F4": F4,
                    "p_phi": PP, "p_X": PX, "c_phi": CP, "c_X": CX},
        "indices": INDICES, "kappa": KAPPA, "H": H,
        "hessian": _factor_matrix(matrix), "source": _factor_matrix(source),
        "quadratic_density": sp.expand((vector.T*matrix*vector)[0]/2),
        "linear_density": sp.expand((source.T*vector)[0]),
        "bare_reduced_density": sp.expand(bare), "curvature_coefficient": P+C*S**2/2,
        "invariants": (h_squared, tr_h**2, vhv*tr_h, vh2v, vhv**2),
        "trace_H": tr_h, "vHv": vhv,
    }


@cache
def cd_substitution():
    """The CD pointwise values; forcing derivatives are NOT specialized here."""
    return {C: (4*P**2-2*P)/S**2, F4: (sp.Rational(1, 2)-2*P**2)/S**4}


@cache
def quadratic():
    """Full CD quadratic form and four exact projective gauge columns."""
    data = literal()
    gauge = sp.zeros(64, 4)
    for a, c in product(range(4), repeat=2):
        gauge[index(a, a, c), c] = 1
    matrix = _factor_matrix(data["hessian"].subs(cd_substitution()))
    source = _factor_matrix(data["source"].subs(cd_substitution()))
    return {"hessian": matrix, "source": source,
            "gauge": sp.ImmutableMatrix(gauge),
            "gauge_hessian_residual": matrix*gauge,
            "gauge_source_residual": gauge.T*source}


@cache
def quotient():
    """Trace gauge kappa^a_{ac}=0, all 60 directions, with exact block factors."""
    full = quadratic()
    omitted = tuple(index(0, 0, c) for c in range(4))
    kept = tuple(i for i in range(64) if i not in omitted)
    embed = sp.zeros(64, 60)
    for column, component in enumerate(kept):
        embed[component, column] = 1
        a, b, c = INDICES[component]
        if a == b:
            embed[index(0, 0, c), column] = -1
    matrix = _factor_matrix(embed.T*full["hessian"]*embed)
    remaining = set(range(60))
    blocks = []
    while remaining:
        queue = [min(remaining)]
        remaining.remove(queue[0])
        component = []
        while queue:
            row = queue.pop()
            component.append(row)
            for column in sorted(remaining):
                if matrix[row, column] != 0:
                    remaining.remove(column)
                    queue.append(column)
        blocks.append(tuple(sorted(component)))
    block_matrices = tuple(matrix.extract(block, block) for block in blocks)
    determinants = tuple(sp.factor(block.det(method="domain-ge")) for block in block_matrices)
    return {"embedding": sp.ImmutableMatrix(embed), "kept": kept,
            "hessian": matrix, "blocks": tuple(blocks),
            "block_matrices": block_matrices, "block_determinants": determinants,
            "determinant": sp.factor(sp.prod(determinants)),
            "source": _factor_matrix(embed.T*full["source"])}


@cache
def eliminate():
    """Solve every Euler equation for p>0, 8p²!=1, with generic derivative forcing.

    The returned rational identity is not an assertion of uniqueness at the
    exceptional value p²=1/8, even if the particular source cancels that pole.
    """
    data = quotient()
    solution = sp.zeros(60, 1)
    for indices, matrix in zip(data["blocks"], data["block_matrices"], strict=True):
        forcing = data["source"].extract(indices, (0,))
        answer = -matrix.inv(method="DM")*forcing
        for i, value in zip(indices, answer, strict=True):
            solution[i] = sp.factor(value)
    full_solution = _factor_matrix(data["embedding"]*solution)
    reduced = sp.factor(literal()["bare_reduced_density"].subs(cd_substitution())
                        + (quadratic()["source"].T*full_solution)[0]/2)
    return {"quotient_solution": sp.ImmutableMatrix(solution),
            "solution": full_solution, "reduced_density": sp.expand(reduced),
            "full_euler_residual": _factor_matrix(quadratic()["hessian"]*full_solution
                                                  + quadratic()["source"])}


@cache
def coefficients():
    """Read all five independent quadratic and both linear invariant coefficients."""
    density = eliminate()["reduced_density"]
    poly = sp.Poly(density, *H_SYMBOLS)
    def coefficient(monomial):
        return sp.factor(poly.coeff_monomial(monomial))

    a1 = coefficient(H[1, 2]**2)/2
    a2 = coefficient(H[1, 1]*H[2, 2])/2
    a3 = (coefficient(H[0, 0]*H[1, 1])+2*a2)/S**2
    a4 = (coefficient(H[0, 1]**2)+2*a1)/S**2
    a5 = (coefficient(H[0, 0]**2)-a1-a2+(a3+a4)*S**2)/S**4
    q1 = coefficient(H[1, 1])
    q2 = (coefficient(H[0, 0])+q1)/S**2
    out = {name: sp.factor(value) for name, value in zip(
        ("P", "Q1", "Q2", "A1", "A2", "A3", "A4", "A5"),
        (coefficient(1), q1, q2, a1, a2, a3, a4, a5), strict=True)}
    reconstructed = out["P"]+out["Q1"]*literal()["trace_H"]+out["Q2"]*literal()["vHv"]
    reconstructed += sum(out[f"A{i+1}"]*item
                         for i, item in enumerate(literal()["invariants"]))
    out["density_reconstruction_residual"] = sp.factor(density-reconstructed)
    return out


@cache
def source_formula_comparison():
    """Compare the literal result with (4.5)--(4.12), correcting ONLY (4.8).

    Both versions of Q2 are retained.  This comparison does not enter the
    construction, inversion, or coefficient extraction above.
    """
    x = -S**2
    f = P-C*x/2
    f_phi = PP-CP*x/2
    delta = 2*P-C*x+2*F4*x**2
    t = PP-F3*x
    expected = {
        "P": 3*x*t**2/delta,
        "Q1": -2*f_phi+4*P*t/delta,
        "Q2": 2*f_phi/x-4*(P-3*x*PX)*t/(x*delta),
        "A1": -C/2-P*(C-2*F4*x)/delta,
        "A2": C/2+P*(C-2*F4*x)/delta,
        "A3": 2*CX+(4*P*F4+(4*PX-C)*(C-2*F4*x))/delta,
        "A4": (-2*CX+2*PX*(3*PX-C)/P
               + PX*x*(PX*C-4*P*CX)/P**2
               + (C**2-4*P*F4-2*C*F4*x)/delta),
        "A5": (-PX*(PX*C-4*P*CX)/P**2
               + 2*PX*(4*P*F4+(3*PX-C)*(C-2*F4*x))/(P*delta)),
    }
    expected = {key: sp.factor(value.subs(cd_substitution()))
                for key, value in expected.items()}
    direct = coefficients()
    printed = (2*f_phi/x-4*(P-3*PX)*t/(x*delta)).subs(cd_substitution())
    return {"expected": expected,
            "residuals": {key: sp.factor(direct[key]-value)
                          for key, value in expected.items()},
            "printed_Q2": sp.factor(printed),
            "printed_Q2_defect": sp.factor(direct["Q2"]-printed),
            "curvature_coefficient": sp.factor(f.subs(cd_substitution()))}


@cache
def palatini_control():
    """Independent pure p(phi,X)R_Gamma control, at arbitrary positive p.

    The trace-gauge connection is obtained from the conformal Levi-Civita
    solution and is checked against all 64 literal Euler equations.
    """
    data = literal()
    specialization = {C: 0, CP: 0, CX: 0, F3: 0, F4: 0}
    gradient = sp.Matrix([S*(PP*int(a == 0)-2*PX*H[a, 0]) for a in range(4)])
    solution = sp.Matrix([
        (int(a == c)*gradient[b]-int(b == c)*SIGN[b]*SIGN[a]*gradient[a])/(2*P)
        for a, b, c in INDICES])
    matrix = data["hessian"].subs(specialization)
    forcing = data["source"].subs(specialization)
    density = sp.factor((forcing.T*solution)[0]/2)
    expected = 3*sum(SIGN[a]*gradient[a]**2 for a in range(4))/(2*P)
    poly = sp.Poly(density, *H_SYMBOLS)
    q1 = poly.coeff_monomial(H[1, 1])
    q2 = sp.factor((poly.coeff_monomial(H[0, 0])+q1)/S**2)
    return {"solution": sp.ImmutableMatrix(solution), "density": density,
            "expected_density": expected, "Q1": q1, "Q2": q2,
            "full_euler_residual": _factor_matrix(matrix*solution+forcing),
            "density_residual": sp.factor(density-expected),
            "Q2_residual": sp.factor(q2-6*PP*PX/P),
            "printed_Q2_defect": sp.factor(q2+6*PP*PX/(P*S**2))}


@cache
def exceptional():
    """Extra algebraic kernel at p²=1/8; not a dynamical ghost diagnosis."""
    data = quotient()
    value = sp.sqrt(2)/4
    kernel = []
    for indices, matrix in zip(data["blocks"], data["block_matrices"], strict=True):
        for local in matrix.subs(P, value).nullspace():
            vector = sp.zeros(60, 1)
            for i, coefficient in zip(indices, local, strict=True):
                vector[i] = coefficient
            kernel.append(sp.ImmutableMatrix(vector))
    kernel_matrix = sp.ImmutableMatrix.hstack(*kernel)
    return {"p": value, "quotient_kernel": kernel_matrix,
            "quotient_nullity": len(kernel), "quotient_rank": 60-len(kernel),
            "kernel_residual": _factor_matrix(data["hessian"].subs(P, value)*kernel_matrix),
            "source_compatibility_residual": _factor_matrix(
                kernel_matrix.T*data["source"].subs(P, value))}


def require_domain(p, s):
    """Validate exact finite real parameters before claiming an invertible quotient."""
    values = []
    for name, raw in (("p", p), ("s", s)):
        if isinstance(raw, (bool, float, str)):
            raise TypeError(f"{name} must be exact, not a bool, binary float, or string")
        value = sp.sympify(raw)
        if not isinstance(value, sp.Expr):
            raise TypeError(f"{name} must be an exact scalar")
        if value.has(sp.Float) or value.free_symbols:
            raise ValueError(f"{name} must be an exact constant")
        if value.is_real is not True or value.is_finite is not True:
            raise ValueError(f"{name} must be finite and real")
        if value.is_positive is not True:
            raise ValueError(f"{name} must be strictly positive")
        values.append(value)
    if sp.simplify(8*values[0]**2-1).is_zero is not False:
        raise ValueError("The quotient must exclude p²=1/8")
    return tuple(values)


@cache
def inverse_bound():
    """Explicit coordinate l-infinity inverse bound on the entire CD tube.

    The proof uses 9/40<=p²<=11/40, hence 9/20<p<11/20 and 8p²-1>=4/5.
    Numerator coefficients are bounded by absolute values; denominator factors
    are retained.  This is the row-sum norm in the displayed trace-gauge
    coordinates in an orthonormal scalar rest frame, not an invariant
    Lorentz-frame norm, kinetic residue, or heavy mass/gap.
    """
    p_min, p_max, d_min = sp.Rational(9, 20), sp.Rational(11, 20), sp.Rational(4, 5)

    def entry_bound(expression):
        if expression == 0:
            return sp.S.Zero
        numerator, denominator = sp.fraction(sp.cancel(expression))
        poly = sp.Poly(numerator, P)
        numerator_bound = sum(abs(value)*p_max**degree[0] for degree, value in poly.terms())
        constant, factors = sp.factor_list(denominator, P)
        denominator_bound = abs(constant)
        for factor, power in factors:
            if factor == P:
                denominator_bound *= p_min**power
            elif factor == 8*P**2-1:
                denominator_bound *= d_min**power
            else:
                raise ArithmeticError(f"Unproved denominator factor: {factor}")
        return numerator_bound/denominator_bound

    row_bounds = []
    inverse_residuals = []
    for matrix in quotient()["block_matrices"]:
        inverse = matrix.inv(method="DM")
        inverse_residuals.append(_factor_matrix(matrix*inverse-sp.eye(matrix.rows)))
        row_bounds.extend(sum(entry_bound(inverse[i, j]) for j in range(matrix.cols))
                          for i in range(matrix.rows))
    return {"row_bounds": tuple(row_bounds), "norm_upper": max(row_bounds),
            "p_lower_used": p_min, "p_upper_used": p_max,
            "exceptional_factor_lower_used": d_min,
            "inverse_residuals": tuple(inverse_residuals)}


@cache
def checks():
    """Exact zero residuals; all finite-dimensional checks replay the literal action."""
    data = quadratic()
    q = quotient()
    raw = literal()
    gauge = data["gauge"]
    result = {
        "generic_projective_hessian": _factor_matrix(raw["hessian"]*gauge),
        "generic_projective_source": _factor_matrix(gauge.T*raw["source"]),
        "cd_projective_hessian": data["gauge_hessian_residual"],
        "cd_projective_source": data["gauge_source_residual"],
        "trace_gauge_embedding": gauge.T*q["embedding"],
        "quotient_determinant": sp.factor(q["determinant"]+2**52*P**72*(8*P**2-1)**3),
        "full_euler": eliminate()["full_euler_residual"],
        "density_reconstruction": coefficients()["density_reconstruction_residual"],
        "palatini_euler": palatini_control()["full_euler_residual"],
        "palatini_density": palatini_control()["density_residual"],
        "palatini_Q2": palatini_control()["Q2_residual"],
        "exceptional_kernel": exceptional()["kernel_residual"],
        "exceptional_source_compatibility": exceptional()["source_compatibility_residual"],
    }
    result.update({f"coefficient_{key}": value
                   for key, value in source_formula_comparison()["residuals"].items()})
    for i, residual in enumerate(inverse_bound()["inverse_residuals"]):
        result[f"block_inverse_{i}"] = residual
    return result


@cache
def calibration():
    """Exact replay metadata, no kinetic or UV conclusion."""
    return {"connection_components": 64, "projective_nullity": 4,
            "quotient_dimension": 60,
            "block_sizes": tuple(len(block) for block in quotient()["blocks"]),
            "block_determinants": quotient()["block_determinants"],
            "determinant": quotient()["determinant"],
            "CD_p_squared_lower": sp.Rational(9, 40),
            "CD_p_squared_upper": sp.Rational(11, 40),
            "CD_exceptional_factor_lower": sp.Rational(4, 5),
            "inverse_norm_upper": inverse_bound()["norm_upper"],
            "exceptional_quotient_rank": exceptional()["quotient_rank"],
            "printed_Q2_defect": source_formula_comparison()["printed_Q2_defect"],
            "algebraic_only": True, "UV_or_kinetic_health_claim": False}
