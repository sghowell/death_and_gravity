"""Universal Laplace-type input, explicit Hodge traces and pole convention."""
from functools import cache
from itertools import combinations

import sympy as sp

R, Ric2, Riem2, boxR = sp.symbols("R Ricci_squared Riemann_squared box_R", real=True)
m2 = sp.Symbol("m0_squared", positive=True)


@cache
def clock_reduction():
    from p8_affine import dictionary
    from p8_aligned_quantum import kernel
    clock_p = sp.simplify(dictionary.lift()["p"].subs(dictionary.x, -1))
    masses = kernel.masses()
    a_clock, b_clock = [sp.factor(masses[name].subs(kernel.geometry.P, clock_p)) for name in ("a", "b")]
    return {"actual_clock_p": clock_p-sp.Rational(1, 2),
            "actual_clock_temporal_mass_one": a_clock-1,
            "actual_clock_spatial_mass_one": b_clock-1,
            "actual_clock_temporal_mass_derivative_zero": sp.diff(a_clock, dictionary.u),
            "actual_clock_spatial_mass_derivative_zero": sp.diff(b_clock, dictionary.u)}


@cache
def cochain_control():
    """Exact finite Hodge complex, including a scalar constant mode.

    The finite contact determinant is not discarded here. Setting its
    continuum trace to zero is specific to the stated dimensional scheme.
    """
    gradient = sp.Matrix([[-1, 1, 0], [0, -1, 1], [0, 0, 0], [0, 0, 0]])
    curl = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 2]])
    delta1 = gradient*gradient.T+curl.T*curl
    delta0 = gradient.T*gradient
    K = curl.T*curl+m2*sp.eye(4)
    factor = sp.eye(4)+gradient*gradient.T/m2
    return {"gradient": gradient, "curl": curl, "K": K, "Delta1": delta1, "Delta0": delta0,
            "d_squared_zero_cochain": curl*gradient,
            "Hodge_Proca_operator_factorization": (K*factor-delta1-m2*sp.eye(4)).applyfunc(sp.expand),
            "determinant_ratio_including_contact_factor": sp.factor(K.det()*(delta0+m2*sp.eye(3)).det()
                                        -m2**3*(delta1+m2*sp.eye(4)).det()),
            "scalar_constant_mode_is_retained": delta0*sp.ones(3, 1)}


@cache
def coefficients():
    # D=-(nabla²+E); Hodge one-forms have E=-Ricci and
    # tr(Omega_mu_nu Omega^mu_nu)=-Riemann², not +Riemann².
    universal = (5*R**2-2*Ric2+2*Riem2+12*boxR)/360
    scalar = {0: sp.Integer(1), 2: R/6, 4: universal}
    vector = {0: sp.Integer(4), 2: 4*R/6-R,
              4: 4*universal+(-60*R**2+180*Ric2-60*boxR-30*Riem2)/360}
    proca = {order: sp.factor(vector[order]-scalar[order]) for order in (0, 2, 4)}
    pole_weight = proca[4]-m2*proca[2]+m2**2*proca[0]/2
    return {"scalar": scalar, "vector": vector, "proca": proca, "pole_weight": sp.expand(pole_weight)}


@cache
def explicit_tensor_traces():
    pairs = list(combinations(range(4), 2))
    values = iter(sp.symbols("curvature0:21", real=True))
    blocks = sp.zeros(6)
    for i in range(6):
        for j in range(i, 6):
            blocks[i, j] = blocks[j, i] = next(values)
    # The unique independent algebraic Bianchi relation in four dimensions.
    blocks[2, 3] = blocks[3, 2] = blocks[1, 4]-blocks[0, 5]
    tensor = sp.MutableDenseNDimArray.zeros(4, 4, 4, 4)
    for i, (a, b) in enumerate(pairs):
        for j, (c, d) in enumerate(pairs):
            value = blocks[i, j]
            tensor[a, b, c, d] = tensor[b, a, d, c] = value
            tensor[b, a, c, d] = tensor[a, b, d, c] = -value
    Ricci = sp.Matrix(4, 4, lambda a, b: sum(tensor[c, a, c, b] for c in range(4)))
    E = -Ricci
    omega = sum(sp.trace(sp.Matrix(4, 4, lambda a, b, mu=mu, nu=nu: tensor[a, b, mu, nu])**2)
                for mu in range(4) for nu in range(4))
    norm_Riemann = sum(tensor[a, b, c, d]**2 for a in range(4) for b in range(4)
                       for c in range(4) for d in range(4))
    norm_Ricci = sum(value**2 for value in Ricci)
    return {"curvature_parameter_count": len(set().union(*(value.free_symbols for value in blocks))),
            "Ricci_symmetry": (Ricci-Ricci.T).applyfunc(sp.expand),
            "trace_E_is_minus_R": sp.expand(sp.trace(E)+sp.trace(Ricci)),
            "trace_E_squared_is_Ricci_norm": sp.expand(sp.trace(E**2)-norm_Ricci),
            "trace_connection_curvature_squared_is_minus_Riemann_norm": sp.expand(omega+norm_Riemann)}


@cache
def checks():
    d = coefficients()
    expected = {0: 3, 2: -R/2, 4: -R**2/8+sp.Rational(29, 60)*Ric2-Riem2/15-boxR/15}
    out = {"Hodge_vector_minus_scalar_coefficient_"+str(order): sp.expand(d["proca"][order]-value)
           for order, value in expected.items()}
    t = sp.Symbol("proper_time", positive=True)
    expanded = sp.series(sp.exp(-m2*t)*sum(d["proca"][order]*t**(order//2) for order in (0, 2, 4)), t, 0, 3).removeO()
    out["mass_weighted_heat_coefficient"] = sp.expand(expanded.coeff(t, 2)-d["pole_weight"])
    out["flat_S6_47_pole_weight"] = sp.expand(d["pole_weight"].subs({R: 0, Ric2: 0, Riem2: 0, boxR: 0})-3*m2**2/2)
    GB = sp.Symbol("Gauss_Bonnet", real=True)
    expected_gb = -GB/15+sp.Rational(13, 60)*Ric2-sp.Rational(7, 120)*R**2-boxR/15
    out["Gauss_Bonnet_basis_with_total_derivative"] = sp.expand(d["proca"][4].subs(Riem2, GB+4*Ric2-R**2)-expected_gb)
    out.update({name: value for name, value in explicit_tensor_traces().items() if name != "curvature_parameter_count"})
    out.update(clock_reduction())
    out.update({name: value for name, value in cochain_control().items()
                if name in ("d_squared_zero_cochain", "Hodge_Proca_operator_factorization",
                            "determinant_ratio_including_contact_factor", "scalar_constant_mode_is_retained")})
    return out
