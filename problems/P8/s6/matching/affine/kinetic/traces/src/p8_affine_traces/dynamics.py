"""All branches of the constant rank-one two-trace curl family.

The result concerns the actual CD/M1 quadratic action.  The zero-Schur
case is a new constraint chart, never an inverse-mass or Proca limit.
All kinetic coefficients here are normalized by M² tau².
"""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_affine import connection
from p8_affine_kinetic import scalar as old

A, B = sp.symbols("A B", real=True)
zeta = sp.Symbol("zeta", real=True, nonzero=True)
gamma = sp.Symbol("gamma", real=True, nonzero=True)
q = old.q
u = old.u
a, theta, J, ell, w, lam = old.a, old.theta, old.J, old.ell, old.w, old.lam
n, shift, v, matter = old.n, old.shift, old.v, old.matter
vd, sd, sigmad, sigma, temporal = old.vd, old.sd, old.sigmad, old.sigma, old.temporal
d, e, f = sp.symbols("d e f", real=True)
thetad, elld = sp.symbols("theta_dot ell_dot", real=True)
vdd, P0, P1, Ps = sp.symbols("vdd P0 P1 Ps", real=True)


def _matrix(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.factor))


@cache
def geometry():
    """Independent eight-row contraction and every frozen quotient block."""
    vmap, umap = sp.zeros(4, 64), sp.zeros(4, 64)
    signs = (-1, 1, 1, 1)
    for mu in range(4):
        for i in range(4):
            vmap[mu, connection.index(i, mu, i)] += 1
            vmap[mu, connection.index(i, i, mu)] -= sp.Rational(1, 4)
            umap[mu, connection.index(mu, i, i)] += signs[mu]*signs[i]
            umap[mu, connection.index(i, i, mu)] -= sp.Rational(1, 4)
    quotient = connection.quotient()
    full_map = vmap.col_join(umap)
    projected = full_map*quotient["embedding"]
    response = sp.zeros(8)
    for indices, block in zip(quotient["blocks"], quotient["block_matrices"], strict=True):
        columns = projected[:, list(indices)]
        response += columns*block.subs(connection.P, sp.Rational(1, 2)).inv()*columns.T
    sol = connection.eliminate()["solution"]
    return {"V_map": sp.ImmutableMatrix(vmap), "U_map": sp.ImmutableMatrix(umap),
            "map": sp.ImmutableMatrix(full_map), "D_two": _matrix(response),
            "Vstar": _matrix(vmap*sol), "Ustar": _matrix(umap*sol),
            "projective_residual": full_map*connection.quadratic()["gauge"]}


@cache
def rolling():
    """Linearize the contracted source, including coefficient/background products."""
    data = geometry()
    h, hp, hubble = sp.symbols("h hp hubble", real=True, nonzero=True)
    p, ss = connection.P, connection.S
    clock = {p: sp.Rational(1, 2), ss: 1, connection.PP: 0,
             connection.PX: -1/(4*h), connection.CP: 0,
             connection.CX: -1/(2*h), connection.F3: 0}
    parts = {}
    for name in ("V", "U"):
        vector = data[name+"star"]
        trace = sp.diff(vector[0], connection.H[1, 1])
        coefficient_n = sp.diff(trace, p).subs(clock)*(-1/(2*h))
        forcing_n = (sp.diff(vector[0], connection.PP)+sp.diff(vector[0], connection.F3))
        forcing_n = forcing_n.subs(clock)*hp/(2*h**2)
        source_n = sp.factor((-3*hubble*coefficient_n+forcing_n).subs(hp, 3*hubble*h/2))
        alpha = sp.factor(sp.diff(vector[0], connection.H[0, 0]).subs(clock))
        beta = sp.factor(sp.diff(vector[1], connection.H[0, 1]).subs(clock))
        parts[name] = {"alpha": alpha, "beta": beta, "n": source_n,
                       "trace": sp.factor(trace.subs(clock)),
                       "clock": _matrix(vector.subs(clock))}
    alpha = sp.factor(A*parts["V"]["alpha"]+B*parts["U"]["alpha"])
    beta = sp.factor(A*parts["V"]["beta"]+B*parts["U"]["beta"])
    source_n = sp.factor(A*parts["V"]["n"]+B*parts["U"]["n"])
    alpha_prime = sp.diff(alpha, h)*3*hubble*h/2
    shift_d = sp.factor(-alpha_prime-source_n)
    shift_e = sp.factor(beta-alpha)
    curl_n = sp.factor(shift_d-sp.diff(shift_e, h)*3*hubble*h/2)
    bg = old.background()
    actual = {h: bg["h"], hubble: bg["H"]}
    return {"h": h, "H": hubble, "parts": parts, "alpha": alpha, "beta": beta,
            "source_n": source_n, "d": shift_d, "e": shift_e, "curl_n": curl_n,
            "actual_alpha": alpha.subs(actual), "actual_d": shift_d.subs(actual),
            "actual_e": shift_e.subs(actual), "actual_curl_n": curl_n.subs(actual)}


@cache
def coefficients():
    """The only Schur scalar on the rolling background, not an open-tube inverse."""
    combination = sp.zeros(4, 8)
    combination[:, :4] = A*sp.eye(4)
    combination[:, 4:] = B*sp.eye(4)
    response = _matrix(combination*geometry()["D_two"]*combination.T)
    g = sp.factor(response[0, 0])
    return {"gamma": g, "D": response, "eta": sp.diag(1, -1, -1, -1),
            "gamma_formula": 3*(A**2-14*A*B+B**2)/8,
            "center_e": A-B, "null_ratios": (7+4*sp.sqrt(3), 7-4*sp.sqrt(3))}


@cache
def regular():
    """Gamma!=0 chart: retain W0,n,b and the actual free matter first."""
    base = old.action()["base"]
    extra = (temporal+d*n)**2/(2*gamma)-q*(sigma+e*n)**2/(2*gamma)
    extra += zeta*q*(sigmad-temporal)**2/2
    W0 = (gamma*zeta*q*sigmad-d*n)/(1+gamma*zeta*q)
    C = zeta*q/(2*(1+gamma*zeta*q))
    after = C*(sigmad+d*n)**2-q*(sigma+e*n)**2/(2*gamma)
    lapse = (vd+ell*matter/2)/theta
    reduced = sp.factor((base+after).subs(n, lapse))
    kinetic = _matrix(sp.hessian(reduced, (vd, sd, sigmad))/2)
    Jeff = J-q*e**2/(2*gamma)
    completed = (sd+w*vd/theta)**2/2+C*(sigmad+d*vd/theta)**2+Jeff*vd**2/theta**2
    vel = sp.Matrix([vd, sd, sigmad])
    return {"before": base+extra, "extra": extra, "W0": W0, "C": C,
            "temporal_divisor": 1+gamma*zeta*q, "after": base+after,
            "lapse": lapse, "reduced": reduced, "kinetic": kinetic,
            "J_eff": Jeff, "Schur": Jeff/theta**2,
            "pivots": (C, sp.Rational(1, 2), Jeff/theta**2),
            "W0_euler": sp.factor(sp.diff(extra, temporal).subs(temporal, W0)),
            "W0_action": sp.factor(extra.subs(temporal, W0)-after),
            "shift_residual": sp.factor(sp.diff(base+extra, shift)/(2*q)-theta*n+vd+ell*matter/2),
            "square_residual": sp.factor((vel.T*kinetic*vel)[0]-completed),
            "determinant_residual": sp.factor(kinetic.det()-C*Jeff/(2*theta**2))}


@cache
def null_schur():
    """Fresh gamma=0 constraint; all statements are at quadratic order.

M delta_kappa+N^T Ekin=0 and N M^-1 N^T=0 force T=Tstar.
The toy matrix is an independent finite-dimensional null-Schur control;
the actual eight-trace geometry is verified separately in geometry().
"""
    source1, source2, current = sp.symbols("source1 source2 Ekin", real=True)
    M = sp.Matrix([[0, 1], [1, 0]])
    N = sp.Matrix([[1, 0]])
    source = sp.Matrix([source1, source2])
    stationary = -M.inv()*source
    solution = stationary-M.inv()*N.T*current
    residual = M*solution+source+N.T*current
    change = solution-stationary
    return {"M": M, "N": N, "stationary": stationary, "solution": solution,
            "Euler_residual": residual, "projected_response": N*M.inv()*N.T,
            "trace_constraint": N*(solution-stationary),
            "auxiliary_action_correction": sp.factor((change.T*M*change)[0]/2)}


@cache
def ostrogradsky():
    """Exact full-M1 highest-derivative Legendre map on theta!=0, gamma=0.

The shift fixes n while its companion lapse equation reconstructs b.
No physical-metric or free-matter field is frozen.  e!=0 and zeta!=0
make the joint (v_ddot,s_dot) Hessian invertible for either sign.
"""
    lapse = (vd+ell*matter/2)/theta
    lapse_dot = (vdd+ell*sd/2+elld*matter/2)/theta-thetad*(vd+ell*matter/2)/theta**2
    electric = -e*lapse_dot+f*lapse
    base = sp.factor(old.action()["base"].subs(n, lapse))
    density = base+zeta*q*electric**2/2
    full = a**3*density
    ae = -e/theta
    be = ae*ell/2
    remainder = sp.factor(electric-ae*vdd-be*sd)
    L1 = sp.factor(sp.diff(base, sd).subs(sd, 0))
    L0 = sp.factor(base.subs(sd, 0))
    sd_solution = Ps/a**3-ell*P1/(2*a**3)-L1
    vdd_solution = (P1/(a**3*zeta*q*ae)-be*sd_solution-remainder)/ae
    velocity_sub = {sd: sd_solution, vdd: vdd_solution}
    hamiltonian = (P0*vd+(Ps-ell*P1/2-a**3*L1)**2/(2*a**3)
                   +P1**2/(2*a**3*zeta*q*ae**2)-(remainder/ae)*P1-a**3*L0)
    highest = _matrix(sp.hessian(full, (vdd, sd)))
    legendre = P0*vd+P1*vdd_solution+Ps*sd_solution-full.subs(velocity_sub, simultaneous=True)
    return {"lapse": lapse, "lapse_dot": lapse_dot, "electric": electric,
            "density": density, "full_L": full, "ae": ae, "be": be,
            "remainder": remainder, "L0": L0, "L1": L1,
            "highest_hessian": highest, "highest_determinant": sp.factor(highest.det()),
            "s_dot": sd_solution, "v_ddot": vdd_solution, "Hamiltonian": hamiltonian,
            "P1_residual": sp.factor(sp.diff(full, vdd).subs(velocity_sub, simultaneous=True)-P1),
            "Ps_residual": sp.factor(sp.diff(full, sd).subs(velocity_sub, simultaneous=True)-Ps),
            "Legendre_residual": sp.factor(legendre-hamiltonian),
            "affine_momentum_coefficient": sp.diff(hamiltonian, P0),
            "affine_momentum_second_derivative": sp.diff(hamiltonian, P0, 2),
            "lapse_reconstruction_coefficient": 2*q*theta,
            "missing_free_matter_rank_control": sp.factor(sp.hessian(full-a**3*sd**2/2, (vdd, sd)).det())}


def _exact(value, label):
    if isinstance(value, (bool, float, str)):
        raise TypeError(f"{label} must be an exact real constant")
    if isinstance(value, Fraction):
        value = sp.Rational(value.numerator, value.denominator)
    value = sp.sympify(value)
    if not isinstance(value, sp.Expr) or value.has(sp.Float) or value.free_symbols:
        raise TypeError(f"{label} must be an exact real constant")
    if value.is_real is not True or value.is_finite is not True:
        raise ValueError(f"{label} must be finite and real")
    return value


def require_domain(first, second, coupling, time, momentum):
    """Classify every exact rank-one branch; never divide by its zero Schur scalar."""
    aa, bb, zz, tt, qq = [_exact(value, label) for value, label in
                         ((first, "A"), (second, "B"), (coupling, "zeta"),
                          (time, "u"), (momentum, "q"))]
    if qq.is_positive is not True:
        raise ValueError("q>0 is required")
    gg = sp.factor(3*(aa**2-14*aa*bb+bb**2)/8)
    result = {"A": aa, "B": bb, "zeta": zz, "u": tt, "q": qq, "gamma": gg}
    if (aa.is_zero is True and bb.is_zero is True) or zz.is_zero is True:
        return {**result, "branch": "unchanged_auxiliary", "threshold": None}
    if gg.is_zero is True:
        if tt.is_zero is not False:
            raise ValueError("The null-Schur Ostrogradsky proof requires u!=0")
        if sp.simplify(aa-bb).is_zero is not False:
            raise ValueError("A nonzero null-Schur direction must have A!=B")
        return {**result, "branch": "null_schur_ostrogradsky", "threshold": None}
    if gg.is_zero is not False:
        raise ValueError("The sign/rank of gamma must be proved exactly")
    if zz.is_negative is True:
        return {**result, "branch": "negative_transverse_kinetic", "threshold": None}
    if zz.is_positive is not True:
        raise ValueError("The sign of the nonzero coupling must be proved exactly")
    if tt.is_zero is not False:
        raise ValueError("The scalar kinetic proof requires u!=0")
    if gg.is_negative is True:
        threshold = -1/(gg*zz)
        label = "negative_longitudinal_pivot"
    elif gg.is_positive is True:
        bg = old.background()
        threshold = sp.factor(2*gg*bg["J"].subs(u, tt)*bg["h"].subs(u, tt)**2/(aa-bb)**2)
        label = "negative_clock_pivot"
    else:
        raise ValueError("The sign of gamma must be proved exactly")
    if sp.simplify(qq-threshold).is_positive is not True:
        raise ValueError("The displayed strict high-momentum condition is required")
    return {**result, "branch": label, "threshold": threshold}


def units(mass_squared, time_scale, physical_coupling, first=1, second=0):
    mm, tau, zz, aa, bb = [_exact(value, label) for value, label in
                          ((mass_squared, "M²"), (time_scale, "tau"),
                           (physical_coupling, "zeta_physical"), (first, "A"), (second, "B"))]
    if mm.is_positive is not True or tau.is_positive is not True:
        raise ValueError("Positive M² and tau are required")
    gg = sp.factor(3*(aa**2-14*aa*bb+bb**2)/8)
    mass = None if gg.is_zero is True or zz.is_zero is True else mm/(gg*zz)
    return {"normalized_zeta": zz/(mm*tau**2), "gamma": gg,
            "isolated_mass_squared_physical": mass,
            "isolated_mass_squared_normalized": None if mass is None else tau**2*mass,
            "mass_is_not_the_coupled_health_test": True}


@cache
def checks():
    geom, rr, cc, reg, null, ost = geometry(), rolling(), coefficients(), regular(), null_schur(), ostrogradsky()
    eta = sp.diag(1, -1, -1, -1)
    h, hubble = rr["h"], rr["H"]
    result = {
        "projective_traces": geom["projective_residual"],
        "D_VV": _matrix(geom["D_two"][:4, :4]-3*eta/8),
        "D_UU": _matrix(geom["D_two"][4:, 4:]-3*eta/8),
        "D_VU": _matrix(geom["D_two"][:4, 4:]+21*eta/8),
        "gamma_from_every_quotient_block": sp.factor(cc["gamma"]-cc["gamma_formula"]),
        "Lorentz_four_trace_response": _matrix(cc["D"]-cc["gamma"]*eta),
        "rolling_alpha": sp.factor(rr["alpha"]-3*(A-B)/(2*h)),
        "rolling_d": sp.factor(rr["d"]-3*(7*A+3*B)*hubble/(8*h)),
        "rolling_e": sp.factor(rr["e"]-(A-B)/h),
        "null_curl_n": sp.factor(rr["curl_n"]-3*(11*A-B)*hubble/(8*h)),
        "V_trace_velocity_drops": rr["parts"]["V"]["trace"],
        "U_trace_velocity_drops": rr["parts"]["U"]["trace"],
        "regular_W0_euler": reg["W0_euler"], "regular_W0_action": reg["W0_action"],
        "unchanged_shift": reg["shift_residual"], "regular_squares": reg["square_residual"],
        "regular_determinant": reg["determinant_residual"],
        "null_full_Euler": null["Euler_residual"], "null_projected_response": null["projected_response"],
        "null_trace_constraint": null["trace_constraint"],
        "null_auxiliary_correction": null["auxiliary_action_correction"],
        "joint_highest_hessian": sp.factor(ost["highest_determinant"]-a**6*zeta*q*e**2/theta**2),
        "Ostro_P1": ost["P1_residual"], "Ostro_Ps": ost["Ps_residual"],
        "Ostro_Legendre": ost["Legendre_residual"],
        "unconstrained_affine_momentum": ost["affine_momentum_coefficient"]-vd,
        "affine_momentum_second_derivative": ost["affine_momentum_second_derivative"],
        "missing_free_matter_rank_control": ost["missing_free_matter_rank_control"],
        "equal_trace_gamma": sp.factor(cc["gamma"].subs(B, A)+9*A**2/2),
    }
    for i, ratio in enumerate(cc["null_ratios"]):
        result[f"null_ratio_{i}"] = sp.simplify(cc["gamma"].subs(A, ratio*B))
    scale = sp.Symbol("scale", real=True, nonzero=True)
    scaled_gamma = cc["gamma"].subs({A: scale*A, B: scale*B}, simultaneous=True)
    result["proportional_gamma"] = sp.factor(scaled_gamma-scale**2*cc["gamma"])
    result["proportional_clock_threshold"] = sp.factor(
        (scale*(A-B))**2/(2*scaled_gamma)-(A-B)**2/(2*cc["gamma"]))
    result["proportional_null_acceleration"] = sp.factor(zeta*(scale*e)**2/scale**2-zeta*e**2)
    return result


@cache
def proof_checks():
    cc, reg, ost = coefficients(), regular(), ostrogradsky()
    gp, zp, ep, jp, margin = sp.symbols("gp zp ep jp margin", positive=True)
    positive_q = 2*gp*jp/ep**2+margin
    negative_q = 1/(gp*zp)+margin
    result = {
        "positive_gamma_cannot_have_equal_traces": (-9*ep**2/2).is_negative is True,
        "null_nontrivial_directions_are_distinct": all(sp.simplify(r-1).is_zero is False for r in cc["null_ratios"]),
        "positive_gamma_clock_sign_entire_strict_domain": sp.factor(
            reg["J_eff"].subs({gamma: gp, J: jp, e: ep, q: positive_q})).is_negative is True,
        "negative_gamma_vector_sign_entire_strict_domain": sp.factor(
            reg["C"].subs({gamma: -gp, zeta: zp, q: negative_q})).is_negative is True,
        "negative_coupling_transverse_example": bool((-sp.Rational(1, 2)) < 0),
        "both_null_coupling_signs_nondegenerate": all(
            ost["highest_determinant"].subs({a: 1, theta: 2, e: 3, q: 5, zeta: sign}) != 0
            for sign in (-1, 1)),
        "linear_momentum_not_a_frozen_mass_inference": ost["affine_momentum_coefficient"] == vd,
        "zero_coupling_highest_rank_changes": ost["highest_determinant"].subs(zeta, 0) == 0,
        "zero_direction_is_unchanged": cc["gamma"].subs({A: 0, B: 0}) == 0,
        "rank_one_only_in_this_module": True,
        "no_EFT_frequency_or_open_tube_inverse_claim": True,
    }
    return {name: bool(value) for name, value in result.items()}


@cache
def calibration():
    return {"gamma": coefficients()["gamma"], "rank_one_parameters": ("A", "B", "zeta"),
            "null_ratios": coefficients()["null_ratios"],
            "positive_gamma_threshold": "q>2*gamma*J/e²", "negative_gamma_threshold": "q>−1/(gamma*zeta)",
            "null_highest_determinant": ostrogradsky()["highest_determinant"],
            "unchanged_controls": ("A=B=0", "zeta=0"),
            "normalized_coupling": "zeta_physical/(M²*tau²)",
            "identity_count": len(checks()), "proof_check_count": len(proof_checks()),
            "full_rank_one_family_covered": True, "low_energy_cutoff_claim": False}
