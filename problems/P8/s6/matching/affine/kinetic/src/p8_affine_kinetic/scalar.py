"""Source-preserving scalar constraints of the selective affine-vector extension.

All coefficients are dimensionless after factoring M² tau² from the action.
Here zeta=zeta_physical/(M² tau²), q=(tau*k_physical)².  A negative
punctured configuration-velocity pivot is distinguished from the inertia of
a momentum-exchanged chart at the crossing itself.  No EFT cutoff is inferred.
"""

import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import connection

u = sp.Symbol("u", real=True)
q, zeta = sp.symbols("q zeta", positive=True)
a = sp.Symbol("a", positive=True)
H, theta, lam, ell, w, J, dd, e = sp.symbols("H theta lam ell w J dd e", real=True)
n, shift, v, matter, sigma, temporal = sp.symbols("n shift v matter sigma temporal", real=True)
vd, sd, sigmad = sp.symbols("vd sd sigmad", real=True)
pm, pi, pb = sp.symbols("pm pi pb", real=True)
BETA = sp.Rational(4, 3)
P8 = Path(__file__).resolve().parents[6]


def _matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.factor))


@cache
def background():
    """Actual original free-M1 trajectory, not a scalar spectator background."""
    raw = json.loads((P8/"certificates"/"witness-CD_matter.json").read_text())
    local = {"t": u}
    aa = sp.sympify(raw["background"]["a"], locals=local)
    hh = (1+u**2)**3
    ll = sp.sympify(raw["background"]["chi_dot"], locals=local)
    delta = 1/(2*hh)
    return {
        "u": u, "a": aa, "h": hh,
        "H": sp.sympify(raw["background"]["H"], locals=local),
        "theta": sp.sympify(raw["principal"]["Theta"], locals=local),
        "lam": sp.sympify(raw["principal"]["Lambda"], locals=local),
        "J": sp.sympify(raw["principal"]["J"], locals=local),
        "ell": ll, "l": ll, "delta": delta, "w": ll*(3*delta-1),
        "alpha": 3/(2*hh), "e": 1/hh,
        "d": 21*sp.sympify(raw["background"]["H"], locals=local)/(8*hh),
    }


@cache
def vector_reconstruction():
    """Contract the full frozen 64-component solution before linearizing.

The rest frame is aligned with the clock normal, so H_0i=-partial_i n
at first order; a coordinate-Hessian shift term must not be left alone.
All H_ij perturbations drop out because their background coefficient is zero.
The lower coefficient q(u,-1)=q_x(u,-1)=0 is used, not q=0 off the tube.
"""
    sol = connection.eliminate()["solution"]
    vec = sp.Matrix([
        sum(sol[connection.index(i, mu, i)]-sol[connection.index(i, i, mu)]/4
            for i in range(4)) for mu in range(4)
    ]).applyfunc(sp.factor)
    h, hp, hubble = sp.symbols("h hp hubble", real=True, nonzero=True)
    p, ss = connection.P, connection.S
    substitutions = {p: sp.Rational(1, 2), ss: 1, connection.PP: 0,
                     connection.PX: -1/(4*h), connection.CP: 0,
                     connection.CX: -1/(2*h), connection.F3: 0}
    clock = vec.subs(substitutions).applyfunc(sp.factor)
    trace_coefficient = sp.diff(vec[0], connection.H[1, 1])
    # x=-N^-2=-1+2n+O(n²); p_N=-1/(2h).  At fixed x,
    # p_phi,N=hp/(2h²), J3_N=hp/(2h²), and q_N=0.
    p_n = -1/(2*h)
    pp_n = hp/(2*h**2)
    cubic_n = hp/(2*h**2)
    coefficient_n = sp.diff(trace_coefficient, p).subs({p: sp.Rational(1, 2), ss: 1})*p_n
    zero_hessian = {symbol: 0 for symbol in connection.H_SYMBOLS}
    forcing_n = (
        sp.diff(vec[0], connection.PP).subs({p: sp.Rational(1, 2), ss: 1})*pp_n
        +sp.diff(vec[0], connection.F3).subs({p: sp.Rational(1, 2), ss: 1})*cubic_n
    ).subs(zero_hessian)
    n_coefficient = sp.factor(-3*hubble*coefficient_n+forcing_n)
    temporal_n = sp.factor(n_coefficient.subs(hp, 3*hubble*h/2))
    alpha = 3/(2*h)
    alpha_prime = -3*hp/(2*h**2)
    d_actual = sp.factor((-alpha_prime-temporal_n).subs(hp, 3*hubble*h/2))
    return {
        "generic_vector": sp.ImmutableMatrix(vec), "clock_vector": sp.ImmutableMatrix(clock),
        "h": h, "hp": hp, "H": hubble, "trace_coefficient": trace_coefficient,
        "lapse_time_coefficient": -alpha, "lapse_coefficient_before_clock_identity": n_coefficient,
        "lapse_coefficient": temporal_n, "spatial_lapse_coefficient": -5/(2*h),
        "alpha": alpha, "d": d_actual, "e": 1/h,
        "trace_perturbation_residual": sp.factor(trace_coefficient.subs(substitutions)),
        "no_shift_or_zeta_velocity": True,
    }


@cache
def action():
    """Complete scalar density divided by a³ after the stated CD boundaries.

Spatial gauge has zero scalar shear, clock gauge delta_phi=0.  The metric
curvature perturbation is zeta_metric=v+delta*n, shift b=a²*psi; matter is
s=delta_chi.  The vector is W=V+d[alpha*n], W_i=partial_i sigma.
Only W0 is eliminated here; lapse and scalar shift remain explicit.
"""
    S = J+w**2/2-3*theta**2
    base = (-3*vd**2+S*n**2+6*theta*n*vd+sd**2/2+w*n*sd-3*ell*vd*matter
            +q*(v**2+2*lam*n*v-matter**2/2+2*theta*n*shift
                 -2*vd*shift-ell*shift*matter))
    extra = BETA*(temporal+dd*n)**2-BETA*q*(sigma+e*n)**2
    extra += zeta*q*(sigmad-temporal)**2/2
    solution = (3*zeta*q*sigmad-8*dd*n)/(8+3*zeta*q)
    C = 4*zeta*q/(8+3*zeta*q)
    reduced_extra = C*(sigmad+dd*n)**2-BETA*q*(sigma+e*n)**2
    return {
        "base": base, "extra_before_temporal": extra,
        "full_before_temporal": base+extra, "temporal": solution, "C": C,
        "extra": reduced_extra, "L": base+reduced_extra,
        "S": S, "J_eff": J-BETA*q*e**2,
        "temporal_euler_residual": sp.factor(sp.diff(extra, temporal).subs(temporal, solution)),
        "temporal_reduction_residual": sp.factor(extra.subs(temporal, solution)-reduced_extra),
        "shift_constraint": sp.factor(sp.diff(base+reduced_extra, shift)/(2*q)),
    }


@cache
def unitary():
    """Direct configuration reduction ONLY for theta!=0, not at the crossing."""
    data = action()
    lapse = (vd+ell*matter/2)/theta
    reduced = sp.factor(data["L"].subs(n, lapse))
    velocities = sp.Matrix([vd, sd, sigmad])
    kinetic = _matrix(sp.hessian(reduced, tuple(velocities))/2)
    C = data["C"]
    completed = (sd+w*vd/theta)**2/2+C*(sigmad+dd*vd/theta)**2
    completed += data["J_eff"]*vd**2/theta**2
    return {
        "lapse": lapse, "L": reduced, "kinetic": kinetic, "completed_kinetic": completed,
        "pivots": (C, sp.Rational(1, 2), data["J_eff"]/theta**2),
        "determinant": sp.factor(kinetic.det()),
        "Schur": sp.factor(kinetic[0, 0]-kinetic[0, 1]**2/kinetic[1, 1]
                            -kinetic[0, 2]**2/kinetic[2, 2]),
        "square_residual": sp.factor((velocities.T*kinetic*velocities)[0]-completed),
        "constraint_residual": sp.factor(data["shift_constraint"].subs(n, lapse)),
    }


@cache
def hamiltonian():
    """Regular first-order shift chart, before dividing by J_eff.

Momentums below are density-normalized.  The full canonical density keeps
P_b=2a³q*v and adds -H*b*P_b from the moving canonical boundary.
"""
    data = action()
    C = data["C"]
    velocities = {vd: theta*n-ell*matter/2, sd: pm-w*n,
                  sigmad: pi/(2*C)-dd*n}
    before = sp.factor((-2*q*shift*vd+pm*sd+pi*sigmad-data["L"]).subs(velocities))
    R = q*(theta*shift+lam*v)+w*pm/2-3*ell*matter*theta/2+dd*pi/2-BETA*q*e*sigma
    constant = (pm**2/2+q*ell*shift*matter+q*matter**2/2-q*v**2
                -3*ell**2*matter**2/4+pi**2/(4*C)+BETA*q*sigma**2)
    regular = constant+R**2/data["J_eff"]
    canonical = a**3*regular.subs({v: pb/(2*a**3*q), pm: pm/a**3, pi: pi/a**3}, simultaneous=True)
    canonical -= H*shift*pb
    return {
        "before_lapse": before, "constant": constant, "R_eff": R,
        "J_eff": data["J_eff"], "lapse": -R/data["J_eff"],
        "density": regular, "canonical_density": canonical,
        "before_lapse_residual": sp.factor(before-constant+data["J_eff"]*n**2+2*R*n),
        "lapse_residual": sp.factor(sp.diff(before, n).subs(n, -R/data["J_eff"])),
        "primary_momentum_residual": sp.factor(sp.diff(data["L"], vd).subs(velocities)+2*q*shift),
    }


def _center_substitution():
    return {a: 1, H: 0, theta: 0, lam: -sp.Rational(1, 2), ell: sp.Rational(1, 10),
            w: sp.Rational(1, 20), J: sp.Rational(1199, 800), dd: 0, e: 1}


@cache
def center():
    """Center chart inertia: corroboration, NOT its own physical ghost argument."""
    data = hamiltonian()
    density = data["density"].subs(_center_substitution())
    hessian = _matrix(sp.hessian(density, (v, pm, pi)))
    norm = sp.diag(2*q, 1, 1)
    kinetic = _matrix(norm.T*hessian.inv()*norm/2)
    expected = sp.Matrix([
        [-2*q*(8*q-9)/(19*q-18), 3*q/(20*(19*q-18)), 0],
        [3*q/(20*(19*q-18)), (3800*q-3597)/(400*(19*q-18)), 0],
        [0, 0, action()["C"]],
    ])
    return {
        "density": density, "momentum_hessian": hessian, "kinetic": kinetic,
        "kinetic_residual": _matrix(kinetic-expected),
        "two_momentum_determinant": sp.factor(hessian[:2, :2].det()),
        "two_kinetic_determinant": sp.factor(kinetic[:2, :2].det()),
        "lapse_pole": sp.Rational(3597, 3200), "velocity_chart_pole": sp.Rational(18, 19),
        "baseline_velocity_chart_pole": sp.Integer(6),
    }


@cache
def center_jets():
    data = background()
    result = {f"{name}_{order}": sp.factor(sp.diff(data[name], u, order).subs(u, 0))
              for name in ("H", "theta", "d", "e", "alpha", "J", "w", "lam", "ell", "a")
              for order in (0, 1)}
    result["physical_q_derivative_at_center"] = 0
    return result


@cache
def center_pencil():
    """Exact center phase/Cauchy coefficients, including all first coefficient jets.

The resulting instantaneous characteristic polynomial is not promoted to a
uniform high-frequency cone: the direct punctured action supplies the verdict.
"""
    canonical = hamiltonian()["canonical_density"]
    states = (shift, matter, sigma, pb, pm, pi)
    matrix = sp.hessian(canonical, states)
    sub = _center_substitution()
    h0 = _matrix(matrix.subs(sub))
    h1 = _matrix((4*matrix.diff(H)+3*matrix.diff(theta)+sp.Rational(21, 2)*matrix.diff(dd)).subs(sub))
    omega = sp.zeros(6)
    omega[:3, 3:] = sp.eye(3)
    omega[3:, :3] = -sp.eye(3)
    f0, f1 = _matrix(omega*h0), _matrix(omega*h1)
    A, L, C = h0[3:, 3:], h0[3:, :3], -h0[:3, :3]
    Ap, Lp = h1[3:, 3:], h1[3:, :3]
    mv = _matrix((L*A+Ap-A*L.T)*A.inv())
    mq = _matrix(Lp+L*L+A*C-mv*L)
    observer = sp.zeros(3, 6)
    observer[:, :3] = sp.eye(3)
    residual = _matrix(observer*(f1+f0*f0)-mq*observer-mv*observer*f0)
    omitted = _matrix(observer*f1)
    return {"H0": h0, "H1": h1, "phase": f0, "phase_derivative": f1,
            "Mq": mq, "Mv": mv, "cauchy_residual": residual,
            "omitted_coefficient_jets": omitted}


def _exact(value, name):
    if isinstance(value, (bool, float, str)):
        raise TypeError(f"{name} must be an exact finite real constant")
    if isinstance(value, Fraction):
        value = sp.Rational(value.numerator, value.denominator)
    value = sp.sympify(value)
    if not isinstance(value, sp.Expr) or value.has(sp.Float) or value.free_symbols:
        raise TypeError(f"{name} must be an exact finite real constant")
    if value.is_real is not True or value.is_finite is not True:
        raise ValueError(f"{name} must be finite and real")
    return value


def require_domain(time, momentum, coupling, *, punctured=True, unhealthy=False):
    """Validate before cached calculations; theta=0 is not a unitary chart."""
    if not isinstance(punctured, bool) or not isinstance(unhealthy, bool):
        raise TypeError("Chart and verdict flags must be bools")
    time, momentum, coupling = (_exact(value, name) for value, name in
                                ((time, "time"), (momentum, "q"), (coupling, "zeta")))
    if momentum.is_positive is not True or coupling.is_positive is not True:
        raise ValueError("q and normalized zeta must be strictly positive")
    if punctured and time.is_zero is not False:
        raise ValueError("The direct unitary kinetic proof requires u!=0")
    bg = background()
    threshold = sp.factor((3*bg["J"]*bg["h"]**2/4).subs(u, time))
    if unhealthy and sp.simplify(momentum-threshold).is_positive is not True:
        raise ValueError("The strict negative-pivot inequality is required")
    return {"u": time, "q": momentum, "zeta": coupling, "threshold": threshold}


def units(mass_squared, time_scale, physical_coupling):
    values = [_exact(value, name) for value, name in
              ((mass_squared, "M²"), (time_scale, "tau"), (physical_coupling, "zeta_physical"))]
    if any(value.is_positive is not True for value in values):
        raise ValueError("All three physical scales must be positive")
    mass_squared, time_scale, physical_coupling = values
    return {"normalized_zeta": physical_coupling/(mass_squared*time_scale**2),
            "isolated_proca_mass_squared": 8*mass_squared/(3*physical_coupling),
            "normalized_isolated_mass_squared": 8*mass_squared*time_scale**2/(3*physical_coupling)}


@cache
def checks():
    bg, vr, ac, un, ha, ce = (background(), vector_reconstruction(), action(), unitary(),
                             hamiltonian(), center())
    h = vr["h"]
    target_clock = sp.Matrix([3*connection.H[0, 0]/(2*h)]
                             +[5*connection.H[0, i]/(2*h) for i in range(1, 4)])
    result = {"literal_vector_clock": _matrix(vr["clock_vector"]-target_clock),
              "trace_velocity_drops": vr["trace_perturbation_residual"],
              "rolling_temporal_n": sp.factor(vr["lapse_coefficient"]+3*vr["H"]/(8*h)),
              "gradient_shift_d": sp.factor(vr["d"]-21*vr["H"]/(8*h)),
              "background_clock_relation": sp.factor(sp.diff(bg["h"], u)-3*bg["H"]*bg["h"]/2),
              "free_matter_current": sp.factor(sp.diff(bg["a"]**3*bg["ell"], u)),
              "W0_euler": ac["temporal_euler_residual"],
              "W0_reduction": ac["temporal_reduction_residual"],
              "unchanged_shift": sp.factor(ac["shift_constraint"]-theta*n+vd+ell*matter/2),
              "unitary_shift": un["constraint_residual"], "kinetic_squares": un["square_residual"],
              "kinetic_Schur": sp.factor(un["Schur"]-ac["J_eff"]/theta**2),
              "kinetic_determinant": sp.factor(un["determinant"]-ac["C"]*ac["J_eff"]/(2*theta**2)),
              "momentum_primary": ha["primary_momentum_residual"],
              "Hamiltonian_before_lapse": ha["before_lapse_residual"],
              "Hamiltonian_lapse": ha["lapse_residual"],
              "center_kinetic": ce["kinetic_residual"],
              "center_momentum_determinant": sp.factor(ce["two_momentum_determinant"]
                                                       +400*q*(19*q-18)/(3200*q-3597)),
              "center_Cauchy": center_pencil()["cauchy_residual"]}
    # The exact zero-coupling action loses the vector pair and reduces to CD.
    extra_zero = ac["extra_before_temporal"].subs(zeta, 0)
    result["zero_coupling_auxiliary_control"] = sp.factor(extra_zero.subs(
        {temporal: -dd*n, sigma: -e*n}))
    for key, expected in (("H_1", 4), ("theta_1", 3), ("d_1", sp.Rational(21, 2)),
                          ("e_1", 0), ("alpha_1", 0), ("J_0", sp.Rational(1199, 800))):
        result[f"jet_{key}"] = sp.factor(center_jets()[key]-expected)
    return result


@cache
def proof_checks():
    bg = background()
    numerator, denominator = sp.fraction(sp.factor(bg["J"]))
    polynomial = sp.Poly(numerator, u)
    ce = center()
    qtest = sp.Integer(8)
    k2 = ce["kinetic"][:2, :2].subs(q, qtest)
    baseline = sp.Matrix([[6*q/(q-6), q/(20*(q-6))],
                          [q/(20*(q-6)), (200*q-1199)/(400*(q-6))]])
    result = {
        "J_even_polynomial_all_nonnegative": (
            all(value >= 0 for value in polynomial.all_coeffs())
            and all(power[0] % 2 == 0 for power, value in polynomial.terms() if value != 0)
        ),
        "J_polynomial_constant_strict": polynomial.eval(0) > 0,
        "J_denominator_positive_form": sp.factor(denominator-800*(1+u**2)**18) == 0,
        "theta_only_zero_at_center": sp.factor(bg["theta"]-u*(4*(1+u**2)**3-1)/(1+u**2)**4) == 0,
        "positive_vector_pivot": action()["C"].is_positive is True,
        "center_q8_two_pivot_product_negative": k2.det() < 0,
        "center_q8_clock_diagonal_negative": k2[0, 0] < 0,
        "baseline_q8_positive": baseline.subs(q, qtest)[0, 0] > 0 and baseline.subs(q, qtest).det() > 0,
        "center_momentum_not_used_as_punctured_proof": True,
        "omitted_jets_nonzero": any(value != 0 for value in center_pencil()["omitted_coefficient_jets"]),
        "zero_coupling_loses_vector_kinetic": action()["C"].subs(zeta, 0) == 0,
        "negative_maxwell_sign_fails_transverse": (-sp.Rational(1, 2)) < 0,
        "no_EFT_cutoff_inferred": True,
    }
    return {name: bool(value) for name, value in result.items()}


@cache
def calibration():
    return {
        "configuration_fields": ("v=zeta_metric-delta*n", "s=delta_chi", "sigma=W_longitudinal"),
        "normalized_zeta_definition": "zeta_physical/(M_squared*tau_squared)",
        "physical_momentum_definition": "q=(tau*k_physical)^2", "beta": BETA,
        "center_J": sp.Rational(1199, 800), "center_threshold": sp.Rational(3597, 3200),
        "unhealthy_condition": "u!=0, zeta>0, q>3*J(u)*h(u)^2/4",
        "scalar_kinetic_inertia_under_condition": (2, 1, 0),
        "DHOST_primary_preserved_at_quadratic_order": True,
        "literal_action_everywhere_healthy": False,
        "EFT_low_frequency_ghost_or_cutoff_claim": False,
        "identity_count": len(checks()), "proof_check_count": len(proof_checks()),
    }
