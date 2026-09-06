"""Regular phase reduction and physical reconstruction, through first beta.

beta=c_C/M^2 in cosmic units; in fixed centre ell0 units beta is replaced by
epsilon=c_C/(M*ell0)^2. These are formal first-order action identities and a
specified finite-dimensional reduced representative, not a resummation theorem.
"""

from functools import cache

import sympy as sp

V, S, P, PS, N, B = sp.symbols("v s p P_s n b", real=True)
Q, H, THETA, LAMBDA, DELTA, L, J = sp.symbols("q H theta Lambda delta l J", real=True)
VDOT, SDOT, PDOT, BDOT = sp.symbols("vdot sdot pdot bdot", real=True)
EP, LDOT = sp.symbols("E_p Lensing_dot", real=True)
QG, QS, PG, PM = sp.symbols("Q_g Q_s P_g P_s_chart", real=True)
COORDINATES = (V, S)
MOMENTA = (P, PS)


def cd(value):
    """The pinned CD algebraic relation, not a generic DHOST identity."""
    return sp.factor(value.subs(LAMBDA, 1-3*DELTA))


@cache
def phase_data():
    w = L*(3*DELTA-1)
    r = -THETA*P/2+LAMBDA*Q*V+w*PS/2-3*L*THETA*S/2
    n0, b0 = -r/J, -P/(2*Q)
    h0 = -Q*V**2+PS**2/2+(Q/2-3*L**2/4)*S**2-L*P*S/2+r**2/J
    ep = PDOT+3*H*P-2*Q*(V+LAMBDA*n0)
    off = (1-DELTA)*n0-V-(PDOT+H*P)/(2*Q)
    lensing = (1-DELTA-LAMBDA)*n0-2*(V+H*b0)
    return {"w": w, "R": r, "n0": n0, "b0": b0, "h0": h0,
            "E_p": ep, "lensing_off": off, "lensing0": lensing,
            "lensing_CD": 2*DELTA*n0-2*(V+H*b0),
            "h1": -sp.Rational(4, 3)*Q**2*lensing**2,
            "v_dot": sp.diff(h0, P), "s_dot": sp.diff(h0, PS),
            "p_dot": -sp.diff(h0, V)-3*H*P, "P_s_dot": -sp.diff(h0, S)-3*H*PS}


@cache
def auxiliary_checks():
    data = phase_data()
    w = data["w"]
    sigma = J-3*THETA**2+w**2/2
    lagrangian = (-3*VDOT**2+sigma*N**2+6*THETA*N*VDOT+SDOT**2/2
                  +w*N*SDOT-3*L*VDOT*S
                  +Q*(V**2+2*LAMBDA*N*V-S**2/2+2*THETA*N*B-2*VDOT*B-L*B*S))
    velocities = {VDOT: THETA*N-L*S/2-Q*B/3-P/6, SDOT: PS-w*N}
    hamiltonian = sp.expand((P*VDOT+PS*SDOT-lagrangian).subs(velocities, simultaneous=True))
    hessian = sp.hessian(hamiltonian, (N, B))
    stationary = {N: data["n0"], B: data["b0"]}
    values = {"partial_Legendre_metric_momentum": sp.expand(sp.diff(lagrangian, VDOT)+6*VDOT-6*THETA*N+3*L*S+2*Q*B),
              "partial_Legendre_matter_momentum": sp.expand(sp.diff(lagrangian, SDOT)-SDOT-w*N),
              "auxiliary_Hessian_nn": sp.expand(hessian[0, 0]+2*J),
              "auxiliary_Hessian_nb": hessian[0, 1],
              "auxiliary_Hessian_bb": sp.expand(hessian[1, 1]+2*Q**2/3),
              "auxiliary_Hessian_determinant": sp.expand(hessian.det()-4*J*Q**2/3),
              "stationary_n": sp.factor(sp.diff(hamiltonian, N).subs(stationary, simultaneous=True)),
              "stationary_b": sp.factor(sp.diff(hamiltonian, B).subs(stationary, simultaneous=True)),
              "stationary_Hamiltonian": sp.factor(hamiltonian.subs(stationary, simultaneous=True)-data["h0"])}
    return values


def field_map_checks():
    data = phase_data()
    lensing = data["lensing0"]
    off = lensing-EP/(2*Q)
    full_map = (EP-4*Q*lensing)/3
    first_order_density = sp.Rational(4, 3)*Q**2*off**2-EP*full_map
    return {"off_shell_lensing_residual": sp.factor(data["lensing_off"]-lensing+data["E_p"]/(2*Q)),
            "CD_lensing_reduction": cd(lensing-data["lensing_CD"]),
            "full_derivative_phase_map_action_identity": sp.expand(first_order_density-sp.Rational(4, 3)*Q**2*lensing**2),
            "lapse_regular_old_velocity": sp.factor(data["v_dot"]-THETA*data["n0"]+L*S/2),
            "matter_regular_old_velocity": sp.factor(data["s_dot"]-PS+data["w"]*data["n0"]),
            "shift_old_flow": sp.factor(-data["p_dot"]/(2*Q)-H*P/Q+V+LAMBDA*data["n0"]+H*data["b0"])}


def reconstruction():
    data = phase_data()
    lensing = data["lensing0"]
    v1 = -sp.Rational(4, 3)*Q*lensing
    n1_before_map = sp.Rational(4, 3)*Q**2*(DELTA-1)*lensing/J
    n1 = sp.factor(n1_before_map+sp.diff(data["n0"], V)*v1)
    b1 = 4*LDOT
    return {"v1": v1, "n1_before_phase_map": n1_before_map,
            "n1": n1, "b1": b1, "zeta0": V+DELTA*data["n0"],
            "zeta1": sp.factor(v1+DELTA*n1),
            "pre_mixed_metric_momentum": P+3*L*S,
            "matter_density_background_and_perturbation": L+PS+3*L*V,
            "matter_density_first_order_map": 3*L*v1}


def reconstruction_checks():
    data, rec = phase_data(), reconstruction()
    lensing = sp.Symbol("Lensing", real=True)
    # qdot=-2Hq and (a^3)'=3Ha^3: the b Euler derivative is exactly
    # -(8/3)q^2*Lensing_dot, not -(8/3)q^2*(Lensing_dot+H*Lensing).
    b_euler = -sp.Rational(8, 3)*Q**2*H*lensing-sp.Rational(8, 3)*Q**2*(LDOT-H*lensing)
    n_euler = sp.Rational(8, 3)*Q**2*(1-DELTA)*lensing
    # Differentiate the branch v map before comparing with the ORIGINAL
    # partial Legendre velocities. The qdot term cannot be omitted.
    v1_dot = sp.Rational(8, 3)*H*Q*data["lensing0"]-sp.Rational(4, 3)*Q*LDOT
    return {"first_order_lapse_stationarity": sp.expand(-2*J*(sp.Rational(4, 3)*Q**2*(DELTA-1)*lensing/J)-n_euler),
            "first_order_shift_stationarity": sp.expand(-2*Q**2/3*rec["b1"]-b_euler),
            "lapse_field_map_pullback": sp.factor(rec["n1"]-sp.Rational(4, 3)*Q**2*(LAMBDA+DELTA-1)*data["lensing0"]/J),
            "CD_physical_lapse_map": cd(rec["n1"]+sp.Rational(8, 3)*DELTA*Q**2*data["lensing0"]/J),
            "physical_spatial_metric_map": sp.factor(rec["zeta1"]-rec["v1"]-DELTA*rec["n1"]),
            "original_matter_Legendre_velocity_through_first_order": sp.factor(sp.diff(data["h1"], PS)+data["w"]*rec["n1"]),
            "original_metric_Legendre_velocity_through_first_order": sp.factor(6*(sp.diff(data["h1"], P)+v1_dot)-6*THETA*rec["n1"]+2*Q*rec["b1"])}


def old_flow_derivative(value):
    """Cosmic-time derivative along the old phase flow, including densities.

    This symbolic route keeps all coefficient jets and qdot=-2Hq. In fixed
    local units the separate compact conversion must include ell drift.
    """
    data = phase_data()
    time_jets = {H: sp.Symbol("H_dot"), THETA: sp.Symbol("theta_dot"),
                 J: sp.Symbol("J_dot"), DELTA: sp.Symbol("delta_dot"),
                 LAMBDA: sp.Symbol("Lambda_dot"), L: -3*H*L, Q: -2*H*Q}
    time_jets.update({V: data["v_dot"], S: data["s_dot"],
                     P: data["p_dot"], PS: data["P_s_dot"]})
    return sum(sp.diff(value, variable)*rate for variable, rate in time_jets.items())


@cache
def chart_data(chart):
    data = phase_data()
    if chart == "unitary":
        mapping = {V: QG, S: QS, P: PG, PS: PM}
        generator = 0
    elif chart == "gamma":
        mapping = {V: PG/(2*Q), S: QS, P: -2*Q*QG, PS: PM}
        generator = -H*QG*PG
    else:
        raise ValueError("Use unitary or gamma")
    lensing = sp.expand(data["lensing_CD"].subs(mapping, simultaneous=True))
    variables = (QG, QS, PG, PM)
    row = sp.Matrix([sp.diff(lensing, item) for item in variables])
    return {"variables": variables, "lensing": lensing, "lensing_row": row,
            "h0": data["h0"].subs(mapping, simultaneous=True)+generator,
            "h1": -sp.Rational(4, 3)*Q**2*lensing**2,
            "time_dependent_gamma_generator": generator, "mapping": mapping}


def rank_one_checks():
    values = {}
    symplectic = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1],
                           [-1, 0, 0, 0], [0, -1, 0, 0]])
    for chart in ("unitary", "gamma"):
        data = chart_data(chart)
        correction = sp.hessian(data["h1"], data["variables"])
        target = -sp.Rational(8, 3)*Q**2*data["lensing_row"]*data["lensing_row"].T
        for index, value in enumerate(correction-target):
            values[f"{chart}_rank_one_Hessian_{index}"] = sp.factor(value)
        generator = symplectic*correction
        for index, value in enumerate(generator.T*symplectic+symplectic*generator):
            values[f"{chart}_Hamiltonian_generator_{index}"] = sp.factor(value)
    return values


def controls():
    data = phase_data()
    lensing = data["lensing0"]
    # A branch-only map is not an off-shell action identity.
    wrong_action = sp.Rational(4, 3)*Q**2*(lensing-EP/(2*Q))**2+sp.Rational(4, 3)*Q*lensing*EP
    bounce = {H: 0, THETA: 0, DELTA: sp.Rational(1, 2), LAMBDA: -sp.Rational(1, 2), L: sp.Rational(1, 10)}
    return {"omitting_E_p_field_map_term": sp.factor(wrong_action-sp.Rational(4, 3)*Q**2*lensing**2),
            "matter_momentum_enters_lensing_at_bounce": sp.diff(data["lensing0"], PS).subs(bounce),
            "omitting_lapse_reconstruction_pullback": sp.factor(sp.diff(data["n0"], V)*(-sp.Rational(4, 3)*Q*lensing)),
            "omitting_expansion_density_in_shift_flow": -3*H*P,
            "omitting_mixed_matter_boundary_changes_metric_momentum": 3*L*S,
            "omitting_gamma_time_generator": H*QG*PG}
