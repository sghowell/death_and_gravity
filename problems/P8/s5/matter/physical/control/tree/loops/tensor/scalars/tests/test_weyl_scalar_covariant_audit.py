"""Independent scalar curvature, gauge, clock and phase-map audits."""

import itertools

import sympy as sp


def test_weyl_scalar_literal_four_index_contraction():
    lapse = sp.symbols("A_tt A_tz A_zz")
    curvature = sp.symbols("zeta_tt zeta_tz zeta_zz")
    shift = sp.symbols("B_ttz B_tzz B_zzz")
    shear = sp.symbols("E_ttzz E_tzzz E_zzzz")
    signs = (1, -1, -1, -1)
    derivative_index = {(0, 0): 0, (0, 3): 1, (3, 0): 1, (3, 3): 2}

    def metric_hessian(i, j, a, b):
        # Conformal physical metric: h00=2A, h0i=-B_,i,
        # hij=-2*zeta*delta_ij-2*E_,ij; all spatial dependence is along z.
        if (a, b) not in derivative_index:
            return sp.Integer(0)
        jet = derivative_index[a, b]
        result = 2*lapse[jet] if i == j == 0 else 0
        if (i, j) in ((0, 3), (3, 0)):
            result -= shift[jet]
        if i == j and i > 0:
            result -= 2*curvature[jet]
        if i == j == 3:
            result -= 2*shear[jet]
        return result

    riemann = {}
    for a, b, c, d in itertools.product(range(4), repeat=4):
        riemann[a, b, c, d] = (
            metric_hessian(a, d, c, b)+metric_hessian(b, c, d, a)
            - metric_hessian(a, c, d, b)-metric_hessian(b, d, c, a))/sp.Integer(2)
    ricci = {(b, d): sum(signs[a]*riemann[a, b, a, d] for a in range(4))
             for b, d in itertools.product(range(4), repeat=2)}
    scalar = sum(signs[a]*ricci[a, a] for a in range(4))
    riemann_squared = sum(sp.prod(signs[i] for i in indices)*value**2
                          for indices, value in riemann.items())
    ricci_squared = sum(signs[b]*signs[d]*value**2 for (b, d), value in ricci.items())
    weyl_squared = sp.factor(riemann_squared-2*ricci_squared+scalar**2/3)
    lensing_zz = lapse[2]-curvature[2]+shift[1]-shear[0]
    assert sp.expand(weyl_squared-sp.Rational(4, 3)*lensing_zz**2) == 0
    assert not weyl_squared.has(lapse[0], lapse[1], curvature[0], curvature[1],
                               shift[0], shift[2], shear[1], shear[2])
    # The scalar perturbation can have nonzero Weyl even though the FLRW
    # background has zero Weyl. Removing the shear term is not gauge safe.
    assert weyl_squared.subs({lapse[2]: 1, curvature[2]: 0,
                             shift[1]: 0, shear[0]: 0}) == sp.Rational(4, 3)


def test_weyl_scalar_lensing_gauge_invariance_and_physical_shift_clock():
    eta = sp.Symbol("eta", real=True)
    lapse, zeta, shift, shear, time_shift, spatial_shift, hc = (
        sp.Function(name)(eta) for name in ("A", "zeta", "B", "E", "T", "L", "Hc"))
    lensing = lapse-zeta+sp.diff(shift, eta)-sp.diff(shear, eta, 2)
    transformed = (lapse-sp.diff(time_shift, eta)-hc*time_shift)-(zeta-hc*time_shift)
    transformed += sp.diff(shift+time_shift-sp.diff(spatial_shift, eta), eta)
    transformed -= sp.diff(shear-spatial_shift, eta, 2)
    assert sp.simplify(transformed-lensing) == 0
    t = sp.Symbol("t", real=True)
    a, physical_shift = sp.Function("a", positive=True)(t), sp.Function("b")(t)
    # B_conformal=b_physical/a, d_eta=a*d_t.
    converted = a*sp.diff(physical_shift/a, t)
    assert sp.simplify(converted-sp.diff(physical_shift, t)
                       + sp.diff(a, t)*physical_shift/a) == 0


def test_weyl_scalar_full_phase_field_map_not_bare_equation_substitution():
    q = sp.Symbol("q", positive=True)
    lensing, momentum_equation = sp.symbols("L0 E_p", real=True)
    off_shell = lensing-momentum_equation/(2*q)
    density1 = sp.Rational(4, 3)*q**2*off_shell**2
    map1 = (momentum_equation-4*q*lensing)/3
    reduced_density = sp.Rational(4, 3)*q**2*lensing**2
    assert sp.expand(density1-momentum_equation*map1-reduced_density) == 0
    branch_only = -4*q*lensing/3
    assert sp.expand(density1-momentum_equation*branch_only-reduced_density) == momentum_equation**2/3


def test_weyl_scalar_rank_one_hamiltonian_normalization():
    q = sp.Symbol("q", positive=True)
    phase = sp.Matrix(sp.symbols("Q1 Q2 P1 P2", real=True))
    row = sp.Matrix(sp.symbols("r1 r2 r3 r4", real=True))
    lensing = (row.T*phase)[0]
    h1 = -sp.Rational(4, 3)*q**2*lensing**2
    hessian = sp.hessian(h1, phase)
    assert hessian == -sp.Rational(8, 3)*q**2*row*row.T
    symplectic = sp.BlockMatrix([[sp.zeros(2), sp.eye(2)], [-sp.eye(2), sp.zeros(2)]]).as_explicit()
    generator = symplectic*hessian
    assert (generator*generator).applyfunc(sp.expand) == sp.zeros(4)
    assert sp.expand(generator.trace()) == 0
    # Nilpotence of this difference generator does not remove time ordering
    # with the baseline oscillator evolution.
    baseline = sp.BlockMatrix([[sp.zeros(2), sp.eye(2)], [-q*sp.eye(2), sp.zeros(2)]]).as_explicit()
    assert (baseline*generator-generator*baseline).applyfunc(sp.expand) != sp.zeros(4)


def test_weyl_scalar_public_phase_and_auxiliary_reconstruction_bridge():
    from p8_m1_weyl_scalar import reduction as r

    data, reconstruction = r.phase_data(), r.reconstruction()
    actual_lensing = (1-r.DELTA)*data["n0"]-r.V+r.BDOT-r.H*data["b0"]
    # b0=-p/(2q), with qdot=-2Hq, before imposing the p equation.
    shift_derivative = -r.PDOT/(2*r.Q)-r.H*r.P/r.Q
    actual_lensing = actual_lensing.subs(r.BDOT, shift_derivative)
    assert sp.factor(actual_lensing-data["lensing_off"]) == 0
    branch_lensing = actual_lensing.subs(r.PDOT, -sp.diff(data["h0"], r.V)-3*r.H*r.P)
    assert sp.factor(branch_lensing-data["lensing0"]) == 0
    assert sp.factor(data["h1"]+sp.Rational(4, 3)*r.Q**2*branch_lensing**2) == 0

    # Derive the auxiliary Euler derivatives from the unreduced Weyl term.
    off = (1-r.DELTA)*r.N-r.V+r.BDOT-r.H*r.B
    density1 = sp.Rational(4, 3)*r.Q**2*off**2
    lapse_euler = sp.diff(density1, r.N)
    shift_jet = sp.diff(density1, r.BDOT)
    lensing_jet = sp.Symbol("calL", real=True)
    shift_jet = sp.Rational(8, 3)*r.Q**2*lensing_jet
    # a^-3*d_t(a^3*q^2*calL)=q^2*(calL_dot-H*calL).
    shift_jet_derivative = (3*r.H*shift_jet
                            + sp.diff(shift_jet, r.Q)*(-2*r.H*r.Q)
                            + sp.diff(shift_jet, lensing_jet)*r.LDOT)
    shift_euler = -sp.Rational(8, 3)*r.Q**2*r.H*lensing_jet-shift_jet_derivative
    assert sp.expand(shift_euler+sp.Rational(8, 3)*r.Q**2*r.LDOT) == 0
    lapse_before_map = -lapse_euler/(2*r.J)
    lapse_before_map = lapse_before_map.subs(r.BDOT, lensing_jet-(1-r.DELTA)*r.N+r.V+r.H*r.B)
    lapse_before_map = lapse_before_map.subs(lensing_jet, data["lensing0"])
    assert sp.factor(lapse_before_map-reconstruction["n1_before_phase_map"]) == 0
    assert sp.factor(-shift_euler/(2*r.Q**2/3)-reconstruction["b1"]) == 0
    lapse_pulled = lapse_before_map+sp.diff(data["n0"], r.V)*reconstruction["v1"]
    assert sp.factor(lapse_pulled-reconstruction["n1"]) == 0
    assert r.cd(lapse_pulled+sp.Rational(8, 3)*r.DELTA*r.Q**2*data["lensing0"]/r.J) == 0
    assert sp.factor(reconstruction["zeta1"]-reconstruction["v1"]-r.DELTA*lapse_pulled) == 0
    # The old gamma coordinate is the zeroth-order physical shift only.
    assert data["b0"].subs(r.chart_data("gamma")["mapping"], simultaneous=True) == r.QG
    assert reconstruction["b1"] != 0


def test_weyl_scalar_public_full_canonical_hamiltonian_bridge():
    from p8_m1_weyl_scalar import normalization as n

    q, eps, a, omega = sp.symbols("q epsilon a Omega", positive=True)
    variables = sp.Matrix(sp.symbols("Y1 Y2 P1 P2"))
    y, p = variables[:2, 0], variables[2:, 0]
    t = sp.Matrix([[2, 0], [3, 5]])
    boundary = sp.Matrix([[7, 11], [11, 13]])
    cq, cp = sp.Matrix([17, 19]), sp.Matrix([23, 29])
    potential = sp.Matrix([[q+31, 37], [37, q+41]])
    connection = sp.Matrix([[0, omega], [-omega, 0]])
    old_q = a**(-sp.Rational(3, 2))*t.inv()*y
    old_p = a**(-sp.Rational(3, 2))*t.T*(p-boundary*y)
    direct_lensing = (cq.T*old_q+cp.T*old_p)[0]
    old_hamiltonian = (p.dot(p)+y.dot(potential*y))/2-p.dot(connection*y)
    full_hamiltonian = old_hamiltonian-sp.Rational(4, 3)*eps*a**3*q**2*direct_lensing**2
    symplectic = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1],
                           [-1, 0, 0, 0], [0, -1, 0, 0]])
    direct_generator = symplectic*sp.hessian(full_hamiltonian, variables)
    ry, rp = n.row_map(cq, cp, t, boundary)
    public = n.canonical_generator(potential, connection, ry, rp, q, eps)["representative"]
    assert (direct_generator-public).applyfunc(sp.simplify) == sp.zeros(4)
    # This catches omission of either the symmetric boundary or connection.
    wrong_ry, wrong_rp = n.row_map(cq, cp, t, sp.zeros(2))
    wrong = n.canonical_generator(potential, sp.zeros(2), wrong_ry, wrong_rp, q, eps)["representative"]
    assert (direct_generator-wrong).applyfunc(sp.simplify) != sp.zeros(4)
