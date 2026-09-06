"""Independent covariant-lapse, minisuperspace and source Schur audits."""

import itertools

import sympy as sp


def test_literal_covariant_FLRW_curvature_and_Einstein_boundary():
    t = sp.Symbol("t", real=True)
    a, lapse = sp.Function("a", positive=True)(t), sp.Function("N", positive=True)(t)
    mass = sp.Symbol("M2", positive=True)
    metric = sp.diag(lapse**2, -a**2, -a**2, -a**2)
    inverse = metric.inv()
    derivative = lambda expression, index: sp.diff(expression, t) if index == 0 else 0
    connection = {}
    for r, m, n in itertools.product(range(4), repeat=3):
        connection[r, m, n] = sp.simplify(sum(inverse[r, s]*(
            derivative(metric[s, n], m)+derivative(metric[s, m], n)
            -derivative(metric[m, n], s)) for s in range(4))/2)
    ricci = {}
    for m, n in itertools.product(range(4), repeat=2):
        ricci[m, n] = sum(derivative(connection[r, m, n], r)
                         -derivative(connection[r, m, r], n)
                         +sum(connection[r, r, s]*connection[s, m, n]
                              -connection[r, n, s]*connection[s, m, r] for s in range(4))
                         for r in range(4))
    curvature = sp.simplify(sum(inverse[m, n]*ricci[m, n]
                               for m, n in itertools.product(range(4), repeat=2)))
    hubble = sp.diff(a, t)/(lapse*a)
    assert sp.simplify(curvature+6*(sp.diff(hubble, t)/lapse+2*hubble**2)) == 0
    covariant_density = -mass*lapse*a**3*curvature/2
    boundary = 3*mass*a**2*sp.diff(a, t)/lapse
    mini = -3*mass*a*sp.diff(a, t)**2/lapse
    assert sp.simplify(covariant_density-sp.diff(boundary, t)-mini) == 0


def test_full_lapse_variation_from_symmetric_polynomial_potential():
    t = sp.Symbol("t", real=True)
    a, b, ng, nf, phi, chi = (sp.Function(name)(t) for name in ("a", "b", "Ng", "Nf", "phi", "chi"))
    mg, mf, nu = sp.symbols("Mg2 Mf2 nu", positive=True)
    z = sp.Symbol("z")
    y, c = b/a, nf/ng
    generating = sp.Poly((1+z*c)*(1+z*y)**3, z)
    betas = (-3, 1, 0, 0, -1)  # Extract beta1 into nu=m4*beta1.
    potential = nu*ng*a**3*sum(betas[index]*generating.nth(index) for index in range(5))
    kinetic = sp.diff(phi, t)**2+sp.diff(chi, t)**2
    field_potential = sp.Function("V")(phi)
    matter = a**3*kinetic/(2*ng)-ng*a**3*field_potential
    lagrangian = (-3*mg*a*sp.diff(a, t)**2/ng-3*mf*b*sp.diff(b, t)**2/nf-potential+matter)
    hg, hf = sp.diff(a, t)/(ng*a), sp.diff(b, t)/(nf*b)
    rho, pressure = kinetic/(2*ng**2)+field_potential, kinetic/(2*ng**2)-field_potential
    g_constraint = 3*mg*hg**2-rho-3*nu*(y-1)
    f_constraint = 3*mf*hf**2-nu*(y**(-3)-1)
    assert sp.simplify(sp.diff(lagrangian, ng)/a**3-g_constraint) == 0
    assert sp.simplify(sp.diff(lagrangian, nf)/b**3-f_constraint) == 0
    euler_a = sp.diff(sp.diff(lagrangian, sp.diff(a, t)), t)-sp.diff(lagrangian, a)
    euler_b = sp.diff(sp.diff(lagrangian, sp.diff(b, t)), t)-sp.diff(lagrangian, b)
    g_pressure = -mg*(2*sp.diff(hg, t)/ng+3*hg**2)-pressure-nu*(3-c-2*y)
    f_pressure = -mf*(2*sp.diff(hf, t)/nf+3*hf**2)-nu*(1-1/(c*y**2))
    assert sp.simplify(euler_a/(3*ng*a**2)-g_pressure) == 0
    assert sp.simplify(euler_b/(3*nf*b**2)-f_pressure) == 0
    assert sp.simplify(rho+pressure-kinetic/ng**2) == 0
    # Fixing a lapse before varying really loses its constraint.
    assert sp.diff(lagrangian.subs(nf, 1).doit(), nf) == 0
    assert sp.simplify(f_constraint.subs(nf, 1).doit()) != 0


def test_covariant_source_balance_selects_the_regular_dynamic_branch():
    t = sp.Symbol("t", real=True)
    a, b, ng, nf, rho, pressure = (sp.Function(name)(t) for name in ("a", "b", "Ng", "Nf", "rho", "p"))
    mg, nu = sp.symbols("Mg2 nu", positive=True)
    y, c = b/a, nf/ng
    hg, hf = sp.diff(a, t)/(ng*a), sp.diff(b, t)/(nf*b)
    constraint = 3*mg*hg**2-rho-3*nu*(y-1)
    acceleration = -2*mg*sp.diff(hg, t)/ng-rho-pressure-nu*(y-c)
    matter_balance = sp.diff(rho, t)/ng+3*hg*(rho+pressure)
    noether = sp.diff(constraint, t)/ng+3*hg*acceleration+matter_balance
    assert sp.simplify(noether+3*nu*c*(y*hf-hg)) == 0
    # c and nu cannot be set to zero in the stated regular branch.
    assert sp.factor(-3*nu*sp.Symbol("c", positive=True)) != 0


def test_bounce_elimination_retains_the_second_metric_acceleration():
    mg, mf, nu, c, y = sp.symbols("Mg2 Mf2 nu c y", positive=True)
    h_dot, density_plus_pressure = sp.symbols("Hdot rho_plus_p", real=True)
    # y>0: the f constraint at Hg=Hf=0 selects y=1, not a free ratio.
    f_constraint = nu*(y**(-3)-1)
    assert sp.factor(-y**3*f_constraint/nu) == (y-1)*(y**2+y+1)
    assert sp.Poly(y**2+y+1, y).all_coeffs() == [1, 1, 1]
    g_acceleration = -2*mg*h_dot-density_plus_pressure-nu*(1-c)
    f_acceleration = -2*mf*h_dot-nu*(c-1)
    obstruction = -2*(mg+mf)*h_dot-density_plus_pressure
    assert sp.expand(g_acceleration+f_acceleration-obstruction) == 0
    solution = sp.solve(f_acceleration, c)
    assert len(solution) == 1
    assert sp.factor(solution[0]-(1-2*mf*h_dot/nu)) == 0
    assert sp.simplify(g_acceleration.subs(c, 1-2*mf*h_dot/nu)-obstruction) == 0
    # The target has positive Hdot; no positive kinetic matter sum can supply it.
    tau, matter_sum = sp.symbols("tau matter_sum", positive=True)
    assert (-obstruction).subs({h_dot: 4/tau**2, density_plus_pressure: matter_sum}).is_positive is True


def test_source_preserving_TT_schur_is_not_a_resummed_four_derivative_parent():
    mg, mf, nu = sp.symbols("Mg2 Mf2 nu", positive=True)
    d = sp.Symbol("D")
    kinetic = sp.Matrix([[mg*d+nu, -nu], [-nu, mf*d+nu]])
    schur = sp.factor(kinetic[0, 0]-kinetic[0, 1]*kinetic[1, 0]/kinetic[1, 1])
    total, mass_squared = mg+mf, nu*(1/mg+1/mf)
    assert sp.factor(schur-mg*d*(d+mass_squared)/(d+nu/mf)) == 0
    physical_propagator = sp.factor(kinetic.inv()[0, 0])
    target_propagator = (1/d+(mf/mg)/(d+mass_squared))/total
    assert sp.factor(physical_propagator-target_propagator) == 0
    assert sp.simplify(sp.limit(d*physical_propagator, d, 0)-1/total) == 0
    assert sp.simplify(sp.limit((d+mass_squared)*physical_propagator, d, -mass_squared)-mf/(mg*total)) == 0
    # S_T=-h K_eff h/8. Thus +cC(Dh)^2/2 corresponds to -4cC D^2 in K_eff.
    c_c = mf**2/(4*nu)
    fourth_order = total*d-4*c_c*d**2
    exact_remainder = sp.factor(schur-fourth_order)
    assert sp.factor(exact_remainder-mf**3*d**3/(nu*(nu+mf*d))) == 0
    assert sp.factor(c_c-total*(mf/mg)/(4*mass_squared)) == 0
    # The artificial extra root of the truncated inverse is not the exact pole.
    false_root = total/(4*c_c)
    assert sp.factor(schur.subs(d, false_root)) != 0


def test_public_lapse_equations_and_source_normalization_match_independent_algebra():
    from p8_bimetric import background, matching, vacuum

    data = background.equations()
    t, a, b, ng, nf = background.t, background.a, background.b, background.Ng, background.Nf
    mg, mf, nu = background.MG2, background.MF2, background.NU
    hg, hf = sp.diff(a, t)/(ng*a), sp.diff(b, t)/(nf*b)
    density = background.KIN/(2*ng**2)+background.POT
    pressure = background.KIN/(2*ng**2)-background.POT
    assert sp.simplify(data["EL_Ng"]-3*mg*hg**2+density+3*nu*(b/a-1)) == 0
    assert sp.simplify(data["EL_Nf"]-3*mf*hf**2+nu*((a/b)**3-1)) == 0
    assert sp.simplify(data["EL_a"]-mg*(2*sp.diff(hg, t)/ng+3*hg**2)
                       -pressure-nu*(3-nf/ng-2*b/a)) == 0
    assert sp.simplify(data["EL_b"]-mf*(2*sp.diff(hf, t)/nf+3*hf**2)
                       -nu*(1-ng*a**2/(nf*b**2))) == 0
    d = vacuum.D
    independent_matrix = sp.Matrix([[mg*d+nu, -nu], [-nu, mf*d+nu]])
    actual = matching.tt_schur()
    assert sp.factor(actual["physical_g_kernel"]-1/independent_matrix.inv()[0, 0]) == 0
    assert sp.factor(actual["c_C"]-mf**2/(4*nu)) == 0
    assert sp.factor(vacuum.spectrum()["m_FP_squared"]-nu*(1/mg+1/mf)) == 0
    expected_null_stress = -8*(mg+mf)/background.TAU**2
    assert sp.factor(background.bounce_equations()["CD_required_null_stress"]-expected_null_stress) == 0


def test_scalar_source_projector_from_lapse_constraint_not_TT_guess():
    from p8_bimetric import matching

    # In the rest frame of any non-null timelike Fourier momentum, conservation
    # makes T_0mu zero. Use the six orthonormal symmetric spatial components.
    # The relative lapse q remains in the full Fierz-Pauli mass, enforcing tr v=0.
    values = sp.Matrix(sp.symbols("v11 v22 v33 sqrt2v12 sqrt2v13 sqrt2v23"))
    source = sp.Matrix(sp.symbols("j11 j22 j33 sqrt2j12 sqrt2j13 sqrt2j23"))
    q, d = sp.symbols("relative_lapse D")
    reduced_mass, spring, total_mass = sp.symbols("mu nu M2", positive=True)
    trace = sum(values[:3, 0])
    quadratic_form = values.dot(values)-trace**2
    massive_action = (-reduced_mass*d*quadratic_form-spring*(quadratic_form+2*q*trace)
                      +2*source.dot(values))/8
    p = matching.projectors()
    solution = p["P2"]*source/(reduced_mass*d+spring)
    substitution = dict(zip(values, solution, strict=True))
    substitution[q] = sum(source[:3, 0])/(3*spring)
    for coordinate in (*values, q):
        assert sp.factor(sp.diff(massive_action, coordinate).subs(substitution, simultaneous=True)) == 0
    # The massless off-shell scalar constraint has a different trace coefficient;
    # its sign is not a negative physical on-shell residue.
    massless_action = (-total_mass*d*quadratic_form+2*source.dot(values))/8
    massless_solution = (p["P2"]-p["P0"]/2)*source/(total_mass*d)
    mapping = dict(zip(values, massless_solution, strict=True))
    for coordinate in values:
        assert sp.factor(sp.diff(massless_action, coordinate).subs(mapping, simultaneous=True)) == 0
    wrong_solution = (p["P2"]-p["P0"]/2)*source/(reduced_mass*d+spring)
    assert sp.factor(sp.diff(massive_action, q).subs(dict(zip(values, wrong_solution, strict=True)))) != 0
