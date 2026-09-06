"""Independent action, common-clock cone, branch and outside-contract audits.

The independent derivations do not call the main cone engine or import
the printed reduced-vector normalization. A separate final test bridges
their results to the public API. The local bounce controls do not assert a subcutoff growth rate,
a global CD reconstruction, or a healthy ultraviolet completion.
"""

import sympy as sp


def test_literal_two_tensor_inertias_and_physical_clock_principal_roots():
    g, f, a, ng, y, c, alpha, beta = sp.symbols("G F a Ng y c alpha beta", positive=True)
    r, s = alpha+beta*y, alpha+beta*c
    b, nf, ae, ne = a*y, ng*c, a*r, ng*s
    inertia = sp.diag(g*a**3/(8*ng), f*b**3/(8*nf))
    gradient = sp.diag(g*ng*a/8, f*nf*b/8)
    # Actual dT=Ne dt and k_physical=k/Ae, not just a field redefinition.
    physical_principal = sp.simplify(ae**2*(ne*inertia).inv()*(gradient/ne))
    expected = sp.diag((r/s)**2, (c*r/(y*s))**2)
    assert all(sp.factor(entry) == 0 for entry in physical_principal-expected)
    speed_squared = sp.Symbol("lambda", real=True)
    assert sp.factor((physical_principal-speed_squared*sp.eye(2)).det()
                     -(speed_squared-(r/s)**2)*(speed_squared-(c*r/(y*s))**2)) == 0


def test_both_subluminal_tensor_inequalities_force_equal_ratio_and_lapse():
    y, c, alpha, beta = sp.symbols("y c alpha beta", positive=True)
    r, s = alpha+beta*y, alpha+beta*c
    cg, cf = r/s, c*r/(y*s)
    assert sp.factor(cg-1-beta*(y-c)/s) == 0
    assert sp.factor(cf-1-alpha*(c-y)/(y*s)) == 0
    gfactor, ffactor = beta*(r+s)/s**2, alpha*(c*r+y*s)/(y**2*s**2)
    assert gfactor.is_positive and ffactor.is_positive
    assert sp.factor(cg**2-1-gfactor*(y-c)) == 0
    assert sp.factor(cf**2-1-ffactor*(c-y)) == 0
    assert sp.factor((cg-1)*(cf-1)+alpha*beta*(c-y)**2/(y*s**2)) == 0
    # Independent exact controls for the two sides of the ratio ordering.
    assert (cg**2).subs({alpha: 1, beta: 1, y: 2, c: 1}) == sp.Rational(9, 4)
    assert (cf**2).subs({alpha: 1, beta: 1, y: 1, c: 2}) == sp.Rational(16, 9)


def test_literal_source_aware_tensor_stiffness_has_pressure_factor_at_cone_coincidence():
    epsilon, zz = sp.symbols("epsilon z", real=True)
    y, c, alpha, beta, m4 = sp.symbols("y c alpha beta m4", positive=True)
    pressure = sp.Symbol("p", real=True)
    betas = sp.symbols("beta0:5", real=True)
    roots = [c, y*sp.exp(epsilon/2), y*sp.exp(-epsilon/2), y]
    generating = sp.Poly(sp.prod(1+zz*root for root in roots), zz)
    interaction = -m4*sum(betas[j]*generating.nth(j) for j in range(5))
    r, s = alpha+beta*y, alpha+beta*c
    matter = s*pressure*sp.prod(alpha+beta*root for root in roots[1:])
    # The literal norm-two diagonal polarization has L2/(Ng a^3)=-mu/4.
    mu = sp.simplify(-2*sp.diff(interaction+matter, epsilon, 2).subs(epsilon, 0))
    p_poly = m4*(betas[1]+2*betas[2]*y+betas[3]*y**2)
    q = p_poly-alpha*beta*r**2*pressure
    remainder = m4*(betas[2]+betas[3]*y)-alpha*beta**2*r*pressure
    assert sp.factor(mu-y*q-y*(c-y)*remainder) == 0
    assert sp.factor(mu.subs(c, y)-y*q) == 0
    xi = sp.Symbol("Xi", positive=True)
    assert sp.factor(((r/s)**2*mu/xi).subs(c, y)-y*q/xi) == 0


def test_literal_shared_source_null_weights_and_full_lapse_monotone_reduction():
    a, b, ng, nf, alpha, beta = sp.symbols("a b Ng Nf alpha beta", positive=True)
    chi_dot, potential = sp.symbols("chi_dot V", real=True)
    ae, ne = alpha*a+beta*b, alpha*ng+beta*nf
    lm = ae**3*chi_dot**2/(2*ne)-ne*ae**3*potential
    rho = chi_dot**2/(2*ne**2)+potential
    pressure = chi_dot**2/(2*ne**2)-potential
    rho_g, pg = -sp.diff(lm, ng)/a**3, sp.diff(lm, a)/(3*ng*a**2)
    rho_f, pf = -sp.diff(lm, nf)/b**3, sp.diff(lm, b)/(3*nf*b**2)
    r, s, y, c = ae/a, ne/ng, b/a, nf/ng
    assert sp.factor(rho_g-alpha*r**3*rho) == 0
    assert sp.factor(pg-alpha*s*r**2*pressure) == 0
    assert sp.factor(rho_f-beta*r**3*rho/y**3) == 0
    assert sp.factor(pf-beta*s*r**2*pressure/(c*y**2)) == 0
    assert sp.factor((rho_g+pg-alpha*r**3*(rho+pressure)).subs(nf, y*ng)) == 0

    t = sp.Symbol("t", real=True)
    av, bv, nv, fv = (sp.Function(name)(t) for name in ("a", "b", "Ng", "Nf"))
    ratio = bv/av
    lapse_ratio = fv/nv
    B = nv*sp.diff(bv, t)-fv*sp.diff(av, t)
    assert sp.factor(sp.diff(ratio, t)-B/(nv*av)-(lapse_ratio-ratio)*sp.diff(av, t)/av) == 0
    # Cone coincidence plus QB=0 with Q>0 sets c=y, B=0 and hence y'=0.
    yc, g = sp.symbols("constant_y G", positive=True)
    null = sp.Symbol("scalar_null", nonnegative=True)
    rc = alpha+beta*yc
    hg = sp.diff(av, t)/(nv*av)
    he = sp.diff(rc*av, t)/((rc*nv)*(rc*av))
    assert sp.factor(he-hg/rc) == 0
    physical_derivative = sp.diff(he, t)/(rc*nv)
    # Retain Ndot before imposing the actual g-null acceleration equation.
    acceleration = (sp.diff(av, t)*sp.diff(nv, t)/nv+sp.diff(av, t)**2/av
                    -nv**2*av*alpha*rc**3*null/(2*g))
    assert sp.factor(physical_derivative.subs(sp.diff(av, t, 2), acceleration)
                     +alpha*rc*null/(2*g)) == 0


def test_arbitrarily_degenerate_sign_change_cannot_follow_from_nonincreasing_hubble():
    physical_time = sp.Symbol("T", real=True)
    g, alpha, r = sp.symbols("G alpha r", positive=True)
    proposed = physical_time**3
    required_null = -2*g*sp.diff(proposed, physical_time)/(alpha*r)
    assert sp.factor(required_null+6*g*physical_time**2/(alpha*r)) == 0
    assert required_null.subs(physical_time, 1).is_negative
    # This exact example is a control. The all-degeneracy theorem follows
    # from monotonicity on the whole connected interval, not from jet order.


def test_allowing_zero_vector_principal_speed_admits_regular_luminal_tensor_bounce():
    y, rho = sp.symbols("y rho", positive=True)
    r = 1+y
    pressure = 2*y/r**2
    null = rho+pressure
    x = sp.sqrt(y**2+r**3*rho/3)
    z = -sp.sqrt(y**-2+r**3*rho/(3*y**3))
    # An autonomous prescription, distinct from both the free and CD
    # frozen backgrounds, makes c=y identically along the pressure branch.
    h = (x+y**2*z)/r**2
    ng, nf = 1/r, y/r
    vy, vrho = y*(y*z-x)/r, -3*h*null

    def derivative(expression):
        return sp.diff(expression, y)*vy+sp.diff(expression, rho)*vrho

    assert sp.simplify(-2*derivative(x)-ng*r**3*null) == 0
    assert sp.simplify(-2*derivative(z)-nf*r**3*null/y**3) == 0
    assert sp.simplify(h-(ng*x+nf*y*z)/r) == 0
    assert sp.simplify(vy-y*(nf*z-ng*x)) == 0
    assert sp.factor(2*y-r**2*pressure) == 0
    assert sp.factor(vrho+3*h*null) == 0
    # Canonical realization: phi'=sqrt(null), V=(rho-pressure)/2;
    # n>0 gives a local analytic inverse field clock and potential.
    scalar = derivative(null)/(2*sp.sqrt(null))+3*h*sp.sqrt(null)
    scalar += derivative((rho-pressure)/2)/sp.sqrt(null)
    assert sp.simplify(scalar) == 0
    initial = {y: 1, rho: sp.Rational(1, 2)}
    assert h.subs(initial) == 0
    assert sp.simplify(derivative(h).subs(initial)) == sp.Rational(1, 6)
    assert vy.subs(initial) == -sp.sqrt(sp.Rational(7, 3))
    assert null.subs(initial) == 1
    # TT cones are exactly luminal, whereas mu=y Q is exactly zero.
    assert nf/ng == y
    mu = y*(2*y-r**2*pressure)
    assert sp.factor(mu) == 0
    xi = 2*y**2*(2*y+r**2*rho+r**2*null/2)/(2*y)
    assert xi.is_positive


def test_asymmetric_local_bounce_can_have_positive_vector_but_one_outside_tensor_cone():
    y, rho, pressure = sp.S(2), sp.Rational(1, 2), sp.Rational(4, 9)
    r = 1+y
    x, z = sp.sqrt(34)/2, -sp.sqrt(13)/4
    c = -x/(y*z)
    ng, nf = 1/(1+c), c/(1+c)
    assert sp.simplify(3*x**2-3*y**2-r**3*rho) == 0
    assert sp.simplify(3*z**2-3/y**2-r**3*rho/y**3) == 0
    assert sp.simplify((ng*x+nf*y*z)/r) == 0
    assert sp.simplify(4-c**2) == sp.Rational(18, 13)
    mu = y*(y-1)*(y-c)/r
    zv = 2*y+r**2*rho+r*(1+c)*y*(rho+pressure)/(c+y)
    xi = 2*y**2*zv/(c+y)
    cg2, cf2 = (r/(1+c))**2, (c*r/(y*(1+c)))**2
    assert mu.is_positive and xi.is_positive
    assert (cg2-1).is_positive and (1-cf2).is_positive
    assert (cg2*mu/xi).is_positive
    # The inherited regular reconstruction ODE can prescribe any analytic
    # h with h(0)=0, h'(0)>0 at these positive-root/lapse data. This is not
    # an all-time health or old-matter-frame matching control.


def test_proportional_de_sitter_satisfies_principal_contract_without_bouncing():
    # G=F=m4=alpha=beta=1, beta1=beta3=1, others zero,
    # g=f, H_g=H_f=sqrt(2), constant canonical field at V=1/4.
    rho, pressure = sp.Rational(1, 4), -sp.Rational(1, 4)
    h2, potential_g, potential_f = sp.S(2), sp.S(4), sp.S(4)
    assert 3*h2-potential_g-8*rho == 0
    assert 3*h2-potential_f-8*rho == 0
    assert 3*h2-potential_g+8*pressure == 0
    assert 3*h2-potential_f+8*pressure == 0
    q = 2-4*pressure
    mu = q
    xi = 2+4*rho+2*(rho+pressure)
    assert q == mu == xi == 3
    assert mu/xi == 1
    assert rho+pressure == 0


def test_public_api_matches_independently_derived_physical_roots_and_pressure_identity():
    from p8_composite_cones import background, cone
    from p8_composite_modes import model as m

    y, c, alpha, beta, m4, xi = sp.symbols("audit_y audit_c audit_alpha audit_beta audit_m4 audit_Xi", positive=True)
    pressure, b1, b2, b3 = sp.symbols("audit_p audit_b1 audit_b2 audit_b3", real=True)
    r, s = alpha+beta*y, alpha+beta*c
    assignment = {m.y: y, m.c: c, m.alpha: alpha, m.beta: beta, m.M4: m4,
                  m.p: pressure, m.betas[1]: b1, m.betas[2]: b2, m.betas[3]: b3}
    expected_mu = y*(m4*(b1+b2*(y+c)+b3*c*y)-alpha*beta*r*s*pressure)
    expected_q = m4*(b1+2*b2*y+b3*y**2)-alpha*beta*r**2*pressure
    d = cone.derive()
    assert sp.factor(d["c_g_squared"].subs(assignment, simultaneous=True)-(r/s)**2) == 0
    assert sp.factor(d["c_f_squared"].subs(assignment, simultaneous=True)-(c*r/(y*s))**2) == 0
    assert sp.factor(d["mu"].subs(assignment, simultaneous=True)-expected_mu) == 0
    assert sp.factor(d["Q"].subs(assignment, simultaneous=True)-expected_q) == 0
    linked_vector = d["c_vector_squared"].subs(d["vector_mu"], d["mu"])
    assert sp.factor(linked_vector.subs(assignment | {d["Xi"]: xi}, simultaneous=True)
                     -(r/s)**2*expected_mu/xi) == 0
    result = background.coincident_pressure_control()
    hsymbols = {symbol.name: symbol for symbol in result["h_function"].free_symbols}
    rho = sp.Symbol("audit_rho", positive=True)
    x = sp.sqrt(y**2+(1+y)**3*rho/3)
    z = -sp.sqrt(y**-2+(1+y)**3*rho/(3*y**3))
    converted = result["h_function"].subs({hsymbols["reconstructed_y"]: y,
                                            hsymbols["reconstructed_rho"]: rho}, simultaneous=True)
    assert sp.simplify(converted-(x+y**2*z)/(1+y)**2) == 0
    assert result["initial_H_prime"] == sp.Rational(1, 6)
