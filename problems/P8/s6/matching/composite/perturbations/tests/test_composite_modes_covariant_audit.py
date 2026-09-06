"""Independent anisotropy, shift-elimination, clock and jet audits."""

import sympy as sp


def test_covariant_diagonal_exponential_tt_stiffness_retains_matter():
    epsilon, zz = sp.symbols("epsilon z", real=True)
    a, ng, y, c, alpha, beta, m4 = sp.symbols("a Ng y c alpha beta m4", positive=True)
    pressure = sp.Symbol("p", real=True)
    betas = sp.symbols("beta0:5", real=True)
    roots = [c, y*sp.exp(epsilon/2), y*sp.exp(-epsilon/2), y]
    polynomial = sp.Poly(sp.prod(1+zz*value for value in roots), zz)
    interaction = -m4*ng*a**3*sum(betas[j]*polynomial.nth(j) for j in range(5))
    effective_volume = a**3*sp.prod(alpha+beta*value for value in roots[1:])
    matter = ng*(alpha+beta*c)*effective_volume*pressure
    second_coefficient = sp.diff(interaction+matter, epsilon, 2).subs(epsilon, 0)/2
    mu = y*(m4*(betas[1]+betas[2]*(c+y)+betas[3]*c*y)
            -alpha*beta*(alpha+beta*c)*(alpha+beta*y)*pressure)
    assert sp.simplify(second_coefficient+ng*a**3*mu/4) == 0
    without_matter = sp.diff(interaction, epsilon, 2).subs(epsilon, 0)/2
    assert sp.factor(without_matter-second_coefficient) != 0


def test_two_shift_schur_complement_including_zero_tensor_stiffness():
    aa, bb, cc = sp.symbols("A B C", positive=True)
    x, z, first, second = sp.symbols("xdot zdot shift_g shift_f", real=True)
    lagrangian = aa*(x-first)**2+bb*(z-second)**2+cc*(first-second)**2
    stationary = sp.solve([sp.diff(lagrangian, first), sp.diff(lagrangian, second)], [first, second])
    reduced = sp.factor(lagrangian.subs(stationary))
    harmonic = 1/(1/aa+1/bb+1/cc)
    assert sp.factor(reduced-harmonic*(x-z)**2) == 0
    # The vector shift stiffness Xi is independent input here; setting
    # the tensor/spatial stiffness mu=0 does not set this cc to zero.
    g, f, a, y, c, k, xi = sp.symbols("G F a y c k Xi", positive=True)
    actual = {aa: k**2/2, bb: f*y**3*k**2/(2*g*c), cc: a**2*xi/(2*g)}
    coefficient = sp.factor(harmonic.subs(actual)/k**2)
    expected = 1/(2*(1+g*c/(f*y**3)+g*k**2/(a**2*xi)))
    assert sp.factor(coefficient-expected) == 0
    assert expected.subs({g: 1, f: 1, a: sp.Rational(1, 2), y: 1, c: 1, xi: 6}).is_positive


def test_composite_physical_vector_speed_clock():
    g, f, a, ng, y, c, k, xi, r, s = sp.symbols("G F a Ng y c k Xi r s", positive=True)
    mu = sp.Symbol("mu", real=True)
    inverse_c = 1+g*c/(f*y**3)+g*k**2/(a**2*xi)
    coordinate_omega_squared = ng**2*mu*inverse_c/g
    coordinate_principal = sp.diff(coordinate_omega_squared, k, 2)/2
    physical_speed = coordinate_principal*(a*r)**2/(ng*s)**2
    assert sp.factor(physical_speed-(r/s)**2*mu/xi) == 0


def test_general_bounce_acceleration_and_independent_parent_scale():
    accel = sp.Symbol("A", real=True)
    root = sp.sqrt(sp.Rational(7, 3))
    # At the common initial data, both null equations give X'=Y'=-2.
    yprime, xprime, zprime = -root, -sp.S(2), -sp.S(2)
    denominator, denominator_prime = 2*root, xprime-yprime*(-root)-zprime
    nf_numerator, nf_numerator_prime = root, xprime-2*accel
    nfprime = (nf_numerator_prime*denominator-nf_numerator*denominator_prime)/denominator**2
    cprime = sp.simplify(4*nfprime)
    quadratic_mu = sp.simplify(yprime*(yprime-cprime)/2)
    assert sp.simplify(cprime+(5+12*accel)/(3*root)) == 0
    assert quadratic_mu == sp.Rational(1, 3)-2*accel
    assert quadratic_mu.subs(accel, sp.Rational(7, 36)) == -sp.Rational(1, 18)
    assert quadratic_mu.subs(accel, 4) == -sp.Rational(23, 3)
    assert quadratic_mu.subs(accel, sp.Rational(1, 25)) > 0
    mtau = sp.Symbol("m_tau", positive=True)
    assert sp.factor(quadratic_mu.subs(accel, 4/mtau**2)) == (mtau**2-24)/(3*mtau**2)


def test_zero_algebraic_stiffness_does_not_determine_canonical_frequency():
    t, rate = sp.symbols("t rate", real=True, nonzero=True)
    inertia = sp.exp(2*rate*t)
    omega_squared_when_mu_zero = -sp.diff(sp.sqrt(inertia), t, 2)/sp.sqrt(inertia)
    assert sp.simplify(omega_squared_when_mu_zero+rate**2) == 0
    assert omega_squared_when_mu_zero != 0
