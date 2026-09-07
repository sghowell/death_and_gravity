"""Separately authored full-action, clock, sign and global-scope audits."""

from fractions import Fraction as Q

import sympy as sp
from p8_trimetric import matter, model
from p8_trimetric_global import background, monotonic


def zero(value):
    entries = list(value) if isinstance(value, sp.MatrixBase) else [value]
    assert all(sp.cancel(entry) == 0 for entry in entries)


def test_full_covariant_Einstein_maps_give_each_lapse_and_all_spatial_equations():
    ne, ae, nv, av, nu, au = sp.symbols("ne ae nv av nu au", positive=True)
    g, f = sp.symbols("G F", positive=True)
    pg, pf, bg, bf, hg, hf, dg, df = sp.symbols("pg pf bg bf Hg Hf DgH DfH", real=True)
    e, v, u = sp.diag(ne, ae, ae, ae), sp.diag(nv, av, av, av), sp.diag(nu, au, au, au)
    eg = sp.diag(3*hg**2/ne**2, *([- (2*dg+3*hg**2)/ae**2]*3))
    ef = sp.diag(3*hf**2/nv**2, *([- (2*df+3*hf**2)/av**2]*3))
    result = model.euler_maps(e, v, u, pg=pg, pf=pf, bg=bg, bf=bf,
                              g=g, f=f, einstein_g=eg, einstein_f=ef)
    for key, ni, ai, gi, pi, bi, hi, dhi in (
        ("E_e", ne, ae, g, pg, bg, hg, dg),
        ("E_v", nv, av, f, pf, bf, hf, df),
    ):
        ratio, speed = au/ai, au*ni/(nu*ai)
        lapse = result[key][0, 0]/ai**3
        zero(lapse-(3*gi*hi**2-2*pi*ratio**3-2*bi))
        for index in (1, 2, 3):
            spatial = result[key][index, index]/(ni*ai**2)
            zero(spatial-(gi*(2*dhi+3*hi**2)-2*pi*ratio**3/speed-2*bi))
            zero((spatial-lapse)/(2*gi)-(dhi-pi*ratio**3*(1/speed-1)/gi))


def test_literal_EH_and_interaction_minisuperspace_variation_independently_checks_signs():
    # Boundary-reduced -G R_B/2 equals -3G a adot²/n in the B convention.
    n, a, nu, au, g = sp.symbols("n a nu au G", positive=True)
    adot, addot, ndot, p, b = sp.symbols("adot addot ndot p b", real=True)
    lag = -3*g*a*adot**2/n-2*p*(n*au**3+3*nu*au**2*a)-2*b*n*a**3
    momentum = sp.diff(lag, adot)
    dt_momentum = sp.diff(momentum, a)*adot+sp.diff(momentum, adot)*addot+sp.diff(momentum, n)*ndot
    lapse, spatial = sp.diff(lag, n)/a**3, (sp.diff(lag, a)-dt_momentum)/(3*n*a**2)
    h = adot/(n*a)
    dh = addot/(n**2*a)-adot*ndot/(n**3*a)-adot**2/(n**2*a**2)
    ratio, speed = au/a, au*n/(nu*a)
    zero(lapse-(3*g*h**2-2*p*ratio**3-2*b))
    zero(spatial-(g*(2*dh+3*h**2)-2*p*ratio**3/speed-2*b))


def test_regular_flat_vacuum_lapse_calibration_also_forces_linked_clock_proportionality():
    n, a, nu, au = sp.symbols("n a nu au", positive=True)
    p = sp.Symbol("p", nonzero=True, real=True)
    e, u = sp.diag(n, a, a, a), sp.diag(nu, au, au, au)
    out = model.euler_maps(e, sp.eye(4), u, pg=p, pf=0, bg=-p*(au/a)**3)
    zero(out["E_e"][0, 0])
    for index in (1, 2, 3):
        zero(out["E_e"][index, index]/(2*p*nu*au**2)-(au*n/(nu*a)-1))


def test_algebraic_negative_link_cap_does_not_assert_a_full_flat_vacuum():
    # The endpoint coefficients calibrate both Einstein lapse equations;
    # choosing the wrong B still violates the full auxiliary equation.
    out = model.euler_maps(sp.eye(4), sp.eye(4), sp.eye(4), pg=-2, pf=1, bg=2, bf=-1, b=0)
    zero(out["E_e"])
    zero(out["E_v"])
    assert out["E_u"] != sp.zeros(4)


def test_vacuum_weight_bound_has_one_positive_factorization_for_both_link_signs():
    r, r0, g = sp.symbols("R R0 G", positive=True)
    p, h = sp.symbols("p H", real=True)
    gap = p/r0-p/r
    divided_difference = r*r0*(r**2+r*r0+r0**2)
    zero(gap-p*(r**3-r0**3)/divided_difference)
    # Use the actual Friedmann equation, not an assumption p>0.
    positive = 3*g*h**2/(2*divided_difference)
    assert positive.is_nonnegative is True
    zero(2*divided_difference*(gap-positive)-(2*p*(r**3-r0**3)-3*g*h**2))


def test_fraction_points_check_both_sign_branches_and_uniform_mixed_denominator_gap():
    pg, pf, rg0, rf0 = Q(-2), Q(1), Q(1), Q(1)
    s0 = pg/rg0+pf/rf0
    assert s0 == -1
    for rg, rf in ((Q(1), Q(1)), (Q(1, 2), Q(3, 2)), (Q(3, 4), Q(2)), (Q(1, 8), Q(7))):
        for p, r, r0 in ((pg, rg, rg0), (pf, rf, rf0)):
            h_squared = 2*p*(r**3-r0**3)/3
            assert h_squared >= 0
            assert p/r <= p/r0
            assert p/r0-p/r == 3*h_squared/(2*r*r0*(r*r+r*r0+r0*r0))
        assert pg/rg+pf/rf <= s0 < 0


def test_arbitrary_vacuum_ratios_use_physical_kinetics_in_the_mass_dictionary():
    g, f, rg, rf = sp.symbols("G F Rg0 Rf0", positive=True)
    pg, pf, denominator = sp.symbols("Pg0 Pf0 S0", nonzero=True, real=True)
    kinetic = sp.diag(g/rg**2, f/rf**2)
    spring = 2*pg*pf/denominator
    mass = spring*sp.Matrix([[1, -1], [-1, 1]])
    eigenvalues = (kinetic.inv()*mass).eigenvals()
    assert 0 in eigenvalues
    zero(sum(eigenvalues)-spring*(rg**2/g+rf**2/f))
    # Concrete mixed-sign non-unit ratios: a positive eigenvalue needs S0<0.
    point = {g: 3, f: 5, rg: 2, rf: 3, pg: -4, pf: 1, denominator: -3}
    assert spring.subs(point) > 0
    assert (spring*(rg**2/g+rf**2/f)).subs(point) == sp.Rational(376, 45)


def test_full_tensor_principal_coefficients_follow_from_measure_and_clock_pullback():
    gi, ni, ai, nu, au = sp.symbols("Gi ni ai nu au", positive=True)
    physical_measure = ni*ai**3/(nu*au**3)
    kinetic = gi*physical_measure*(nu/ni)**2
    gradient = gi*physical_measure*(au/ai)**2
    ratio, speed = au/ai, au*ni/(nu*ai)
    zero(kinetic-gi/(speed*ratio**2))
    zero(gradient-gi*speed/ratio**2)
    zero(gradient/kinetic-speed**2)


def test_algebraic_relative_spring_does_not_change_nondegenerate_full_principal_cones():
    kg, kf, cg, cf = sp.symbols("Kg Kf cg cf", positive=True)
    frequency, momentum, scale, spring = sp.symbols("omega k scale spring", real=True)
    kinetic, gradient = sp.diag(kg, kf), sp.diag(kg*cg**2, kf*cf**2)
    matrix = scale**2*(frequency**2*kinetic-momentum**2*gradient)-spring*sp.Matrix([[1, -1], [-1, 1]])
    leading = sp.Poly(matrix.det(), scale).coeff_monomial(scale**4)
    zero(leading-kg*kf*(frequency**2-cg**2*momentum**2)*(frequency**2-cf**2*momentum**2))


def test_mixed_zero_denominator_really_changes_auxiliary_constraint_instead_of_an_inverse():
    p = sp.Symbol("P", nonzero=True, real=True)
    ge, gf, gu = sp.symbols("gamma_e gamma_f gamma_u", real=True)
    lag = -p*((ge-gu)**2-(gf-gu)**2)/4
    zero(sp.diff(lag, gu, 2))
    zero(sp.diff(lag, gu)-p*(ge-gf)/2)
    assert sp.diff(lag, gu) != 0


def test_finite_three_slice_scale_test_is_strict_and_does_not_identify_two_clocks():
    r = Q(1, 4)  # Actual middle ratio to the fixed algebraic cap.
    assert r*(1+Q(1)**2)**2 == 1  # Equality alone is not the strict test.
    assert r*(1+Q(2)**2)**2 > 1
    left, middle, right = Q(-5, 7), Q(1, 11), Q(9, 10)
    left_weight, right_weight = (right-middle)/(right-left), (middle-left)/(right-left)
    assert left_weight > 0 and right_weight > 0
    assert left_weight+right_weight == 1
    # Both Einstein-scale log endpoint gaps >0 force the chord above f0.
    assert left_weight*Q(1, 7)+right_weight*Q(2, 9) > 0


def test_ordered_finite_Einstein_proper_interval_does_not_evade_the_concavity_argument():
    s = sp.Symbol("s", real=True)
    physical_time = s/(1-s**2)
    scale = 1/(1-s**2)
    zero(sp.diff(physical_time, s)-(1+s**2)/(1-s**2)**2)
    log_scale_second = sp.diff(sp.log(scale), s, 2)
    zero(log_scale_second-2*(1+s**2)/(1-s**2)**2)
    # s in (-1,1) covers all physical times with positive lapse. A scale
    # growing at both ends violates concavity even in this finite clock.
    assert log_scale_second.subs(s, 0) == 2


def test_nonnegative_weights_and_subluminal_cones_cannot_support_positive_null_stress():
    for pe, pf, ce, cf in ((Q(1), Q(2), Q(1), Q(1)),
                           (Q(1), Q(0), Q(1), Q(1, 2)),
                           (Q(0), Q(3), Q(1, 4), Q(1)),
                           (Q(2), Q(3), Q(1, 2), Q(2, 3))):
        null_stress = 2*(pe*(ce-1)+pf*(cf-1))
        assert null_stress <= 0
        if null_stress == 0:
            assert pe == 0 or ce == 1
            assert pf == 0 or cf == 1


def test_negative_link_with_nonpositive_endpoint_has_no_real_regular_Friedmann_point():
    magnitude, ratio, g = sp.symbols("magnitude ratio G", positive=True)
    minus_endpoint = sp.Symbol("minus_endpoint", nonnegative=True)
    h_squared = -2*(magnitude*ratio**3+minus_endpoint)/(3*g)
    assert h_squared.is_negative is True


def test_no_Hubble_division_is_needed_at_a_zero_Hubble_slice():
    p, r, c, g = sp.symbols("p R c G", nonzero=True, real=True)
    b, dh = sp.symbols("b DH", real=True)
    lapse_at_h0 = -2*p*r**3-2*b
    space_at_h0 = 2*g*dh-2*p*r**3/c-2*b
    zero((space_at_h0-lapse_at_h0)/(2*g)-(dh-p*r**3*(1/c-1)/g))


def test_all_zero_links_have_no_auxiliary_tensor_inverse():
    g, f, u = sp.symbols("gamma_g gamma_f gamma_u", real=True)
    pg, pf = sp.symbols("Pg Pf", real=True)
    lag = -(pg*(g-u)**2+pf*(f-u)**2)/4
    zero(sp.diff(lag, u, 2).subs({pg: 0, pf: 0}))


def test_noncommuting_traceful_coframe_jets_recover_full_vacuum_FP_invariant():
    z = sp.Symbol("z", real=True)
    ejet = sp.Matrix([[2, 0, 0, 0], [0, 1, 2, 0], [0, 2, 3, 1], [0, 0, 1, -1]])
    vjet = sp.Matrix([[1, 0, 0, 0], [0, 2, 0, 1], [0, 0, -2, 1], [0, 1, 1, 3]])
    assert ejet*vjet != vjet*ejet
    rg, rf, pg, pf = sp.Integer(2), sp.Integer(3), sp.Integer(2), sp.Integer(-6)
    pe, pv = pg/rg, pf/rf
    total = pe+pv
    e, v = (sp.eye(4)+z*ejet/2)/rg, (sp.eye(4)+z*vjet/2)/rf
    u = (pg*e+pf*v)/total
    # On the zero-source stationary branch this is the literal source-
    # eliminated potential, including BOTH endpoint terms.
    lag = sp.cancel(model.potential_density(e, v, u, b=-3*total, pg=pg, pf=pf,
                                             bg=-pg*rg**3, bf=-pf*rf**3))
    coefficient = sp.diff(lag, z, 2).subs(z, 0)/2
    delta = ejet-vjet
    expected = pe*pv*(sp.trace(delta)**2-sp.trace(delta**2))/(4*total)
    zero(coefficient-expected)
    assert sp.trace(delta) == 1 and sp.trace(delta**2) == 53
    assert coefficient == -26


def test_primary_background_module_bridges_to_literal_covariant_dictionary():
    d = background.derive()
    g, r, c, h, dh, p, b = (d[key] for key in ("G", "R", "c", "H", "D_H", "p", "b"))
    zero(d["lapse_residual"]-(3*g*h**2-2*p*r**3-2*b))
    zero(d["spatial_residual"]-(g*(2*dh+3*h**2)-2*p*r**3/c-2*b))
    zero(d["physical_tensor_kinetic"]-g/(c*r**2))
    zero(d["physical_tensor_gradient"]-g*c/r**2)


def test_actual_Einstein_Bianchi_identity_factors_without_dividing_a_Hubble_rate():
    p, b, hi, hu = sp.symbols("p b Hi Hu", real=True)
    ratio, speed = sp.symbols("R c", positive=True)
    density = 2*p*ratio**3+2*b
    pressure = -2*p*ratio**3/speed-2*b
    # D_i R includes the actual lapse conversion n_u/n_i=R/c.
    di_ratio = ratio*(ratio*hu/speed-hi)
    continuity = sp.diff(density, ratio)*di_ratio+3*hi*(density+pressure)
    zero(continuity-6*p*ratio**3*(ratio*hu-hi)/speed)
    zero(continuity.subs(p, 0))  # A disconnected metric supplies no clock lock.
    zero(continuity.subs(hi, ratio*hu))
    assert continuity.subs({hi: 0, hu: 1, p: -1, ratio: 2, speed: 3}) != 0


def test_Bianchi_clock_lock_and_raychaudhuri_give_the_undivided_physical_equation():
    g, ratio, speed = sp.symbols("G R c", positive=True)
    p, hu, hup = sp.symbols("p Hu Hup", real=True)
    rp = ratio*(1-speed)*hu
    di_hi = ratio*(rp*hu+ratio*hup)/speed
    ray_residual = di_hi-p*ratio**3*(1/speed-1)/g
    physical = g*(hup+(1-speed)*hu**2)/ratio**2-p*(1-speed)/ratio
    zero(g*speed*ray_residual/ratio**4-physical)


def test_sum_of_exact_Einstein_and_u_residuals_is_the_monotonicity_equation():
    kg, kf, cg, cf = sp.symbols("Kg Kf cg cf", positive=True)
    pg, pf, h, hp, nh = sp.symbols("Pg Pf Hu Hup nh", real=True)
    eg = kg*(hp+(1-cg)*h**2)-pg*(1-cg)
    ef = kf*(hp+(1-cf)*h**2)-pf*(1-cf)
    auxiliary = pg*(cg-1)+pf*(cf-1)-nh/2
    k, kp = kg+kf, -2*h*(kg*(1-cg)+kf*(1-cf))
    zero(2*(eg+ef-auxiliary)-(2*k*hp-h*kp+nh))
    assert k.is_positive is True


def test_normalized_Hubble_monotonicity_has_the_actual_physical_time_factor():
    k = sp.Symbol("K", positive=True)
    h, hp, kp, nh = sp.symbols("H Hp Kp nh", real=True)
    normalized_derivative = sp.diff(h/sp.sqrt(k), h)*hp+sp.diff(h/sp.sqrt(k), k)*kp
    zero(normalized_derivative-(2*k*hp-h*kp)/(2*k**sp.Rational(3, 2)))
    zero(normalized_derivative.subs(hp, (h*kp-nh)/(2*k))+nh/(2*k**sp.Rational(3, 2)))
    assert normalized_derivative.subs({h: 0, hp: -3, kp: 17, k: 4}) == -sp.Rational(3, 2)


def test_positive_bounce_acceleration_has_a_strict_source_residual_for_any_positive_K():
    k, acceleration = sp.symbols("K acceleration", positive=True)
    nh = sp.Symbol("nh", nonnegative=True)
    # At the physical bounce H=0, exact equations would set this to zero.
    required_residual = 2*k*acceleration+nh
    assert required_residual.is_positive is True
    # This is a necessary correction size, not an omitted-operator estimate.
    assert required_residual.subs({k: 3, acceleration: 4, nh: sp.Rational(1, 100)}) == sp.Rational(2401, 100)


def test_primary_monotonicity_module_handles_each_active_link_count_with_actual_clock():
    for active_count in (1, 2):
        d = monotonic.derive(active_count)
        h, hp, nh = d["H_u"], d["H_u_prime"], d["n_h"]
        k = sum(link["G"]/link["R"]**2 for link in d["links"])
        kp = sum(-2*link["G"]*(1-link["c"])*h/link["R"]**2 for link in d["links"])
        zero(d["K"]-k)
        zero(d["K_prime"]-kp)
        zero(d["combined_residual"]-(2*k*hp-h*kp+nh))
        zero(d["normalized_Hubble_prime"]-hp/sp.sqrt(k)+h*kp/(2*k**sp.Rational(3, 2)))
        for link in d["links"]:
            zero(link["interaction_balance"]-6*link["p"]*link["R"]**3*(link["R"]*h-link["H_i"])/link["c"])
        assert k.is_positive is True


def test_primary_bounce_residual_bound_distinguishes_strict_obstruction_and_no_verdict():
    bounded = monotonic.bounce_residual_budget(3, 4, sp.Rational(1, 100), sp.Rational(1, 2), 7)
    assert bounded["combined_residual_lower_bound"] == sp.Rational(2051, 100)
    assert bounded["status"] == "POSITIVE_ACTUAL_EQUATION_RESIDUAL"
    uncontrolled = monotonic.bounce_residual_budget(3, 4, 0, 4, 7)
    assert uncontrolled["combined_residual_lower_bound"] == -4
    assert uncontrolled["status"] == "INCONCLUSIVE"


def test_conditional_endpoint_bound_is_derived_from_the_identified_CD_physical_curve():
    t = sp.Symbol("T", real=True)
    tau, length = sp.symbols("tau L", positive=True)
    scale = (1+(t/tau)**2)**2
    h = sp.diff(scale, t)/scale
    zero(h.subs(t, length*tau)-4*length/(tau*(1+length**2)))
    zero(h.subs(t, -length*tau)+4*length/(tau*(1+length**2)))
    for selected in (sp.Rational(1, 7), sp.Rational(1, 2), sp.Integer(1), sp.Integer(3)):
        threshold = sp.cancel(tau*h.subs(t, selected*tau))
        assert monotonic.cd_endpoint_test(selected, 0)["necessary_normalized_endpoint_error"] == threshold
        assert monotonic.cd_endpoint_test(selected, threshold)["status"] == "INCONCLUSIVE"


def test_negative_link_expanding_control_solves_every_actual_coframe_component():
    # An independent q<0 point of the local analytic family, not a claim of
    # stable flat massive spectrum: 5/6<A<1, q(A-1)>0, positive lapse/source.
    a, q = sp.Rational(9, 10), sp.Integer(-1)
    n = a/(6*a-5)
    rho = 12*q*(a-1)/a
    hr2 = 2*q*(a**3-1)/3
    hrprime = -6*q*a**3*(a-1)/(6*a-5)
    u = sp.diag(n, a, a, a)
    source = matter.canonical_source(u, [n*sp.sqrt(2*rho), 0, 0, 0])
    einstein = sp.diag(3*hr2, *([-2*hrprime-3*hr2]*3))
    out = model.euler_maps(sp.eye(4), sp.eye(4), u, pg=q, pf=q, b=-6*q,
                           bg=-q, bf=-q, g=1, f=1, epsilon=1,
                           matter_gradient=source["J_u"], einstein_g=einstein, einstein_f=einstein)
    assert n == sp.Rational(9, 4) and rho == sp.Rational(4, 3) and hr2 > 0
    for key in ("E_e", "E_v", "E_u"):
        zero(out[key])
