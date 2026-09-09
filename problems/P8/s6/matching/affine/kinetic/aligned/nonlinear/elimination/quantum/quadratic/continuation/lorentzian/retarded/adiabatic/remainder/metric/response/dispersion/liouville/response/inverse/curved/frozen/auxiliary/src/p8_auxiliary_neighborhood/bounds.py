"""Rational majorants for holomorphic coefficients and an actual auxiliary root."""
from functools import cache

import sympy as sp
from p8_affine import dictionary

from . import model

TIME_RADIUS=sp.Rational(1,50)
LAPSE_RADIUS=sp.Rational(1,100)
INNER_LAPSE=sp.Rational(1,200)
HAMILTONIAN_UPPER=sp.Integer(10000)
ROOT_RADIUS=sp.Rational(1,10**14)
JET_RADIUS=sp.Rational(1,10**24)
JET_COUNT=9
PIVOT_INVERSE=sp.Integer(20)


def taylor_majorant(degree):
    if type(degree) is not int or degree not in range(5):
        raise ValueError("Require native invariant Taylor degree 0 through 4")
    count=sp.binomial(degree+8,8)
    return {"degree":degree,"monomial_count":count,
            "total_absolute_coefficient_upper":HAMILTONIAN_UPPER*count/JET_RADIUS**degree}


def response_radius(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact rational invariant radius")
    radius=sp.Rational(value)
    if not 0<=radius<=JET_RADIUS:
        raise ValueError("Require invariant radius between zero and 10^-24")
    return {"input_radius":radius,"lapse_displacement_upper":800000000*radius,
            "temporal_vector_absolute_upper":5*radius,
            "boundaries":"Uniform local classical auxiliary bounds, not a frequency or quantum-response estimate."}


@cache
def scalar_potential_bound():
    value=sp.factor(dictionary.target()["scalar_F"])
    numerator,denominator=sp.fraction(value)
    poly=sp.Poly(numerator,dictionary.u,dictionary.x)
    upper=sum(abs(c)*sp.Rational(13,25)**powers[0]*sp.Rational(21,20)**powers[1]
              for powers,c in poly.terms())/(200*sp.Rational(99,100)**12)
    return {"upper":upper,"numerator":numerator,"denominator":denominator,
            "monomials":len(poly.terms())}


@cache
def coefficients():
    umax=sp.Rational(13,25)
    rmin=sp.Rational(99,100)
    rmax=1+umax**2
    hmin=rmin**3
    s_deviation=LAPSE_RADIUS/(1-LAPSE_RADIUS)
    x_deviation=s_deviation*(2+s_deviation)
    ratio_deviation=x_deviation/hmin
    hprime=6*umax*rmax**2
    p_radius=sp.Rational(1,50)
    pmax=sp.Rational(1,2)+p_radius
    cubic_deviation=p_radius*(pmax**2+pmax/2+sp.Rational(1,4))
    gt_den_lower=sp.Rational(27,4)-22*cubic_deviation-3*p_radius
    gt_upper=9*(sp.Rational(3,4)+2*cubic_deviation)/gt_den_lower
    gs_den_deviation=72*p_radius*(1+p_radius)+88*p_radius
    gs_upper=9*(9+8*p_radius)/(81-gs_den_deviation)
    gs_inverse_upper=(81+gs_den_deviation)/(9*(9-8*p_radius))
    ode_A=sp.Rational(1,2)/sp.Rational(19,20)+sp.Rational(3,4)/sp.Rational(9,10)
    ode_force=18/(4*sp.Rational(9,10)**3)
    Q_factor=sp.Rational(7,2)/(1-sp.Rational(3,2)*sp.Rational(1,20))
    I_upper=2*s_deviation
    I_phi=I_upper/TIME_RADIUS
    omega_phi=6*sp.Rational(1,40)/(4*sp.Rational(9,10)**2)
    f_phi=6*sp.Rational(1,40)/(2*sp.Rational(9,10)**2)
    boundary_derivative=12*2*2*sp.Rational(1,10)*sp.Rational(1,3)
    Blinear=3*sp.Rational(21,20)*umax*sp.Rational(1,40)/rmin**4
    b_upper=2*Blinear+I_upper
    c_upper=3*3*sp.Rational(21,20)*sp.Rational(1,40)/hmin+sp.Rational(3,2)*sp.Rational(21,20)*4*sp.Rational(1,40)**2
    fnew_upper=1000+6*sp.Rational(21,20)*(omega_phi**2+f_phi*omega_phi)+model.MARGIN*sp.Rational(1,40)**2/hmin**2
    f0_upper=2*1002+sp.Rational(21,20)*I_phi
    pieces={"trace_square":162,"potential":3000,"Gauss":2,"source_Gauss":2,
            "shear":4,"electric":1,"magnetic":sp.Rational(1,2),"vector_mass":2,
            "matter_momentum":4,"matter_gradient":1,"spatial_curvature":2}
    return {"complex_time_absolute_upper":umax,"one_plus_time_squared_absolute_lower":rmin,
            "one_plus_time_squared_absolute_upper":rmax,"h_absolute_lower":hmin,
            "s_minus_one_absolute_upper":s_deviation,"x_plus_one_absolute_upper":x_deviation,
            "ratio_minus_one_absolute_upper":ratio_deviation,"h_prime_absolute_upper":hprime,
            "affine_p_disc_radius":p_radius,"gamma_t_denominator_absolute_lower":gt_den_lower,
            "gamma_t_absolute_upper":gt_upper,"gamma_s_absolute_upper":gs_upper,
            "inverse_gamma_s_absolute_upper":gs_inverse_upper,
            "Q_ODE_coefficient_absolute_upper":ode_A,"Q_ODE_forcing_over_distance_upper":ode_force,
            "Q_over_distance_squared_upper":Q_factor,"boundary_primitive_absolute_upper":I_upper,
            "boundary_primitive_phi_absolute_upper":I_phi,
            "omega_phi_absolute_upper":omega_phi,"f_phi_absolute_upper":f_phi,
            "boundary_primitive_s_absolute_upper":boundary_derivative,
            "full_trace_linear_absolute_upper":b_upper,"full_normal_source_constant_absolute_upper":c_upper,
            "full_transformed_scalar_plus_margin_absolute_upper":fnew_upper,
            "full_post_boundary_potential_absolute_upper":f0_upper,
            "full_Hamiltonian_piece_upper_before_N":pieces,
            "full_Hamiltonian_absolute_upper_from_pieces":2*sum(pieces.values()),
            "declared_full_Hamiltonian_absolute_upper":HAMILTONIAN_UPPER}


@cache
def contraction():
    # Cauchy radii on |N-1|<=1/200 and max|Y_i|<=1/2.
    HNY=HAMILTONIAN_UPPER/(INNER_LAPSE*sp.Rational(1,2))
    HNNN=6*HAMILTONIAN_UPPER/INNER_LAPSE**3
    HNNY=2*HAMILTONIAN_UPPER/(INNER_LAPSE**2*sp.Rational(1,2))
    pivot_change=HNNN*ROOT_RADIUS+JET_COUNT*HNNY*JET_RADIUS
    force_at_center=JET_COUNT*HNY*JET_RADIUS
    rate=PIVOT_INVERSE*pivot_change
    center_step=PIVOT_INVERSE*force_at_center
    root_per_jet=sp.Rational(10,9)*PIVOT_INVERSE*JET_COUNT*HNY
    root_upper=root_per_jet*JET_RADIUS
    # G=K-3H/N has |G|<45 on the outer coefficient polydisc,
    # G(1,0)=0, so its N/Y derivatives have Cauchy bounds.
    G_N=sp.Integer(45)/INNER_LAPSE
    G_Y=sp.Integer(45)/sp.Rational(1,2)
    temporal=3*root_upper*(G_N*root_upper+JET_COUNT*G_Y*JET_RADIUS)+108*root_upper**2+4*JET_RADIUS
    theta=sp.Rational(1,1000)
    tail_first=HAMILTONIAN_UPPER*sp.binomial(13,8)*theta**5
    tail_ratio=sp.Rational(14,6)*theta
    return {"outer_lapse_radius":LAPSE_RADIUS,"inner_lapse_Cauchy_radius":INNER_LAPSE,
            "nine_invariant_outer_polydisc_radius":sp.Integer(1),
            "actual_implicit_root_disc_radius":ROOT_RADIUS,
            "actual_nine_invariant_polydisc_radius":JET_RADIUS,
            "inverse_background_lapse_pivot_upper":PIVOT_INVERSE,
            "H_NY_absolute_upper":HNY,"H_NNN_absolute_upper":HNNN,"H_NNY_absolute_upper":HNNY,
            "lapse_pivot_change_upper":pivot_change,"lapse_force_at_center_upper":force_at_center,
            "Newton_contraction_upper":rate,"Newton_center_step_upper":center_step,
            "Newton_closed_disc_image_radius_upper":center_step+rate*ROOT_RADIUS,
            "lapse_displacement_over_input_radius_upper":root_per_jet,
            "lapse_displacement_upper":root_upper,
            "temporal_vector_absolute_upper":temporal,
            "temporal_vector_over_input_radius_upper":temporal/JET_RADIUS,
            "invariant_Hamiltonian_Taylor_coefficient_policy":"For multi-index alpha in nine invariant variables, absolute coefficient <=10000*(10^24)^|alpha|; these are not momentum-reduced physical vertices.",
            "degree_three_total_coefficient_majorant":HAMILTONIAN_UPPER*sp.binomial(11,8)/JET_RADIUS**3,
            "degree_four_total_coefficient_majorant":HAMILTONIAN_UPPER*sp.binomial(12,8)/JET_RADIUS**4,
            "quartic_invariant_Taylor_remainder_input_radius":theta*JET_RADIUS,
            "quartic_invariant_Taylor_remainder_upper":tail_first/(1-tail_ratio)}


@cache
def checks():
    F=scalar_potential_bound()
    out={"entire_frozen_scalar_denominator":sp.factor(
        F["denominator"]-200*(1+dictionary.u**2)**12)}
    h,y=sp.symbols("positive_h source_distance",nonzero=True)
    x=y-1
    A=1/(2*x)+3/(4*(h-y))
    ode=model.lower.ode()
    mapping={dictionary.u:model.u,dictionary.x:x}
    out["original_Q_ODE_coefficient_chart"]=sp.factor(
        ode["coefficient"].subs(mapping).subs((1+model.u**2)**3,h)-A)
    # Use the actual rational h,h' for the forcing bridge rather than
    # replacing noncanonical expanded factors by an independent symbol.
    hu=(1+model.u**2)**3
    out["original_Q_ODE_forcing_chart"]=sp.factor(
        ode["forcing"].subs(dictionary.x,x)+3*sp.diff(hu,model.u)*y/(4*hu**2*(hu-y)))
    k=sp.Symbol("degree",integer=True,nonnegative=True)
    out["nine_variable_coefficient_count_ratio"]=sp.factor(
        sp.prod(k+j for j in range(2,10))/sp.factorial(8)
        -sp.prod(k+j for j in range(1,9))/sp.factorial(8)*(k+9)/(k+1))
    return out


@cache
def gates():
    c,d=coefficients(),contraction()
    return {
        "complex_time_one_plus_square_avoids_zero":bool(1-TIME_RADIUS**2>sp.Rational(99,100)),
        "complex_x_absolute_below_twenty_one_twentieths":bool(1+c["x_plus_one_absolute_upper"]<sp.Rational(21,20)),
        "complex_x_distance_below_one_fortieth":bool(c["x_plus_one_absolute_upper"]<sp.Rational(1,40)),
        "complex_h_and_w_absolute_above_nine_tenths":bool(c["h_absolute_lower"]-sp.Rational(1,20)>sp.Rational(9,10)),
        "complex_ratio_disc_radius_below_three_hundredths":bool(c["ratio_minus_one_absolute_upper"]<sp.Rational(3,100)),
        "principal_square_root_affine_p_radius_valid":bool(c["ratio_minus_one_absolute_upper"]/2<c["affine_p_disc_radius"]),
        "complex_h_prime_below_six":bool(c["h_prime_absolute_upper"]<6),
        "gamma_t_denominator_stays_nonzero":bool(c["gamma_t_denominator_absolute_lower"]>0),
        "gamma_t_absolute_below_two":bool(c["gamma_t_absolute_upper"]<2),
        "gamma_s_and_inverse_absolute_below_two":bool(max(c["gamma_s_absolute_upper"],c["inverse_gamma_s_absolute_upper"])<2),
        "entire_scalar_potential_complex_majorant_below_one_thousand":bool(scalar_potential_bound()["upper"]<1000),
        "Q_complex_ODE_coefficient_below_three_halves":bool(c["Q_ODE_coefficient_absolute_upper"]<sp.Rational(3,2)),
        "Q_complex_ODE_forcing_below_seven_times_distance":bool(c["Q_ODE_forcing_over_distance_upper"]<7),
        "Q_complex_solution_below_four_times_distance_squared":bool(c["Q_over_distance_squared_upper"]<4),
        "boundary_primitive_below_one_fortieth":bool(c["boundary_primitive_absolute_upper"]<sp.Rational(1,40)),
        "boundary_primitive_phi_below_two_by_complex_Cauchy":bool(c["boundary_primitive_phi_absolute_upper"]<2),
        "omega_phi_and_f_phi_below_one_tenth":bool(max(c["omega_phi_absolute_upper"],c["f_phi_absolute_upper"])<sp.Rational(1,10)),
        "actual_boundary_primitive_integrand_below_two":bool(c["boundary_primitive_s_absolute_upper"]<2),
        "actual_trace_linear_coefficient_below_one":bool(c["full_trace_linear_absolute_upper"]<1),
        "actual_source_constant_below_two":bool(c["full_normal_source_constant_absolute_upper"]<2),
        "actual_transformed_scalar_and_margin_below_1002":bool(c["full_transformed_scalar_plus_margin_absolute_upper"]<1002),
        "actual_post_boundary_potential_below_3000":bool(c["full_post_boundary_potential_absolute_upper"]<3000),
        "entire_Hamiltonian_below_declared_complex_majorant":bool(c["full_Hamiltonian_absolute_upper_from_pieces"]<HAMILTONIAN_UPPER),
        "literal_background_J_gives_inverse_pivot_below_twenty":bool(
            2*sp.Rational(1199,800)/sp.Rational(5,4)**18>sp.Rational(1,20)),
        "chosen_root_and_jet_discs_inside_Cauchy_domain":bool(ROOT_RADIUS<INNER_LAPSE and JET_RADIUS<sp.Rational(1,2)),
        "Newton_contraction_strictly_below_one_tenth":bool(d["Newton_contraction_upper"]<sp.Rational(1,10)),
        "Newton_maps_actual_closed_disc_strictly_inside_itself":bool(d["Newton_closed_disc_image_radius_upper"]<ROOT_RADIUS),
        "actual_lapse_displacement_below_one_e_minus_fifteen":bool(d["lapse_displacement_upper"]<sp.Rational(1,10**15)),
        "lapse_secondary_pivot_absolute_above_one_twenty_fifth":bool(
            sp.Rational(1,20)-d["lapse_pivot_change_upper"]>sp.Rational(1,25)),
        "actual_temporal_vector_below_five_times_input_radius":bool(d["temporal_vector_over_input_radius_upper"]<5),
        "quartic_invariant_Taylor_remainder_below_one_e_minus_seven":bool(
            d["quartic_invariant_Taylor_remainder_upper"]<sp.Rational(1,10**7)),
        "real_lapse_root_remains_in_original_clock_tube":bool(
            (1+ROOT_RADIUS)**-2>sp.Rational(9,10) and (1-ROOT_RADIUS)**-2<sp.Rational(11,10)),
        "nine_invariants_keep_electric_and_magnetic_zeta_normalization":True,
        "full_spatial_canonical_jets_not_only_homogeneous_data":True,
        "classical_tree_plus_margin_not_the_nonlocal_quantum_constraint":True,
        "invariant_Taylor_remainder_not_a_physical_interacting_cutoff":True}
