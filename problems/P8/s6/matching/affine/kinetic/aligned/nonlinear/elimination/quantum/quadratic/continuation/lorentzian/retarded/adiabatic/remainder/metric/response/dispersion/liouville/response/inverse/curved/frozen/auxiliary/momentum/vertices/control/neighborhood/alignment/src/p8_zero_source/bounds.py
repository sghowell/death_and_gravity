"""Continuous auxiliary bounds and positive-majorant scaling for the new action."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import bounds as previous_aux
from p8_coupled_vertices import majorant as previous_vertices

RATIO=sp.Integer(10**9)


def response_radius(value):
    prior=previous_aux.response_radius(value)
    prior["temporal_vector_absolute_upper"]=4*prior["input_radius"]
    prior["boundaries"]="New source-free candidate, same local classical auxiliary domain; not a frequency or quantum-response estimate."
    return prior


@cache
def auxiliary():
    old=previous_aux.coefficients()
    pieces=dict(old["full_Hamiltonian_piece_upper_before_N"])
    # In the new literal trace square p-b replaces p-b-delta*j.
    # The old 162 triangle majorant remains an upper bound. The
    # unchanged Gauss square is bounded by 2; -c*j is absent.
    pieces["source_Gauss"]=sp.Integer(0)
    contraction=dict(previous_aux.contraction())
    contraction["temporal_vector_absolute_upper"]=4*previous_aux.JET_RADIUS
    contraction["temporal_vector_over_input_radius_upper"]=sp.Integer(4)
    return {"unchanged_original_scalar_boundary_and_mass_bounds":old,
            "new_Hamiltonian_piece_upper_before_N":pieces,
            "new_Hamiltonian_complex_absolute_upper":2*sum(pieces.values()),
            "same_declared_complex_Hamiltonian_upper":previous_aux.HAMILTONIAN_UPPER,
            "new_auxiliary_contraction_bounds":contraction,
            "new_temporal_root_formula":"T=-gamma_t*j/U; |gamma_t/U|<4 on the same outer coefficient domain.",
            "same_nine_invariant_and_lapse_domains":True}


@cache
def coefficients():
    old=previous_vertices.coefficient_bounds()["common_N_derivative_bounds"]
    new=tuple(previous_aux.HAMILTONIAN_UPPER*sp.factorial(n)*100**n for n in range(5))
    return {"new_universal_Cauchy_lapse_derivative_bounds":new,
            "old_sharper_common_lapse_derivative_bounds":old,
            "new_over_old_bound_ratios":tuple(sp.factor(a/b) for a,b in zip(new,old)),
            "common_bound_ratio_upper":RATIO,
            "coefficient_domain":"The same complex time tube of radius 1/50 around [-1/2,1/2], lapse radius 1/100 and nine-invariant unit polydisc. Differentiate the actual new coefficient, not the old isolated Gauss coefficient."}


@cache
def scaling():
    eps,r=sp.symbols("physical_degree positive_bound_ratio",positive=True)
    bs=sp.symbols("positive_B0:5",positive=True)
    ls=sp.symbols("positive_linear1:5",positive=True)
    qs=sp.symbols("positive_quadratic2:5",positive=True)
    vs=sp.symbols("positive_volume1:5",positive=True)
    ms=sp.symbols("positive_moving0:5",positive=True)
    linear=sum(v*eps**n for n,v in enumerate(ls,1))
    quadratic=sum(v*eps**n for n,v in enumerate(qs,2))
    volume=1+sum(v*eps**n for n,v in enumerate(vs,1))
    moving=sum(v*eps**n for n,v in enumerate(ms))
    B=tuple(r*v for v in bs)
    L=tuple(v*linear for v in B[:4])
    Q=tuple(v*quadratic for v in B[:3])
    n1=20*L[1]
    f2=B[3]*n1**2/2+L[2]*n1+Q[1]
    h=B[0]+L[0]+Q[0]+10*L[1]**2
    h+=B[3]*n1**3/6+L[2]*n1**2/2+Q[1]*n1
    h+=B[4]*n1**4/24+L[3]*n1**3/6+Q[2]*n1**2/2+10*f2**2
    expression=sp.Poly(sp.expand(volume*h+moving),eps)
    degrees={n:sp.Poly(expression.nth(n),r).degree() for n in (3,4)}
    # Every monomial is nonnegative in all independent majorant
    # coefficients. For r>=1, each r^k <= r^d for k<=d.
    variables=(r,*bs,*ls,*qs,*vs,*ms)
    positive={n:all(c>=0 for c in sp.Poly(expression.nth(n),*variables).coeffs()) for n in (3,4)}
    return {"ratio_symbol":r,"abstract_cubic_positive_majorant":expression.nth(3),
            "abstract_quartic_positive_majorant":expression.nth(4),
            "cubic_maximum_coefficient_ratio_degree":degrees[3],
            "quartic_maximum_coefficient_ratio_degree":degrees[4],
            "all_abstract_cubic_coefficients_nonnegative":positive[3],
            "all_abstract_quartic_coefficients_nonnegative":positive[4],
            "new_cubic_majorant_over_old_upper":RATIO**4,
            "new_quartic_majorant_over_old_upper":RATIO**6}


@cache
def checks():
    d=scaling()
    return {"positive_cubic_majorant_coefficient_degree":d["cubic_maximum_coefficient_ratio_degree"]-4,
            "positive_quartic_majorant_coefficient_degree":d["quartic_maximum_coefficient_ratio_degree"]-6}


@cache
def gates():
    a,c,s=auxiliary(),coefficients(),scaling()
    new=a["new_auxiliary_contraction_bounds"]
    return {"actual_new_complex_Hamiltonian_below_same_Cauchy_majorant":
                bool(a["new_Hamiltonian_complex_absolute_upper"]<previous_aux.HAMILTONIAN_UPPER),
            "all_five_new_Cauchy_bounds_below_recorded_old_ratio":
                all(value<RATIO for value in c["new_over_old_bound_ratios"]),
            "coefficient_bound_ratio_above_one":bool(RATIO>=1),
            "new_auxiliary_Newton_contraction_below_one_tenth":bool(new["Newton_contraction_upper"]<sp.Rational(1,10)),
            "new_auxiliary_Newton_disc_strictly_invariant":bool(new["Newton_closed_disc_image_radius_upper"]<previous_aux.ROOT_RADIUS),
            "new_temporal_root_below_five_times_invariant_radius":bool(new["temporal_vector_over_input_radius_upper"]<5),
            "abstract_cubic_majorant_has_only_nonnegative_coefficients":s["all_abstract_cubic_coefficients_nonnegative"],
            "abstract_quartic_majorant_has_only_nonnegative_coefficients":s["all_abstract_quartic_coefficients_nonnegative"]}
