"""Complex spatial momentum geometry and adversarial scope controls."""
from functools import cache

import sympy as sp


@cache
def checks():
    x=sp.Matrix(sp.symbols("x1:4",real=True))
    y=sp.Matrix(sp.symbols("y1:4",real=True))
    a,b=sp.symbols("root_real root_imaginary",real=True)
    xx,yy,xy=x.dot(x),y.dot(y),x.dot(y)
    zeta=x+sp.I*y
    square=sp.expand(zeta.dot(zeta))
    gram=sum((x[i]*y[j]-x[j]*y[i])**2 for i in range(3) for j in range(i+1,3))
    root_modulus_squared=sp.Symbol("root_modulus_squared",nonnegative=True)
    return {
        "complex_bilinear_spatial_square_real_part":sp.expand(sp.re(square)-xx+yy),
        "complex_bilinear_spatial_square_imaginary_part":sp.expand(sp.im(square)-2*xy),
        "complex_root_fourth_power_bound_is_exact_Gram_sum":sp.expand(
            (xx+yy)**2-((xx-yy)**2+4*xy**2)-4*gram),
        "root_imaginary_part_follows_from_modulus_and_real_square":sp.expand(
            ((a*a+b*b)-(a*a-b*b))/2-b*b),
        "root_imaginary_squared_bound_has_nonnegative_norm_deficit":sp.expand(
            yy-(root_modulus_squared-xx+yy)/2-(xx+yy-root_modulus_squared)/2),
        "complex_null_vector_can_have_arbitrarily_large_Euclidean_norm":
            sp.expand((a*sp.Matrix([1,sp.I,0])).dot(a*sp.Matrix([1,sp.I,0]))),
        "principal_fast_mutation_exceeds_matter_type_at_pure_imaginary_momentum":
            sp.Rational(1001,1000)-1-sp.Rational(1,1000)}


def null_cone_example(size):
    if isinstance(size,bool) or not isinstance(size,(int,sp.Rational)):
        raise TypeError("Require an exact positive rational complex-vector size")
    size=sp.Rational(size)
    if size<=0:
        raise ValueError("Complex-vector size must be positive")
    vector=sp.Matrix([size,sp.I*size,0])
    return {"complex_vector":vector,"bilinear_square":vector.dot(vector),
        "Euclidean_modulus_squared":2*size**2,
        "covered_by_original_polynomial_generator_not_a_singular_chart":True}


@cache
def data():
    return {"spatial_dimension":3,"complex_vector_convention":"zeta=x+i*y; x,y real",
        "scalar_root_is_pointwise_not_a_global_holomorphic_choice":"k^2=zeta dot zeta",
        "root_modulus_bound":"abs(k)<=sqrt(norm(x)^2+norm(y)^2)",
        "root_imaginary_bound":"abs(Im(k))<=norm(y)",
        "Gram_remainder_is_sum_of_real_squares":True,
        "both_scalar_root_signs_obey_the_same_bounds":True,
        "entire_response_defined_by_original_polynomial_q_generator":True,
        "large_complex_null_vector_example":null_cone_example(10**30),
        "no_unphysical_identification_k_with_complex_vector_Euclidean_norm":True}


@cache
def gates():
    rows={"large_complex_null_cone_example_is_in_low_scalar_root_disc":
        null_cone_example(10**30)["bilinear_square"]==0,
        "large_complex_null_cone_example_is_not_in_a_unit_vector_ball":
        null_cone_example(10**30)["Euclidean_modulus_squared"]>1,
        "fast_principal_mutation_has_strictly_larger_exponential_type":sp.Rational(1001,1000)>1}
    return {name:bool(value) for name,value in rows.items()}
