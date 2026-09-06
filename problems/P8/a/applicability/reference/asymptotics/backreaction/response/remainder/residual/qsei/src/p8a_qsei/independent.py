"""Fraction-only jet recurrence and Fourier/proper-time constant assembly."""

from fractions import Fraction as F
from math import comb


def replay(prior_a7_report):
    if prior_a7_report.get("claim") != "P8-A.7":
        raise ValueError("The independent mode replay requires the pinned A.7 input")
    data = prior_a7_report["derived_finite_geometry"]
    b = list(map(F, data["U_derivative_bounds_through_5_per_delta"][:4]))
    cap = F(data["delta_bar"])
    duration = F(3)
    # Polynomial coefficient arrays in delta. Differentiate U*b by the
    # Leibniz convolution before evaluating a bound at the amplitude cap.
    derivatives = [[F(2)]]
    products = []
    for order in range(4):
        product = [F(0)]*(order+2)
        for potential_order in range(order+1):
            for degree, coefficient in enumerate(derivatives[order-potential_order]):
                product[degree+1] += comb(order, potential_order)*b[potential_order]*coefficient
        products.append(product)
        if order < 3:
            derivatives.append([duration*coefficient for coefficient in product])
    per_delta = lambda polynomial: sum(coefficient*cap**(degree-1)
                                      for degree, coefficient in enumerate(polynomial) if degree > 0)
    naive = [per_delta(polynomial) for polynomial in derivatives[1:]]
    product_norms = list(map(per_delta, products))
    inverse_frequency = [2*duration*b[0]]
    inverse_frequency += [(product_norms[j-1]+duration*product_norms[j])/2 for j in range(1, 4)]
    h = [F(2), 4+cap*b[0], cap*b[1]+4*(4+cap*b[0])]
    error = [inverse_frequency[j]+inverse_frequency[j+1]
             +sum(comb(j, ell)*h[ell]*inverse_frequency[j-ell] for ell in range(j+1))
             for j in range(3)]
    infrared = F(33)*b[0]
    p1, p0 = F(5, 6), F(25, 36)
    physical_h, hdot = F(8, 31), F(132, 961)
    f = [p0/2, p1+F(3, 2)*physical_h*p0,
         1+2*physical_h*p1+(F(3, 4)*physical_h**2+F(3, 2)*hdot)*p0]
    weighted = physical_h*p1+(hdot+physical_h**2/2)*p0
    root0 = f[2]+F(3, 2)*weighted
    error_root = 2*(error[0]*f[2]+2*error[1]*f[1]+error[2]*f[0])+2*infrared*f[0]
    delta = F(1, 10**14)
    coefficient = (root0+delta*error_root)**2
    if not 0 < delta <= cap or cap*b[0]*duration**2 > 1 or coefficient >= 5:
        raise ValueError("The independent finite QSEI coefficient check failed")
    if any(value < 0 for polynomial in derivatives+products for value in polynomial):
        raise ValueError("A purported monotone mode majorant has a negative coefficient")
    return {"naive_b_derivatives_per_delta": list(map(str, naive)),
            "product_Ub_derivatives_per_delta": list(map(str, product_norms)),
            "b_inverse_frequency_constants": list(map(str, inverse_frequency)),
            "conformal_Hc_derivative_caps": list(map(str, h)),
            "ultraviolet_error_derivatives_per_delta": list(map(str, error)),
            "infrared_error_per_delta": str(infrared),
            "sampler_norms_per_proper_second_derivative": list(map(str, f)),
            "weighted_first_sampler_norm": str(weighted),
            "auxiliary_flat_kernel_root_coefficient": str(root0),
            "mode_error_root_coefficient_per_delta": str(error_root),
            "maximum_QSEI_coefficient": str(coefficient),
            "strict_rounding_margin": str(5-coefficient),
            "nonnegative_mode_majorant_polynomials": [[str(c) for c in polynomial]
                                                       for polynomial in derivatives+products],
            "UV_frequency_integral": "1/3",
            "IR_frequency_measure": "1/2",
            "IR_Parseval_is_full_complex_not_half": True}
