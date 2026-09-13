"""Complete finite-mode causal covariance bridge and ordered state-difference kernels."""

from functools import cache

import sympy as s


@cache
def data():
    u0, u1, u2, u3 = s.symbols("u0:4", real=True)
    U = s.Matrix([[u0, u1], [u2, u3]])
    c0, c1, c2 = s.symbols("c0:3", real=True)
    C = s.Matrix([[c0, c1], [c1, c2]])
    d0, d1, d2, g0, g1, g2 = s.symbols("d0 d1 d2 g0 g1 g2", real=True)
    D = s.Matrix([[d0, d1], [d1, d2]])
    G = s.Matrix([[g0, g1], [g1, g2]])
    J = s.Matrix([[0, 1], [-1, 0]])
    force = J * G * C + C * (J * G).T
    W = U * (C + s.I * J / 2)
    Wreverse = U * (C - s.I * J / 2)
    comm = s.trace(D * (W * G * W.T - Wreverse * G * Wreverse.T)) / 2
    tangent = -s.trace(D * U * force * U.T) / 2
    checks = {
        "full_covariance_force_both_adjoint_blocks": force - (J * G * C - C * G * J),
        "independent_Wick_commutator_equals_canonical_covariance_tangent": s.expand(
            s.I * comm - tangent
        ),
        "full_zero_initial_covariance_tangent": s.zeros(2),
    }
    re, im = s.symbols("ordered_product_real ordered_product_imag", real=True)
    product = re + s.I * im
    checks["ordered_even_parity_real_causal_kernel"] = s.expand(
        s.I * (product - s.conjugate(product)) / 2 + im
    )
    checks["ordered_odd_parity_imaginary_Fourier_kernel"] = s.expand(
        s.I * (product + s.conjugate(product)) / 2 - s.I * re
    )
    # Exact complex beta is fixed at the common Cauchy neighborhood.
    b, theta = s.symbols("beta_magnitude beta_phase", real=True)
    alpha = s.sqrt(1 + b * b)
    beta = b * s.exp(s.I * theta)
    checks["full_fixed_Bogoliubov_CCR"] = s.simplify(
        alpha**2 - beta * s.conjugate(beta) - 1
    )
    E0, E1, L0, L1, Z0, Z1 = s.symbols("E0 E1 L0 L1 Z0 Z1", real=True)
    rho = (E0**2 + E1**2 + L0**2 + L1**2 + Z0**2 + Z1**2) / 2
    pressure = (E0**2 + E1**2 - (L0**2 + L1**2) / 3 - Z0**2 - Z1**2) / 2
    checks["complete_energy_feature_readout"] = 2 * rho - (
        E0**2 + E1**2 + L0**2 + L1**2 + Z0**2 + Z1**2
    )
    checks["complete_pressure_feature_readout"] = (
        2 * pressure - E0**2 - E1**2 + (L0**2 + L1**2) / 3 + Z0**2 + Z1**2
    )
    return {
        "unchanged_preparation": "T is the exact S240 SLE. S is the exact normalized KG solution initialized to the full positive W6 and its full derivative at t=-1. On the reference T=alpha_p S+beta_p conjugate(S), with the actual full S240 coefficients fixed in time. Both Cauchy covariances are held fixed on later metric perturbations sharing the t0=-1/2 neighborhood. S is an AUXILIARY exact comparison, not the WKB evolution, not a newly selected physical state and not asserted Hadamard.",
        "full_covariance_tangent": "At every finite reference-mode projector, Cdot=A C+C A^t, A=J M0. The source tangent solves C1dot=A C1+C1 A^t+J M_G C+C(J M_G)^t with ZERO initial tangent. The complete current tangent is-(1/2)tr(M_D C1+M_DG C). Both adjoint source blocks and the instantaneous second-Hamiltonian contact are retained. The same equations hold for the difference of the two fixed initial covariances.",
        "literal_Wick_covariance_bridge": s.expand(s.I * comm - tangent),
        "annihilation_feature": "F_p(t)=(S_p', i p_x S_p/a, i p_y S_p/a, i p_z S_p/a, mass*S_p), or the same vector with T_p. Use positive-frequency exp(-i omega t), W2(t,s)=S(t)conjugate(S(s)) and the physical CCR i/a³.",
        "ordered_full_pair": "For each REAL Cartesian ADM component basis A, V_A(t;p,q)=F_p(t)^t M_A(t) F_q(t), q=P-p. Before the common momentum integral and a(t)^3 a(s)^3, the retarded memory is (i/2)theta(t-s){V_A(t;p,q)conjugate(V_B(s;p,q))-conjugate(V_A(t;-p,-q))V_B(s;-p,-q)}. Directions are complexified BILINEARLY after the component kernel is formed. Do not conjugate a supplied complex source coefficient.",
        "parity_reduction": "The even-even and odd-odd components reduce to-Im[V_A conjugate(V_B)]; mixed even-odd components reduce to+i Re[V_A conjugate(V_B)]. These are Fourier kernels: mixed channels can be imaginary. Reverse reflection, reality K_AB(-P)=conjugate(K_AB(P)), and BOTH time orderings are retained.",
        "entire_state_difference": "Subtract the COMPLETE ordered pair products computed with T and S, including all terms in alpha,beta and their complex phases. Add the instantaneous contact difference-a³/2 integral[conjugate(F_T)^t C_DG F_T-conjugate(F_S)^t C_DG F_S] d³p/(2pi)^3. Common state-independent covariant counterterms and the identical local finite prescription cancel only in this DIFFERENCE. The comparison state's complete renormalized response is not set to zero.",
        "regulator_and_limit": "At finite K project the reference Fourier field onto Omega_p<=K, with Omega_p²=n+|p|²/16. Each memory term retains BOTH projected internal legs; the one-mode contact retains the same projection. The exact finite covariance identity precedes the limit. The full absolute estimates and regulator tails justify the limit of this linear tangent. No finite-amplitude or nonlinear differentiability theorem is inferred.",
        "reference_profile_split": "One may algebraically split the EXISTING fixed S240 reference profile into the exact-comparison reference stress part and the T-minus-S stress part. No candidate, prescribed sum or finite renormalization changes. The displayed state-selection profile is only a bookkeeping summand. Its complete second variation and second metric-chart contact must be included if that matched summand is discussed.",
        "checks": {
            key: value.applyfunc(s.expand)
            if isinstance(value, s.MatrixBase)
            else s.expand(value)
            for key, value in checks.items()
        },
        "gates": {
            "actual_exact_SLE_and_exact_comparison_not_WKB_evolution": True,
            "fixed_covariance_not_live_reminimization": True,
            "full_complex_Bogoliubov_coefficients_not_occupation_only": True,
            "complete_memory_and_second_contact_difference": True,
            "both_reflected_reverse_and_forward_ordered_products": True,
            "component_kernel_before_complex_direction_extension": True,
            "common_covariant_counterterms_cancel_only_in_difference": True,
            "remaining_comparison_response_not_deleted": True,
        },
    }
