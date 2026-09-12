"""Ordered variable-coefficient coordinate inverses and curved reference bound."""

from functools import cache

import sympy as s
from p8_vacuum_affine_flat_scalar_quotient_inverse import resolvent as flat

from . import geometry as g

C0 = flat.C
CMAX = s.Integer(8)
CPRIME = s.Integer(35)
FORWARD = s.Integer(20)
ADJOINT = s.Integer(108)
JMAX = flat.JBOUND


@cache
def data():
    source = s.Matrix([g.wg, g.cg])
    detector = s.Matrix([g.wd, g.cd])
    primal = g.L * source.diff(g.eta, 2) + g.Q * source + g.C * source.diff(g.eta)
    adjoint = (
        g.L.T * detector.diff(g.eta, 2)
        + g.Q.T * detector
        - (g.C.T * detector).diff(g.eta)
    )
    boundary = (
        detector.T * g.L * source.diff(g.eta)
        - detector.diff(g.eta).T * g.L * source
        + detector.T * g.C * source
    )[0]
    identity = s.expand(
        (detector.T * primal - adjoint.T * source)[0] - s.diff(boundary, g.eta)
    )
    delta, T, u, v = s.symbols("elapsed window u v", positive=True)
    simplex = s.integrate(s.integrate(delta - u, (v, 0, u)), (u, 0, delta))
    kernel = C0 * C0 * JMAX * simplex * s.exp(ADJOINT * delta)
    operator = C0 * C0 * JMAX * T**4 * s.exp(ADJOINT * T) / 24
    kap = g.bridge.KAPPA
    checks = {
        "complete_formal_adjoint_in_conformal_density": identity,
        "retained_derivative_of_transposed_coefficient": adjoint
        - (
            g.L.T * detector.diff(g.eta, 2)
            + g.Q.T * detector
            - g.C.T * detector.diff(g.eta)
            - g.C.diff(g.eta).T * detector
        ),
        "primal_Volterra_majorant": C0 * CMAX - FORWARD,
        "adjoint_Volterra_strict_margin": ADJOINT
        - C0 * (CMAX + CPRIME)
        - s.Rational(1, 2),
        "ordered_triangle_with_left_linear_kernel": simplex - delta**3 / 6,
        "full_curved_reference_kernel_bound": kernel
        - s.Rational(375, 16) * delta**3 * s.exp(ADJOINT * delta),
        "full_curved_reference_operator_bound": operator
        - s.Rational(375, 64) * T**4 * s.exp(ADJOINT * T),
        "physical_force_density_bound_keeps_kappa_and_a4": 6
        * 64
        * s.pi**2
        * kap
        * operator
        - 2250 * s.pi**2 * kap * T**4 * s.exp(ADJOINT * T),
    }
    return {
        "complete_formal_adjoint_action": adjoint,
        "exact_integration_by_parts_boundary": boundary,
        "primal_velocity_equation": "u=R0'*f-R0'*(C u), y=Iu. Its Volterra kernel is bounded by20; both causal coordinate inverse identities follow from this equation and the unchanged initial boundary.",
        "adjoint_velocity_equation": "v=(R0')^T*f+(R0')^T*(C^T v)+(R0')^T*((C')^T Iv), z=Iv. Its reordered Volterra kernel is below108 forT<=1. Bc*=B0^T-D(C^T), not B0^T-C^T D.",
        "coordinate_inverse_bounds": "||partial_t Y(t,s)||<=(5/2)exp(20(t-s)), ||Y(t,s)||<=(5/2)(t-s)exp(20(t-s)); ||partial_t Z(t,s)||<=(5/2)exp(108(t-s)), ||Z(t,s)||<=(5/2)(t-s)exp(108(t-s)). Y=Bc^-1 andZ=(Bc*)^-1 are retarded kernels in their actual order.",
        "reference_definition": "Aref=Bc* diag(Ftrace(D_eta^2+q),(8/3)F2(D_eta^2+q)) Bc. The scalar factors remain the fixed flat reference factors, not an asserted actual curved state/mass substitution.",
        "ordinary_composite_kernel": "Ecurv(t,s)=integral_s^t du integral_s^u dv Y(t,u)Jdiag(u-v)partial_v Z(v,s). Z(s,s)=0 permits moving the middle causal derivative onto Z. This is the ordered inverseY Kdiag Z, not a commuting convolution ansatz.",
        "uniform_kernel_bound": kernel,
        "uniform_operator_bound": operator,
        "spaces": "For every realr, the normalized curvature-adapted reference inverse is bounded onC_tH^r by375T^4 exp(108T)/64, all comoving momenta, withT<=1. A scalar dominating difference kernel gives the sameL2_tH^r bound. These are amplitude and dual-source norms; the metric embedding costs the retained factor14 when used.",
        "graph": "The causal forward graph includes the complete initial boundary. Both inverse products hold there. Smooth coefficient Volterra inverses and the explicit scalar forward distributions define the ordered composition; no initial atom or shear pole is deleted.",
        "physical_units": "Define the conformal-density reference force by Qbar_ref=M_(1/(64pi^2 kappa a^4)) Aref. Its inverse is Ecurv composed on the RIGHT withM_(64pi^2 kappa a^4), with norm below2250pi^2 kappa T^4 exp(108T). The time-dependent density is not commuted through the kernel.",
        "scope": "The actual finite local Hessian's explicit lower-order remainder is retained separately. Conformal mass-scale/state/contact/tree/matter differences and the full S222 auxiliary/clock bridge are not proved bounded perturbations here.",
        "checks": checks,
        "gates": {
            "primal_Volterra_constant_finite": FORWARD > 0,
            "adjoint_Volterra_constant_strict_enclosure": C0 * (CMAX + CPRIME)
            < ADJOINT,
            "uniform_reference_constant_below_six": s.Rational(375, 64) < 6,
            "actual_nonzero_coefficient_derivative_retained": g.C.diff(g.eta)
            != s.zeros(2),
            "physical_force_bound_not_kappa_small": 2250 * s.pi**2 * kap > 1,
            "original_shifted_shear_primitive_bound_retained": JMAX
            == s.Rational(45, 2),
        },
    }
