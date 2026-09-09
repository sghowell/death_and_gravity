"""Constrained physical Proca symbol and actual nearby canonical block."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_offclock_scalar import principal, spatial
from p8_zero_source import neighborhood as previous


@cache
def symbol():
    frequency,k1,k2,k3=sp.symbols("physical_frequency physical_k1 physical_k2 physical_k3",real=True)
    m2=sp.Symbol("positive_Proca_mass_squared",positive=True)
    metric=sp.diag(-1,1,1,1)
    covector=sp.Matrix([frequency,k1,k2,k3])
    raised=covector.T*metric
    norm=(raised*covector)[0]
    operator=(norm+m2)*sp.eye(4)-covector*raised
    inverse=(sp.eye(4)+covector*raised/m2)/(norm+m2)
    return {"physical_cotangent_norm":norm,"mass_squared":m2,
            "covector":covector,"raised_row":raised,
            "full_massive_plane_wave_operator":operator,
            "full_massive_plane_wave_inverse":inverse,
            "constrained_wave_operator":(norm+m2)*sp.eye(4),
            "constraint":"delta_g W=0 follows from delta_g(delta_g*d+m^2)W=m^2*delta_g W=0. The equivalent wave equation is (delta_g*d+d*delta_g+m^2)W=0 with this constraint.",
            "principal_boundary":"The Proca operator itself is not normally hyperbolic. Its constrained wave system has principal cone g^{mu nu} xi_mu xi_nu=0. This is a vector finite-propagation result on admitted smooth globally hyperbolic physical metrics, not full coupled scalar or nonlinear semiclassical stability."}


@cache
def neighborhood():
    old=previous.data()
    d=old["unchanged_constraint_compatible_light_family"]
    t,q=spatial.t,principal.q
    A=sp.MutableDenseMatrix(old["new_momentum_Hessian"])
    C=sp.MutableDenseMatrix(old["new_coordinate_Hessian"])
    gt,gs=(spatial.at_lapse(d[name]) for name in ("gamma_t","gamma_s"))
    A[2,2]+=(1-gt)/t
    C[2,2]+=t**3*q*(1-1/gs)
    return {"constant_mass_momentum_Hessian":sp.ImmutableMatrix(A.applyfunc(principal.clean)),
            "constant_mass_coordinate_Hessian":sp.ImmutableMatrix(C.applyfunc(principal.clean)),
            "unchanged_zero_scalar_vector_cross_block":old["new_momentum_coordinate_block"],
            "new_temporal_auxiliary_pivot":-model.N**sp.Rational(5,2),
            "physical_squared_spatial_momentum":q/t**2,
            "physical_lapse":t**2,
            "scope":"The exact same light constraint family is retained. The raw central vector Hamiltonian has positive kinetic and mass/gradient coefficients. Background time jets are not removed from the vector Euler equation; the covariant constrained-wave argument fixes its characteristic cone."}


@cache
def checks():
    d=symbol()
    P,R=d["full_massive_plane_wave_operator"],d["full_massive_plane_wave_inverse"]
    m2,norm=d["mass_squared"],d["physical_cotangent_norm"]
    xi,raised=d["covector"],d["raised_row"]
    n=neighborhood()
    A,C=n["constant_mass_momentum_Hessian"],n["constant_mass_coordinate_Hessian"]
    t,q=spatial.t,principal.q
    return {"full_Proca_operator_implies_divergence_constraint":(raised*P-m2*raised).applyfunc(sp.factor),
            "constraint_completion_is_physical_Klein_Gordon_symbol":(P+xi*raised-d["constrained_wave_operator"]).applyfunc(sp.factor),
            "full_massive_Proca_inverse_left":(P*R-sp.eye(4)).applyfunc(sp.factor),
            "full_massive_Proca_inverse_right":(R*P-sp.eye(4)).applyfunc(sp.factor),
            "three_massive_polarizations_in_full_determinant":sp.factor(P.det()-m2*(norm+m2)**3),
            "unconstrained_principal_Proca_symbol_is_singular":sp.factor((P-m2*sp.eye(4)).det()),
            "new_nearby_vector_kinetic_is_positive_physical_Proca":principal.clean(A[2,2]-1/t-10**6*t/q),
            "new_nearby_vector_gradient_is_physical_Proca":principal.clean(C[2,2]-t**3*q),
            "physical_lapse_and_metric_give_the_correct_central_dispersion":principal.clean(
                A[2,2]*C[2,2]-t**4*(q/t**2+10**6)),
            "new_all_nearby_scalar_vector_cross_block_zero":n["unchanged_zero_scalar_vector_cross_block"]}
