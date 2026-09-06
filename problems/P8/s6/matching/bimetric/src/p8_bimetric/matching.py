"""Source-aware flat-vacuum matching, with the actual massive mode retained.

The mode action convention is -h^T K h/8+j*h_g/4, so K*h=(j,0).
The source j is normalized; for a +h:T/2 matter variation j=2*T. D is the
flat Lorentzian wave operator, not spatial k^2 or a cosmological frequency.
"""

import itertools

import sympy as sp

from . import vacuum
from .background import MF2, MG2, NU

D = vacuum.D


def tt_schur():
    spectrum = vacuum.spectrum()
    total, mass = spectrum["M_squared"], spectrum["m_FP_squared"]
    kernel = sp.Matrix([[MG2*D+NU, -NU], [-NU, MF2*D+NU]])
    schur = sp.factor(kernel[0, 0]-kernel[0, 1]*kernel[1, 0]/kernel[1, 1])
    response = 1/total/D+MF2/(MG2*total)/(D+mass)
    local_kernel = total*D-MF2**2*D**2/NU
    return {"kernel": kernel, "physical_g_kernel": schur,
            "physical_g_response": response,
            "four_derivative_kernel": local_kernel,
            "kernel_remainder": MF2**3*D**3/(NU*(NU+MF2*D)),
            "c_C": MF2**2/(4*NU), "c_R_flat_quadratic_tree": sp.Integer(0),
            "beta_C": MF2**2/(4*NU*total),
            "heavy_source_contact": MF2**2/(NU*total**2),
            "inverse_f_kernel_pole": -NU/MF2,
            "physical_massive_pole": -mass}


def source_matching_checks():
    data, spectrum = tt_schur(), vacuum.spectrum()
    total, mass = spectrum["M_squared"], spectrum["m_FP_squared"]
    h0, hm, source = vacuum.H0, vacuum.HM, vacuum.SOURCE
    reduced = spectrum["relative_mode_kinetic_coefficient"]
    # True massive-mode integration: source kept in physical g=h0-F/M^2*hm.
    massive_action = -reduced*(D+mass)*hm**2/8-source*MF2*hm/(4*total)
    hm_stationary = -source*MF2/(total*reduced*(D+mass))
    source_contact = sp.factor(massive_action.subs(hm, hm_stationary))
    expected_contact = source**2*MF2/(8*MG2*total*(D+mass))
    g_from_modes = h0-MF2*hm/total
    h0_stationary = source/(total*D)
    physical_response = g_from_modes.subs({h0: h0_stationary, hm: hm_stationary}, simultaneous=True)/source
    return {"Schur_matches_full_matrix_response": sp.factor(1/data["physical_g_kernel"]-data["kernel"].inv()[0, 0]),
            "physical_g_positive_pole_decomposition": sp.factor(1/data["physical_g_kernel"]-data["physical_g_response"]),
            "massive_mode_own_equation": sp.factor(sp.diff(massive_action, hm).subs(hm, hm_stationary)),
            "massive_mode_source_source_functional": sp.factor(source_contact-expected_contact),
            "physical_metric_observable_after_true_massive_integration": sp.factor(physical_response-data["physical_g_response"]),
            "exact_local_kernel_remainder": sp.factor(data["physical_g_kernel"]-data["four_derivative_kernel"]-data["kernel_remainder"]),
            "Weyl_TT_normalization": sp.factor(-sp.Rational(1, 8)*(data["four_derivative_kernel"]-total*D)-data["c_C"]*D**2/2),
            "Weyl_mass_gap_dictionary": sp.factor(data["c_C"]-total*spectrum["alpha_squared"]/(4*mass)),
            "heavy_contact_Weyl_dictionary": sp.factor(data["heavy_source_contact"]-4*data["c_C"]/total**2)}


def projectors():
    """Conserved non-null source space in an orthonormal symmetric basis.

    The first three entries are diagonal components; the other three are
    sqrt(2) times off-diagonal components. This is a projector algebra test,
    not a claim that the off-shell scalar constraint has positive norm.
    """
    scalar = sp.zeros(6)
    scalar[:3, :3] = sp.ones(3)/3
    spin2 = sp.eye(6)-scalar
    return {"P2": spin2, "P0": scalar}


def projector_checks():
    p = projectors()
    p2, p0 = p["P2"], p["P0"]
    data, spectrum = tt_schur(), vacuum.spectrum()
    total, mass = spectrum["M_squared"], spectrum["m_FP_squared"]
    response = data["physical_g_response"]*p2-p0/(2*total*D)
    kernel = data["physical_g_kernel"]*p2-2*total*D*p0
    local = total*D*(p2-2*p0)-4*data["c_C"]*D**2*p2
    values = {"spin2_rank_five": p2.trace()-5, "constraint_scalar_rank_one": p0.trace()-1}
    for name, matrix in {"spin2_projector": p2*p2-p2, "scalar_projector": p0*p0-p0,
                         "orthogonal_projectors": p2*p0,
                         "full_conserved_source_inverse": kernel*response-sp.eye(6),
                         "covariant_curvature_squared_dictionary": kernel-local-data["kernel_remainder"]*p2}.items():
        values.update({f"{name}_{index}": sp.factor(value) for index, value in enumerate(matrix)})
    # I_sym-1/3*trace for massive exchange versus I_sym-1/2*trace
    # for massless exchange. This is the scalar/source distinction TT misses.
    v = sp.Matrix(sp.symbols("T11 T22 T33 sqrt2T12 sqrt2T13 sqrt2T23", real=True))
    norm, trace = v.dot(v), sum(v[:3, 0])
    values["massive_source_trace_coefficient_one_third"] = sp.expand((v.T*p2*v)[0]-norm+trace**2/3)
    values["massless_source_trace_coefficient_one_half"] = sp.expand((v.T*(p2-p0/2)*v)[0]-norm+trace**2/2)
    values["massive_source_functional_after_mass_basis_integration"] = sp.factor(
        (v.T*(response-(p2-p0/2)/(total*D))*v)[0]
        -MF2/(MG2*total)*(norm-trace**2/3)/(D+mass))
    return values


def physical_residue_checks():
    """Positive on-pole contractions; do not count off-shell P0 as a ghost."""
    x, y, z, xy, xz, yz = sp.symbols("T11 T22 T33 T12 T13 T23", real=True)
    signs = (1, -1, -1, -1)
    # Null p=(1,0,0,1), T_{0 nu}+T_{3 nu}=0; six independent entries.
    source = sp.Matrix([[z, -xz, -yz, -z], [-xz, x, xy, xz],
                        [-yz, xy, y, yz], [-z, xz, yz, z]])
    null_conservation = source*sp.Matrix([1, 0, 0, 1])
    contraction = sum(signs[i]*signs[j]*source[i, j]**2 for i in range(4) for j in range(4))
    trace = sum(signs[i]*source[i, i] for i in range(4))
    massive_norm = x*x+y*y+z*z+2*(xy*xy+xz*xz+yz*yz)-(x+y+z)**2/3
    massive_squares = ((x-y)**2+(y-z)**2+(z-x)**2)/3+2*(xy*xy+xz*xz+yz*yz)
    return {**{f"null_source_conservation_{index}": value for index, value in enumerate(null_conservation)},
            "massless_two_helicity_positive_contraction": sp.expand(contraction-trace**2/2-(x-y)**2/2-2*xy**2),
            "massive_five_polarization_positive_contraction": sp.expand(massive_norm-massive_squares)}


def fierz_pauli_source_checks():
    """Derive the trace projectors from all ten metric components at rest.

    At a non-null timelike momentum a conserved source has j_0mu=0.
    The FP lapse enforces the spatial trace constraint; its vector equations
    eliminate h_0i. This prevents a TT-only assumption about scalar exchange.
    The massless temporal components can instead be gauge fixed to zero.
    """
    kinetic, mass = sp.symbols("kinetic mass_squared", positive=True)
    lapse = sp.Symbol("h00", real=True)
    shift = sp.Matrix(sp.symbols("h01 h02 h03", real=True))
    spatial = sp.Matrix(sp.symbols("h11 h22 h33 sqrt2h12 sqrt2h13 sqrt2h23", real=True))
    source = sp.Matrix(sp.symbols("j11 j22 j33 sqrt2j12 sqrt2j13 sqrt2j23", real=True))
    trace, source_trace = sum(spatial[:3, 0]), sum(source[:3, 0])
    einstein = -kinetic*D*(spatial.dot(spatial)-trace**2)/8
    fp = -kinetic*mass*(spatial.dot(spatial)-trace**2+2*lapse*trace-2*shift.dot(shift))/8
    action = einstein+fp+source.dot(spatial)/4
    p2, p0 = projectors()["P2"], projectors()["P0"]
    solution = dict(zip(spatial, p2*source/(kinetic*(D+mass)), strict=True))
    solution.update(dict.fromkeys(shift, 0))
    solution[lapse] = source_trace/(3*kinetic*mass)
    checks = {"FP_lapse_enforces_trace": sp.expand(sp.diff(action, lapse)+kinetic*mass*trace/4)}
    for index, component in enumerate([lapse, *shift, *spatial]):
        checks[f"all_ten_FP_stationary_equations_{index}"] = sp.factor(
            sp.diff(action, component).subs(solution, simultaneous=True))
    checks["FP_full_stationary_source_functional"] = sp.factor(
        action.subs(solution, simultaneous=True)
        -(source.dot(source)-source_trace**2/3)/(8*kinetic*(D+mass)))
    massless_action = einstein+source.dot(spatial)/4
    massless_solution = dict(zip(spatial, (p2-p0/2)*source/(kinetic*D), strict=True))
    for index, component in enumerate(spatial):
        checks[f"massless_gauge_fixed_stationary_equations_{index}"] = sp.factor(
            sp.diff(massless_action, component).subs(massless_solution, simultaneous=True))
    checks["massless_full_stationary_source_functional"] = sp.factor(
        massless_action.subs(massless_solution, simultaneous=True)
        -(source.dot(source)-source_trace**2/2)/(8*kinetic*D))
    return checks


def flat_curvature_checks():
    """Spin-projector coefficient check using independently stated invariants.

    For a transverse metric h: R_lin^2=D^2*(tr h)^2 and, modulo the
    quadratic Euler boundary, C_lin^2=D^2*(h:h-(tr h)^2/3)/2.
    These determine both curvature-squared coefficients at flat quadratic
    order; no nonlinear or rolling-CD matching is inferred.
    """
    h = sp.Matrix(sp.symbols("h11 h22 h33 sqrt2h12 sqrt2h13 sqrt2h23", real=True))
    p = projectors()
    trace = sum(h[:3, 0])
    weyl = D**2*(h.dot(h)-trace**2/3)/2
    ricci_scalar_squared = D**2*trace**2
    return {"Weyl_only_spin2_curvature_structure": sp.expand(weyl-D**2*(h.T*p["P2"]*h)[0]/2),
            "R2_scalar_curvature_structure": sp.expand(ricci_scalar_squared-3*D**2*(h.T*p["P0"]*h)[0])}


def direct_curvature_checks():
    """Four-index contraction of arbitrary transverse timelike Fourier jets.

    Six arbitrary spatial metric entries are retained. Covariant polynomial
    identities extend the open timelike-momentum calculation; the null-pole
    source-residue test is separate and does not divide by p^2.
    """
    x, y, z, xy, xz, yz = sp.symbols("h11 h22 h33 h12 h13 h23", real=True)
    h = sp.Matrix([[x, xy, xz], [xy, y, yz], [xz, yz, z]])
    signs = (1, -1, -1, -1)

    def hessian(i, j, c, d):
        return D*h[i-1, j-1] if i > 0 and j > 0 and c == d == 0 else 0

    riemann = {(a, b, c, d): (hessian(a, d, c, b)+hessian(b, c, d, a)
                             -hessian(a, c, d, b)-hessian(b, d, c, a))/2
               for a, b, c, d in itertools.product(range(4), repeat=4)}
    ricci = {(b, d): sum(signs[a]*riemann[a, b, a, d] for a in range(4))
             for b, d in itertools.product(range(4), repeat=2)}
    scalar = sum(signs[a]*ricci[a, a] for a in range(4))
    riemann2 = sum(sp.prod(signs[i] for i in indices)*value**2 for indices, value in riemann.items())
    ricci2 = sum(signs[b]*signs[d]*value**2 for (b, d), value in ricci.items())
    weyl = sp.expand(riemann2-2*ricci2+scalar**2/3)
    norm = sum(value**2 for value in h)
    trace = sp.trace(h)
    return {"literal_scalar_curvature_squared": sp.expand(scalar**2-D**2*trace**2),
            "literal_Weyl_spin2_structure": sp.expand(weyl-D**2*(norm-trace**2/3)/2),
            "literal_Euler_boundary_for_time_only_mode": sp.expand(riemann2-4*ricci2+scalar**2)}


def remainder_record():
    """Exact rational-kernel remainder, not a cosmological cutoff bound."""
    eta, magnitude = sp.symbols("eta abs_D", nonnegative=True)
    return {"domain": "0<=eta<1 and |MF2*D|<=eta*NU; away from poles, with a specified vacuum Green prescription",
            "kernel_absolute_remainder": MF2**3*magnitude**3/(NU**2*(1-eta)),
            "kernel_relative_remainder": eta**2/(1-eta),
            "warning": "D is a Lorentz-invariant wave operator; this alone bounds neither arbitrary spatial frequency, interaction strength nor rolling-time error"}


def remainder_checks():
    ratio = sp.Symbol("ratio", real=True)
    data = tt_schur()
    relative = data["kernel_remainder"]/((MG2+MF2)*D)
    # The absolute-value inequalities in the proof follow from this identity,
    # |1+ratio|>=1-|ratio|, and 0<MF2/(MG2+MF2)<1.
    return {"relative_kernel_remainder_in_disk_variable": sp.factor(
                relative.subs(D, NU*ratio/MF2)-MF2*ratio**2/((MG2+MF2)*(1+ratio))),
            "strict_weight_less_than_one_margin": sp.factor(1-MF2/(MG2+MF2)-MG2/(MG2+MF2))}


def controls():
    data, spectrum = tt_schur(), vacuum.spectrum()
    # Wrong-field elimination, with sources set to zero, is not heavy-mode
    # integration. Its changed action need not describe the same solutions.
    f_from_g_equation = 1+MG2*D/NU
    wrong_kernel = MG2*D+NU-2*NU*f_from_g_equation+(MF2*D+NU)*f_from_g_equation**2
    wrong_d2 = sp.expand(wrong_kernel).coeff(D, 2)
    right_d2 = sp.expand(data["four_derivative_kernel"]).coeff(D, 2)
    return {"discarding_massive_source_contact_loses_exchange": spectrum["massive_source_residue"]/(D+spectrum["m_FP_squared"]),
            "wrong_field_equation_changes_four_derivative_kernel": sp.factor(wrong_d2-right_d2),
            "massless_and_massive_trace_projectors_differ": sp.Rational(1, 6),
            "physical_g_not_the_massless_eigenmode": -MF2*vacuum.HM/(MG2+MF2),
            "zero_curvature_R2_tree_coefficient_not_radiative_closure": sp.Integer(0),
            "auxiliary_inverse_pole_not_physical_massive_pole": sp.factor(data["inverse_f_kernel_pole"]-data["physical_massive_pole"])}
