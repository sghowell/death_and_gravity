"""Finite derivative-coordinate lift and exact auxiliary constraint density."""

from functools import cache

import sympy as s

from . import split


def omega(n):
    return s.zeros(n).row_join(s.eye(n)).col_join((-s.eye(n)).row_join(s.zeros(n)))


@cache
def point_lift():
    q, v, z1, z2 = s.symbols(
        "base_position independent_velocity auxiliary1 auxiliary2", real=True
    )
    Pq, Pv, P1, P2 = s.symbols(
        "new_base_momentum new_velocity_momentum new_auxiliary_momentum1 new_auxiliary_momentum2",
        real=True,
    )
    newq = s.Matrix([q, v, z1, z2])
    newp = s.Matrix([Pq, Pv, P1, P2])
    shift = s.Matrix([q * v + q * q, v * v - q])
    oldq = s.Matrix([q, v, z1 + shift[0], z2 + shift[1]])
    J = oldq.jacobian(newq)
    oldp = J.inv().T * newp
    phase = oldq.col_join(oldp)
    jac = phase.jacobian(newq.col_join(newp))
    return {
        "whole_finite_jet_configuration_variables": newq,
        "whole_finite_jet_momenta": newp,
        "whole_velocity_dependent_algebraic_shift": shift,
        "whole_jet_point_coordinates": oldq,
        "whole_jet_point_Jacobian": J,
        "whole_jet_cotangent_momenta": oldp,
        "whole_eight_phase_cotangent_Jacobian": jac,
        "whole_jet_lift_argument": "Promote every finite derivative in the exact connection shift to an independent jet coordinate, retaining its defining kinematic multiplier constraint. On this extended configuration space the invertible full affine split is a point transformation, with momenta pulled back by its entire Jacobian transpose. The auxiliary momentum constraints then remove the induced base/jet momentum shifts involving those momenta. This is not an assertion that a velocity-dependent field redefinition was point canonical on the original unextended space.",
        "whole_actual_parent_binding": "The complete S174 source-centered connection depends on clock Hessian and metric/clock derivative data. Those finite jets, the original shifted retained vector W=T-B du, the original time/spatial boundaries and all remaining kinematic/metric constraints are retained. The full 64-coordinate inverse and non-unit configuration Jacobian are in split.projective. The general cotangent theorem, not this small diagnostic alone, applies to every one of those directions.",
        "checks": {
            "whole_four_coordinate_canonical_one_form": s.simplify(J.T * oldp - newp),
            "whole_eight_phase_symplectic_identity": s.simplify(
                jac.T * omega(4) * jac - omega(4)
            ),
            "whole_eight_phase_Liouville_Jacobian": s.factor(jac.det() - 1),
            "whole_jet_point_shift_Jacobian": s.factor(J.det() - 1),
            "whole_primary_auxiliary_momenta_unchanged": oldp[2:, :] - newp[2:, :],
        },
        "gates": {
            "shift_really_depends_on_independent_velocity": shift.diff(v)
            != s.zeros(2, 1),
            "full_jet_momentum_contact_is_nonzero": oldp[:2, :] != newp[:2, :],
            "all_eight_diagnostic_phase_coordinates_retained": jac.shape == (8, 8),
        },
    }


@cache
def primary_secondary():
    q, v, z1, z2, p, Pv, P1, P2 = s.symbols("q v z1 z2 p P_v p_z1 p_z2", real=True)
    positions = s.Matrix([q, v, z1, z2])
    momenta = s.Matrix([p, Pv, P1, P2])
    variables = positions.col_join(momenta)
    z = s.Matrix([z1, z2])
    A = s.Matrix([[2 + q * q, q * v], [q * v, -3 - v * v]])
    Lred = (v * v - q * q) / 2
    H = p * v - Lred - (z.T * A * z)[0] / 2
    Cz = s.Matrix([s.diff(H, item) for item in z])
    Cv = s.diff(H, v)
    constraints = s.Matrix([P1, P2, *Cz, Pv, Cv])
    derivatives = constraints.jacobian(variables)
    bracket = derivatives * omega(4) * derivatives.T
    aux_zero = {z1: 0, z2: 0, P1: 0, P2: 0}
    surface = s.simplify(bracket.subs(aux_zero, simultaneous=True))
    remaining = surface[4:, 4:]
    D = -A
    wanted = s.zeros(2).row_join(-D).col_join(D.T.row_join(s.zeros(2)))
    surface_derivatives = derivatives.subs(aux_zero, simultaneous=True)
    pq_bracket = s.simplify(
        omega(4)
        - omega(4)
        * surface_derivatives.T
        * surface.inv()
        * surface_derivatives
        * omega(4)
    )
    reduced = s.factor(H.subs(aux_zero, simultaneous=True).subs(v, p))
    f = s.Matrix([q * v + q * q, v * v - q])
    y = s.Matrix(s.symbols("original_auxiliary1 original_auxiliary2", real=True))
    original = Lred + ((y - f).T * A * (y - f))[0] / 2
    velocity_hessian = s.simplify(
        s.diff(original, v, 2).subs(dict(zip(y, f, strict=True)), simultaneous=True)
    )
    singular_square = (s.sqrt(10) - 3) / 2
    diagnostic = s.factor(velocity_hessian.subs(q, 0))
    return {
        "whole_finite_jet_base_auxiliary_variables": variables,
        "whole_finite_indefinite_auxiliary_Hessian": A,
        "whole_exact_first_order_jet_Hamiltonian": H,
        "whole_six_primary_and_secondary_constraints": constraints,
        "whole_six_constraint_bracket": bracket,
        "whole_six_constraint_bracket_on_auxiliary_surface": surface,
        "whole_retained_kinematic_constraint_bracket": remaining,
        "whole_reduced_physical_Hamiltonian": reduced,
        "whole_original_unextended_velocity_Hessian_at_auxiliary_solution": velocity_hessian,
        "whole_nonempty_unextended_Legendre_singularity": singular_square,
        "whole_general_56_constraint_argument": "After the full finite-jet point lift the exact action is the retained first-order jet action plus lambda z^T A z/2, with A nonsingular. The 56 auxiliary primary momenta and 56 secondary equations have invertible bracket [[0,-D],[D^T,E]], D=partial_z C. Its determinant is (det D)^2 and its inverse lower-right block is zero. Retained functions commute with the auxiliary primaries, so their Dirac bracket is their original canonical bracket. The positive density |det D| cancels the delta-function Jacobian on the one regular auxiliary root. Remaining kinematic, lapse, spatial and other constraints are not discarded or assumed nonsingular; they are pulled back to the root with the induced Dirac bracket.",
        "whole_singular_Legendre_boundary": "The small diagnostic has a regular reduced oscillator and a nonsingular auxiliary A everywhere, but its original unextended velocity Hessian vanishes at q=0 and v^2=(sqrt(10)-3)/2. Thus an argument that inverts that original velocity Hessian globally would fail. The finite-jet constraint reduction does not require that inverse. It also does not prove a nonlinear quantum ordering/time-slicing prescription or a physical evolution theorem for the complete parent.",
        "checks": {
            "whole_exact_two_auxiliary_Dirac_block": surface[:4, :4] - wanted,
            "whole_auxiliary_remaining_constraint_cross_on_zero_source_root": surface[
                :4, 4:
            ],
            "whole_complete_six_constraint_density_factorization": s.factor(
                surface.det() - A.det() ** 2 * remaining.det()
            ),
            "whole_remaining_kinematic_constraint_pair": remaining
            - s.Matrix([[0, 1], [-1, 0]]),
            "whole_retained_physical_canonical_bracket": s.simplify(
                pq_bracket[0, 4] - 1
            ),
            "whole_reduced_physical_oscillator": s.factor(
                reduced - (p * p + q * q) / 2
            ),
            "whole_original_velocity_Hessian_not_assumed_regular": s.simplify(
                diagnostic.subs(v, s.sqrt(singular_square))
            ),
        },
        "gates": {
            "full_auxiliary_A_inertia_indefinite_not_ghost_count": A.det().subs(
                {q: 0, v: 0}
            )
            < 0,
            "all_remaining_kinematic_constraints_retained": constraints.rows == 6,
            "remaining_reduced_velocity_is_regular": s.diff(Lred, v, 2) == 1,
            "original_velocity_singularity_is_real": singular_square > 0,
        },
    }


@cache
def source_dependent_jet():
    q, v = s.symbols("q independent_velocity", real=True)
    z = s.Matrix(s.symbols("z1 z2", real=True))
    J = s.Matrix(
        s.symbols("original_connection_source1 original_connection_source2", real=True)
    )
    A = s.Matrix([[2 + q * q, q * v], [q * v, -3 - v * v]])
    f = s.Matrix([q * v + q * q, v * v - q])
    Lred = (v * v - q * q) / 2
    L = Lred + (z.T * A * z)[0] / 2 + J.dot(z + f)
    solution = -A.inv() * J
    reduced = s.factor(L.subs(dict(zip(z, solution, strict=True)), simultaneous=True))
    expected = Lred + J.dot(f) - (J.T * A.inv() * J)[0] / 2
    zero = {J[0]: 0, J[1]: 0}
    return {
        "whole_jet_source_dependent_auxiliary_solution": solution,
        "whole_original_source_pulled_back_Lagrangian": L,
        "whole_source_dependent_retained_jet_Lagrangian": reduced,
        "whole_retained_source_dependent_momentum_equation": s.diff(reduced, v),
        "whole_source_dependent_jet_rule": "An original source J dot y becomes J dot (z+f(q,v)), not J dot z alone. Eliminate the source-dependent auxiliary equations, retaining -J^T A^-1 J/2 and J dot f. The remaining kinematic momentum equation is differentiated from this entire reduced jet Lagrangian. Its source variations include the nonlinear f derivatives and the auxiliary inverse contact. The same source pullback rule applies to the whole actual 64-connection map.",
        "checks": {
            "whole_source_dependent_jet_auxiliary_Euler": split.clean(A * solution + J),
            "whole_full_source_dependent_jet_action": s.factor(reduced - expected),
            "whole_source_zero_recovers_retained_jet_action": s.factor(
                reduced.subs(zero) - Lred
            ),
            "whole_original_first_source_pullback": split.clean(
                s.Matrix([s.diff(reduced, item).subs(zero) for item in J]) - f
            ),
            "whole_original_second_source_contact": split.clean(
                s.hessian(reduced, J) + A.inv()
            ),
            "whole_nonzero_second_source_velocity_contact": s.factor(
                s.diff(reduced, v, 2, J[1]).subs(zero) - 2
            ),
        },
        "gates": {
            "whole_source_contact_nonzero": A.inv() != s.zeros(2),
            "source_changes_remaining_momentum_constraint": s.diff(
                reduced, v, J[1]
            ).subs(zero)
            != 0,
            "classical_zero_source_root_not_used_at_nonzero_source": True,
        },
    }


@cache
def ordering_boundary():
    Q, hbar = s.symbols("new_coordinate hbar", real=True)
    wave = s.Function("test_wavefunction")(Q)
    B = s.exp(-Q)

    def transformed_momentum(f):
        return -s.I * hbar * (B * s.diff(f, Q) + s.diff(B, Q) * f / 2)

    unitary_H = s.simplify(transformed_momentum(transformed_momentum(wave)) / 2)
    a = B * B
    weyl_H = (
        -hbar
        * hbar
        * (
            a * s.diff(wave, Q, 2)
            + s.diff(a, Q) * s.diff(wave, Q)
            + s.diff(a, Q, 2) * wave / 4
        )
        / 2
    )
    contact = hbar * hbar * s.exp(-2 * Q) * wave / 8
    return {
        "independent_nonlinear_point_map": "x=exp(Q) on x>0, with wavefunctions transformed as half-densities sqrt(dx/dQ) psi(x(Q)).",
        "whole_unitarily_pulled_back_free_Hamiltonian": unitary_H,
        "incorrect_unchanged_naive_Weyl_symbol_quantization": weyl_H,
        "whole_nonzero_ordering_contact": contact,
        "quantum_boundary": "A unit Liouville Jacobian and classical cotangent identity do not by themselves imply covariance of an unchanged Weyl symbol under nonlinear coordinates. This exact free-particle counterexample has contact hbar^2 exp(-2Q)/8. For the present parent, retain the transported finite time-slicing/operator prescription and state/boundaries. The conditional complement insertion and finite constraint-density identities do not separately construct the full interacting continuum quantum theory.",
        "checks": {
            "whole_exact_nonlinear_ordering_contact": s.simplify(
                unitary_H - weyl_H - contact
            )
        },
        "gates": {
            "ordering_contact_is_nonzero": contact != 0,
            "finite_constraint_density_not_full_quantum_equivalence_by_itself": True,
        },
    }
