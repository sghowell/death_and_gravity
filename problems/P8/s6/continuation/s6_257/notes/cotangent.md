# Extended point map and the complete primitive boundary shift

Do not use a position-only measure Jacobian as a canonical phase-space
Jacobian. Before imposing primary constraints, retain all 16 positions:
six spatial metric coordinates, three Wi, three Ni, T,N, and two matter
spectators. The physical old positions are

    h_phys=C gamma, Wi_old=Wi, Ni_old=Ni,
    W0_old=N T+Ni Wi, N_old=N, matter_old=matter.

Their position Jacobian is C^6 N. Write the six dual metric momenta in a
consistent independent-coordinate convention, including the off-diagonal
tensor multiplicities. The full cotangent lift gives

    Pi_hat=C Pi_phys,
    pT=N pW0,
    pWi_new=pWi_old+Ni pW0,
    pNi_new=pNi_old+Wi pW0,
    pN_new=pN_old+C_N Pi_phys dot gamma+T pW0.

Both matter momenta are unchanged. The inverse is explicit in the report.
The identity p_old dot dF=P_new dot dq proves symplecticity on each clock
slice. Position and momentum Jacobians multiply to one, not C^6 N. Time
dependence gives

    H_point=H_old pulled back-(C_u/C) Pi_hat dot gamma.

The original covariant C^9 metric-position result is consistent with
C^6 in ADM metric positions: the ADM-to-ten-metric-coordinate Jacobian
has absolute value 2N det(gamma), and det(C gamma)=C^3 det(gamma).
Its ratio times C^6 gives C^9. The independent W0-to-T factor is N.

## Primitive boundary is a second canonical transformation

Within the original S174 ADM variational boundary convention, let
F_boundary=V I(u,N). The remaining linear lapse-velocity term is

    V I_s(s_dot-Ni partial_i s)
      =-N V I K-V I_u+partial_u(V I)-partial_i(V Ni I).

This follows from V_dot-partial_i(V Ni)=N V K; the packet keeps the
entire shift divergence. Thus the final Lagrangian equals the point-
transformed Lagrangian minus the displayed boundary, and

    P_final=P_point-partial_q F_boundary,
    H_final=H_point+partial_u F_boundary.

The complete inverse momenta substitute P_point=P_final+partial_q F_boundary
in the point inverse. In final momenta the combined additive time term is

    -(C_u/C)(Pi_final+partial_gamma F_boundary) dot gamma
      +partial_u F_boundary.

Omitting either the lapse-momentum shift, the six metric shifts, or the
extra time contact changes this canonical bookkeeping. The boundary's
symmetric full Hessian makes its momentum translation symplectic and
its phase Jacobian one. The report independently differentiates the
actual V I along a general 16-coordinate path and checks the combined
inverse. Finite-parameter tests inspect the complete 32-dimensional
phase Jacobian, not just a selected coordinate block.

Keep both endpoint phases and the spatial flux under the same boundary
conditions. On a closed-time contour transport the same initial density
and its two branch phases. Unit Liouville volume is not permission to
reset the quantum state or erase a physical boundary contribution.

