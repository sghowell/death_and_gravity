# Compensate an analysis variable, not the physical solution

Let V=<v^2> denote only the selected positive covariance
addition. In the original mean system write

    c=-9H=(9/2)p0,
    d=delta_p-c V=delta_p+9H V.

This is an invertible change of the forced linear analysis
variable. It is not an action change, counterterm, independent
state or new trace observable. In particular

    delta_p=d-9H V

must be restored. Since H(0)=0, d and delta_p have the same
zero anchor. Existence/uniqueness then identifies the same
S6.106 physical solution and the same center jets.

First transform each natural scalar source Hessian with
B^T K B into the old canonical density phase. Let
V_mat=diag(2,0,0,0), so V=Tr(V_mat Sigma)/2, and

    Vdot_mat=M_old^T V_mat+V_mat M_old.

The compensated source matrices are

    N_c=K_N+alpha c V_mat,
    X_c=K_X-c V_mat/2,
    P_c=K_P-(3Hc+c')V_mat-c Vdot_mat,
    C_c=K_psi.

The variance derivative is transported by the ACTUAL coupled
generator, not set to zero. Independent symbolic substitution
checks all four original equations against

    n=(alpha d-3ell beta xi+<N_c>)/(2J),
    xi'=-d/2+alpha n/3+<X_c>,
    d'=-3H d+ell beta n-3ell^2 xi+<P_c>,
    delta psi'=beta n-3ell xi+<C_c>.

The homogeneous mean matrix is exactly the already bounded
S6.104 matrix. No scalar source has been replaced by a
conserved-fluid or minimally coupled source.

For Y=(xi,100t^3 d), the forcing and lapse/field kernels,
expressed in W coordinates, are

    alpha N_c/(6J)+X_c,
    100t^3[ell beta N_c/(2J)+P_c],
    N_c/(2J),
    C_c.

A source entry with m momentum indices gains t^(-m/2)
when paired with the extra weighted covariance. Its bound
therefore has effective power r+m/2. The half entry sums
of the four actual half-line integrals are respectively
less than 257, 2682, 158 and 11. All entries are reconstructed
and coefficientwise bounded, retaining off-diagonal terms.

Two exact asymptotic controls show why this change is useful.
The unshifted scale diagonal times u tends to 396; its
weighted trace diagonal divided by u^4 tends to -79200.
The corresponding compensated limits both vanish. A
nonintegrable unshifted forcing estimate is NOT a proof
that the actual physical solution diverges. The compensating
terms restore precisely the omitted canonical-volume chains.
