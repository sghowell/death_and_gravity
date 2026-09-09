# Actual leading mean equations, charge and Ward completion

Write A=hat_a, b=log A and p=2 pi_trace/(3 A^3).
The homogeneous symplectic term is 3 A^3 p b', so its
canonical density momentum is P_b=3 A^3 p. The free matter
charge is Q_psi and ell=Q_psi/A^3. The actual scalar
Hamiltonian has H=A^3 Hcal(N,p,ell,u), including the
original boundary primitive and the new margin.

At fixed Q_psi, Hamilton's equations imply

    b'=Hcal_p/3,
    p'=-Hcal+ell Hcal_ell.

These follow by differentiating at fixed P_b, not by
treating p as a canonical density momentum. The additional
physical Proca Hamiltonian has no P_b dependence and
its b derivative is -3 A^3 s at the background.
Thus it adds s to p', not minus its energy density.

The exact existing lapse jets give, on the full original
clock, Hcal_N=0, Hcal_NN=-2J, Hcal_p=3H,
Hcal_pp=-3/2, Hcal_ell=ell, Hcal_ell,ell=1,
Hcal_p,ell=0 and Hcal_N,ell=beta. In particular

    alpha=Hcal_Np=3u[4(1+u^2)^3-1]/(1+u^2)^4,
    beta=ell[1-3/(2h)].

The literal retuned Hamiltonian, not an old nonconstant
Proca model, supplies J. The native reconstruction also
checks the original trace equation and conserved charge.

Let xi=delta b and dp=delta p. Keeping Q_psi fixed gives
delta ell=-3ell xi. Varying the lapse constraint and the
two evolution equations gives exactly FORMULATION.md.
The two-by-two generator has trace -3H, the required
expanding-density normalization. Its determinant need
not have a fixed sign to give a well-defined forced
linear Cauchy problem.

Initial xi=dp=0 do NOT mean n=0: the constraint fixes
n(0)=80[r(0)-3s(0)/2]/243. The initial scalar field
response is zero and its derivative is fixed by
delta psi'=beta n-3ell xi. This is the actual free
matter equation, not an independent energy assignment.

With B=xi+n/(2h), its normal velocity variation is
delta psi'-ell n=-3ell B. Hence delta rho_m=delta p_m
=-3ell^2 B=-6rho_m,0 B. Since rho_m,0'=-6H rho_m,0,
combining this with vector conservation gives

    (delta rho_m+r)'
       +3H(2 delta rho_m+r+s)
       +6rho_m,0 B'=0.

The last term is the connection/volume variation acting
on the background source. Omitting it would not give a
conserved physical mean correction. This is the complete
homogeneous identity for this relative state variation,
not the full covariant Ward identity for all seven modes.

At leading covariance order the changed mode functions
and all source contractions use the same classical
background. Their further response to the induced mean
is higher order. The unchanged fluctuation sectors and
state-independent counterterms cancel from the difference
of the two leading equations. This cancellation neither
computes nor sets to zero the common vacuum tadpole.
