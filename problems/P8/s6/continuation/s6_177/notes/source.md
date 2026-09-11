# The complete canonical affine source and all removable quotients

Use the literal S6.174 source, not a newly defined effective source.
With Q=Box Phi, Z=Phi^mu Phi_(mu nu) Phi^nu and R=1+r/kappa,
the exact canonical substitution gives

    S_mu = Phi_mu r/kappa (Q/Y-Z/Y^2)
           -3 Phi_mu r Href/kappa^(3/2)
           +3 Phi_mu r/(kappa^2 R)
                [r_Phi/4+r_Y Z/(2Y)].

This includes the reference-H term and both derivatives of the full
r. Setting q=0 was not used: its cancellation follows from the full
retained shift B. For r=Y^2 a, write the manifestly regular vectors

    V_mu=Phi_mu a(YQ-Z),
    B_mu=-3 Phi_mu r Href,
    C_mu=3 Phi_mu r[r_Phi/4+(2a+Y a_Y)Z/2],
    U_mu=V_mu+B_mu/sqrt(kappa)+C_mu/(kappa R).

Then S=U/kappa, including at Y=0. The null value is
S_mu=-Phi_mu a Z/kappa and generally is NOT zero.

The EXACT canonical retained vector action is

    -F(A)^2/4 + A^2/(2zeta) - A.J + kappa S^2/2,
    A=sqrt(kappa*zeta)W, J=U/sqrt(kappa*zeta).

The squared source contact is not dropped. For Euclidean component
norms of vectors at eta, let v,b,c bound V,B,C on a fixed compact
canonical jet set. Since R>1/2 and kappa>=kappa0,

    ||U|| <= Ustar = v+b/sqrt(kappa0)+2c/kappa0,
    ||J|| <= Ustar/sqrt(kappa*zeta),
    |kappa S^2/2| <= Ustar^2/(2kappa).

The last is an ABSOLUTE Lorentzian contraction bound, not positivity
of that contraction. For bounded A, the source coupling is bounded
by ||A|| Ustar/sqrt(kappa*zeta). Near eta the bounded inverse metric
and volume merely change the constants. All functions are complete;
no Fourier support or leading vacuum-source replacement is required.

Likewise the two dependent canonical scalar coefficients satisfy

    |Delta A4| <= 7r_Y^2/(2kappa),
    |Delta A5| <= 2|Y|(2a+Y a_Y)^2/kappa.

Every fixed finite derivative of1/R is bounded uniformly on compact
canonical sets. Differentiating1/R introduces derivatives of r divided
by kappa; the zeroth inverse is at most2. Product rules therefore give
finite C_j/kappa bounds for the dependent Ia terms and the contact,
and C_j/sqrt(kappa) for J. Their constants depend on the full compact
set and derivative order, not on kappa.

The source-free normalized Proca field remains in the limit as a
free spectator. A classical source coefficient tending to zero does
not justify commuting the limit with a real-time inverse near its
mass shell, a functional determinant, regulator removal, state
preparation or an unbounded contour integral. No such quantum or
momentum-uniform interchange is claimed.
