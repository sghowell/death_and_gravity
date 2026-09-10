# Full action, exact elimination and generated vertices

Let S[Phi,H] be the actual local polynomial action, with mass-one light
field, quadratic heavy field and interaction G H Phi^2/2 plus
lambda4 Phi^4/24 in the potential. Set Phi=F(Psi)=Psi+R(Psi), where
R is exactly the cubic, at-most-four-derivative expression of S6.111.
This is a formal near-identity map in field degree. Its inverse begins
Psi=Phi-R(Phi)+O(Phi^5). No global differential inverse is required or
asserted.

F does not depend on H. At the regulated formal-integral level,
Gaussian integration in H therefore commutes with this substitution.
The resulting reduced action is S2[F]+S4[F], where S4 retains the
complete inverse heavy kinetic operator. No expansion in loop
momentum over heavy mass is permitted. The field-independent heavy
determinant is the same before and after substitution.

Write D_R for directional variation with an independent placeholder
R, so its repeated use below does not differentiate R itself. The
complete generated field-degree pieces are

    S2'  = S2,
    S4'  = S4 + D_R S2,
    S6'  = D_R S4 + (D_R^2 S2)/2,
    S8'  = (D_R^2 S4)/2,
    S10' = (D_R^3 S4)/6,
    S12' = (D_R^4 S4)/24.

The Lorentzian Hessian of S2 is -(Box+1), so its S6 term has the
corresponding minus sign. The native two-site calculation leaves
the free quadratic matrix general and checks all six degrees with
an arbitrary symmetric full heavy inverse. It is a polarization
check, not a discretized loop-momentum approximation.

For connected diagrams with E external light lines and vertices V_n,
the half-edge and Euler identities give

    L = 1 + [sum_n (n-2)V_n - E]/2.

At E=4,L=1 the only interaction topologies are two quartic vertices
or one sextic vertex. Thus the generated S6 tadpole cannot be dropped.
The S8 and higher terms start at two loops for four external lines.
They nevertheless remain part of the defining full action.

Transform all parent counterterms, with their existing fixed finite
parts. A one-loop quadratic counterterm generates a quartic tree
insertion D_R Sct2 in addition to Sct4. Generated sextic counterterms
in a tadpole would carry a second loop order. All such higher-order
terms are defined but not used as one-loop contributions. None of
the resulting higher-derivative insertions is resummed into a new
free propagator.
