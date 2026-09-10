# A fixed contact from the literal full Hessian

Use the same polynomial potential

    V=Phi^2/2+M H^2/2+G H Phi^2/2+L Phi^4/24.

Differentiate its full two-field Hessian before
putting H=H_*=-G Phi^2/(2M), and use g=G^2.
The light Schur block at radial y is

    y+1+F(y) Phi^2,
    F(y)=A-g/(y+M), A=(L-g/M)/2.

This identity is checked literally and against
the frozen S6.110 kernel. At the actual parameters,
L-3g/M>0, so 0<F(y)<A at every finite y>=0.

The quartic Taylor derivative of
(1/2)ln[1+F(y)Phi^2/(y+1)] is
-6F(y)^2/(y+1)^2. S6.110 fixes its cancellation.
The later on-shell mass condition changes a
quadratic term, not that quartic condition.
A one-loop change in the stationary heavy
solution affects the potential value only at
higher order because the tree potential is
stationary there.

At zero external momenta the full vertex is

    V0(y)=C0+2g/(y+M)=-2F(y),
    C0=-L+g/M=-2A.

The S6.113 local reference already subtracts
the high-momentum A^2 part with the entire I0.
The remaining finite potential contact is fixed as

    delta L_fin
      =6/(16pi^2) integral_0^infinity
                  y[F(y)^2-A^2]/(y+1)^2 dy.

It is strictly negative. Indeed F^2-A^2 is
negative everywhere and the remaining radial
weight is positive for y>0. It is not a free
choice made after inspecting b2.

The exact decomposition

    F^2-A^2=C0 g/(y+M)+g^2/(y+M)^2

reduces the contact to the two convergent
integrals in notes/radial.md. Vacuum and
quadratic counterterms do not change this
quartic Taylor cancellation.
