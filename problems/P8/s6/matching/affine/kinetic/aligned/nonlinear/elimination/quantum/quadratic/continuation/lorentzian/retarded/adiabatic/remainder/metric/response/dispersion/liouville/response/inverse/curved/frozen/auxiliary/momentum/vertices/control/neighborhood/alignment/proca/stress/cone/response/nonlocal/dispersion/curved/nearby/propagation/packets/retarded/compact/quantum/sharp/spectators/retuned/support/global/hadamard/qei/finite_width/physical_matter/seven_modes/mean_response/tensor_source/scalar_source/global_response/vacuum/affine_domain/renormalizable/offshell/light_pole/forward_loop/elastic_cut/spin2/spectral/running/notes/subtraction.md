# Anchored reference change and full-channel compensation

For r=nu^2>0, the common regulated radial reference bubble
has the finite difference

I_r-I_1=integral_0^infinity y[(y+r)^(-2)-(y+1)^(-2)]dy
       =-ln(r)=-2ln(nu).

The written primitive is
ln(y+r)+r/(y+r)-ln(y+1)-1/(y+1).
It tends to zero at infinity and equals ln(r) at zero.
Its derivative, endpoints, and the independent differentiated
integral -2 integral y/(y+r)^3 dy=-1/r are checked.
Only the convergent difference is used; no divergent bubble
is assigned a finite unregulated value.

The actual tree amplitude is
A0=-L+g sum_z 1/(M-z), z=s,t,u.
The full-model subtraction in each channel is
-I_r C(z)^2/(32pi^2), C(z)=-L+g/(M-z).
With ell=ln(nu)/(16pi^2), the explicit derivative of the
renormalized one-loop term is sum_z C(z)^2.

Changing the tree parameters by

dL/dell=3L^2, dg/dell=2Lg, dM/dell=g

gives dA0/dell=-sum_z C(z)^2. Every channel and the
heavy-mass derivative are retained. Omitting dM leaves a
nonzero g^2 sum_z (M-z)^(-2) mismatch.

This is compensation through one-loop order. Differentiating
the couplings inside a one-loop term generates two-loop
contributions, which are not set to zero or bounded here.
The closed solution of these beta equations must not be
described as exact scale invariance of the full quantum amplitude.

The additional finite, momentum-independent potential contact
is fixed once at nu=1 by S6.113. It is transported as that
same fixed-order coefficient. Reimposing zero one-loop
potential quartic at every reference would instead change
the scheme and the flow just derived. Its coupling dependence
under this transport starts at the next loop order.

The light pole remains at mass one with unit LSZ residue
under the same on-shell prescription. Local quadratic and
one-point terms retain their original physical conditions.
Heavy reference parameters do not become exact stable-heavy
on-shell observables. No additional external-leg rescaling
or finite b2 counterterm is introduced.
