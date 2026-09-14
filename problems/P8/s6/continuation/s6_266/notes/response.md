# Full implicit derivatives and nonlinear chart contacts

All derivatives below are evaluated at the actual entire-source
root (Nstar(z),z). Write d=C_N, b_i=C_Ni, c_i=C_i and e=C_NN.
Then

    N_i=-c_i/d,
    N_ij=-(C_ij+b_i N_j+b_j N_i+e N_i N_j)/d.

The packet retains all12 first and144 mixed second components and
checks the complete differentiated constraint and symmetry.

In a fixed parameter direction v, put
c=C_z[v], b=C_Nz[v], a=C_zz[v,v], f=C_NNN,
g=C_NNz[v], k=C_Nzz[v,v], r=C_zzz[v,v,v].
The full third derivative is

    N'''=-(r+3k N'+3g(N')^2+f(N')^3
             +3b N''+3e N' N'')/d.

This is checked by substituting the resulting cubic root into the
complete cubic Taylor polynomial of C. Polarization determines
every mixed third component of the symmetric tensor.
The literal constraint has degree at most2 in the twelve stated
invariants, so r=0 in that chart ONLY. A nonlinear canonical
chart need not have C_zzz=0; the general formula retains it.

Outward full-source bounds give
|C_Nz|<100, |C_NNz|<1000, |C_zz|<10 and |C_Nzz|<100,
with only seven nonzero entries in C_zz before implicit mixing.
Together with |C_N|>3 and the previous lapse-derivative bounds,
the first derivative is<2, the complete second derivative is
bounded by1470<2000, and the complete third by9353600/3<10^7.
All norms in these three statements are individual-component
ceilings. For arbitrary directions retain the corresponding
l1 norm factors; no dimension-free Euclidean tensor bound is
silently inferred.

For actual canonical field/spatial-jet coordinates w^A, write the
entire invariant map z^i(w). Then

    N_A=N_i z^i_A,
    N_AB=N_ij z^i_A z^j_B+N_i z^i_AB.

The second term includes the density, gradient-energy, shear and
vector quadratic contacts. It cannot be discarded by regarding
the invariant variables as new independent Gaussian fields.
The reusable pullback function keeps both terms. Independent
nonlinear finite-dimensional charts check it by direct second
differentiation. The third-order formula is likewise a full chain
rule, not permission to erase the nonlinear map's third jets.

These written identities and interval bounds are classical.
No canonical Gaussian state, Weyl ordering, interaction mean or
continuum product is assigned to the invariant chart.
