# Small action-normalized basis and both endpoint factors

Write A_c=sqrt(K_c omega_c), A_m=sqrt(K_m omega_m).
The same four-mode basis as S6.90, in order (+c,+m,-c,-m), is
S=(1/sqrt(2)) times
[ 1/A_c,         0,  1/A_c,          0;
 -ell/A_c,   1/A_m, -ell/A_c,    1/A_m;
 i A_c, i ell A_m, -i A_c, -i ell A_m;
 0,          i A_m,        0,    -i A_m ].
The native calculation checks this matrix and its inverse against
the full parent basis, not just its eigenvalues or principal cone.

On the actual complex disc, |N|,|e|,|R| are between .99 and 1.01;
the same is true of |c_clock|. The inherited pivot and M annuli
give 2<|A_c|<4. Since A_m=e/sqrt(R), .9<|A_m|<1.1.
The matter density has |ell|<1/8. These are complex moduli, with
branches fixed as in jets.md.

Use the full derivative
S'=S_{A_c} A_c (A_c'/A_c)
  +S_{A_m} A_m (A_m'/A_m)+S_ell ell'.
The exact Laurent entry majorants and row sums give
norm(S)<=331 sqrt(2)/80<6,
norm(S^-1)<=73 sqrt(2)/32<4,
norm(S')<=200407 sqrt(2)/20000000<1/50.
No derivative of the moving basis is dropped. The logarithmic
jet bounds apply on the whole complex disc, allowing the
subsequent Cauchy estimate, not just on its real diameter.

The scalar response takes the packet matrix's (chi,pi_chi)
entry, times P(s,t)/k, where
P=(R_s/R_t)^(3/2) N_s e_s^3.
For ordered real endpoints its square is at most
(101/99)^3*(101/100)^8<4, so 0<P<2.
The leading front amplitudes are P ell_s ell_t/(A_cs A_ct)
and P/(A_ms A_mt). Hence each is less than 3. Their strict
positive lower bounds remain the old S6.91 bounds.

These estimates preserve both endpoint volumes and both endpoint
basis normalizations. They do not replace the physical source
factor with one, transfer a vacuum canonical normalization,
or change the overall positive action prefactor of S6.93.
