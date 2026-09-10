# Full two-field tensor normalization

Label the canonical fields by 0=Phi and 1=H. In the convention
V contains lambda_abcd phi_a phi_b phi_c phi_d/4! and
h_abc phi_a phi_b phi_c/3!, the nonzero interactions are

lambda_0000=L,
h_001=h_010=h_100=G,
m^2_ab=diag(m_light^2,M).

There are no gauge or Yukawa couplings. In units where the
one-loop beta coefficient is multiplied by 1/(16pi^2), the
pure-scalar tensor expressions are

B4_abcd=sum_ef [lambda_abef lambda_efcd
              +lambda_acef lambda_efbd
              +lambda_adef lambda_efbc],

B3_abc=sum_ef [lambda_abef h_efc
             +lambda_acef h_efb
             +lambda_bcef h_efa],

B2_ab=sum_ef [m^2_ef lambda_abef+h_aef h_bef].

The code contracts every internal index and checks all 16
quartic, eight cubic and four mass components. The only
nonzero results are

B4_0000=3L^2,
B3_001=B3_010=B3_100=LG,
B2_00=L m_light^2+2G^2,
B2_11=G^2.

Consequently beta_g=2G beta_G has coefficient 2Lg, not Lg.
The heavy mass coefficient is g. These match twice the
existing reference-bubble counterterm coefficients of S6.113,
because its radial bubble changes by -2 ln(nu).
The matching is checked against the actual imported ancestor
counterterm expressions.

The mass tensor is shown in full to avoid silently omitting
the light UV divergence. Its actual pole and residue remain
fixed by their existing on-shell subtractions. This diagnostic
tensor contraction is not a conversion of the whole scheme
to MS or an import of higher-loop coefficients.

At this loop order the cubic self-energy derivative is UV
finite and the quartic tadpole has no momentum dependence.
There is no one-loop UV wave-function anomalous dimension
from these interactions. The pre-existing finite on-shell
field normalization is retained once; changes of its couplings
contribute only at the uncomputed next loop order.
