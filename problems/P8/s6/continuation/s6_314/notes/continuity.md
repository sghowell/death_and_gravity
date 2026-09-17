# Total-energy continuity of the full finite angular conversion

Use S301's fixed domain t in[0,1],n on S2 and transverse projected
current T_sigma(t,n). All norms here are Euclidean nuclear norms.
They bound both Frobenius norm and absolute trace. The comparison
state0 is elastic at the same incoming energy and outgoing pair-rest
direction. All momenta in the current use their proper signed
incoming/outgoing convention.

Only the two outgoing massive legs move. Each has spatial change<=3R
and energy change<=R; its Doppler denominators are>=1/4. Orthogonal
projection does not increase vector norm. With massive spatial
norms<2, the rank-one numerator difference is<=12R. Dividing by
the Doppler denominator costs<=48R. The inverse denominator change
is<=(R+3R)*16=64R, and multiplying by a numerator<=4 costs<=256R.
Thus each massive current changes by<=304R. All null marked currents
have total nuclear norm<=2R, independently of their multiplicity,
radial location or relative directions. Therefore
 ||T_sigma-T0||_*<=610R<640R.

The existing norms ||T_sigma||_*<=65/4 and ||T0||_*<=16 give, for
H(t)=mean[tr(T^2)-tr(T)^2/2],U(t)=mean[tr(T)^2],
 |delta H(t)|<=(3/2)*(65/4+16)*640R=30960R<32000R,
 |delta U(t)|<=(65/4+16)*640R=20640R<21000R.
These bounds hold before the angular average as well.

The separate S301 radial Holder bounds hold for both states:
|H(t)-H(1)|<32000(1-t)^(1/4),
|U(t)-U(1)|<21000(1-t)^(1/4).
Combining the two available estimates, the STATE difference of the
H radial difference is<64000*min(R,(1-t)^(1/4)); for U it is
<42000 times the same minimum. This joint estimate is essential;
a constant O(R) bound alone cannot be integrated with1/(1-t).

Writing tau=1-t and splitting at tau=R^4 gives exactly
 int_0^1 min(R,tau^(1/4))/tau dtau=4R(1-lnR).
The K1 formula retains BOTH its trace boundary term and radial
difference. Consequently
 |delta K1|<10500R+256000R(1-lnR)<300000R(1-lnR).
Together with |delta K0|<1440R,|p'(0)|<4 and pi>3,
 |Delta_sigma-Delta0|
 <(300000+4*1440)R(1-lnR)/(72kappa)
 <5000R(1-lnR)/kappa.

For the entire finite-e continuation use the exact S301 radial formula.
Its c(e)<=e/2,A(e)<=2 and e<=1/8 imply
 |delta K_e-delta K0|/e
 <[10500+2*256000+168000/8]R(1-lnR)
 =543500R(1-lnR).
The phase obeys0<p(e)<=1 and |p(e)-1|<=4e.
Thus for b_e=[p(e)K_e-K0]/(8pi^2*kappa*e),
 |delta b_e|<549260R(1-lnR)/(72kappa)
            <8000R(1-lnR)/kappa.
The modulus R(1-lnR) is continuous with value0 at R=0. No exchange
of angular differentiation and integration is based on boundedness
alone, and no finite-e positivity of the continued contraction is used.

Thirty-six exact original recoil-current configurations calibrate the
component and quadratic estimates, including a collinear bounded
representative. Three physical boundary-current differences from the
elastic current are explicitly nonzero. These finite tests do not
replace the nuclear-norm and joint Holder proof.
