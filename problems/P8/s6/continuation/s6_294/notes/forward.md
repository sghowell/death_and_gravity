# Entire finite mixed amplitude and explicit forward majorant

The known loop representative is conditional on the inherited S6.278 analytic soft division.
The offshell/graph probe passed25named30scalar, including D4/5/6 full tensor
contractions and twelve labelled boxes with each of six ordered channels twice.
S288, not this packet, already proves physical-forward Bsoft is pure imaginary.

## Exact IR subtraction of the entire mixed box

Let q_a=n-a, A_b(v)=mu-b*v*(1-v), C_b(u,v)=u^2 A_b+n(1-u),
Delta=u^2 A_b+h[n-a(1-u-h)] on0<u,v<1,0<h<1-u.
All are positive on the subthreshold neighbourhood used below.
D0(a,b)=Gamma(2-e)(4pi nu^2)^(-e) integral u Delta^(e-2) dh du dv.
Integration by parts in h using Delta_h=q_a+a(u+2h) gives EXACTLY

D0(a,b)=Gamma(1-e)(4pi nu^2)^(-e)/q_a *
[J_b(e)/(2e)-I_b(e)-(1-e)a K_ab(e)],

J_b(e)=integral_0^1 A_b^(e-1)dv,
I_b(e)=integral_0^1du dv u C_b^(e-1),
K_ab(e)=integral_0^1du dv integral_0^(1-u)dh
        u(u+2h)Delta^(e-2).
The complete boundary h=1-u is C_b and is independent of a.
K gains one soft power; its value and required a/b/e derivatives are
integrable in a fixed compact positive-gap neighbourhood. At e=0,
T(b)=-I_b(0), exactly the active-light matter triangle, not a new master.
No finite term is guessed from a pole coefficient.

## Complete finite nonendpoint bracket

Define Lll(a)=integral ln[mu-a v(1-v)]dv,
Lmix=integral ln[nv+mu(1-v)^2]dv,
Lh(a)=integral ln[v(n-a(1-v))]dv.
The full offshell U(a) at e=0 is
-integral_0^1dr dv /[(n-a)v+r(a v+mu(1-v)^2)].
The full T(a) at e=0 is
-integral_0^1dz dv (1-z)/[nz+(1-z)^2(mu-a v(1-v))].
J0(a)=integral1/A_a, J1(a)=integral ln(A_a)/A_a.
V(a)=(a-2mu)^2-2mu^2, V_e'(a)|0=2mu^2.
The noncusp Gamma coefficient simplifies, before expansion, to

K_a(e)=-2(a-2mu)Mll(a,e)/q_a+4n Mmix(e)/q_a
 -2mu^(1+e)/[(1+e)q_a]
 +{4[-a^2/(2+2e)+2na-n^2]Mh(a,e)
   -2a n^(1+e)/(1+e)}/q_a^2.

M's are whole bubble parameter powers with their displayed polynomials.
K_a(0)=2mu/q_a. Its complete first coefficient is

K1(a)=-2(a-2mu)Lll(a)/q_a+4n Lmix/q_a
 -2mu(ln(mu)-1)/q_a
 +{2a^2+4[-a^2/2+2na-n^2]Lh(a)-2an(ln(n)-1)}/q_a^2.

All proper, selfenergy, four external residues and contact-contact terms
are included in this K, not an ultraviolet-only subtraction.
The complete cusp+box pole is
[sum_a 1/q_a]*[sum_b V(b)J0(b)-2mu]/e
=2*Htree*Bsoft/e. Therefore the S278 division subtracts exactly the whole
heavy Born amplitude's soft pole.

Put c=EulerGamma-ln(4pi E^2), and channels(s,t,u),s+t+u=4mu.
After that stated conditional analytic soft division, the ENTIRE
nonendpoint minimal g^2/kappa finite amplitude is
g^2/(16pi^2 kappa) times F, where

F=Htree*{sum_b[V(b)J1(b)+2mu^2 J0(b)]
        +c*[sum_b V(b)J0(b)-2mu]}
  -sum_a K1(a)
  +sum_a{4(a-2mu)T(a)-4[a^2-2nmu]U(a)/q_a}
  +2sum_(a!=b) V(b)[T(b)-a K_ab(0)]/q_a.

This formula is obtained from whole Feynman graphs, not reconstruction
from forward cuts or a high-energy contour assumption.
Add the whole S290 g^2/kappa matter-graviton endpoint, including both
mass assignments and H-metric bubble. Its finite R H anchor is unassigned.
The covariant D-dimensional R H subtraction must precede D->4;
subtracting only its D4 contracted scalar leaves an extra mu/3.
There is no choice of local/higher-EFT finite coefficients here.

## Quantitative bound at original mu=1, n>=8, |v|<=1/2

Use s=2+v,t=0,u=2-v, E^2 in[1/4,1], L=ln n.
All varying channels lie[3/2,5/2], all channels[0,5/2].
A>=3/8, q>=n/2, |ln A|<1, |c|<3.
Bounds for derivatives0,1,2 in the channel variable:
J0 <=(3,2,3), |J1| <=(3,4,6);
|V| <=2, |V'|<=4,V''=2;
|Lll|<=1, |Lll'|<=2/3, |Lll''|<=4/9.
The Htree product in F has second v derivative <=2016/n.

0<=Lmix<=L. Also0<=Lh<=L, |Lh'|<=1/n,
|Lh''|<=2/n^2. Differentiating the four displayed K1 groups gives:
summed first group<=38/n;
second<=16L/n; third<=1/n;
fourth<=(56L+71)/n.
Thus |sum K1|'' <=(110+72L)/n.

The T integral split at z=1/2 gives |T|<=(L+4)/n:
on the first half Delta>=nz+1/16; on the second Delta>=n/2.
|T'|<=(2/3)|T|, |T''|<=(8/9)|T|.
Hence the summed4(a-2)T term second v derivative<=20(L+4)/n.

For U split v at1/2. On the first half denominator>=nv/2+r/4;
on the second it is>=n/4. Thus |U|<=2(L+3)/n.
Its a derivative ratios are <=2/n, so
|U'|<=2|U|/n, |U''|<=8|U|/n^2.
The coefficient cU=-4(a^2-2n)/(n-a) obeys
|cU|<=23,|cU'|<=125/n,|cU''|<=60/n.
The summed cU*U second v derivative<=73(L+3)/n.

For K_ab, Delta>=(3/8)u^2+(n/2)h.
Integrate the numerator u^2 over h up to infinity; its bound is16/(3n).
For the2uh term keep h<=1, integrate h exactly and then u:
it is<=4(L+3)/n^2. Thus K<=(L+14)/(2n).
Along any ordered pair of changing channels, |Delta_v|/Delta
<=2/3+2/n<1 and Delta_vv=0, giving
|K_v|<=2K,|K_vv|<=6K.
For d_ab=2a V(b)/q_a:
|d|<=20/n,|d'|<=56/n,|d''|<=92/n.
The six d*K second derivatives sum<=164(L+14)/n.
For d0_ab=2V(b)/q_a:
|d0|<=8/n,|d0'|<=18/n,|d0''|<=17/n.
The six d0*T second derivatives sum<=37(L+4)/n.

TOTAL |F''(v)|<=(4869+366L)/n.
At the exact original n=10^200/512+2, ln n<500, so
|F''|<200000/n. The forward coefficient is F''(0)/2:
|b20_nonendpoint|<2000*g^2/(kappa*n).
Adding S290's known endpoint bound gives
|b20_known_minimal_g2_gravity|<2001*g^2/(kappa*n)<10^-1001.
The S293 C sector adds less than g^2/(kappa*n) if combining the two.
These are conservative bounds for selected known one-loop representatives,
not the unmatched R H/heavy residue, higher operators, all-loop corrections,
a physical IR-safe observable, fixed-transfer Regge remainder or P8 closure.

The exact derivative identities and26 nonnegative rational safety margins are checked independently. The final original-parameter exponent is an exact integer/rational inequality, not a floating-point quadrature estimate.

## Differentiation and regulator boundary

The box integration-by-parts sign follows directly from
-d_h Delta^(e-1)=(1-e)[q_a+a(u+2h)]Delta^(e-2);
the aK term therefore enters with a minus sign. The rejected plus sign
has the explicit nonzero defect2(1-e)a(u+2h). The science tests check
the actual implemented K coefficient as well as this identity.

The positive compact gap bounds hold on a neighbourhood of v=0,
not only at the center. Each invariant derivative of K has a bounded
Delta-derivative ratio; epsilon derivatives add logarithms but keep
the gained soft power. Split the soft corner as in the displayed
u/h majorant. Integrable powers times finitely many logarithms dominate
these derivatives for0<=e<=1/8, with a finite n-dependent factor at
large Delta. The same argument applies to the integrable offshell U
corner and the gapped T/bubble moments. Thus the finite limit and
two forward derivatives commute after the exact subtraction.

This proves the known coefficient and its stated magnitude. It does
not give an evaluated physical detector-error bound, interchange a
heavy-resonance or physical cut limit, or supply a fixed-transfer
complex-energy Regge estimate.
