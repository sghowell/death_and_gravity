# Full energy identity and finite chart conversions

For a real solution define Eenergy=(ydot^T K ydot+q y^T G y)/2. Differentiate before freezing any time coefficient and use the complete Euler equation. Since Omega is skew,
Eenergy'=-ydot^T K' ydot/2-3H ydot^T K ydot-ydot^T R y
         +q y^T(G'-2HG)y/2.
R need not be symmetric. The exact identity is checked both symbolically for arbitrary matrices and independently along a nontrivial curve with its computed forcing. Complex Fourier data are handled by summing the identities for real and imaginary parts.

Let m=1e-8, M=2e18. Positivity gives |ydot|²<=2Eenergy/m and q|y|²<=2Eenergy/m. Since q>=1, the complete coefficient bounds imply
|Eenergy'| <= [(M+2M+M+4M)/m+12] Eenergy <1e28 Eenergy.
The cross estimate has harmless extra slack. Thus sqrt(Eenergy(t))<=exp(1e28 |t-s|)sqrt(Eenergy(s)); the root-energy rate is safely twice the rate directly needed.

The phase conversion is explicit. Put lambda=1+|P|² and use a>=1, hence q<=|P|². From the original lapse formula and Jc>1/100,
|n| <= [135+100q+3/20000]|Z| <200lambda|Z|.
In the outer chart vdot=Theta n-ell sigma/2 and sigmadot=ps-wn. In the central chart b=-pv/(2q) and
bdot=-(1+9A/(2q))v-(E+3Tcorr/(2q))n-Hb,
with the same sigmadot. Including both velocity components, the stated coefficient box gives |ydot|<=500lambda|Z| and |y|<=2|Z| in either chart. Since K,G<1e4 I and sqrt(q)<=lambda, sqrt(Eenergy)<=50200lambda|Z|<1e12lambda|Z|.

Conversely, in the outer chart
n=(vdot+ell sigma/2)/Theta,
ps=sigmadot+wn,
pv=[2Jc n+w ps-3Theta ell sigma+(2Eq+3Tcorr)v]/Theta.
The absolute coefficient sums give |Z|<=5000(|ydot|+lambda|y|). For example the ydot coefficient sum is below3202 and the y coefficient sum below165+8q. Positivity gives |ydot|+|y|<=30000sqrt(Eenergy), so the inverse bound is below1e12lambda sqrt(Eenergy).

Centrally, x=(v,ps)=M^-T K[ydot+M^-1 S y] and pv=-2qb. The exact four entries of M^-1 S have the safe absolute bounds
(2+300q,45,20q,3).
The first includes the tiny Tcorr remainder using q>=4096. Their sum is below1000(1+q)<=1000lambda. Since ||M^-T||<=1 and ||K||<=1e4,
|Z|<=2e7(|ydot|+lambda|y|)<=6e11lambda sqrt(Eenergy).
An independent finite-q stationary solve and exact original-phase reconstruction check the complete formulas.

Therefore both directions of phase/energy conversion use the safe common A=1e12:
sqrt(Eenergy)<=A lambda|Z| and |Z|<=A lambda sqrt(Eenergy).
These are finite comparison estimates, not a coercive claim about the original Hamiltonian or an optimized stability constant.
