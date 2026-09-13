# Symmetric on-shell value, strict loop sign and the new finite contact

Use the point s=t=u=s0=4/3 and external virtualities1. It is on shell and crossing symmetric, but subthreshold rather than real scattering kinematics. All parameter integrals below are convergent positive real integrals.

Write n=D+2 and A=C+g²/(n-s0). The exact values are

A=-2g²(3D²-2)/[D²(3D+2)]<0,
A_tree=C+3g²/(n-s0)=4g²/[D²(3D+2)]>0.

The entire first-loop value is

A1_sym=3[A²B/2+2Ag²C_sym+2g4D_sym]/(16pi²).

All six boxes have the same value here because their two arguments coincide, not because their mass order can be interchanged generally.

## Complete positive parameter bounds

With L(z)=2/3+D z, exact polynomial decompositions give

Q_C-L=s0(1-z)²(xi-1/2)²+(2/3)z+(2/3)z²,
Q_D-L=s0(1-z)²(xi-1/2)²+s0 z²(eta-1/2)²+(2/3)z+z²/3.

Thus both denominators are bounded below on the FULL parameter cubes. The bubble integrand has0<=s0x(1-x)<=1/3, so0<B<1/2 using -ln(1-u)<=u/(1-u). Integrating the complete triangle and box majorants yields

0<C_sym<ln(1+3D/2)/D,
0<D_sym<ln(1+3D/2)/D².

For the box, integral z/L²=[ln(1+3D/2)+(2/3)/(D+2/3)-1]/D² is strictly smaller than the logarithmic bound. The actual1+3D/2<10^198 and ln10<7/3 give logarithm below462. Also abs(A)<3g²/D.

The triangle has a useful independent LOWER bound: Q_C<1+(n-1)z almost everywhere, so C_sym>B21(n). Since n>10^197 and ln10>2, ln n>394, and B21(n)>393/n. The actual n<(101/100)D therefore gives C_sym>39300/(101D)>375/D. The elementary ln10>2 follows from e<3 and e²<9<10; the factorial-series geometric tail proves e<3.

Finally3D²-38D-40>0 at the actual D, so abs(A)>(19/10)g²/D. Combining the positive upper and negative lower contributions without dropping a term,

A²B/2+2Ag²C_sym+2g4D_sym
 <[9/4-2(19/10)375+924]g4/D²
 =-(1995/4)g4/D².

Hence A1_sym<-5985g4/(64pi²D²)<0. With pi²<10, its magnitude divided by the exact positive tree exceeds5985g²(3D+2)/2560>10^190.

## Explicit new condition and classical comparison

The separately named V2S-T1-OS4 prescription adds deltaC_fin=-A1_sym>0. It fixes the one-loop value at this point. No finite derivative or full angular condition is implied, and the S234 base scheme is not rewritten.

The complete triangle inequality gives

abs(A1_sym)<3g4[9/4+8*462]/(16pi²D²)
          <14793g4/(192D²).

Consequently0<deltaC_fin/(24q)<14793g²(D+2)/[768(D-1)]<10^-6. In a CONTACT-ONLY classical comparison, C to C+deltaC_fin changes the full heavy-square remaining quartic from q to q-deltaC_fin/24>(1-10^-6)q. The whole comparison potential remains nonnegative and coercive with its classical origin minimum. This does not claim that all finite source/mass counterterms keep the old bare minimum, or that the full quantum effective potential is controlled.

The unadjusted loop is large relative to a finely cancelled observable, not a proof of strong coupling or no UV parent. The new condition is essential value matching at this order; omitted-loop errors and the rest of the amplitude still require calculation.

## One-dimensional evaluation and retained box enclosure

Let a=(1-z)²+nz, b=s0(1-z)², A0=a-b/4 and r=sqrt(b)/(2sqrt(A0)). Complete angular integration gives

J1=2atan(r)/sqrt(A0 b),
J2=1/[2A0(A0+b/4)]+atan(r)/[A0 sqrt(A0 b)].

Their removable b0 limits are1/A0 and1/A0². Then C_sym=integral(1-z)J1 dz and Dbar(s0,0)=integral z(1-z)J2 dz. The omitted heavy-angle fraction in the latter is bounded by s0/(4n), giving the EXACT enclosure

Dbar(s0,0)<=D_sym<Dbar(s0,0)/(1-s0/(4n))².

D_sym is retained inside this enclosure; it is not replaced by its lower endpoint. For independent numerical diagnostics, the logarithmic coordinate w=ln[1+(n-1)z] resolves the heavy endpoint layer. Those values check the formulas but do not replace the preceding continuum proof.
