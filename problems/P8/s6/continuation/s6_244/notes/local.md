# Entire finite local current and unchanged profile

The convergent state-minus-adiabatic bound is only one term. Here the FULL same-prescription finite heat action and fixed profile are included, without relying on a cancellation to hide either term.

## Complex jet neighborhood for the full density

Put J=log h=2log(a)I+Q. The homogeneous heat density is an autonomous holomorphic function L(J0,J1,J2), with Jj the raw time jets and v=exp(tr J0/2). There are six independent complex symmetric matrix coordinates per jet. Allow each of the J0,J1,J2 coordinates to vary in a disk of radius1/12. In addition allow a complex source amplitude of radius1/100 along a normalized C12 direction.

Throughout these disks, conservative operator-norm bounds are

||J0||<4, ||J1||<6, ||J2||<10.

For the external jets in the Euler rule, ||J3||<18 and ||J4||<50. To see the latter bounds directly, |H''|<=8 follows from3x<=1+3x²<=(1+x²)³ for x=|t|<=1/2. Also |H'''|<=24 because |1-6t²+t⁴|<=1 on this interval and its denominator is(1+t²)^4. Add the history/source jet budgets to2H'',2H'''. The lower jet bounds have still larger margins after all six coordinate disks.

Use the coordinate mixed extrinsic matrix K=.5exp(-J0)Dexp(J0)[J1]. Its coordinate derivative is denoted A=K'. The curvature trace polynomials are unchanged by conjugating to a Fermi frame, so their bounds can be taken directly in these coordinates, with NO additional uncontrolled frame norm.

Complete matrix-exponential Frechet bounds give

||K||<3exp(8)<10^5,
||A||<41exp(8)<10^6.

The second estimate retains both differentiated exponential factors and J2: .5exp(8)(2*6²+10)=41exp(8). The complete invariant polynomials then imply |R|<2*10^11, |Ric²|,|Riem²|<10^22. The full volume obeys |v|<exp(6)<1000.

Retain394<ell<462 and n>1. The ENTIRE finite scalar density, not a selected leading term, has bound

1000*[462+154*(2*10^11)+(462/36)*(2*10^11)²+(462/90)*(2*10^22)]n²
<10^28 n².

Dropping the overall denominator64pi² only enlarges this positive bound.

## Complete Euler current and source derivatives

For each independent symmetric component i, the FULL Euler derivative is

E_i=L_(0i)-sum_(j,b)L_(1i,jb)J_(j+1,b)
+sum_(j,b)L_(2i,jb)J_(j+2,b)
+sum_(j,k,b,c)L_(2i,jb,kc)J_(j+1,b)J_(k+1,c).

Here j,k=0,1,2 and b,c run over all six symmetric entries. This is obtained by differentiating L_(0i)-d_t L_(1i)+d_t² L_(2i), retaining every product and contact.

Cauchy bounds first, second and third jet derivatives by12,288,10368 times the whole density bound. There are18,18,324 terms in the three sums. Using raw J1..J3<=18 and J2..J4<=50 gives the exact coefficient

12+18*288*18+18*288*50+324*10368*18²
=1088743692 <2*10^9.

Contract all six detector components and then take source derivatives0,1,2 on the separate radius1/100 disk. The largest factor is2!*100². Hence EVERY complete heat-current derivative considered here is below

6*1088743692*2*100²*10^28 n² ||D||F <10^44 n²||D||F.

The same argument retains all finite logarithmic, curvature, Euler and volume terms. It leaves no detector time derivative in the current bound.

## Entire fixed profile and the final sum

At N=1,X=1 the ENTIRE unchanged S240 profile has density -v P_H_ref(t). This is its full value on all homogeneous spatial histories, not a linearized profile. Its reference pressure is held fixed under the source variation, with |P_H_ref|<10^400 from S240.

The current is -v P_H_ref tr D/2. Its first/second source derivatives multiply by tr G/2 and(tr G/2)², respectively. With |v|<128 on the larger real background interval, all three full bounds are below432*10^400||D||F<10^406||D||F for normalized sources. The profile and its contacts are NOT zero individually.

Together with the all-momentum state-minus-adiabatic term, the entire current and first two directional derivatives have bound

[10^51/n+10^44 n²+10^406]/kappa <10^45 n²/kappa <10^-350.

At the unchanged reference Q=0, the full profile matches the complete reference spatial mean. This statement does not cancel the response, and it is not true of a finitely projected mean. All later clock-chart contacts must retain that distinction.

The first response therefore has the stated detector time L2 times prepared source H13 bound. The second directional derivative supplies a finite C2 control on PRESCRIBED homogeneous histories; polarization can be used for mixed directions with an explicit factor. Neither this finite neighborhood nor the numerical smallness proves a nonlinear feedback solution or compatible inverse.
