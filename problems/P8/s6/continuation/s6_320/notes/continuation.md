# Multiplicity-independent complex recoil and hard cores

## Full hard source in a relative epsilon=1e-12 tube

For the full physical recoil map, use epsilon=1e-12.
The perturbation of total energy and spatial momentum is at most
epsilon*W, even at arbitrarily large N. These are smaller than the
S313 bounds2cW used to prove the recoil square-root tube. Its
massive-momentum bounds |dp0|<5epsilon W,
||dpvec||<13epsilon W, unchanged branches, component bounds<3
(light) or<5 (heavy), and |sqrt(rho)|<2 therefore still hold.

For any light-scalar subset, let W_S be its REAL center energy.
The real gap is>=3W_S/8. The perturbation of its kinetic denominator
is<100epsilon W_S: expand2p.Q+Q^2 and use
||dp||<14epsilon W, ||Q||<=sqrt(2)W_S,
||dQ||<=sqrt(2)epsilon W_S, ||p||<=4 and W<=1/8.
Thus every complex light inverse is<4/W_S.
This statement includes compound subsets, not only singleton axes.

The same hard-channel argument retains
 |Dmixed|>(tau+W^2)/600000,
 |Dtimelike|>=45/16, |Dheavy|>=n/2.
A hard vertex has maximum component<=sqrt(tau)+5W, so its
Euclidean momentum square over the next hard denominator is
<4*50*600000=120000000. The one initial hard inverse keeps
600000/(delta_Born+W^2). Use the positive REAL center Born pieces;
they are not analytically redefined or replaced by absolute sectors.

The S319 maximal-soft/core bijection and energy charge ownership
are combinatorial and unchanged. Charge4/W_block at a scalar
vertex, and add1/W_block for every uncharged block. The complex
pure-current bound above retains the needed W_block^(m-1) factor.

The new core functions, with x marking a complete pure-soft block, are

 L=(1-1024x)/(1-5120x),
 D=1/(1-4x),
 H=(1-1024x)/(1-3072x),
 E=1024/(1-1024x)^2,
 V=120000000*1024*(2/(1-32x)^3-2),
 Fm=L^4*(D+(3/2)D^2H),
 FG=(3*600000/8)*L^4 E^2/(1-V).

Their sequence/set coefficients are nonnegative. Compose x=d(z).
At cap1e-17,radius5e-18 exact arithmetic gives
 radius+819200*((1-416cap)^(-2)-1-832cap)<cap,
 V(cap)<1/1000, Fm(cap)<3, FG(cap)<3e11.
The original3n^2 exceeds3e11. Therefore the complete continued
finite tree obeys, throughout the closed relative polydisc,

 |M_N(z)|/A0
 <=3n^2*N!*(2e17)^N/[kappa^(N/2)*product_i(wi)].

Here wi are positive real center energies, not the complex zi.
The equality of original and complete temporal representatives
continues algebraically away from all internal poles. No source
type guard is weakened: a successor can own a separate exact
complex evaluator while preserving every frozen real point API.
