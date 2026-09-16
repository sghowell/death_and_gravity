# Original recoil and propagator bounds

Use the S300 exact recoil against Q=(W,v), |v|<=W<=1/8,
5/4<=E<=2, and a fixed unit outgoing rest-frame direction.
All massive energies are at most2; each outgoing energy loss is
W/2 +/- rprime(u.v)/(2Eprime) >= W/16, since rprime/Eprime<7/8.
The spatial recoil shift from the associated elastic Born is at most3W.

For an incoming/outgoing scalar pair, let tau_j be its positive elastic
transfer magnitude and D its actual pair square. On |p|<=sqrt(3), the
massive energy gradient is bounded bysqrt(3)/2, so
-D>=|Delta p|^2/4 and -D>=(Delta E)^2/3>=W^2/768.
If W<=sqrt(tau_j)/6 then -D>=tau_j/16; otherwise tau_j<36W^2.
Both imply -D>(tau_j+W^2)/30000. The timelike pair invariants are at
least25/4 (incoming) and45/8 (outgoing).

For an external scalar emitting the entire Q,
|2pi.Q+Q^2|>=3W/8: incoming pi=-ki gives
2ki.Q-Q^2>=W/2-W^2>=3W/8; outgoing pi gives at leastW/2.
Every scalar component is at most2, every Q component at most1/8.
The emitted and shifted-hard scalar stress components are each below27.
Every trace-reversed internal graviton numerator has components at most9.

## Complete matter hard-current budget

The original |C|<4g^2/n, n>=128, 1/(n-a_channel)<2/n, and the
internal heavy stress has components below n+96<=2n. Thus each hard
current component, with1/sqrt(kappa) factored out, is bounded by:

- C external:4*72/W*4g^2/n=1152g^2/(nW);
- heavy external:3*4*72/W*2g^2/n=1728g^2/(nW);
- C density, heavy internal and heavy density:at most
  (4+3*8+3*2*2)g^2/n=40g^2/n.

At W<=1/8 the sum is at most2885g^2/(nW), below4096g^2/(nW).
The spatial Frobenius norm is at most three times a component bound.
Using Am>4g^2/n^3 gives

||Um_spatial||F/Am <3072n^2/(sqrt(kappa)*W).

At the actual n=10^200/512+2,kappa=10^800, the coefficient is just
above3/256 and strictly below1/64. No cancellation of large intermediate
matter graphs or floating-point subtraction is used in this bound.

## Complete Einstein hard-current budget

Let delta=min(1,min(tau_t,tau_u)), d=delta+W^2<=65/64.
Every internal trace-reversed metric field has component bound270000/d.
The12 external-emission graphs give the coefficient

12*72*16*27*270000 =100776960000

multiplying1/(kappa^(3/2)*W*d).
For a general unit-entry root basis tensor E and internal field bound B:
|trE|<=4, |trH|<=4B, |inner(E,H)|<=16B, |d2|<=12B.
The scalar seagull is bounded by
4(17*12+128+128+256+256)B=3888B.
All six give6298560000/d before the same coupling factor.

The cubic Einstein component estimate, with maximum momentum component L,
is138240 L^2 ABC. It does NOT assume the root tensor is TT.
For either mixed channel, L<=sqrt(tau_j)+3W and hence
L^2<=18(tau_j+W^2); its two internal fields each have bound
270000/(tau_j+W^2). Both mixed cubic graphs therefore contribute at most

2*138240*18*270000^2/d =362797056000000000/d.

The timelike channel has L<=4 and both fields below2, giving8847360.
Multiplying by the harmless upper bound65/64 for d gives8985600/d.
Since W<=1, combining all these into the W*d denominator is valid.
The total coefficient362797163084505600 is below10^18. Consequently

||UG_spatial||F <3*10^18/(kappa^(3/2)*W*(delta+W^2)).
AG>8/(kappa*delta) then implies
||UG_spatial||F/AG <(3/8)*10^18/(sqrt(kappa)*W) <1/(64W).

## Uniform original full-Born result and its limit

Both Am and AG are positive. Weighting the two current bounds by their
positive Born fractions gives

||U_spatial||F/(Am+AG) <1/(64W).

Together with the exact pair-vertex coefficient budget this proves

|M_pair|/(Am+AG)
 <(530/64) W/(sqrt(kappa)*a*b)
 <9*10^-400*(1/a+1/b).

The bound is uniform in emitted directions, energy sharing and all
nonforward elastic hard angles in the stated compact window.
It is deliberately conservative and applies to the isolated47-graph
class only. It is not the full nonleading two-real correction.
Its explicit1/a+1/b behavior is NOT square-integrable at both soft
endpoints without the correct overlap subtraction/virtual prescription.
No integrated finite detector error, all-N theorem, uniform full
amplitude remainder, matching value or original P8 closure follows.
