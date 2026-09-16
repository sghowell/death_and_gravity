# Explicit finite real-minus-soft error

Work in mass-one units and the unchanged S295 recoil coordinates:
5/4<=E<=2, s=4E^2, all emitted and hard rest-frame directions.
Let tau=min(-t,-u),delta=min(1,tau), and0<omega<=x<=delta/192.
The original recoil proof gives final energy shifts at mostomega,
spatial-vector shifts at most3omega,0<rho<=1 and|rho-1|<=2omega.
Every scalar component is at most2. For any pair transfer, the pair
four-vector changes by at most3omega componentwise: incoming pairs do not
change, mixed pairs contain one recoiling scalar, and the outgoing pair
changes by exactly-k. Both old and new pair components are at most4.
Therefore|D-D0|<=4*3*8 omega=96omega. Every original channel has|D0|>=delta;
hence|D|>=delta/2. This is why the resolution condition is necessary for
this proof.

Use entrywise absolute matrix bounds, and a physical TT polarization with
unit Frobenius norm, so every epsilon entry is at most1 and tr(epsilon)=0.
The old scalar stress is bounded by25 per entry. The trace-reversed
propagator numerator Hnum=p q+q p+eta is bounded by9.
A shifted emitting hard momentum differs from its Born value by at most
4omega per component; its partner differs by at most3omega.
Thus the hard stress difference is at most84omega, and the opposite
Hnum difference at most24omega. The contracted hard numerator satisfies
|N-N0|<=16(84*9+25*24)omega=21696omega and|N0|<=3600.
Consequently a channel quotient changes by at most
43392omega/delta+691200omega/delta^2
<=734592omega/delta^2.
All three hard channels change by at most2203776omega/delta^2,
strictly below2300000omega/delta^2.

For a physical massive scalar, |p.k|>=omega/4 and|p.epsilon.p|<=4.
Thus sum_i|J_i|<=64/omega. The two final soft currents differ from their
Born counterparts by at most608, using
|delta J|<=12*4+4*4*16=304 on each final leg.
The external hard-difference contribution is bounded by
64*2300000/delta^2, plus608 times G0.

Each internal field has entries at most18/delta.
For a seagull, |tr H|<=72/delta, |tr(epsilon H)|<=288/delta,
|d2|<=144/delta, |p.q+1|<=17 and |p.epsilon.q|<=4.
Each p epsilon H q term has at most64 products of size2*1*18*2/delta.
The six seagulls are therefore bounded by
6*4*(17*144+72*4/2+2*64*2*18*2)/delta=283392/delta.

For the cubic Einstein vertex all momenta have components at most4.
If field entries are bounded by A,B,C, then
|C1(B)|<=6B, |C2(A;B)|<=24AB, |M1(A)|<=3A.
A general bilinear has at most512 terms, or128 when its metric is diagonal.
The six permutations and overall4 give
4*6*(512*3*6*6+2*128*24*6)ABC=2211840ABC.
With two18/delta internal fields and one unit external field, all three
cubic graphs are bounded by2149908480/delta^2.

The unchanged positive gravity Born obeys G0=kappa A_G>8/delta:
near either endpoint this follows from S303's pole-residue lower bound
257/32>8; in the interior S297 gives G0>257/18>8.
Dividing the complete external, seagull and cubic budget by G0 gives
|M_G5-A_G S0|/A_G
< [608+(64*2300000+283392+2149908480)/(8delta)]/sqrt(kappa)
<=287174592/(sqrt(kappa)delta)
< B/(sqrt(kappa)delta), B=300000000.
S295 gives the matter analogue330000/sqrt(kappa). Since Am,AG>0 and
delta<=1, their positive Born-weighted combination satisfies the same
B/(sqrt(kappa)delta) bound for the full47-graph tree. This uses amplitudes
before squaring, so it retains every interference.

For each of two polarizations |S0|<=64/(sqrt(kappa)omega).
Write M5/A0=S0+R. Using rho<=1 and|rho-1|<=2omega gives
|rho sum|S0+R|^2-sum|S0|^2|
<= [256B/(omega delta)+2B^2/delta^2+16384/omega]/kappa.
The angle-integrated graviton measure is omega d_omega/(4pi^2).
The lower endpoint is integrable after this subtraction, with bound
[256B x/delta+B^2 x^2/delta^2+16384x]/(4pi^2 kappa).
At x<=delta/192 the numerator is at most7325418750256/3.
Since pi>3, the full result is below1831354687564/(27kappa)
<1e11/kappa=1e-789 at the unchanged original point.

This is uniform in emitted direction and valid at each nonzero hard
transfer, with angle-dependent resolution. It does not assert the same
bound when x is held fixed and delta tends to zero.

## Uniform extension to every fixed0<x<=1/8

The adaptive estimate above is not extrapolated outside its domain.
A direct physical gap bound controls the remaining radiation.

The final scalar energies never exceed E<=2. On the convex ball|p|<=sqrt(3),
the massive energy sqrt(1+|p|^2) has gradient norm at mostsqrt(3)/2.
For a mixed pair between an incoming and outgoing scalar,
-D=|Delta p|^2-(Delta E)^2>=|Delta p|^2/4 and
-D>=(Delta E)^2/3.
The exact recoil energy loss is
omega[1+/-beta_prime cos(theta)]/2>=omega/16, because
beta_prime<=sqrt(3)/2<7/8.
Consequently-D>=omega^2/768. For the corresponding Born transfer tau_j,
|Delta p|>=sqrt(tau_j)-3omega.
Ifomega<=sqrt(tau_j)/6, then-D>=tau_j/16.
Otherwise tau_j<36omega^2. These two cases prove
-D>(tau_j+omega^2)/30000; the worst denominator is768*37=28416.
The incoming and outgoing timelike pair invariants are at least25/4 and45/8.
These are physical lower bounds for every emitted direction.

Put y=sqrt(delta) and split the energy integral at a=y/192 if x>a.
Foromega<=a a mixed transfer satisfies
|D-D0|<=6sqrt(tau_j)omega+10omega^2<=7sqrt(tau_j)omega.
Its gap is at leasttau_j/2. Thus its quotient difference is bounded by
43392omega/tau_j+50400omega/tau_j^(3/2).
The timelike quotient difference is at most100992omega.
All three hard channels differ by less than310000omega/delta^(3/2).

The general cubic bound, before fixing its momentum maximum L, is
138240L^2 ABC. In a mixed channel,
L<=sqrt(tau_j)+3omega<2sqrt(tau_j), with internal field maxima18/tau_j.
Each mixed cubic is therefore at most179159040/tau_j; the timelike one
is at most179159040 because its field maxima are below9 and L<=4.
All cubic graphs are bounded by537477120/delta.
Combining external hard differences, seagulls and soft-current recoil,
and dividing by G0>8/delta, gives
|M_G5/A_G-S0|<69700672/(sqrt(kappa)sqrt(delta))
<1e8/(sqrt(kappa)sqrt(delta)).
The positive full-Born weighting also covers the330000 matter remainder.
The earlier finite-subtraction integration, now with B=1e8 and y instead
ofdelta, bounds the entire low-energy band by
610651562692/(81kappa)<1e10/kappa.

Foromega>=a use the direct full amplitude rather than a soft Taylor series.
The global physical gaps imply internal field maxima270000/(tau_j+omega^2).
The12 external graphs sum to at most
29859840000/[omega(delta+omega^2)].
The six seagulls sum to at most4250880000/(delta+omega^2).
In a mixed cubic, L<=sqrt(tau_j)+3omega and
L^2<=10(tau_j+omega^2), so the two mixed cubics sum to at most
201553920000000000/(delta+omega^2).
The timelike cubic179159040 is below358318080/(delta+omega^2),
since delta+omega^2<=65/64<2.
Divide by G0>8/delta and useomega<=1/8 to obtain
|M_G5|/A_G <D delta/[sqrt(kappa)omega(delta+omega^2)], D=1e16.
The unsimplified constant is3149283804498720<D.

The Born soft current also improves near forward scattering. Pair outgoing
and incoming momenta at the nearer Bose endpoint; their spatial difference
has normsqrt(tau). The numerator changes by at most4sqrt(tau), the
eikonal denominator by at mostsqrt(tau), and each denominator is at
least1/4 after removingomega. Each paired current is at most
80sqrt(tau)/omega, so|S0|<=160sqrt(delta)/(sqrt(kappa)omega);
the original64/omega bound suffices when delta=1.

Let wm=Am/(Am+AG),wg=AG/(Am+AG), both positive and summing to1.
Convexity gives|wm X+wg Y|^2<=wm|X|^2+wg|Y|^2.
With|M_m5/Am|<=|S0|+Bm/sqrt(kappa),Bm=330000, andrho<=1,
the high-band absolute real-minus-soft integrand per polarization is at
most3|S0|^2+2Bm^2/kappa+|M_G5/A_G|^2.
This bounds all mixed interferences, rather than discarding them.
After summing the two polarizations, its coefficient is
6*160^2 delta/omega^2+4Bm^2+
2D^2 delta^2/[omega^2(delta+omega^2)^2], all divided bykappa.

The needed integrals obey
delta log(x/a)<=delta log(24/sqrt(delta))<9/2,
using log24<4 and-delta log(delta)<=1/e<1.
After z=omega/sqrt(delta), the gravity integral is bounded by
integral_(1/192)^infinity dz/[z(1+z^2)^2]
<log192+1/4<25/4.
Split this at z=1 and bound the two pieces by1/z and1/z^5.
Finite positive Taylor sums prove exp4>24 and exp6>192.
The constant matter term contributes at mostBm^2/32.
Using4pi^2>36, the high-band coefficient is at most
312500000000000000000000850954050/9.
Adding the low-band coefficient gives
2812500000000000000000618310149142/81 <1e32.

Thus for every fixed0<x<=1/8 the complete selected47-graph finite
real-minus-leading-soft normalized rate error is below1e32/kappa,
uniformly over all physical nonforward hard angles.
At the original point this is1e-768.
The theorem concerns the normalized finite correction, not the divergent
forward Born cross section, an uncomputed hard matching term, or all-N
and higher-loop radiation.
