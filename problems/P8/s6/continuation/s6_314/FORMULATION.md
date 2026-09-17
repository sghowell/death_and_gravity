# State-transport connector and matched real-tree sectors

Keep the unchanged original parameters n=10^200/512+2,g=1/8192,
kappa=10^800, tuned contact and full positive Born A0=Am+AG.
Use the S300 recoil at5/4<=E<=2, every nonforward hard Born direction,
all physical emitted directions and total marked energy R<=x<=1/8.
Retain the S301 dimensional angular/phase convention, unit physical
TT polarizations and the S313 phase-weighted F1,F2.

Let P_sigma(y)=exp(Delta_sigma)*exp(-gamma_E*a_sigma)
*y^a_sigma/Gamma(1+a_sigma), and P0 its elastic counterpart.
These are known soft reference factors, not complete hard amplitudes.
Let r1=density(|F1|^2-K0) be the unchanged S307 signed one-real
measure and D1 its unchanged S309 same-state leading-soft dressing.

For each two-real polarization pair put
 U=J_a(sigma_b)*F1(b)/sqrt(kappa),
 V=J_b(sigma_a)*F1(a)/sqrt(kappa), O=G00/(ab),
 r2_density=|F2|^2-|U|^2-|V|^2+|O|^2.
The measure dR2 includes both graviton phase measures, the physical
polarization sum and1/2! for identical gravitons. Define
 D2(x)=integral_(a+b<x) P_sigma_ab(x-a-b)dR2.
The full two-marked state sigma_ab is kept fixed within its additional
leading-soft sum. It is not sequentially changed after each unmarked
emission.

Let dB1 be the nonnegative elastic one-soft Born measure and define
 E1(x)=integral_(0<b<x) [P_sigma_b(x-b)-P0(x-b)]dB1.
The difference is taken before the seed integral. The individual
unsubtracted seed integrals are not assigned separate finite values.

The new theorems are:
 |Delta_sigma-Delta0|<5000R(1-lnR)/kappa;
 |b_e(sigma)-b_e(0)|<8000R(1-lnR)/kappa,
 b_e=(a_e-a)/(2e),0<e<=1/8;
 TV(R2;x)<2*10^-725*x+2*10^-652*x^2;
 |E1|/P0<10^5*x*(2-ln x)/kappa^2;
 |D2|/P0<4*10^-725*x+4*10^-652*x^2.
A uniform finite-regulator majorant justifies the E1 limit at fixed
positive x without assuming a_e>=0.

The defined expression P_match2=P0+D1+E1+D2 exactly matches the
zero/one/two-real physical D4 TREE densities, in the radiation-count
sense specified in notes/matching.md. Its relative difference from
P0 is bounded by
 2*min(10^32,2*10^14*x+2*10^37*x^2)/kappa
 +10^5*x*(2-ln x)/kappa^2
 +4*10^-725*x+4*10^-652*x^2
 <10^-653,
and tends uniformly to zero as x tends to zero. This implies
pointwise positivity of this named reference only.

No complete finite hard real-virtual or evanescent matching, all-N
nonleading error, positive event measure, threshold monotonicity,
unitarity, interacting quantum state, absolute complex Regge estimate,
common-parent bounce or original V/G/B/P8 closure is asserted.
