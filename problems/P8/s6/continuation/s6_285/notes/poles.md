# Full prescribed subkernel and unmatched light coefficients

Put p=lambda^2+|q|^2 and kappa_HP=kappa+delta_kappa_HP.
The full H/Proca conserved sectors are

O2=p kappa_HP/kappa+p^2 A2(p)/(16pi^2 kappa),
O0=-2p kappa_HP/kappa+p^2 H0(p)/(192pi^2 kappa).

A2=(log n+2)/60 plus the H and Proca integrals
integral_0^1 p W2(v)/[4nu+p(1-v^2)]dv.
H0=2(log n+2) plus the corresponding W0 integrals.

For a scalar W2=v^6/30 and W0=v^2(3-v^2)^2.
For Proca W2=v^2(30-20v^2+3v^4)/30 and
W0=v^2(3-2v^2+3v^4). The absolute factors64pi^2 and768pi^2
relative to the current densities are checked independently. These
integrals retain the entire thresholds and exact fixed finite terms.

For the actual hierarchy, n>3>e and n<10^198. The series for e^3
already exceeds10 through degree3, so log10<3 and1<log n<600.
Using pi^2>9 gives

0<delta_kappa_HP<10^198,
0<delta_kappa_HP/kappa<10^-602,
0<1-kappa/kappa_HP<10^-602.

No arbitrary finite Newton subtraction is imposed.

On |p|<=1 the local TT correction to O2/p is below
(602/(960*9))/kappa<1/(10kappa). Its nonlocal part is bounded by
8(I2_H+I2_A)/kappa<1/(1000kappa).
For the trace quotient the corresponding bounds are7/(10kappa)
and1/(100kappa). The exact coarse arithmetic margins are checked.
Both total quotient errors are therefore<1/kappa. Since the fixed
constant terms are kappa_HP/kappa>1 and-2kappa_HP/kappa<-2,
neither quotient vanishes in this entire complex disk.
No additional H/Proca response pole is there.

The H/Proca TT response divided by the original Einstein response
has error at most(delta_kappa_HP+1)/(kappa-1)<4*10^-602, including
the continuous p0 ratio. This is a low-complex-momentum result for a
specified subkernel, not a physical cutoff or high-energy stability
theorem.

## Full three-cut kernel without choosing the missing light polynomial

Add the Phi nonlocal scalar integrals at nu1, and retain its unmatched
c_Phi_R,c_Phi_W,c_Phi_R2. In the full kernel they contribute

Delta O2=2c_Phi_R p/kappa-4c_Phi_W p^2/kappa
         +p^2 A_Phi,nonlocal/(16pi^2 kappa),

Delta O0=-4c_Phi_R p/kappa-24c_Phi_R2 p^2/kappa
         +p^2 H_Phi,nonlocal/(192pi^2 kappa).

Thus the complete parameterized Gaussian massless residue is

kappa/[kappa_HP+2c_Phi_R].

Its sign and a full inverse bound have not been established. The cut
cannot fix these finite coefficients. Setting them to the generic H
formula at nu1 would be an additional matching prescription, not a
consequence of the fixed Phi vacuum constant.

All three cosmological volume terms nevertheless cancel before metric
variation. For an insertion c0+delta p+c2 p^2+c3 p^3, the first-loop
inverse coefficient is-c0/p^2-delta/p-c2-c3 p. The source fixes c0=0;
the finite Newton simple pole remains. This proves absence of that
gapped cosmological double pole only, not removal of every branch or
pole-looking term in the complete scalar amplitude.
