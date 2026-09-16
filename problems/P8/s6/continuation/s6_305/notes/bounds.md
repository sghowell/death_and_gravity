# Uniform known-reference bound, sign and endpoint

Use the original mu=1,25/4<=s<=16,t,u<0. Thus every |a|<=16.
The invariant formulas give |N0|<=108 and
|N2|<=2+32+256+54=344. The four nonnegative radial weights integrate to

int w_H2=1/210, int w_H0=68/35,
int w_V2=3/14, int w_V0=36/35.

For every complex |p|<=16, the denominator modulus is at least4nu-16.
Using n>=10^12 and m=10^6 gives

|A_non|<=16[(1/210)/(4*10^12-16)+(3/14)/(4*10^6-16)]<1e-6,
|H_non|<=16[(68/35)/(4*10^12-16)+(36/35)/(4*10^6-16)]<1e-5.

The inherited exact hierarchy and e^3>10 imply log(n)<600. Therefore
|A|<11 and |H|<1205. With pi>3,

|MHV_nonNewton| <
3[344*11/(16*9)+108*1205/(768*9)]/kappa^2
=25981/(192*kappa^2)<136/kappa^2.

S297 proves AG>257/(18kappa)>14/kappa throughout this compact set,
and the complete tuned matter amplitude Am is positive. In particular
the normalized gapped non-Newton interference divided by AG is<20/kappa.

For M1 put tau=min(-t,-u). Then0<tau<=6, while the other spacelike
magnitude is at least(s-4)/2>=9/8 and at most12. Its logarithm and the
physical timelike logarithm each have modulus<7. The remaining log has
modulus|log(tau)|. Since |4N2+10N0|<=2456, a convenient uniform bound is

|MM1|<2[1+|log(tau)|]/kappa^2.

Indeed3*2456*7/(3840*9)=2149/1440<2 bounds all three logs by
7[1+|log(tau)|]. When tau<=1, use AG>8/(kappa*tau) and
tau[1-log(tau)]<=1, whose derivative is-log(tau)>=0 on(0,1].
When tau>=1, use log(tau)<=log6<2 and AG>14/kappa.
Both cases give |2Re MM1/AG|<1/kappa. This is a uniform normalized
interference bound, not a bounded unnormalized massless amplitude.

Now let wG=AG/(Am+AG), which is strictly between0 and1. For the specified
known spectator reference, including the FIXED Newton coefficient,

E_spec=2Re(MHV+MM1)/(Am+AG)
      =wG[-2dk/kappa+e], with |e|<21/kappa.

S285's unchanged
dk=[n(log(n)-1)+5m]/(96pi^2)
satisfies0<dk<10^198. More specifically n>e implies
dk>5*10^6/(96pi^2)>11 using pi<4. Thus the bracket is strictly
negative. Its magnitude is below(2*10^198+21)/10^800<3e-602.
Consequently the entire compact nonforward domain has

-3e-602 < E_spec < 0.

This sign belongs ONLY to the fixed known reference. Unassigned physical
finite matching can change the full rate. It is not a positivity theorem
for a coefficient after removing a graviton pole.

At either forward endpoint, AG has positive simple-pole residue
(s^2-4s+2)/kappa; Am remains finite. Hence wG tends1. The gapped
non-Newton insertion stays finite, while M1 grows at worstlog(tau).
Both normalized remainders vanish, giving

lim E_spec = -2dk/kappa.

The unknown remaining Newton coordinate would add its own endpoint term.
The complex phase and transfer branch are not removed from the Regge
problem, and the forward Born cross section is still divergent.

## Separate Gaussian-only complex disk

For |p|<=16, the normalized TT and trace quotient remainders beyond
1+dk/kappa have moduli at most
16*11/(16pi^2*kappa)<6/kappa and
16*1205/(384pi^2*kappa)<6/kappa respectively.
Thus the fixed H/Proca quotient cannot vanish in this disk: its real
constant part is greater than1 and the remainder is less than1.
This extends the earlier Gaussian-only disk, not the physical cutoff.

Let r=(10^198+6)/kappa<1/2. For either quotient increment z, |z|<r and

|1/(1+z)-(1-z)|=|z^2/(1+z)|<r^2/(1-r)<5e-1204.

This controls repeated insertions of the SAME fixed Gaussian kernels
only. It says nothing about independent interacting higher-loop graphs,
unknown physical counterterms, the full quantum state or UV completion.
