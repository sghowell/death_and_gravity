# Exact coefficient and finite-amplitude majorants

Assume mu=nu=1,25/4<=s<=16,t,u<0,s+t+u=4,
0<delta<=1,min(-t,-u)>=delta. PutL=|ln delta|.
For any channela and otherb,c, use

delta<=|a|<=16,2<|a-4|<=20,
|z=1+2b/(a-4)|<=17,|H=a/(a-4)|<=8.

These deliberately loose inequalities hold for all three rows.
|V_a|<=322, V_a^2<110000, and4|V_a|<1300.
The compact massive-triangle factor is
h=-2a+6-2/a+bc/a, so|h|<=296/delta.
Consequently |cmm|<=190624/delta<200000/delta,
and|cmm1|<=4(296+322)/delta=2472/delta<2500/delta.

For the complete massless angular polynomials use the weighted
coefficient l1 norm with|H|<=8,|Z|<=17 atN=2. It gives

|c00|<=207864500<210000000,
|b00|<=6718985/3<2300000.

This bounds every polynomial term; no physical angle is sampled.
For bmm use its entire S284 rational expression. The intermediate
quantities obey|p|<=644,|b|<=100/delta,|a_coeff|<=140/delta.
The exact resulting bound is84471200/(3delta^2)
<30000000/delta^2.

## Master, completion and leg bounds

From notes/masters.md,|J|<=200,|Q|<=1800,|L_a|<=18.
Here0<ell=ln4pi-Gamma_E<3. The elementary estimates
3<pi<22/7,0<Gamma_E<1,ln100<5,ln16<3,ln12<3 andln10<3
suffice. The logarithmic upper bounds follow from lower partial sums
of the positive exponential series at3 and5; the displayed pi and
Euler-constant intervals are their standard integral bounds.

Thus|Cmm_finite|<=1200,|Bmm_finite|<=21,
|ell_a|<10+L,|B00_finite|<12+L, and
|C00|<10(1+L)/sqrt(delta).

The whole Newton numerator satisfies|N_a|<=290, hence
|T0|<=870/delta and|T1|<=6/delta. Also
|B_soft|<=3*322*200/2+1=96601<100000.

After replacing1/delta by1/delta^2 and constants by their common
larger denominator, the four physical-completion basis magnitudes
are at most(3,3,768,3072)/delta^2. Therefore

|R0|<=19878/(5delta^2)<4000/delta^2,
|R1|<=822059/(75delta^2)<11000/delta^2,
|E_compact|<=1260607/(450delta^2)<5000/delta^2.

The complete finite four-leg contribution is bounded by
[(6*3+14)*870+6*6]/delta=27876/delta<28000/delta.
The completion is at most23000/delta^2.
Each old Gram finite term is at most
12*5*256/4+(5*5+2)*256/2=7296, so their sum is below22000.
Finally|2*T1*B_soft|<1200000/delta.

## One common bound

Every term is now bounded by its listed constant times
(1+L)/delta^2, using0<delta<=1. This also dominates
the triangle factor(1+L)/sqrt(delta), all1/delta terms,
and10+L<=10(1+L),12+L<=12(1+L).
The exact budget is

6*(110000*200*10+1300*200)
+3*(210000000*10+200000*1200+2500*100
    +2300000*12+30000000*21)
+5000+28000+23000+22000+1200000
=10316388000 <10^11.

Since16pi^2>144 and10^11<144*10^9, this proves the advertised
known amplitude bound. The three matching terms are excluded,
not presumed smaller than the known part.

## Unchanged original parameter comparison

n=10^200/512+2<10^200/256 andg=1/8192 imply
A_m>=4g^2/n^3>10^-600 in the chosen physical domain.
S297 establishesA_G>0, soA_m+A_G>=A_m.
Choose the explicit comparison windowdelta=10^-204.
Then1+L=1+204ln10<613<1000 andkappa=10^800 gives

|A_known|<10^(9+3+408-1600)=10^-1180,
|A_known/(A_m+A_G)|<10^-580,
|2Re[A_known/(A_m+A_G)]|<2*10^-580<10^-579.

There is no loop square in this one-loop interference statement.
The estimate deteriorates asdelta approaches zero. It is neither
a forward limit nor a high-energy/Regge remainder estimate.
