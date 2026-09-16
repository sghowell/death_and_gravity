# Uniform complex neighborhood and mixed derivative

All norms in complex continuation are absolute Euclidean/Frobenius
norms; dot products defining the action and momenta stay bilinear.
Only a,b are continued. Directions, TT tensors, E,u, elastic Born
A0, positive transfers tau_t,tau_u and delta=min(1,tau_t,tau_u)
are fixed real-center parameters. No absolute-value operation is
inserted into the analytic amplitude.

## Recoil neighborhood

For every real center a,b>0,W=a+b<=1/8, take the closed polydisc
 |a'-a|,|b'-b|<=cW, c=10^-12.
At its center h=E-W/2, Eprime^2=h^2-v.v/4>=45/32 and
rprime^2>=13/32. Continue both square roots from these positive values.
The perturbations obey |dh|<=cW,||dv||<=2cW and
|d(Eprime^2)|<=4cW+cW^2+2c^2W^2<5cW.
Factoring each square-root difference, or first applying its
nonzero-disk bound and then factoring, gives
|dEprime|<5cW,|drprime|<10cW, |Eprime|>1,|rprime|>1/2,
|h+Eprime|>2, |rprime/Eprime|<1.

Put beta=+/-rprime*(u.v)/[4Eprime(h+Eprime)]-1/2.
For t=rprime/Eprime, |dt|<20cW. The factor
f=t/[4(h+Eprime)] has |f|<1/8 and |df|<3cW.
Consequently |dbeta|<(1/4+3/8)cW<cW and |beta|<1.
In p0=h-/+t*(u.v)/2 and pvec=+/-rprime*u+beta*v this gives
 |dp0|<5cW, ||dpvec||<13cW.
The real recoil-to-Born spatial displacement is<=3W (S300);
its complex extension is<4W. All shifted scalar vertex momenta
have components<3, and heavy ones<5.

For Lipschitz bounds use the convex spatial path between the real
future massive momenta, perturbed by at most13cW. Its real norm is
<=sqrt(3). Throughout this convex tube |pvec|<7/4 and |p0|>7/8.
Indeed the squared-energy perturbation is at most
2sqrt(3)*epsilon+epsilon^2; continuation from p0>=1 bounds the
square-root difference by4epsilon. Thus the energy gradient is<2.
The Doppler gap at the real points is>1/4, and its perturbation
is below65cW, so |p.n|>1/8 on the tube and its gradient is<3.
For S(p)=p.eps.p/(p.n), unit TT gives a numerator bound4
and a gradient bound4; hence
 ||grad S||<4*8+4*3*64=800<1024.
Incoming signed momenta are handled by S(-p)=-S(p) on the same
future-massive tube. This is a complex Lipschitz estimate, not
an unexamined import of a real mean-value inequality.

## Hard, light and collinear denominators

Every hard momentum changes componentwise by<100cW. At a mixed
cut its real maximum component is<=sqrt(tau)+4W and S312 gives
-D>(tau+W^2)/300000, for every ray assignment. Consequently
 |dD|<=800cW(sqrt(tau)+4W)+40000c^2W^2
      <3601c(tau+W^2)<(tau+W^2)/600000.
The complex denominator therefore has magnitude>(tau+W^2)/600000.
Timelike gaps retain45/16 and heavy kinetic gaps retain n/2.

Individual light propagators are2a'p.n1 or2b'p.n2. Their a' or b'
zero is removable after scaling by a'b': the topology inventory
has no repeated single-ray light denominator on a double-external
path. A total-Q light denominator has real magnitude>=3W/8.
Its perturbation is<100cW, so its magnitude is>W/4.
Thus |p.Q'+q1'.q2'|>W/8. Also |a'b'|<W^2,|a'|,|b'|<2W
and W/2<|a'+b'|<2W.

The47 pair graphs must be combined before taking an angular bound.
S311's conserved-current algebra cancels the apparent relative-angle
pole exactly. Its normalized coefficients are polynomials in
z=a'/(a'+b'), with max degree4 and |z|<2. Every real-angle monomial
is r^j/(1+r^2)^m with j<=2m. Absolute polynomial coefficient sums
weighted by2^degree are2642,432,352,1764, each<10^4.
This algebraic identity extends to the complex energy neighborhood.
It is not an absolute bound on uncombined singular pair diagrams.

All remaining denominators and square-root radicands are nonzero
in a neighborhood of the closed polydisc. The scaled full amplitude
G=a'b'*sqrt(rho2)*M6/A0 is therefore holomorphic there, including
removable single-ray axes.

## Regular graphs with at most one light propagator

There are247 such graphs. Multiplication by a'b' makes the light
propagator factor at most8W: for a single-ray propagator this is
4|other energy|; for a total-Q one it is4|a'b'|/W; with no light
propagator it is |a'b'|. Use the canonical S312 vertex bounds1024^r
and the heavy chain bound8/n. On a mixed hard path
L<=sqrt(tau)+5W, L^2<=50(tau+W^2), L^2/|D|<30000000.
The first trace-reversed hard propagator is<1800000/(delta+W^2).
Each subsequent cubic or quartic step costs, respectively,
 R3=3*138240*30000000=12441600000000,
 R4=3*10616832*30000000=955514880000000.
Overcounting by the full regular inventories is safe. The Einstein
coefficient is
 177*1024^4*8*1800000*max(1,R3,R3^2,R4)
 =433798508189475293550870528000000000000000000000<10^49,
multiplying W/[kappa^2*(delta+W^2)].
The matter coefficient is (33*4+177*8)*1024^2*8 multiplying
g^2 W/(n*kappa). These estimates use absolute component sums.

## Double-external groups and the elastic limit

In the ordering identity |Na|,|Nb|<4, |Xa|<5|a'|,|Yb|<5|b'|,
|A/a'|,|B/b'|>1/8, |Z|<=2|a'b'| and |A+B+Z|>W/8.
After a'b' scaling the three correction terms are bounded by
(5120+5120+65536)W=75776W<10^5 W.
The scaled leading individual product is<=1024.

For components<3, each scalar stress has entries<55; its trace
reverse has entries<165, so its two-stress contraction |N|<2*10^5.
Assigned radiation shifts each scalar component by<2W. The exact
bilinear identity bounds the stress change by
72W+24W^2<=75W<100W. Recoil-to-Born component shifts<4W give
144W+96W^2<=156W<200W. Therefore assigned-core |dN|<528000W<10^6 W,
and Born-core |dN|<1056000W<10^7 W.

In a mixed channel use scaled paired currents S_L,S_R. The complex
gradient proved above gives
 |S_L^Born|<=1024sqrt(tau),
 |S_L-S_L^Born|<=4096W,
 |S_L|<10^4(sqrt(tau)+W),
and the same on R. A product therefore obeys
 |Sprod|<2*10^8(tau+W^2),
 |Sprod-Sprod0|<45154304W(sqrt(tau)+W)<10^8W(sqrt(tau)+W).
The hard denominator difference from D0=-tau satisfies
|D-D0|<200W(sqrt(tau)+W). Telescope N*Sprod/D-N0*Sprod0/D0
in its numerator, soft product and inverse denominator.
Use delta<=min(1,tau),
 delta*(sqrt(tau)+W)/(tau+W^2)<=2,
and AG>8/(kappa*delta). The coefficient of W/kappa for one
mixed distribution is at most
 10^7*2*10^8*600000/8
 +2*10^5*10^8*600000*2/8
 +2*10^5*1024^2*200*600000*2/8.
For a timelike distribution use |S_L|,|S_R|<=128,
changes<=8192W and denominators>=45/16, with Born>=25/4.
All8 mixed plus4 timelike distributions cost<10^25 W/kappa.
The60 ordering/core-shift corrections cost at most
 60*(10^5*2*10^5+1024*10^6)*600000/8 *W/kappa
 <10^20 W/kappa.

For matter, |C|<4g^2/n and each heavy propagator<2/n.
|D-D0|<200W gives a heavy inverse difference<800W/n^2<7W/n.
The full four-leg scaled current is<=256 with change<=16384W.
Contact leading variation costs4*(2*256*16384);
heavy leading variation costs12*(2*2*128*8192+7*128^2);
ordering errors cost80*10^5*4. Including the non-double matter
coefficient, their total is<3*10^10 multiplying g^2W/(n*kappa).
Using Am>4g^2/n^3 gives<10^10 n^2 W/kappa.

At zero radiation the grouped leading sum is the complete elastic
(Am+AG)*S_a(Born)*S_b(Born)/kappa, before division by A0.
No approximation to the tuned contact or replacement of the full
positive Born normalization is made.

## Complete pair class and phase

For the off-shell47 hard current, the complex component budgets are
 4*1024*4*4+12*1024*4*2+(16+12288+48)/8
 =165384<10^6
multiplying g^2/(n*sqrt(kappa)*W), and
 12*1024^2*4*10^8+6*1024^2*10^8
 +2*138240*50*10^16+138240*25*64^2*(65/64)
 =138240005662324776960000<10^26
multiplying1/[kappa^(3/2)*W*(delta+W^2)].
These count all external, seagull, heavy internal/density and
Einstein cubic current terms. A spatial Frobenius norm is at most
three times its component bound. Combining the conserved pair
coefficient<10^4, |W'|^2<4W^2 and the positive Born bounds gives
 |a'b'*Mpair/A0|<(3*10^10 n^2+10^32)W/kappa.

The phase is analytic: rho=E*rprime/(r0*Eprime).
For real Eprime, d ln(sqrt(rho))/dEprime
=1/[2Eprime(Eprime^2-1)], and |d sqrt(rho)/dEprime|<2
on the compact interval. Since E-Eprime<W, the real difference
from1 is<2W. In the complex tube |d rho|<40cW and real rho>1/4,
so |d sqrt(rho)|<80cW. In particular |sqrt(rho)|<2 and
|sqrt(rho)-1|<16W. The latter multiplies
|G00|<=48^2/kappa=2304/kappa.

Summing all budgets, multiplying by the phase and adding its elastic
difference gives, everywhere in the closed polydisc,
 |G-G00|<B W, B=(10^40*n^2+10^60)/kappa<10^-350.
The exact arithmetic records check the separate n^2 and constant
coefficients against10^40 and10^60. They are not the proof of
holomorphy; that follows from the denominator and topology arguments.

Two-variable Cauchy on radius cW now gives
 |partial_a partial_b G|<B/(c^2 W)<10^-326/W.
This is the needed uniform derivative theorem. The bounded real
amplitude envelope of S312 alone would not imply it.
