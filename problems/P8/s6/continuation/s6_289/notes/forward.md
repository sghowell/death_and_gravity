# Exact light coefficient and uniform fixed-heavy forward bounds

All coefficients here are formed after consistent subtraction
of known Newton/light poles in all three channels. Set
t=0,s=2mu+v,u=2mu-v. Each nonlocal curvature function has
zero value at t=0, so that channel contributes no v^2 term.
The two crossed gapped terms yield one second derivative.
The circle|a|<=3mu is below4n for n>=mu>0, so their defining
integrals have a positive uniform denominator margin near
the crossing center. Dominated analytic differentiation is
valid before taking these derivatives. No massless physical
IR limit is interchanged.

Let h=1-x^2 and d=4n-2mu h. The exact derivatives are

J2=d_a^2[a(a-4mu)^2/(4n-ah)] at a=2mu
  =-16mu(4n^2-2nmu h-mu^2h^2)/d^3,
J0=d_a^2[a(a+2mu)^2/(4n-ah)] at a=2mu
  =16mu(20n^2-4nmu h+mu^2h^2)/d^3.

The full nonlocal coefficient for one species is

b20_nonlocal=integral_0^1[w2 J2/96+w0 J0/2304]dx
             /(pi^2 kappa^2).

The factors follow from the complete crossed amplitude, including
its trace sector and both s/u contributions. The t remainder is
zero, not a neglected finite curvature constant.

## A full finite-mass bound, not an asymptotic estimate

For0<=h<=1,n>=mu, d>=2n. The uncombined derivative identities are

J2=-4mu/d-8h mu^2/d^2+16h^2 mu^3/d^3,
J0=20mu/d+64h mu^2/d^2+64h^2 mu^3/d^3.

Their exact equality to the rational jets is checked. Taking
absolute values and mu/n<=1 gives|J2|<=6mu/n and
0<J0<=34mu/n throughout the full interval. The complete
weights are nonnegative. For the vector this is explicit as

w2=x^2[13+14(1-x^2)+3(1-x^2)^2]/30,
w0=x^2[8/3+3(x^2-1/3)^2].

Integrating the whole weights gives

|b20_scalar,nonlocal|<=73mu/(2520n pi^2 kappa^2),
|b20_vector,nonlocal|<=mu/(35n pi^2 kappa^2).

For the actual mu1,nH=10^200/512+2,M_A^2=10^6, their sum is
strictly below1/(10^7 pi^2 kappa^2). The arithmetic margin is
exact, with the entire hierarchy retained as integers/rationals.

The fixed local H/Proca coefficient is-L/(320pi^2 kappa^2),
L=log(nH)+2. Hence the COMPLETE selected H/Proca Gaussian
coefficient lies in the closed interval centered there with
radius[73/(2520nH)+1/(35M_A^2)]/(pi^2 kappa^2).
Since L>2, its upper endpoint is strictly negative.

For an absolute bound nH<10^200. The positive exponential
series at5/2 through degree4 exceeds10 by329/384, so
log10<5/2 and L<502. Also pi=4 integral_0^1dx/(1+x^2)>2.
Thus |b20_HP|<(502/320+10^-7)/(4kappa^2)<1/kappa^2,
with final rational margin24312499/40000000. At actual
kappa10^800 this is10^-1600. This is a proved finite-mass
bound on this particular sector, not a full-source sign
test or a permissible isolated finite-gravity positivity bound.

## Exact light nonlocal integral

For a scalar n=mu the dimensionless integrand is

x^2[3x^8-12x^6+22x^4-140x^2+255]/[1920(1+x^2)^3].

The numerator bracket equals
128+120(1-x^2)+4(1-x^2)^2+3(1-x^2)^4,
so the integral is strictly positive. Exact rational integration,
checked again by differentiating its primitive, gives

b20_Phi,nonlocal=(319/4800-23pi/1280)/(pi^2 kappa^2).

The complete light Gaussian contribution still includes
4(r_Phi+4w_Phi/3)/kappa^2, with both local terms unmatched.
Its known-pole subtraction also retains the independent
light Newton coefficient. This does not fix S288's two
full local anchors, and the metric loop is not added twice.

For calibration only, the heavy expansion of
pi^2 kappa^2 b20_nonlocal is
scalar:mu/(240n)+151mu^2/(60480n^2)+O(n^-3);
vector:0/n+mu^2/(20160n^2)+O(n^-3).
The vector leading cancellation is exact. These displayed
series are not used as remainder bounds; the whole-interval
inequalities above supply the actual finite-mass guarantee.
