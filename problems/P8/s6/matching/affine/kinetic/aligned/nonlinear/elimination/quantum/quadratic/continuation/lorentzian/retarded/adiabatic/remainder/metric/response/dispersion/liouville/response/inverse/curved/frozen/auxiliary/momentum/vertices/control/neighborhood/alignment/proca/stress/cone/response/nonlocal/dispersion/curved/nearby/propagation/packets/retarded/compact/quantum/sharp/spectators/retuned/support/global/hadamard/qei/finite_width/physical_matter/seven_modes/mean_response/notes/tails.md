# A global weighted estimate and both time tails

The exact polynomial 800(1+u^2)^18[J-1/(5(1+u^2))]
has nonnegative coefficients and a strictly positive
constant. Thus J>=1/[5(1+u^2)] for every real u.
No asymptotic fit or finite scan supplies this pivot bound.

For u>=0 put t=1+u^2. The vector inequality r+s>=2r/3
and H>=0 give (r a^2)'<=0, hence
r(u)<=r(0)/t^4. Reflection of the ADDED covariance
extends this estimate to negative u. This does not
assume reflection symmetry of the reference vacuum.

Globally |alpha|<=12u/t, |beta|<=ell=1/(10t^6).
For Y=(xi,100t^3 dp), its second diagonal coefficient
contains -6u/t. Drop this negative damping in a Dini
norm bound; its remaining absolute diagonal correction
and the first diagonal correction are each at most
3u/(10t^12). The off-diagonal coefficients are bounded by

    241/(200t^3),  3/t^9+3/(40t^20).

Their total entry-sum integral is below two. This is
a bound on the growth of |Y|, not replacement of the
original system by a Hamiltonian of another sign.

For completeness the elementary integral constants are
I_n=integral_0^infinity (1+u^2)^(-n)du.
Integration by parts gives I_n=(2n-3)I_(n-1)/(2n-2),
I_1=pi/2. The positive integral
integral_0^1 u^4(1-u)^4/(1+u^2)du=22/7-pi
gives I_n<(11/7) binomial(2n-2,n-1)/4^(n-1).
The report retains exact rational upper bounds, including

    C_integral=888061455051153/529139970867200 <2.

The weighted source L1 norm, divided by r(0), is bounded by

    25/6+100 I_1+(25/4) I_12
      <7176850075/44040192 <170.

Here the first component uses integral u/t^4=1/6;
the second retains both direct pressure and lapse-mediated
forcing. Gronwall yields exp(C_integral)<exp(2)<9;
exp(1)<3 follows from the elementary factorial series.
Consequently

    |Y|_infinity <9*170*1500004 eta<2.4*10^9 eta.

This is uniform for u>=0. Uniqueness of the forced
system gives xi,n even and dp odd, extending the bound
to all real u. The scalar field response is odd.

Let C=2.4*10^9 and r_star=1500004. Reconstruction gives

    |n(u)|/eta <=(3C/10)|u|/t^3
                 +(3C/40)/t^11+(25r_star/4)/t^3.

Using |u|/t^3<=1/2 yields |n|<=5.5*10^8 eta.
Then |B|<=2.7*10^9 eta. Integrating
|delta psi'|<=(|n|/10+3|xi|/10)/t^6 gives
|delta psi|<4*10^8 eta. The induced scalar density
has bound 8.1*10^7 eta/t^12 and its ratio to the
background scalar density is at most 1.62*10^10 eta.

In particular n tends to zero. The integrable first
mean equation gives a finite limit xi_infinity.
For u>=1,

    |xi_infinity-xi(u)|/eta
       <=(3C/220)/t^11+(241C/1000)/u^5
          +(25r_star/6)/t^3.

The negative tail follows by reflection. These are
global first-order bounds, not a claim of uniform
control over all nonlinear or higher-loop terms.
