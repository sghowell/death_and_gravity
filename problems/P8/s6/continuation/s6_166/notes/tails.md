# All derivative orders of the actual switching profile

For s(x)=x/(1+x^8)^(1/8) on the positive real root, write

    d_x^n s=P_n(x)/(1+x^8)^(n+1/8).

Direct differentiation proves for every nonnegative integer n

    P_0=x,
    P_(n+1)=(1+x^8)P_n'-(8n+1)x^7 P_n.

For n>=1, induction gives degree P_n=7(n-1), with leading
coefficient (-1)^(n-1)(9)_(n-1), where the last factor is
a rising factorial. The induction step multiplies the
highest coefficient by 7(n-1)-(8n+1)=-(n+8).
There is no hidden loss of an all-order derivative condition
from checking only a finite number of these polynomials.

On real |x|<=2 the denominator is at least one. Thus each
derivative has a finite explicit compact bound
sum_j |[x^j]P_n| 2^j. Add one at n=0 when subtracting
an asymptotic value. The fixed positive time scale rescales
these compact bounds by finite powers of tau^-1.

For the quantitative tail, take real x>=2 and the complex
disk |z-x|<=x/4. There |z|>=3x/4 and

    |z^-8|<=(2/3)^8<1/2.

Use h(w)=(1+w)^(-1/8), analytic on |w|<1/2, and the
branch s(z)=h(z^-8) matching the positive real profile.
It is not an arbitrary eighth-root choice. On this disk

    |h'(w)| <= 2^(9/8)/8 < 1/2,
    |s(z)-1| <= |z^-8|/2
              <= (4/3)^8 x^-8/2 < 5x^-8.

Cauchy's estimate on a disk of radius x/4 gives, for EVERY n,

    |d_x^n(s(x)-1)| <= 5 n!4^n x^(-8-n).

Oddness gives the negative tail with the same absolute cap.
For M(t)=m +/- Delta s(t/tau), this becomes

    |d_t^n(M(t)-M_asym)|
        <= 5 Delta tau^8 n!4^n |t|^(-8-n),
           |t|>=2tau.

Compact control and these tails establish symbol order -8
for all time derivatives. Spatial derivatives vanish.
The constants are for the fixed mass profile; no uniform
limit as tau tends to zero is claimed.

Independent tests compare recurrence derivatives with direct
high-precision differentiation through orders beyond the
report's finite table, including both tails. Complex-circle
diagnostics check the branch and Cauchy estimate. The written
induction and analytic disk argument, not sampling, cover
all derivative orders and all real times.
