# Complete finite masters and a uniform physical Feynman contour

Write A_a(x)=mu-a*x(1-x)-i0 and define
J_a=int dx/A_a,
Q_a=int ln(A_a/mu)dx/A_a,
L_a=int ln(A_a/mu)dx.
Let ell=ln(4pi nu^2/mu)-Gamma_E and
ell_a=ell-ln[(-a-i0)/mu]. C_a is the finite two-massless triangle.

The literal S283 raw convention gives

D0(a,b)=J_b/a[-1/epsilon+ell_a]+O(epsilon),
C0mumu(a)=-J_a/(2epsilon)+(ell J_a-Q_a)/2+O(epsilon),
B00(a)=-1/epsilon+2+ell_a+O(epsilon),
Bmm(a)=-1/epsilon+ell-L_a+O(epsilon).

The S283 uniform box derivation is essential: no additional finite
box term is discarded. The minus signs useD=4+2epsilon, not the opposite
dimensional convention. Every logarithm keeps the Feynman boundary.

## Massive physical contour, mu=1

For s in[25/4,16] deform the x path to

x(r)=r+i*r(1-r)(1-2r)/4,0<=r<=1.

It passes above the left Feynman root and below the right root. The
homotopy r+i*lambda*r(1-r)(1-2r),0<=lambda<=1/4, respects those sides.
Away from the limiting real roots there is no crossed singularity.
Put v=r(1-r),y=v(1-2r)/4. Then

Re A=1-s(v+y^2), Im A=-s*v*(1-4v)/4<=0.

At r=1/2 the negative-real value is approached on the lower branch,
so its logarithm is ln|A|-i*pi, not ln|A|+i*pi.
This specifies the branch continuously over the path.

If |Re A|<1/4, the inequalities
3/4<s*v[1+v(1-4v)/16]<5/4 and
1<=1+v(1-4v)/16<=65/64 imply
3/65<v<1/5. Hence1-4v>1/5 and

|Im A|> (25/4)(3/65)(1/5)/4=3/208>1/100.

Otherwise|Re A|>=1/4. Thus|A|>=1/100 on the whole contour.
Also|A|<=1+16(1/4+1/256+1/16)=97/16<8.
Since1-6r+6r^2 lies in[-1/2,1],|x'|<2.
The lower-branch log has|ln A|<ln100+pi<9.
Therefore |J|<=200,|Q|<=1800,|L|<=18.

For negative channels -12<=a<0 the real parameter denominator is in[1,4],
giving smaller bounds directly. These estimates apply to the complete
boundary values, including their imaginary parts, without taking an
absolute value of a singular real-axis integral.

## Two-massless triangle

The radial change b=r/(1+r) gives

C_a=-int_0^1 dx int_0^infinity dr/[(1+r)(mu*r^2-a*x(1-x)-i0)].

At mu=1, for positive A the exact primitive is
[ln(1+r)-ln(r^2+A)/2+atan(r/sqrt(A))/sqrt(A)]/(1+A).
The endpoint values give, for a=-tau<0 and z=tau*x(1-x),

C_-tau=-int[ln z+pi/sqrt(z)]/[2(1+z)]dx.

Continue A to-z-i0 for positive a=s and z=s*x(1-x):

C_s=-int[ln z+i*pi(1/sqrt(z)-1)]/[2(1-z)]dx.

Atz=1 the real quotient tends to-1, and the imaginary quotient is
1/[sqrt(z)(1+sqrt(z))]. Thus this apparent interior pole is removable.

For positive s>=25/4, |ln z/(1-z)|<=1+(-ln z)_+
<=1-ln x-ln(1-x). Its real integrated contribution is at most3/2.
The imaginary magnitude is at mostpi^2/(2sqrt(s))<2.
Consequently |C_s|<5. The positive-channel threshold value
C_4=ln2-i*pi/2 is retained and independently calibrated.

For delta<=tau<=12 put L=|ln(delta)|.
Since|ln tau|<=L+ln12<L+3 and int(-ln x-ln(1-x))dx=2,

|C_-tau|<5/sqrt(delta)+(L+5)/2
         <10(1+L)/sqrt(delta).

These are analytic integral inequalities, not quadrature error bounds.
High-precision tests only calibrate the formulas and their physical sheet.
