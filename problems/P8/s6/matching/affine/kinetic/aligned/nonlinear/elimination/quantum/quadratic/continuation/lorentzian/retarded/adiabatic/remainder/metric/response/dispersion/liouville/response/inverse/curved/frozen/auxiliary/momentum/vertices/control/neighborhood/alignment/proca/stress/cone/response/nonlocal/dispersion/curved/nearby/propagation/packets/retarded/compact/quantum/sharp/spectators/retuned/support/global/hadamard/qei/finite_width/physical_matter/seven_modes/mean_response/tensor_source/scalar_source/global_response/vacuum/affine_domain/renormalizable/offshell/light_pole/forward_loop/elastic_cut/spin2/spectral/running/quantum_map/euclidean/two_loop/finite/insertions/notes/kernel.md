# The same on-shell kernel, not separately bounded mass terms

Use z=1-s_line and the frozen S6.112 parameters

    b(x)=x(1-x)/[xM+(1-x)^2],  pref=g/(16pi^2).

On 0<=x<=1, 0<=b<1/M except for harmless zero endpoints.
The actual relative inverse correction is

    Q(z)=Pi_R(1-z)/z
        =pref integral [b-ln(1+bz)/z] dx.

The pole convention is D(s)=s-1+Pi_R(s). The corresponding
Euclidean covariance expansion has relative insertion Q.
The same fixed asymptotic multiplier is

    alpha=pref integral b dx,

the finite light-residue counterterm magnitude already
computed independently in S6.112, S6.115 and S6.119.
Write Q=alpha+R, with

    R(z)=-pref integral b integral_0^1 dt/(1+btz) dx.

The anchored logarithmic primitive proves this split on
the branch continued from z=0. At the light shell Q(0)=0
and R(0)=-alpha. The exact universal limit is calculated
before substituting the complicated actual parameter weight.

The auxiliary representation also controls complex shifts.
For u with Re(u)>=-1/2, |u|<=|1+u|; this follows directly
from |1+u|^2-|u|^2=1+2Re(u). Thus |Q|<=alpha in the
corresponding strip. The actual routed arguments have a
much stronger denominator margin and never approach the
logarithm's cut.

After the outer Feynman shift, the argument is
z=(k+w)^2+1 with k real and Hermitian |w|^2<=3. Writing
w=a+ib gives

    Re z >= k^2/2-2.

Indeed the difference from k^2/2+1-|a|^2-|b|^2 is
|k+2a|^2/2. For y=k^2, M>32 and 0<=t<=1,

    Re(1+btz) >= 1-2bt+bty/2 >= (1+bty)/2.

Consequently

    |R(z)| <= 2 pref integral b integral dt/(1+bty)
             <= 2 pref ln(1+y/M)/y.

The last expression has continuous value 2 pref/M at y=0.
It decays sufficiently at infinity for the outer integral.
There is also the uniform bound |R|<=2alpha, which gives
the separately returned rational pointwise upper enclosure.
Neither bound assigns a finite value to the unsubtracted
constant-multiplier bubble.
