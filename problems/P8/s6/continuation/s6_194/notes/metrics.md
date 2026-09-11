# Absolute metric jets and a product-compatible tensor norm

Write the old-signature metric as g=diag(-1,a^2 exp(gamma)),
with a=(1+t^2)^2, gamma=epsilon Gamma. All inverse and exponential
products remain ordered. No simultaneous diagonalization is used.

For time order j<=4 and amplitude order b<=3, the ordered Frechet
formula, its time product rule and ||epsilon Gamma^(r)||<=.01
give the absolute exponential bound

    E(j,b)=2 sum(l=0..j) binomial(j,l) b^(j-l) Bell_l(.01).

Bell_l(.01) is the coefficient derivative of
exp(.01(exp(z)-1)) at zero. Equivalently E(j,b) is twice the
j-th derivative of exp(.01(exp(z)-1)+b z) at zero.
The factor2 bounds exp(.01) strictly. In particular E(0,0)=2,
not the relative zero-jet value1 used in a differently normalized
earlier estimate. Both exp(gamma) and exp(-gamma) obey this bound.

The polynomial a^2=(1+t^2)^4 has nonnegative coefficients.
Absolute derivative coefficients evaluated at |t|=.5 give its
exact displayed time majorants. For a^-2 use the already proved
S190 scale_bound(2,j). Product differentiation gives raw metric
operator bounds M(j,b) and inverse bounds M_inv(j,b).

For a tensor T define N(q,A;T) as the sum over 0<=j<=q,0<=b<=A
of sup over the admitted (t,epsilon) domain of its component-l1
time/amplitude derivative, divided by j!b!. This is a seminorm
on jets, with ordinary norm at order0. The factorial weights
make tensor products submultiplicative. Permutations preserve
the component sum; contracting any indices cannot increase it
beyond the uncontracted product. No hidden dimension factor is
omitted by a contraction.

For a three-dimensional matrix, entry-l1<=3 sqrt(3)||T||op
<6||T||op. Thus for r=0..4,

    N(4-r,3; dt^r g)
      <= 6 sum(j=0..4-r,b=0..3) M(j+r,b)/(j!b!)+1_(r=0),

and identically for the inverse. Exact rational evaluation of all
ten expressions is recorded in metrics.constants(). The largest
is below111807. Consequently G=200000 strictly bounds each one.
G=100000 would fail for the inverse third shifted jet and is
explicitly rejected by an independent negative control.

All later contractions can restrict these norms to fewer time or
amplitude derivatives without increasing them. This finite jet
algebra is not a bound on arbitrarily high derivatives.
