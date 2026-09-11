# Uniform mixed frame derivatives

Consider gamma(t,epsilon)=epsilon Gamma(t), with smooth
compact tracefree symmetric Gamma, zero near the same
initial Cauchy surface, and ||Gamma^(j)||op<=1 for
j=0..12. Work at |epsilon|<=delta=1/100, on the unchanged
unit CD slab. Every time jet of gamma is then bounded
by delta. No state, parent, mass, profile or prescription
is changed.

For a mixed derivative with j time derivatives and
a amplitude derivatives, a<=2, the ordered matrix
exponential formula and the set-partition product bound
are majorized by the scalar exponential

    exp[delta(exp(t)-1)+sigma exp(t)].

Its derivative at t=sigma=0 is
B_(j,a)=sum_(l=0..j)binom(j,l)a^(j-l)Bell_l(delta).
Indeed differentiating sigma a times gives exp(a t);
differentiating time then gives the displayed binomial
sum. A block containing no parameter derivative costs
delta, one parameter derivative costs1, and two costs0.
The zero-order gamma itself is bounded outside this
generating function by the exponential norm. The same
ordered-simplex bound works for noncommuting time
directions; it is not a scalarization of the matrices.

After balancing by exp(-gamma/2), use e_(0,0)=1 and
e_(j,a)=2B_(j,a) otherwise. For a=0, this is exactly the
S6.190 bound. For derivatives of K, the constraint
kk^T/(a^3 m^2) contributes only when a=0. The scale
derivatives therefore combine with these exponential
bounds exactly as recorded in mixed.py. The normalized
omega-squared derivatives have the analogous scale
power2 convolution.

All root, inverse-root, frequency and reciprocal
recurrences are the S6.190 equations with a two-component
multiindex beta=(j,a) and product-binomial coefficient
binom(j,n)binom(a,b). They are evaluated in increasing
total order. The positive-root Sylvester multiplier
b_l/(b_i+b_l) still lies in (0,1); no condition number
enters at any mixed order.

Leibniz gives the mixed derivatives of L=B^-1 B_t and
rho=omega_t/(2omega), and hence those of R and S.
The first parameter derivatives of omega have relative
bounds1 and2. All zero-parameter constants reproduce
S6.190 exactly. Required time jets do not exceed12;
parameter differentiation does not silently request
additional time derivatives of Gamma.
