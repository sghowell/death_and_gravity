# Complete mixed cutoff derivatives through order 196

Set m = 96 and R = 10^20. The Gaussian-frame operator estimate requires
every multiindex with alpha_j <= 2, hence total order at most 192.
The heat remainder adds four phase derivatives. We therefore prove bounds
through total order 196, not just the four derivatives used in S6.270.

## Smooth bump on the entire real line

For b(y) = exp(-1/y) when y > 0 and zero otherwise, put x = 1/y.
On y > 0 its nth derivative is exp(-x) P_n(x), with

    P_0 = 1,    P_(n+1) = x^2(P_n - P_n'),    degree P_n <= 2n.

For B_n = sum_k |p_n,k| k!, the contribution of a monomial x^k to the
next weighted norm is bounded by

    (k+2)! + k(k+1)! = 2(k+1)^2 k!.

Thus B_(n+1) <= 2(2n+1)^2 B_n <= 8(n+1)^2 B_n, and induction gives

    B_n <= 8^n (n!)^2.

For x >= 0, exp(-x)x^k <= k! follows from the exponential series.
These B_n therefore bound the actual derivatives uniformly. The flat
extension through y = 0 and the zero side are included.
The exact first nine weighted polynomial norms are evaluated independently,
and the first five are matched to the unchanged parent values.

## Full transition quotient

Let h(y) = b(y)+b(1-y). On the transition interval [0,1], at least one
argument is >= 1/2, so h >= exp(-2) > 1/9. The transition
theta = b(1-y)/h satisfies 0 <= theta <= 1 and is extended by its exact
constant values outside [0,1].

Differentiating h theta = b(1-y) gives the complete recurrence, n >= 1,

    C_n <= 9 [ B_n + 2 sum_(k=1)^n binom(n,k) B_k C_(n-k) ].

All product-rule terms and both derivatives of the denominator are included.
With C_0 <= 1 and candidate C_n <= 256^n(n!)^2, divide by
256^n(n!)^2. The factorial ratio in the sum is 1/binom(n,k) <= 1.
Writing r = 8/256 = 1/32, the resulting upper bound is no larger than

    9 r^n + 18 sum_(k=1)^n r^k
       <= 9r + 18r/(1-r) = 855/992 < 1.

This proves C_n <= 256^n(n!)^2 for every n. It concerns phase cutoffs,
not any unproved high real-time smoothness of the cosmological profiles.

## All mixed radial partitions

For y_c(w) = (|w|^2-R^2)/(c R^2), c = 1 or 2, the only nonzero inner
derivatives have orders one and two. Wherever theta derivatives matter,
|w| <= 2R, so their multilinear norms are <= 4/R and <= 2/R^2.

In an nth derivative, let j pairs of arguments enter second derivatives of
y_c. The number of full multilinear partitions is

    n! / [(n-2j)! j! 2^j],   0 <= j <= floor(n/2).

After multiplication by the inner derivative bounds, the entire nth
derivative of chi_c is bounded by

    R^-n sum_j
      n! 4^(n-2j) C_(n-j) / [(n-2j)! j!].

This includes repeated coordinates and every mixed coordinate contact.
With the preceding transition estimate, each summand is at most
1024^n(n!)^3. There are at most 2^n terms, also for n = 0. Hence

    ||partial^alpha chi_c||_infinity <= 2048^n(n!)^3/R^n.

The code additionally evaluates all 197 actual finite sums using the
transition majorant and checks the declared bound for each order 0..196.
It does not replace this finite verification with a floating asymptotic.

## Full holomorphic amplitude and Leibniz contacts

Use the full parent complex ball 4R and the simultaneous R/8 polydisc proved
in notes/source.md. For either centered amplitude h = g-g0 or F-F0,
its bound there is 2A, with A = H = 10^1000 or E = 10^-255. Cauchy gives

    ||partial^alpha h|| <= 2A n! 8^n/R^n.

Apply the complete Leibniz rule to chi_c h. For each total order n,

    ||partial^alpha(chi_c h)|| <=
      2A R^-n sum_(k=0)^n
      binom(n,k) 2048^k (k!)^3 (n-k)! 8^(n-k)
      <= 2A 2056^n (n!)^3/R^n = J_n(A).

The inequality uses (k!)^3 (n-k)! <= (n!)^3 and the binomial theorem.
All derivatives of the full implicit g and F, radial cutoffs, spatial
reconstruction and original cross covariance are covered. No separated
vertex or diagonal covariance is substituted for these full amplitudes.
The retained scalar center is annihilated by positive-order derivatives.

For 0 <= n < 196 the exact successive ratio is

    J_(n+1)/J_n = 2056 (n+1)^3/R < 10^-9.

Every one of the 196 rational inequalities is checked. Consequently all
J_n from order 4 through 196 are bounded by J_4. There is no assertion
that the factorial bound decreases at arbitrarily high orders.
