# Matching the actual all-order state, not resetting it

Gamma vanishes in a neighborhood of the initial Cauchy
surface. The actual canonical T/T/L data there are exactly
the source-pinned S6.55 all-order frequency and slope.
In the isotropic frame L=-diag(d_T,d_T,d_L), where
d_T=H/2, d_L=H(1/2+z), lambda=omega'/omega=-Hz and
z=1-m^2/omega^2. Direct differentiation verifies
U=d'+d^2 for both original mode equations. Their scalar
Riccati squeeze is d+lambda/2.

Retain all first five WKB coefficients P_j from the
unchanged recurrence. Let W10=omega w with
w=1+sum_(j=1..5)P_j x^(2j), x=omega^-1, and give W10
its actual differentiated slope. The original all-order
cutoffs are fully on for these terms whenever the
initial nu exceeds twice both fifth thresholds. The
larger threshold is 61205538081189019/32<2e15.
Our high-frequency analysis band nu_minus>=1e16
therefore lies beyond both full turn-ons.

No assertion is made that the remaining terms vanish.
For each active j>=6, the original threshold rule gives
both coefficient and slope bounds at most2^-j nu.
Consequently each frequency or slope contribution is
at most2^-j nu^(2-2j). Their locally finite sum is bounded
by the full positive geometric series,

    sum_(j>=6)2^-j nu^(2-2j)
      <=nu^-10/[64(1-1/(2nu^2))] <nu^-10/32.

Cutoff weights lie in [0,1], so the same bound covers
every partially active higher term. This proof does not
need to compute an infinite number of coefficients.

For real physical u,z fixed in the original coefficient
box, regard x as an independent complex variable.
The exact first-five coefficient and slope bounds give
on |x|<=1/1000:

    |w-1|<1/10000,
    |x(d+lambda/2+w'/(2w))|<1/100.

Here w' includes x'=-lambda x; it is not a derivative
at fixed x. The scalar graph of the W10 data is

    rW=(omega-W10+i(W10'/(2W10)+d))
        /(omega+W10-i(W10'/(2W10)+d)).

It is analytic on that closed x disc and has modulus
below1/10. Holding the physical z fixed when applying
Cauchy's estimate is legitimate: the resulting bound
is uniform over all real z in [0,1], including the actual
z=1-m^2/nu^2.

The exact algebraic WKB-to-Riccati residual identity is

    F(rW)=-i R_W(1+rW)^2/(2omega),
    R_W=omega^2-(d'+d^2)-W^2-W''/(2W)+3W'^2/(4W^2).

The unchanged WKB recurrence through P5 has R_W=O(x^10).
Equating formal powers in this identity uniquely matches
the first ten coefficients of rW to the matrix recurrence
at the isotropic initial surface. This is a finite formal
identity, not convergence of the all-order WKB series.
Independent tenth-order T/L coefficient comparisons
supplement it.

Cauchy's estimate gives a graph Taylor remainder below
1000^11 nu^-11/5 for nu>=2000. To compare the actual
state to W10, both frequencies are at least nu/2 and
|W10'|<=3nu. The frequency/slope tail above then implies
a half-log-rate difference below nu^-11. Using
r=2nu/(nu+W-i v)-1 bounds the graph difference by
(nu^-10/32+nu^-11)/nu<nu^-11.
Together these prove

    ||r_actual(t0)-rhat(t0)|| <1e33 nu^-11.

Both transverse modes and the constrained longitudinal
mode are included. Rotating their initial projectors
does not change the operator norm. The initial positive
frequency phase has no effect on this graph.
