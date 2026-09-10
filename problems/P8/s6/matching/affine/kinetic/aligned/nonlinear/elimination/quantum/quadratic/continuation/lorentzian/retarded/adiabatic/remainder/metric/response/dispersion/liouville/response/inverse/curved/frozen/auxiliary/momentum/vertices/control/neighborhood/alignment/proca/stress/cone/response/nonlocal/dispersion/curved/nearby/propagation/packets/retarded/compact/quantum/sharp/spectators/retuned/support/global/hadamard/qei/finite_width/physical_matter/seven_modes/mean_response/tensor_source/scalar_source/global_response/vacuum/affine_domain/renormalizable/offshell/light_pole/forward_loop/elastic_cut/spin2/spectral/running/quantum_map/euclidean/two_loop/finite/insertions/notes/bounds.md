# Full radial remainder and actual grouped error

Use both full-vertex bounds |V1 V2|<9lambda4^2 and the
two decaying-kernel insertions. The four-dimensional radial
measure is (16pi^2)^-1 y dy. The 1/y in the kernel majorant
cancels this radial factor. Summing all three channels,
including their factor 1/2, gives

    sup_disc |A_remainder|
      <= 54 lambda4^2 g/(16pi^2)^2
           integral_0^infinity ln(1+y/M)/(y+1/4)^2 dy.

The factor 54 is checked as 3 times (1/2) times 9 times
2 line insertions times 2 from the kernel strip estimate.
Both loop measures are retained. No loop-momentum cutoff
appears.

For positive delta and M different from delta, an anchored
primitive is

    -ln(1+y/M)/(y+delta)
      +[ln(y+delta)-ln(y+M)]/(M-delta).

Its derivative is the required integrand, it vanishes at
infinity, and its value at zero gives

    integral_0^infinity ln(1+y/M)/(y+delta)^2 dy
       = ln(M/delta)/(M-delta).

All three statements are checked exactly. With delta=1/4,
the result is ln(4M)/(M-1/4). Absolute domination and the
complex routing argument justify the integrated analytic
remainder on a neighborhood of the unit disc.
Cauchy's coefficient estimate at radius one bounds its
second forward coefficient by the same supremum.
It does not add another factor 1/2 to the whole-amplitude
bound.

For the actual parameters 4M<10^200, so ln(4M)<600 by
the same elementary exp(3)>10 argument of S6.121.
Using pi>3 gives

    E_remainder =
      54 lambda4^2 g times 600/[144^2 (M-1/4)].

The constant part retains the complete S6.113 renormalized
one-loop error bound E1 and the fixed S6.119 upper bound
A=g/(288M) for alpha:

    E_constant=2 A E1,
    E_insertions=E_constant+E_remainder.

Exact rational comparisons give

    E_insertions < 10^-614,
    E_insertions/(4lambda) < 2 times 10^-15.

The tiny constant term is nonzero and below 10^-813.
A diagnostic evaluation of the remainder upper bound is
about 6.2450 times 10^-615; no decimal is used as proof.

Combining this group with the 88 finite refinements of
S6.121 represents 152 bare refinements, together with
this group's specified counterterms. Their total error
upper bound is still below 3 times 10^-607 and below
10^-7 of the tree coefficient. The other 40 raw
refinements and complete two-loop normalization remain
outside that partial sum.
