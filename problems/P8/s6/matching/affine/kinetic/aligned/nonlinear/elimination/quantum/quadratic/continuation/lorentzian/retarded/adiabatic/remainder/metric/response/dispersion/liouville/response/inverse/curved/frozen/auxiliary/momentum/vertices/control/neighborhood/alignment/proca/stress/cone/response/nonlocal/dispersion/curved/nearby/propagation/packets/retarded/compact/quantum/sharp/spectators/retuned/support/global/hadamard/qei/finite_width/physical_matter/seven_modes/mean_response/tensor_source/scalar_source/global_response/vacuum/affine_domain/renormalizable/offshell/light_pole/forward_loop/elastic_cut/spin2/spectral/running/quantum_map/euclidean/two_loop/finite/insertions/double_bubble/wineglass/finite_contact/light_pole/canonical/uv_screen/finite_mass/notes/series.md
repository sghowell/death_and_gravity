# Uniform all-degree finite-mass box tail

Expand each of four propagators in its shift:

    S(q+P)=sum_{j>=0} (-S0(q) i gamma.P)^j S0(q).

At total insertion degree n there are binomial(n+3,3)
ordered weak compositions among the four propagators.
The actual zero-shift propagator makes this a safe
overcount. The closed form is valid at every n; enumerated
degrees zero through twelve are independent diagnostics,
not the proof of the infinite tail.

The spin trace is bounded by four times the product of
operator norms. Six cyclic labelled orders contribute.
Two active flavors have the same y^2 and combine with
tr(Ta Tb)=delta_ab/2 to give exactly delta_ab.
The Yukawa-inert fermions do not enter this transition.
For n>0 the needed radial integral is

    integral d^4q/(2pi)^4 (q^2+m^2)^(-(4+n)/2)
      = 4/[Q n(n+2) m^n].

For t=q^2/m^2, a primitive of t/(1+t)^(2+n/2)
normalized to vanish at zero is

    4/[n(n+2)] - 2(1+t)^(-n/2)/n
      + 2(1+t)^(-1-n/2)/(n+2).

Its derivative, zero anchor and infinity limit are checked.
Both decaying terms vanish for every n>0. This proves
the radial value without treating an odd-degree norm
majorant as an odd integrand cancellation.

Let x=12/m<=1/2. For one color-diagonal amplitude and
each transverse polarization pair, the remainder after
degree two obeys

    |R_box| <= (96 a Y/Q) sum_{n>=3}
                    binomial(n+3,3) x^n/[n(n+2)].

The constant is 4 spin times 6 orders times 4 radial.
For n>=3 the weight is
(n+3)(n+1)/(6n)<=(n+5)/6, with nonnegative difference
(n-3)/(6n). Summing this majorant gives

    x^3(8-7x)/[6(1-x)^2] <= (16/3)x^3 < 6x^3.

The ratio after dividing by x^3 is increasing on [0,1/2]:
its derivative is (9-7x)/[6(1-x)^3]>0.
Consequently, for positive a,Y and finite m,

    |R_box| < 995328 a Y/(Q m^3)
            < 10^6 a Y/(Q m^3).

This is uniform in all routes, cut angles and continued
on-shell configurations in the written domain.
It bounds the entire finite-mass momentum remainder of
the one-loop transition, not a field-amplitude series and
not all loop orders. Ward subtraction is completed before
using any four-dimensional norm bound on this tail.
