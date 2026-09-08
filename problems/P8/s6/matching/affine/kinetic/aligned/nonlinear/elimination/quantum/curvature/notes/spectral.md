# Independent sphere spectral normalization and error control

This is a diagnostic unit four-sphere, not the cosmological metric.
For transverse one-forms the rough-Laplacian eigenvalue is
ell*(ell+3)-1 and multiplicity ell*(ell+3)*(2ell+3)/2, ell>=1.
Adding the Ricci endomorphism +3 gives the Hodge eigenvalue
(ell+1)*(ell+2). These are the d=4 spin-one entries of
[Kluth and Litim, Table IV](https://arxiv.org/abs/1910.00543).
Scalar gradients give the nonconstant scalar Hodge spectrum. Thus
Tr(exp(-s*D1))-Tr(exp(-s*D0)) equals the coexact one-form trace
minus the scalar constant mode, rather than the coexact trace alone.

Write x=ell+3/2, so the polynomial multiplicity is x^3-9x/4 and
the Hodge eigenvalue is x^2-1/4. Extending x to j+1/2, j>=0,
adds exactly -1 at x=1/2 and zero at x=3/2. The desired ratio trace is

    K_ratio(s)=exp(s/4) sum_{j>=0} f_s(j+1/2),
    f_s(x)=(x^3-9x/4)*exp(-s*x^2).

Both extension identities are checked exactly. In particular the
extra -1 is intentional and cannot be dropped as a putative vector
mode. The scalar in this determinant ratio includes its constant.

## Euler-Maclaurin derivation, not a formal divergent series

Apply the Euler-Maclaurin identity through the fifth endpoint
derivative to f_s(x), x from 1/2 to infinity. This follows by shifting
the argument in [DLMF 2.10.1](https://dlmf.nist.gov/2.10.E1) and
integrating its constant B6 term once into the endpoint derivative.
All Gaussian derivatives are absolutely integrable and vanish at
infinity for s>0. The remainder obeys

    |R_EM| <= 1/30240 * integral_{1/2}^infinity |f_s^(6)(x)| dx.

For completeness, B6(z)=1/42-y^2/2-y^3 with y=z(1-z) in [0,1/4]
on one period. Hence -31/1344<=B6(z)<=1/42 and |B6|<=1/42;
dividing by 6! gives the displayed constant. The polynomial identity
and both rational endpoint comparisons are checked independently.

The integral contribution is exactly
exp(-s/4)*[1/(2s^2)-1/s]. Evaluating the three endpoint derivatives
and multiplying by exp(s/4) yields the exact approximant

    1/(2s^2)-1/s-11/30
      -5s/63-197s^2/10080-13s^3/5040
      +11s^4/12096-s^5/30240.

The absolute coefficient sum of the displayed positive powers is
59/576; on 0<s<=1 their total magnitude is therefore <=59s/576.

For the remainder put y=sqrt(s)*x. The integral of the absolute
sixth derivative of x^3 exp(-s*x^2) is at most C3*s, and that
of x exp(-s*x^2) at most C1*s^2. Expanding the corresponding
sixth-derivative Gaussian polynomials and integrating each absolute
monomial over [0,infinity) gives C3=14016 and C1=2124. Each
odd-power Gaussian moment is Gamma((j+1)/2)/2, an exact rational
here. Extending the lower endpoint to zero only increases the bound.
Consequently |R_EM|<=(179/288)*s on this interval. Also
exp(s/4)<=4/3 by the exponential-series/geometric-series comparison.
Combining with the endpoint polynomial proves the continuous bound

    |K_ratio(s)-[1/(2s^2)-1/s-11/30]| <= (1609/1728)*s,
    0<s<=1.

## Coefficient comparison and limits of the control

The unit S4 has R=12, Ricci^2=36, Riemann^2=24 and boxR=0.
Its volume divided by (4*pi)^2 is 1/6. Substitution into the
independently contracted local Proca coefficients gives exactly
1/(2s^2)-1/s-11/30. If the scalar constant mode is lost, the
constant becomes 19/30 instead. That negative control is retained.

This bounded expansion verifies three local coefficients, including
the constant term, from the spectrum. It is not an all-orders
asymptotic assumption. The explicit bound belongs to this unit sphere
and this determinant ratio. The actual FLRW counterterms are computed
from its own metric in the local note. No finite sphere error, global
zero-mode normalization or Euclidean-state property is silently
transferred to the noncompact, time-dependent Lorentzian background.
