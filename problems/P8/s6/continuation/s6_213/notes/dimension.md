# Uniform complex-dimension limit of the complete subtracted comparison

Only the singular unit-W8 comparison is continued. The convergent actual-state correction and finite time remainder remain in physical dimension three, as in the original S193 prescription. No physical state, Borel momentum analyticity or noninteger-dimensional Hilbert space is introduced.

Write d=3+eta with |eta|<=1/4. The canonical general-d dilution rates and the U polynomial are exactly those of S212. The WKB recurrence in dimension.py retains U(d) through all four orders, rather than substituting d3 in its higher coefficients. For each rational coefficient, its denominator is exactly c(1+t^2)^b with c>0 independent of d. Expanding its numerator in t,z,eta and applying

    |t|<=.51, |z|<=1.01, |eta|<=.25, |1+t^2|>.7

gives a rational coefficient majorant. All eight exact denominator reconstructions and physical reductions are checked, as are the first two full dimensional coefficients against the original implementation. Summing the four corrections at omega>=.99m gives relative W8 defects below .001 on the entire disk; their actual displayed maxima are about1.62e-5 and1.51e-5. The far scaled defects are below4.03e-9. These are explicit bounds, not inferred from sampling or a continuity assertion.

The underlying frequency omega is independent of d. S186's original complex-clock domain and S208's joint complex inverse-radius domain therefore remain valid without alteration. For every d in the closed disk, the same positive analytic W branches obey Re W>.5nu on real momenta, and Re What>.5 on the normalized far domain. Their sums have the corresponding strict positive real part. The canonical dilution rates are polynomial in d and bounded on this disk. Clock derivatives follow by Cauchy on the unchanged inner disc. Both opposite-phase branches use the same analytic d parameter; conjugating d in a physical-adjoint formula would destroy holomorphy.

For full general-d tensor contractions, integrate the finite transverse polynomial moments first, using S212's active/passive formula and fixed-source invariant reconstruction. The result is a finite linear combination of rational functions of d times functions of the longitudinal angle u. All denominators d,d-1 and the finite rising factorials remain nonzero. The remaining normalized measure is

    Gamma(d/2)/[sqrt(pi)Gamma((d-1)/2)]
      (1-u^2)^((d-3)/2) du.

Its absolute value is bounded by a finite constant times (1-u^2)^(-1/8). It is integrable at both endpoints. This is an analytic continuation of already evaluated tensor contractions; it does not use a positive measure in a negative-dimensional transverse subspace. At zero individual momentum the full polarization sum has its continuous massive-covariance extension. A direction may be chosen arbitrarily at that measure-zero point; no division by its radius is needed to bound the full sum.

Set U=m+|P|, rho=.01/U and L=2/rho=200U. On the compact scaled domain |x|<=rho, normalized momenta are n and -n+xP. The independent variables xU, P/U, m/U and real unit directions range over a closed bounded set including their high-transfer limits. The complete continued frequencies, amplitudes and all four inverse phases have the preceding strict denominator margins on this set. Consequently each of the five source-time endpoint rows is holomorphic in x and bounded uniformly in d. Cauchy after removing all degrees0,...,4-j gives the full endpoint-j remainder

    |E_j-U_j| <= C_j U^(5-j) r^-4, r>=L,

with every source derivative through orderj retained. Constants C_j are finite and independent of P,r,d. The dimensional radial factor changes r^2 only by r^Re(eta). Its worst far integral is therefore

    C_j U^(5-j) integral_L^infinity r^(-2+1/4)dr
      = (4/3) C_j U^(5-j)L^(-3/4).

Transfer growth is at most U^(4+1/4-j), not a fixed-P assertion. The physical numerical constants remain those of S208; no new numeric constant is deduced merely from compactness.
