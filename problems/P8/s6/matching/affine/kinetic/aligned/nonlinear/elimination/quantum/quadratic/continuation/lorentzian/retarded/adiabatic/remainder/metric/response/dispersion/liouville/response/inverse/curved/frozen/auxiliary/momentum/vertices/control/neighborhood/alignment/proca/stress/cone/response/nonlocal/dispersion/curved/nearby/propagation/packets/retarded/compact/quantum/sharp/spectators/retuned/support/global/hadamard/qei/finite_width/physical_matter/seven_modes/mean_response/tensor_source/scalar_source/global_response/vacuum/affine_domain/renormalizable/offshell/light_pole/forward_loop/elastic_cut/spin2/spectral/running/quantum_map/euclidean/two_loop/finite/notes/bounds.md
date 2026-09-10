# Actual integrated coefficient enclosure

For a graph with h heavy edges, combine the previous estimate
with its absolute coupling lambda4^(3-h)g^h. The frozen
parameters satisfy g/M<lambda4/3, so g/B<lambda4/6 for B=2M.
Also pi>3. The graph's second forward coefficient is bounded by

    lambda4^3 [17(h+2)! 16 (2/3)^h / 41472]
        [C0+C1 ln B+C2(ln B)^2].

Here 41472=2(16 times 3^2)^2 includes the Taylor factor and
both loop measures. The exact source verifies the fourfold
Schwinger rescaling, mass-power cancellation and coupling
powers independently.

Weight and sum the 88 graph bounds using the frozen family
weights and all choice tuples. The resulting polynomial is

    lambda4^3 [
        13688587/17496
       +(515831/972) ln B
       +(59789/972) (ln B)^2 ].

Every coefficient is rebuilt from the exhaustive sector
calculation, rather than inserted as a presumed integral value.

For the actual parameters 64<B<10^200. The positive partial
sum 1+3+3^2/2!+3^3/3!=13 proves exp(3)>10, hence
ln(10)<3 and ln B<600. This avoids floating evaluation of
the very large mass hierarchy. Positivity of every polynomial
coefficient then gives the rational majorant

    Efinite=(393017383387/17496) lambda4^3.

Exact rational comparisons establish

    Efinite < 3 times 10^-607,
    Efinite/(4 lambda) < 10^-7.

A floating diagnostic is approximately 2.69344 times 10^-607,
or 6.73359 times 10^-8 of the positive tree coefficient.
Those decimals are not used by the proof.

The signed finite-subsector correction lies between minus
and plus Efinite. The bound is not a sign determination,
not the complete two-loop coefficient and not a bound for
the excluded ultraviolet-subtraction graphs or LSZ terms.
