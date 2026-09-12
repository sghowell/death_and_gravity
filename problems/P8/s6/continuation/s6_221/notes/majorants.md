# Uniform denominator, positivity and all48 complete coefficient bounds

In the central region take z<=1/4096, |E|<=1/2, |E|>1/4, J<=100, ell<=1/10, |A|,|Tcorr|<=1e-6. If e=1e-6 and z0=1/4096, the difference between Delta and4E² is bounded by
perturb=6e z0+9e²z0²+(2z0+9e z0²)(200+1/400).
The exact rational interval is
2555883493623977491/16777216000000000000
< Delta <
18415636506376022509/16777216000000000000.
It lies strictly inside(1/8,2). This proves XX positive definite for all certified central times and transfers, not just test points.

The complete central K has determinant8J/Delta>1/25. Its diagonal formulas are
K11=4(2J+w²)/Delta,
K22=[(2E+3Tcorr z)²-4Jz-18AJz²]/Delta.
Using positivity already proved, drop the negative -4Jz only for an upper bound. Their sum is below1e4. Outer detK=2J/Theta²>1/200 and its trace is below(200+1/100)*16+1<1e4.

For both G, the determinant is2F divided by E² or Theta². Since F>1/100, |Theta|<=2 and |E|<=1/2 in their domains, both determinants exceed1/200. Their lower-right entry is1, so they are positive definite. Since C=F+w²/2<Jbare+w²/2<100+1/200, both traces are below2(100+1/200)*16+1<1e4.

For a positive two-by-two symmetric matrix, lambda_max<=trace and lambda_min=det/lambda_max>=det/trace. Consequently all four actual K,G obey
1e-8 I < K,G <1e4 I.
The determinant argument is essential; a bound on entries alone would not establish positivity.

It remains to control complete lower-order terms and time variation uniformly, not only positive principal matrices. There are six two-by-two matrices in each chart: K,K',G,G',Omega,R. Substitute q=1/z in the central rational expressions. Multiplication by E^4 Delta² makes every one of the24 entries a polynomial in z and the coefficient jets. Outer multiplication by Theta^4 does the same for all24 entries, with no q dependence remaining.

The certificate constructs every polynomial exactly. It evaluates its absolute coefficient polynomial on:
|Theta|<=2, |E|<=1/2 centrally or1 outside, ell<=1/10, J<=100, |A|,|Tcorr|<=1e-6, |H|<=2, z<=1/4096 and each first/second coefficient jet<=1e6.
Divide back with |E|^-4 Delta^-2<4^4*8² or |Theta|^-4<4^4. Each resulting rational bound is recorded. All48 are below1e18; the largest is below9.88e16. No time, momentum or coefficient sampling is used in these bounds. The independent dynamics tests supplement rather than replace them.

In particular, every corresponding matrix operator norm is at most2e18. The rough common derivative box is intentionally much larger than the known stress jets. It uses only the previously established five stress derivatives, not an invented thirteen-jet smallness estimate.
