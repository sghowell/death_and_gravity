# Exactly two nonreal first-sheet zeros of the complete denominator

Use the ORIGINAL analytic factor

A2(p)=1/30+integral_0^1 p W2(y)/[4m^2+p(1-y^2)]dy,
W2=y^2(30-20y^2+3y^4)/30>0

on the first sheet slit along p<=-4m^2. Its value at p0 is removable. Let D_C(p)=C+pA2(p), where first consider any real C>2m^2/15.

## No real zero and fixed upper-cut image

On(-4m^2,infinity), A2 is strictly increasing. For negativep it is at most1/30, hence pA2>=p/30>-2m^2/15. For positivep, both factors are positive. Therefore D_C>0 on the entire real interval above the cut.

For p=-tau+i0 on the open upper cut, the original factor has A2=D_A+i pi U, with U>0 from S223. Consequently

Im D_C=-pi tau U<0.

The threshold value is C+688m^2/225>0. There is no cut-bank zero. The factor itself has no pole anywhere in either open half-plane.

## Complete argument count, including infinity

Take a sufficiently large upper half-disc, approaching the real boundary from above and making a vanishing indentation at the threshold if needed. Along the lower boundary from-R to-4m^2, the D_C image remains strictly in the lower half-plane. Its continuous argument lies between-pi and0 and ends at0. The remainder of the real boundary maps to the positive axis.

On the full large semicircle p=R exp(i theta), 0<=theta<=pi, the S226 stable first-sheet expansion, now with the explicit error in notes/bounds.md, is

A2(p)=(13/60)Log(p/m^2)-52/225+o(1)

uniformly, including both limiting banks. Thus D_C/p=A2+C/p lies in a small-angle neighborhood of the positive real axis on the arc. Its argument has a continuous small lift; at theta0 it is0 and at theta pi it is a strictly positive angle epsilon_R, because the upper-cut imaginary part of A2 is positive. The D_C argument along the outer arc changes from0 to pi+epsilon_R.

At the starting point-R on the lower boundary, the principal D_C argument is-pi+epsilon_R. The lower boundary therefore contributes pi-epsilon_R, and the outer arc pi+epsilon_R. The TOTAL change is2pi. The threshold indentation has vanishing contribution because its D_C limit is positive. A small positive-height approximation yields ordinary analytic contours with the same winding.

The argument principle gives exactly one upper-half-plane zero COUNTED WITH MULTIPLICITY. It is therefore simple. Schwarz reflection gives exactly one simple lower zero, its conjugate. Combined with the real and bank exclusions, there are exactly two first-sheet zeros and no others. A sampled winding diagnostic is not the proof of this count.

For p=x+ib with b>0, the exact identity

Im[p^2/(M+c p)]
 =b[2xM+c(x^2+b^2)]/[(M+cx)^2+c^2b^2]

shows Im[pA2]>0 when x>=0, using the additional positive b/30 term. Hence the upper zero z has Re z<0, Im z>0. The lower one is its conjugate.

The actual C in notes/bridge.md satisfies the count condition by a very large margin. The isolated S223 A2 zero p=-r m^2 is NOT a zero of D_C: D_C=C there. Its absence from the complete reciprocal is exact composition with the original tree and finite lower term, not manual pole deletion.

For each spatial q=P^2>=0, lambda^2=z-q has a root lambda_plus with positive real part and positive imaginary part; the conjugate zero gives lambda_minus=conj(lambda_plus). Both poles are on the original causal Laplace sheet. A real-axis Fourier prescription that excludes them is not substituted for the original Bromwich inverse.
