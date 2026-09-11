# Complex contour and comparison bounds

At each fixed physical epsilon,t,k select the contour |zeta|=omega/m,
m=1000. The S6.191 ordered coefficient majorants give

    ||rhat||,||rsharp|| <= A0/m < .01,
    ||partial_epsilon rhat||,||partial_epsilon rsharp|| <= A1/m,
    ||partial_epsilon^2 rhat||,||partial_epsilon^2 rsharp|| <= A2/m.

The constants A1 are about125.062 and A2 about6492.621. Derivatives
are taken at fixed zeta before selecting this contour at the
evaluation point. No derivative of an amplitude-dependent contour
is omitted. This is a pointwise Cauchy estimate for each differentiated
holomorphic function.

Both graph factors have norm below .01, so the inverse is bounded
by 1/(1-.01^2) without any positivity claim at complex zeta. Also
||F|| and ||Fsharp|| are bounded by sqrt(1+.01^2). Ordered inverse
and product differentiation, valid for the two independently
bounded graph factors, gives the same ample first and second
Frechet constants16 and128 used in S6.191. Coefficientwise
conjugation does not enlarge these uniform circle bounds.
Consequently the complex covariance and its first two physical
derivatives have norms below

    2, 16 A1/m <3, 16 A2/m+128(A1/m)^2 <128.

The actual vertex and frequency are marker-independent. Combining
their complete physical derivatives with the trace factor6/2 gives
contour current bounds (12,180,6000)omega. In particular the
unrounded second coefficient is5478, not the result obtained by
differentiating only the covariance.

Define J_ad4 solely as the marker Taylor polynomial of this current
at orders0,2,4. It is a mathematical comparison. This definition
does not assert that it is the original covariant subtraction and
does not select new finite counterterms.

Cauchy coefficients, together with the proved even parity, give for
a=0,1,2 and omega>=nu_minus>=K=1e16,

    |partial_epsilon^a(J_ref10-J_ad4)|
      <= B_a omega (m/omega)^6/[1-(m/omega)^2]
      < 2 B_a m^6 omega^-5.

Here B=(12,180,6000). The finite reference is not claimed to be a
convergent infinite adiabatic expansion. The convergent object on
this contour is its finite-graph covariance inverse.

For the finite complement nu_minus<K, each of the three retained
Taylor coefficients is bounded by B_a omega because omega/m>=1.
Thus |partial_epsilon^a J_ad4|<=3 B_a omega<=9 B_a nu_minus.
