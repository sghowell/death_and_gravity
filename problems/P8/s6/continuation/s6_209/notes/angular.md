# Actual angular coefficient filtration and leading shapes

All coefficients are first calculated componentwise in a real tracefree tensor basis. Their complex Fourier action follows by sesquilinear extension. The normalized full projector contractions do not require a global individual polarization frame.

For k=n, l=-n+xP and n.n=1, the 00, TL, LT and LL geometric numerators have angular degree at most four before inverse-radius expansion. The 01 and 10 contractions simplify using P_k C_k^T=C_k^T and C_l P_l=C_l. The 11 contraction uses C_k C_k^T=k^2 I-kk^T and its l counterpart. Their angular degree is also at most four.

Every additional external P in the regular denominator/frequency expansions costs at least one x. The normalized frequency, inverse phase, longitudinal constraint and all W8 coefficients are analytic at x0; their degree-d coefficients have angular degree at most d. Source time differentiation holds comoving momenta fixed. The full degree-d current coefficient therefore has angular polynomial degree at most d+4.

Choose phat=e3. Azimuthal averaging sends odd X or Y monomials to zero and gives

average X^(2a)Y^(2b)
=(1-u^2)^(a+b)(2a)!(2b)!/[4^(a+b)a!b!(a+b)!].

Consequently each Fbar_jd(u) is a polynomial of degree at most d+4<=8. This finite angular filtration is not a claim that the original moving-band integral is a spatial polynomial.

For T=tr(DG), V=(Dphat).(Gphat), W=(phat.D.phat)(phat.G.phat), the complete leading Proca coefficient is

Fbar0=[T(3/4+5u^2/2+3u^4/4)
+V(5/2-3u^2-15u^4/2)
+W(9/8-45u^2/4+105u^4/8)]/(8a).

It includes the fixed-mass longitudinal mode. Its four shapes are

pi(18T-12V-W)/(64a),
-pi(89T-176V)/(420a),
pi(58T-228V+135W)/(1536a),
pi(21T-32V-32W)/(840a).

The sign and Fourier factor are minus K^(4-h)p^h/(2pi)^3. The first shape reproduces S6.206. The last shape is the finite artifact of j0/d0 alone, not the complete conversion; other UV slots cancel it.

For the uniform tail, a degree-at-most-eight polynomial bounded by S on [-1,1] has Chebyshev coefficients bounded by S for the constant and 2S for the others, directly from cosine orthogonality. The first and second derivatives of each T_n have nonnegative Chebyshev coefficients, verified exactly for n1,...,8. Their endpoint sums are n^2 and n^2(n^2-1)/3. Thus the first derivative is bounded by408S and the second by5712S. The sum of the three norms needed below is at most6121S. The argument also applies to complex polynomial coefficients by the same absolute coefficient bounds.
