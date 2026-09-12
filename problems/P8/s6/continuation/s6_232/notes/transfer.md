# Complete original transfer data and positive endpoint derivative

Write s=2-t/2+v and u=2-t/2-v in the unchanged original amplitude

A_tree=2lambda[(s-2)²+(t-2)²+(u-2)²]+3gamma stu-8gamma.

Its v² coefficient is exactly4lambda-3gamma*t, and its v4 coefficient is zero. The potential remains in the full amplitude even though it does not contribute these derivatives. The original higher independent vertices do not change connected four-scalar tree scattering, as established by the source-pinned S177 and S231 ownership arguments.

The normalized identical channels remain

A=(32pi/beta) sum_even(2l+1)t_l P_l(x),
x=1+2t/(S-4), beta=sqrt(1-4/S).

Exact unitarity of an existing full channel matrix implies Im t_l>=0. Rodrigues differentiation gives P_l'(1)=l(l+1)/2, so whenever the stated endpoint series and derivative converge,

rho_t=(32pi/beta) sum_even(2l+1)Im t_l*l(l+1)/(S-4)>=0.

This does not assume that the full amplitude has only l0,l2. All exact partial waves and inelastic channels are allowed. It is a conditional continuum implication, not a claim that checking finitely many Legendre polynomials proves convergence.

At the first original elastic order only the original tree waves contribute. With Q0=a+b/3 and the exact a,b from S231,

rho_first(S,t)=(32pi/beta)[t0_tree²+5t2_tree² P2(1+2t/(S-4))].

Its forward value is the COMPLETE S231 first elastic coefficient. Its derivative is

rho_first,t(S,0)=beta b(S)²/[60pi(S-4)].

Since b=(S-4)²(lambda-3gamma S/4), this is threshold-regular. No massless or potential-free approximation is used.

An independent angular-convolution test obtains the same expression directly from the two-body cut. Take y=n_in.n and z=n_in.n_out. The intermediate azimuth average of(n_out.n)² is z²y²+(1-z²)(1-y²)/2. The cut is therefore

beta/(64pi) integral_-1^1
(a+b y²){a+b[z²y²+(1-z²)(1-y²)/2]} dy.

Its evaluation equals the partial-wave expression. The nonforward cut is a convolution of two amplitudes, not the square of the amplitude at one scattering angle. A negative control detects that error.

The crossing shift t/2 in the dispersive denominator is essential. Differentiating
2rho(S,t)/[pi(S+t/2-2)^3]
gives the negative inverse-fourth term -3rho/[pi(S-2)^4]. Omitting the shift loses this term. By crossing evenness, the first b21 jet evaluated at t0 happens to agree between fixed-s and fixed-v differentiation; no false claim to the contrary is used.
