# From the full Schwinger integral to a convergent box majorant

The light mass squared is one; M is the heavy mass squared.
For a selected two-loop graph, integration of both Euclidean
loop momenta with measures d^4k/(2pi)^4 gives

    (16pi^2)^-2 integral U^-2 exp(-F/U) d^N alpha.

Completing the square in each pair of loop components gives
pi/sqrt(U); its fourth power and the eight factors of 2pi
give precisely this normalization. The full masses and
external invariants remain in F. Every propagator power is
one, so its Schwinger representation introduces no additional
Gamma factor.

The crossing variable is v=s-2, with u=2-v and t=0.
F is affine in v and dF/dv=Pu-Ps. Its second Taylor
coefficient therefore has absolute integrand bounded by

    [2(16pi^2)^2]^-1 U^-2
        [(Pu-Ps)/U]^2 exp(-Re F/U).

The factor one half is the Taylor factor, not a symmetry
factor of a graph. The already independent S6.120 bounds give

    |Pu-Ps|/U <= sum alpha,
    Re F/U >= sum alpha/4+(M-16)Halpha.

These hold throughout the compact unit forward disc.
After alpha=4 beta, homogeneity contributes 4^(N-2)=4^(h+2),
and the exponential is bounded by

    exp[-Lbeta-(4M-63)Hbeta]
       <= exp[-Lbeta-B Hbeta],  B=2M,

since the actual M exceeds 32. No denominator is replaced
by a finite series before integration.

It remains to bound

    I(B)=integral U(beta)^-2 (sum beta)^2
                    exp[-Lbeta-B Hbeta] d^N beta.

Let q=max({light beta_e},{B times heavy beta_e}). The
exponential is at most exp(-q), and the positive layer
identity is exp(-q)=integral_q^infinity exp(-r) dr.
Tonelli's theorem applies to this nonnegative bound, even
before convergence has been proved.

The region q<=r is r times the anisotropic box of
notes/sectors.md. The integrand without its exponential
is homogeneous of degree -2. Consequently

    I(B) <= (h+2)! integral_box U^-2 (sum beta)^2 d^N beta.

The factorial is the exact integral of r^(N-2)exp(-r).
On the box, sum beta<=4+h/B<=4+3/64 and its square is
strictly below 17. The sector theorem therefore gives

    I(B) <= 17(h+2)! B^-h
                [C0+C1 ln B+C2(ln B)^2].

This finite majorant controls the entire unbounded Schwinger
domain. It is not just an interior or finite-box sample.

The same reasoning without one or both powers of sum beta
bounds the zeroth and first derivatives as well: the layer
factor becomes Gamma(h+n+1) for derivative order n=0,1,2.
Thus the full graph integral and the required derivatives
are absolutely controlled, not only a formal second derivative
of a divergent integral.

To specify the first sheet, start at zero external momenta
in the ordinary massive Euclidean representation and scale
the final complex external momenta by sqrt(tau), 0<=tau<=1.
Their squared masses and pair invariants all scale by tau.
The negative kinematic terms in Re F never exceed the same
coarse bounds along this path. The positive majorant gives
dominated continuation to the stated disc. Its strict gap
also permits a small neighborhood; an increment of modulus
at most 1/8 reduces the gap by at most sum alpha/8.
The same sector argument with a weaker positive mass gap
establishes local holomorphy there.

This argument applies only to the selected individually
UV-finite graphs. For the excluded 104 graphs the small
sector exponents need not be positive; subtraction work
cannot be replaced by these estimates.
