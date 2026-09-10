# Regulated local sunset slope

Set d=4-2epsilon, Q=16 pi^2, ell=log(mu^2) in light-mass
units. The Euclidean connected two-vertex Phi^4 sunset has
positive L^2/6. Its contribution to the Euclidean inverse is
negative. With D(s)=-Gamma_E(-s)=s-1+Pi(s), the contribution
to Pi is therefore positive L^2/6, as in S6.126.

On the unit simplex let U=xy+xz+yz, P=xyz/U and w=xyz/U^3.
The dimensional Schwinger integral gives

    Pi(s)=L^2/(6Q^2) e^(2gamma epsilon+2ell epsilon)
      Gamma(-1+2epsilon) integral U^(-2+epsilon)
                                      (1-sP)^(1-2epsilon).

First establish this in the common convergence region
1/2<Re epsilon<1, then continue meromorphically. Differentiating
there before Gamma recurrence gives

    Pi'(s)=L^2/(6Q^2) e^(2gamma epsilon+2ell epsilon)
      Gamma(2epsilon) integral w U^epsilon (1-sP)^(-2epsilon).

The differentiated parameter integral is holomorphic near
epsilon=0: the sector density below is O(t r^Re epsilon).
All proper local Phi^4 cycles contract both external vertices,
so their residual tadpoles have zero momentum derivative.

Let T0=integral w, T1=integral w log U and
C=-integral w log(1-P). The six ordered sectors (1,r,rt),
normalized by a=1+r(1+t), have v=1+t+rt and

    w dx dy=t/v^3 dr dt, U=r v/a^2, P=r t/(a v).

The r primitive -1/(2v^2) yields T0=1/2 exactly.
The Laurent slope has pole L^2/(24Q^2 epsilon), and its
finite MS value at s=1 is

    L^2/(6Q^2) [ell/2+T1/2+C].

The e^(2gamma epsilon) factor removes the Gamma constant.
In particular, neither T1 nor C is a tunable finite field
condition. The slope is recorded BEFORE outer OS. Numerics
are only cross-checks; no guessed closed form for T1/2+C
is used in the enclosure.
