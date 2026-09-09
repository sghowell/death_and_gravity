# Actual new massive block, including contacts

The minimal physical canonical factors and readouts are the new ones, not the old lapse-dependent-mass functions. Independently differentiating the physical Hamiltonian gives pair vectors

    P_T=(0,2(1-z)), P_L=(0,2(1+z)),
    2 P_T P_T^T + P_L P_L^T = diag(0,4(3-2z+3z²)).

The second fixed-canonical Hamiltonian derivative is retained as a seagull contact. Contact plus zero-frequency bubble equals the static vacuum Hessian in both sectors. All sixteen flat adiabatic/frequency-Taylor comparisons through subtraction orders two and four are recomputed with the actual mass jets zero. The exact three-subtraction denominator and radial normalization are checked separately.

The new finite local physical fourth coefficient is diag(0,-4), independently taken from S6.84's actual curved local Euler calculation. Thus the normalized nonzero block is F(p)=-H(p), with

    H(p)=4+integral_0^1 W(y) p/[4m²+p(1-y²)] dy,
    W(y)=y²(3-2y²+3y^4), p=s².

This is the exact massive reference bubble plus its fixed finite local fourth term, divided by s^4/(64 pi²). It is not a massless logarithmic approximation. The integral defines the analytic continuation on the p-plane cut along (-infinity,-4m²].

Writing d=1+4m²/p, exact polynomial division of the three radial moments gives

    H=16/15+d-3d²+sqrt(d)(3-2d+3d²) atanh(1/sqrt(d)).

The branch is fixed first for p>0 and then by the integral's analytic continuation. The apparent p=0 singularities of this closed expression are removable. Independent radial integration yields H(0)=4, H(-4m²)=16/15 and H'(0)=9/(35m²).

At large p on the first sheet,

    H(p)=2 log(p/m²)-14/15+O((m²/p) log(p/m²)).

The explicit logarithm identity at d=1+epsilon makes this expansion uniform on large cut-plane circles. Its leading expression alone has a spurious zero at log(p/m²)=7/15. The exact H is at least four on the positive axis, so that truncated root is not a pole of this massive block.

Multiplying diag(0,F) by diag(0,1/F) yields the active projector diag(0,1), not the identity. No pseudoinverse is silently promoted to a complete two-source inverse.

