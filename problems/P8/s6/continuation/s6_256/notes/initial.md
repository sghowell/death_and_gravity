# Full-current initial root and quantitative persistent margins

All full R,F bindings and their derivatives are applied before setting
u=0,N=1+10^-6. At that point set a_hat=1,Hhat=hbar=hbar_dot=0.
The full R obeys R_u(0,N)=0, but its R_uu and lapse derivatives are
retained. Define

    U=R^-3/4, D=R^1/4/N, Z=U/N,
    Tprimitive=3 U R_uu R_N/(4 R N),
    G=(N U F)_N-Tprimitive,
    mc^2=-2G/Z_N.

This is the whole S6.255 initial lapse constraint, not a tree-only
constraint. The primitive identity I_N=3 U R_u R_N/(4 R N) removes
I_uN exactly. Its integration constant is not changed.

Before evaluating the full jets, derive

    Cnn=[(N U F)_NN-(Tprimitive)_N+Z_NN mc^2/2]/2,
    J=Cnn-Z_N^2 mc^2/(2Z),
    Hhat_dot=[-U R_uu/(2N)-N U F-Z mc^2/2]/(2D),
    Theta_dot=-D_N Hhat_dot-(U R_uu/(2N))_N/2-Tprimitive/2,
    L2=4(N R^3/4/2)_N.

At u=0 the exact tree jets are generated from

    Rtree=N^-2,
    Rtree,uu=-6(N^-2-1),
    Ftree=-(624 N^-4+753 N^-2+224)/200.

For EACH full jet, use an interval of radius 10^-300 around the
exact rational value. This includes the actual profiles and constants.
Exact interval evaluation proves the 17 displayed enclosures, including

    .009 < mc^2 < .011,
    1.51 < Cnn,J < 1.53,
    -3*10^-6 < r=R-1 < -10^-6,
    -3 < R_N < -1,
    -1.001 < L2 < -.999,
    .999 < R,D,Z,Y,Yv < 1.001,
    .499 < C < .501,
    .999*10^-6 < K < 1.001*10^-6,
    .099 < c_M1 < .101, .049 < w_M1 < .051.

The separate initial diagnostics give
3.99<Hhat_dot<4.01 and 2.99<Theta_dot<3.01.
They are NOT asserted as persistent derivative boxes.

For all 14 persistent coefficients, the code evaluates the actual
distance of the full initial enclosure from BOTH box endpoints.
Every distance exceeds 2*10^-20. The subsequent complete first-time
derivative bound 10^40 on an interval 10^-60 therefore cannot exit
these boxes. This quantitative margin check replaces an unspecified
continuity radius.

Theta initially vanishes because R_u(0,N) and all its lapse
derivatives vanish. Hhat initially vanishes by the datum. The positive
M1 square root lies strictly inside [.09,.11]. The datum satisfies the
whole constraint and has strict J,Cnn and Gamma margins. This makes
it eligible for S6.255's unforced local Euler construction at the
specific, fixed epsilon, without claiming it is the original quantum
mean or a physical-scale bounce.
