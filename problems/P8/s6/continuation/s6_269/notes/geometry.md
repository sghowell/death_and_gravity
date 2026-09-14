# Complete nonlinear geometry and mean-zero cotangent lift

Use the physical A_s Fourier norms from notes/reference.md.
The weights obey (1+|k+l|)^s <= (1+|k|)^s(1+|l|)^s
for every s>=0. Scalar, matrix-operator and vector-Euclidean
norms therefore give the required product estimates in A2.
Finite starting modes are smooth. All nonlinear products,
fixed points and inverse operators below are UNTRUNCATED.

The S267 shape proof uses only this Banach product bound
and the bounded projectors P(k), P(0)=I. It therefore applies
in A2, without asserting that the high-P family is small
in its older A8 norm. Set t=10^-274 and assume
||tau||A2<=t, ||f||A2<=t. The full fixed-point map has

    a=2t,
    ||map(f)||A2 <=3a^2+3a^3<t,
    Lip(map)<=6a+9a^2<1/10.

Here trace(Bf)=(2I+Pi0)f, including the generated mean.
Thus Q=I+tau+Bf has ||Q-I||A2<=2t, detQ=1,
divQ=0 and Q>0. Smooth input gives smooth output by
translation-equivariance and the smooth analytic fixed-point
map in A2, not by a same-regularity smooth Diff assumption.
In particular the finite input's full nonlinear harmonics
and homogeneous tracefree shape are retained.

At u=0, gamma=exp(2v)Q^-1 and g=gamma^-1=exp(-2v)Q.
Put v0=10^-329. The Neumann and exponential Banach series give

    ||Q^-1||A2<=1/(1-2t),
    ||exp(+-2v)-1||A2<=2v0/(1-2v0),
    ||gamma-I||A2, ||g-I||A2<e=10^-272.

The exact intermediate rational bounds are in geometry.bounds.
The norms of gamma,g,c=exp(2v),V=exp(3v) and V^-1
are all below2. First and second spatial derivatives of
gamma-I and g-I obey the same e bound because the
physical derivative multipliers are covered by A2.

For the entire original ghost block M=A V_Lie, the
flat nonzero symbol and inverse on the fixed L=1 torus are

    M0(k)=|k|^2 I+kk^T/3,
    M0(k)^-1=|k|^-2(I-kk^T/(4|k|^2)).

All nonzero |k|>=1, so ||M0^-1||A0->A1<=2,
and its first and second derivative multipliers have
operator norm at most1. Constants are not inverted.

The COMPLETE adjoint on the exact divergence-free slice is

    (M* lambda)_l =
      -Q_jk partial_j partial_k lambda_l
      -(1/3)Q_ij partial_l partial_j lambda_i
      -(1/3)(partial_l Q_ij) partial_j lambda_i.

The derivative-Q contact is mandatory. The code first
differentiates the generic off-slice adjoint, including
all divQ and derivative-divQ terms, and only then proves
this exact slice identity. The terms with a second
derivative have at most nine component summands each,
and the derivative-Q term has the same finite component
sum. Bounds by component absolute values, conversion
back to Euclidean vector norm, and all factors1/3 are
covered conservatively by

    ||(M*-M0)M0^-1||A0->A0 <=100||Q-I||A2 <1/10.

Both M* and the perturbation preserve the mean-zero
subspace: their pairing with constants vanishes because
M kills constants on the slice. The right-composed
Neumann series consequently gives a genuine full inverse

    (M*)^-1=M0^-1[I+(M*-M0)M0^-1]^-1,

with A0->A1 norm below4. This is a full convolution
operator result, not the inverse of a Fourier compression,
and not an inverse on the three residual translations.

For curvature, each complete Christoffel component has
A1 norm at most9e: the three inverse-metric summands,
three metric derivative terms, factor1/2 and inverse
norm2 are included. Each Ricci component has A0 norm
at most6(9e)+18(9e)^2. All nine scalar contractions,
each with inverse-metric norm below2, give

    ||R_gamma||A0 <=18[6(9e)+18(9e)^2]<10^-268.

No derivative or quadratic connection term is discarded.

Use the exact S267 reduced-covector extension and full
original generator, with all matter/vector momenta:

    pi0=(Pi_v/6)g + DQ* (P_TT Pi_tau),
    D0=V_Lie* alpha0,
    alpha=alpha0-A*(M*)^-1 Pmean0 D0.

The complete dual is
DQ*[B]=c[(g/3)tr(Bg)-gBg]. Its A1 operator norm is
at most16: use c<=2, g<=2 and |tr(Bg)|<=3||B||||g||.
The Fourier TT projector has operator-norm bound2
including its five-component zero sector. Every one
of the five canonical momentum fluctuations has
A1 norm below10^-200 by the original state bounds, so

    ||pi0||A1 <=(1/3+32)10^-200<10^-198.

The entire original generator has the conservative bound

    ||D0||A0 <=180||pi0||A1
       +12||W||A1 ||Pi_W||A1
       +2[||Pi_M||A0 ||grad M1||A0
          +||Pi_H||A0 ||grad H||A0] <10^-194.

The metric divergence and metric derivative terms are
both included in180. The vector coefficient includes
-W_i div Pi_W, not just its derivative-along-W term.
Pi_M retains the original1/10 background and has norm
below1. The full source identity and its boundary flux
are checked before using these bounds.

For lambda=(M*)^-1 Pmean0 D0, ||lambda||A1<=4||D0||A0.
The full cotangent correction has A0 norm at most
16||lambda||A1, since div* has norm at most1 from
vector A1 to symmetric-matrix A0. Hence

    ||pi||A0<10^-192, ||pi_TF||A0<10^-190.

The exact generic non-diagonal identity
tr(gamma DQ*[B])=0 for EVERY B is important:
Pi_v=2 pi:gamma is unchanged by this correction.
It does not acquire the coarse shear-lift error.
All inverse-metric and trace terms enter the direct
symbolic and independent rational-metric checks.

Finally |tr(pi_TF gamma pi_TF gamma)|/V^2 is bounded
by48||pi_TF||^2<10^-378, including trace dimension3,
both gamma factors and V^-2. Global translations are
still retained, not set to zero on arbitrary coherent
displacements. The S268 symmetry reduction applies
to the operators and same seed.
