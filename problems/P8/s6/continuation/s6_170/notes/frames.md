# Full curved transition and both tails

After the physical mass rotation the off-diagonal axis is sigma2.
The exact two-level rotations used in S6.167 give

    e_(j+1)=sqrt(e_j^2+g_j^2),
    g_(j+1)=(-1)^j partial_t atan(g_j/e_j)/2.

The axes alternate2,1,2,1. Every connection is retained.
The construction is the exact unitary iteration discussed in
[Barbero et al., section IV](https://arxiv.org/pdf/1805.05107);
the new variable-disk constants and tail integrals are derived here.

Remove the final diagonal phase. The remaining exact generator has
norm|g_N|. Unitarity and Duhamel's formula bound the full transition
amplitude by integral_R|g_N|, including every quadratic transition
Dyson order. This is not an all-Feynman-loop assertion.
The same integral also bounds either half-line transition from its
physical asymptotic state to the N-frame instantaneous projector.

Write K=K_N and Ec=sqrt(p^2/16+m0^2).
For |t|<=1, q<=p, E>=Ec, L>=tau. The mass term integrates using
integral w(t/tau)dt<=2tau. For the geometric term use length2,
D>=1 and L^(1-N)<=tau^(1-N). Thus the compact contribution is at most

    C_N=32K^N p(Delta+m tau)/(tau^N Ec^(N+2)).

For |t|>=1, L>=|t|, w(t/tau)<=tau^9/|t|^9 and
L/D<=2/|t|. Since Delta<=.003m and tau<=1, the initial
bracket is at most3m/|t|. Also

    p/(4t^4)<=q<=p/t^4.

The sum of both tails is therefore bounded by

    96mp K^N integral_1^infinity
        t^(-N-5)[m0^2+p^2/(16t^8)]^(-(N+2)/2)dt
      =96mK^N(4/p)^(N/4) integral_0^(p/4)
        y^(N/4)/(y^2+m0^2)^((N+2)/2)dy.

The equality uses y=p/(4t^4), including its Jacobian.
The full y integral converges for every N>=1. It scales as
m0^(-3N/4-1). Therefore the exact full-line and half-line
transitions are O(p^(-N/4)), for any fixed N at high enough p.
The t~(p/m)^(1/4) tail is why four frames would give only p^-1
on this geometry, not the flat p^-5 estimate.

For N20 the full y integral is1/(720m0^16).
A second, low-momentum estimate simply uses E>=m0 before
the t integration. Their minimum gives

    tail<=min(Ap,B/p^5),
    A=4mK^20/m0^22, B=2048mK^20/(15m0^16).

Use its continuous value zero at p=0. Both bounds hold at every
positive p; neither extrapolates a UV asymptotic into the infrared.

For a convenient uniform bound, p/Ec^22<=4/m0^21 and
splitting the tail at p=m0 gives

    |beta_p|<=128K^20(Delta+m tau)/(tau^20 m0^21)
              +2048mK^20/(15m0^21)<1e-1890.

All 42 copies are covered using the common amplitude upper bound.
In particular, the Delta=0 bound stays nonzero because the geometry
can produce pairs of constant-mass particles.
