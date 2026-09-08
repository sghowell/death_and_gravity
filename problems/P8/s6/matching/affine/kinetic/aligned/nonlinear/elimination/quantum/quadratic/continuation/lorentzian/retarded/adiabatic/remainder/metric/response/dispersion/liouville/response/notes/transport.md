# Covariance transport with the actual prepared mode clocks

Let f=d sigma/du>0 and Y=(S,C,P), with S=|v_sigma|^2,
C=Re(p_sigma v_sigma*) and P=|p_sigma|^2. S6.70 gives

    Y=J X, J=diag(f,1,f^-1),
    Y'=M(b_k,D)Y,
    M=[[2b_k,2,0],[-D,0,1],[0,-2D,-2b_k]],
    D=k^2+U, b_k=C_k'/C_k.

Primes here denote acoustic derivatives. The corresponding
physical-time matrix has half-rate d and frequency Omega:
b_k=(d+fdot/(2f))/f, D=Omega^2/f^2. Direct matrix conjugation,
including Jdot, verifies the transformation. The invariant
SP-C^2 is unchanged. Its linearized invariant is also preserved
under the inhomogeneous covariance response.

## Prepared response at a fixed value of sigma

Write r=delta log f and fix sigma_e(u0)=sigma_0(u0). Then

    xi(u)=delta sigma(u)=integral_{u0}^u f r du,
    delta Y|_sigma=J delta X|_u+delta J X-xi Y'.

For W=delta Y|_u and Z=delta Y|_sigma, W=Z+xi M Y.
Differentiating this identity and using xi_dot=f r gives

    Z'=M Z+delta M|_sigma Y,
    delta M|_sigma=delta M|_u-xi M'.

This identity uses matrix order exactly; no stationary or
commuting-coefficient assumption is needed.

For the original clock f=1/a in both sectors, but their variations
are different:

    r_T=n-zeta,
    r_L=[1+(beta-alpha)/2]n-zeta.

The shifts xi_T and xi_L therefore remain distinct. On the
original clock U=a^2m^2 in both sectors; their variations are

    delta U_T=U(2zeta+beta n),
    delta U_L=U(2zeta+alpha n).

The transverse b_T is exactly zero in D=3. Longitudinally,

    b_L=b-U'/[2(k^2+U)],
    delta b_L|_u=delta b-delta U'/[2(k^2+U)]
                     +U' delta U/[2(k^2+U)^2].

Here the varied U' and b already include S6.70's varied derivative.
Then subtract xi b_L' and xi U' to form fixed-sigma sources.
All source coefficients independently reproduce the original
physical moving canonical map and frequency variation.

The potential variation computed from these covariance sources,

    delta V|_sigma=delta U|_sigma
                  -(delta b_k|_sigma)'-2b_k delta b_k|_sigma,

equals S6.70's full longitudinal potential variation, including
its history term, and the transverse delta U_T|_sigma.

The preparation is zero on an initial neighborhood, so r,xi and
all local source jets vanish there. Thus Z(u0)=0 is the same
zero varied covariance as in the original physical problem.
The original covariance is still that of the selected all-order
state; its nonzero mixing is not erased by setting Z(u0)=0.
The retarded initial-value problem is equivalent in both clocks.

## Exact mode Green function

For a normalized actual solution v with
v v*'-v' v*=i, define, for sigma>=s,

    G(sigma,s)=i[v(sigma)v*(s)-v*(sigma)v(s)].

It vanishes on the diagonal, has unit first-derivative jump,
and solves the homogeneous mode equation away from the diagonal.
Consequently, for prepared delta V at fixed acoustic time,

    delta v(sigma)=-integral_{sigma0}^sigma G(sigma,s)
                                      delta V(s)v(s) ds.

The variance response is

    delta S(sigma)=2 integral Im[(v(sigma)v*(s))^2] delta V(s) ds.

In flat vacuum this kernel is -sin(2omega(sigma-s))/(2omega^2),
in agreement with the previously fixed Kubo convention. The
formula uses the actual selected modes, not a new flat state.
For each finite k, ordinary smooth ODE parameter dependence
justifies the derivative and equivalence to covariance transport.

A source compact in physical time need not give xi zero at the
final endpoint. Its later constant coordinate offset is retained
in both the potential source and the physical readout pullback.
No final boundary or covariance reset is introduced.
