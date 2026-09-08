# Literal mass-tensor modes and equivalent causal variations

Use physical lapse N>0, physical scale a>0,
q=k_com²/a², positive temporal/spatial masses a_m,b_m,
and m²>0. The transverse quadratic coefficients are

    g_T²=a/N,
    Omega_bare,T²=N²(q+m² b_m).

For the longitudinal potential sigma, eliminate only W0
from the complete physical action:

    L=a³[q(sigma'-W0)²/(2N)+m²a_m W0²/(2N)
          -N m²b_m q sigma²/2],
    W0=q sigma'/(q+m²a_m).

This gives

    g_L²=a³m²a_m q/[N(q+m²a_m)],
    Omega_bare,L²=N²[b_m q/a_m+m²b_m].

At the old clock N=a_m=b_m=1, set
delta log(a)=zeta, delta N=n, delta a_m=alpha*n,
delta b_m=beta*n, z=q/(q+m²), omega²=q+m².
The actual alpha=4/(9h),beta=28/(81h) are inherited.
Exact directional differentiation gives

    r_T=delta log(g_T)=(zeta-n)/2,
    r_L=[(1+2z)zeta+(alpha*z-1)n]/2,
    delta Omega_bare,T²=-2q*zeta+(2omega²+m²beta)n,
    delta Omega_bare,L²=-2q*zeta+
                         [(2+beta)omega²-alpha*q]n.

All background time dependence of z,alpha,beta must remain
when differentiating r. The old source translation has zero
first variation on the clock, so it does not add a connected
Gaussian kernel term at this order.

## Derivative-free physical canonical form

For either physical mode w with conjugate momentum p_w,
let g² be its kinetic coefficient. Before the oscillator
boundary, the Hamiltonian and canonical matrix are

    H=p_w²/(2g²)+g² Omega_bare² w²/2,
    M=[0,g^-2;-g² Omega_bare²,0].

Its first variation is

    delta M=[0,-2r/g²;
             -g²(delta Omega_bare²+2r Omega_bare²),0].

This form contains n,zeta but no time derivatives of them.
For the exact fundamental matrix Phi, the first variation is

    delta y(u)=Phi(u,u0)delta y(u0)
                +integral_u0^u Phi(u,s)delta M(s)y(s) ds.

It is causal and includes the initial-state term. No momentum
integration, UV subtraction or state choice has been hidden
inside this finite-mode identity.

## Oscillator form and the moving map

The normalized oscillator uses v=g*w and
v'=d*g*w+p_w/g, d=g'/g. Thus

    X=C y, C=[g,0;d*g,1/g], det(C)=1.
    delta C*C^-1=[r,0;r'+2d*r,-r].

The exact identity
delta M_X=S'+S M_X-M_X S+C delta M C^-1
gives the oscillator frequency variation

    delta Omega²=delta Omega_bare²-r''-2d*r'.

The code checks the complete symplectic map and this identity.
It prevents artificial derivative loss from being confused
with the derivative-free physical phase form.

With W(f,f*)=i the real retarded scalar kernel is

    G(u,s)=[f(s)f*(u)-f*(s)f(u)]/i, u>=s.

It has G(s,s)=0,partial_u G(s,s)=1 and the positive flat
limit sin[omega(u-s)]/omega. The oscillator variation
contains its homogeneous initial-data solution and
-integral G delta Omega² f. Dropping the homogeneous term
would choose a state prescription that still needs justification.
