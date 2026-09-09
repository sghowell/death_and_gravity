# Actual leading short-sampling coefficients

At the comoving bounce a=1, kappa_c=243/20,
c0=sqrt(1199/1215), ell=1/10, r0=1/(1215*c0).
The action-normalized mode derivation gives the leading covariance

    W_O(delta_u,k)=hbar/(2*kappa*k) *
                     [r0*exp(-i*c0*k*delta_u)+exp(-i*k*delta_u)].

It is the LEADING homogeneous singular part of the actual global
state, not a new stationary replacement for its finite-frequency
evolution. For a generic positive c and relative weight r, angular
integration with d^3k/(2*pi)^3 gives the radial kernels

    K_field = hbar/(4*pi^2*kappa) *
       integral_0^infinity k [r*exp(-i*c*k*d)+exp(-i*k*d)] dk,
    K_rho = hbar/(8*pi^2*kappa) *
       integral_0^infinity k^3
          [r*(1+c^2)*exp(-i*c*k*d)+2*exp(-i*k*d)] dk.

The second line retains both the actual clock time frequency c*k
and its physical spatial frequency k. They cannot both be replaced
by the free matter frequency.

Define g_eta(u)=eta^-1/2*g(u/eta), real g smooth compactly supported.
The exact reference functional is the one in notes/inequalities.md.
Rescale u=eta*s, v=eta*t, radial k=K/eta and auxiliary Fourier
frequency alpha=A/eta. The leading scalar and derivative kernels
give eta^-2 and eta^-4 respectively. The phase tends to
-K*c0*(s-t), or -K*(s-t), and the principal amplitudes tend to the
displayed bounce values.

Here this passage to the limit is controlled by the actual
all-orders symbol construction, not by assuming a stationary
state. On a fixed compact s,t support and for small eta, all times
remain in the initial gamma slab. The mode frequencies are bounded
below by a positive constant. In the reference Fourier direction
the time phase derivative has magnitude at least A+c_min*K.
Repeated integration by parts against the smooth sampling functions
gives arbitrarily many inverse powers of 1+A+K. The coefficient
symbols, their time derivatives and their rescaled K derivatives
have the required polynomial bounds from the complete asymptotic
mode construction. Choose more integrations than that polynomial
degree plus the radial and A integration dimensions. The resulting
majorant is integrable uniformly away from K=0.

Near K=0 the principal radial powers are K or K^3, hence locally
integrable. Split first at K=delta and at the rescaled fixed
low-frequency cutoff eta*K0. The genuine compact-k part of the
state is smooth and contributes a bounded functional before the
eta^2 or eta^4 normalization. Intermediate frequencies are bounded
by the same classical symbol orders; their normalized contribution
is bounded by a constant times delta^2 or delta^4, with terms
vanishing with eta. Taking eta to zero and then delta to zero
completes the dominated-limit argument. Lower symbol orders and
the smoothing difference are negligible in the normalized limit.
Time variation of the coefficients also vanishes in this scaling.
Thus the limits depend on the leading actual coefficients only,
not on a chosen infrared completion or Borel smoothing choice.

With Fourier convention g_hat(theta)=integral g(u) exp(-i*theta*u) du,
the reference functional integrates |g_hat(alpha+c*k)|^2.
Swapping the two positive integrals uses exactly

    integral_0^(theta/c) k^p dk
       = theta^(p+1)/[(p+1)*c^(p+1)], p=1 or 3.

For real g, integral_0^infinity theta^(2j)|g_hat|^2 dtheta
=pi*integral |g^(j)|^2. The native moments and prefactors yield

    C_field=(1+r/c^2)/(8*pi^2),
    C_time=(1+r/c^2)/(16*pi^2),
    C_spatial=(1+r/c^4)/(16*pi^2),
    C_rho=[1+r*(1+c^2)/(2*c^4)]/(16*pi^2).

Each coefficient multiplies hbar/kappa. C_field multiplies
integral |g'|^2 and the others integral |g''|^2. Substituting the
ACTUAL c0,r0 gives the report's exact algebraic values.
The test-energy enhancement is

    1+10863*sqrt(17985)/1723683599,

strictly between 1.00084 and 1.00085. Decoupling r=0 reproduces
the scalar reference-functional controls 1/(8*pi^2) and 1/(16*pi^2).

No sharpness or attained optimal state is claimed. The remainder
is little-o in this sampling limit, not a computed error at a
fixed cosmological sampling width. One cannot discard it at a
specified finite scale on the evidence in this checkpoint.
