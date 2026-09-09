# An actual normalized continuum mode and its squeezed Hadamard state

## Mass-shell coordinates and the precise null profile

Set c=1 and first fix the affine sampling length to one.
For q>0 and transverse p in R^2 define

    omega=(|p|^2+q^2)/(2q),
    k_z=(|p|^2-q^2)/(2q).

Then omega^2=k_z^2+|p|^2, omega-k_z=q, and the positive
mass-shell measure is exactly

    dmu=d^3k/[(2pi)^3 2omega]
         =dq d^2p/(16pi^3 q).                            (1)

These relations and the absolute Jacobian are checked natively.

Choose epsilon=1/1000 and width=1/10000. Fix any even nonnegative
C-infinity bump b supported in (-width,width), with integral b=1.
For example normalize exp[-1/(1-(z/width)^2)] inside that interval
and set it to zero outside. Put

    C1=4-i/epsilon, C2=-2+i/epsilon,
    F(q)=C1 b(q-epsilon)+C2 b(q-2epsilon).

The supports are disjoint and strictly away from q=0. If

    B2=integral b(z)^2 dz, Bq2=integral z^2 b(z)^2 dz,

then the ACTUAL frequency moments are

    I0=integral |F|^2=(20+2/epsilon^2)B2,
    Iq=integral q|F|^2=(24epsilon+3/epsilon)B2,
    I2=integral q^2|F|^2
       =(32epsilon^2+5)B2+(20+2/epsilon^2)Bq2.             (2)

All are positive finite, and 0<=Bq2<=width^2 B2. No delta functions
or discontinuous momentum cutoffs are the physical state data.

Fix a nonnegative nonzero smooth transverse bump chi supported
in the annulus 1<|p|<2. Define

    J1=integral chi(p) d^2p, J2=integral chi(p)^2 d^2p,
    Jp=integral |p|^2 chi(p)^2 d^2p,
    C=Iq J2/(16pi^3)>0.

For every finite R>=1 use the one-particle mode

    h_R(q,p)=q F(q) chi(p/R)/(R sqrt(C)).                  (3)

Substituting (1) gives integral |h_R|^2 dmu=1 EXACTLY.
Its support in ordinary three-momentum is compact and smooth
for each finite R, because q is bounded away from zero.
No uniform compact momentum support is claimed as R grows.

The positive-frequency solution

    u_R(X)=sqrt(hbar) integral h_R(k) exp(-i k.X) dmu(k)

is therefore a genuine smooth solution with this normalized
one-particle mode, not a prescribed solution germ.
On gamma(s)=(s,0,0,s), (1),(3) give

    u_R(gamma(s))=A_R v(s),
    A_R=sqrt(hbar) R gamma0,
    gamma0=J1/(16pi^3 sqrt(C))>0,
    v(s)=B(s)[C1 exp(-i epsilon s)+C2 exp(-2i epsilon s)],
    B(s)=integral b(z) exp(-i z s) dz.

The transverse integral is R^2 J1 whereas normalization is
proportional to R, leaving the unbounded amplitude R.
Its phase along the null line depends only on q, not p.
This exact mass-shell degeneracy is the mechanism.

## The full quantum state

Decompose the one-particle Hilbert space into span(h_R) and
its orthogonal complement. Keep the vacuum on the complement
and put the chosen oscillator in the pure squeezed vacuum with

    n=<a_dagger a>=1/3, m=<aa>=-2/3.

The quadrature covariance is diag(1/6,3/2), positive with
determinant 1/4. Equivalently n(n+1)=|m|^2. A finite one-mode
Bogoliubov squeeze is implemented by a unitary on Fock space.
This constructs a normalized positive pure quasifree vector state
with the same CCR, not independent assignments of stress and
field strength.

Its vacuum-subtracted covariance is exactly

    Delta W(X,Y)
       =n[u_R(X)conj(u_R(Y))+conj(u_R(X))u_R(Y)]
          +m[u_R(X)u_R(Y)+conj(u_R(X))conj(u_R(Y))].       (4)

Every derivative of u_R exists by its compact smooth momentum
data. Thus (4) is a smooth bisolution; adding it to the vacuum
preserves the Hadamard property. For all null-line calculations
below, first form the SMOOTH renormalized stress and Wick-square
differences in spacetime, then restrict those smooth functions.
There is no claim that the unsmeared vacuum two-point distribution
itself has a pullback to this null line.

## Finite total energy and its unavoidable growth

The actual total energy above the vacuum is

    Energy_R=hbar*n integral omega |h_R|^2 dmu
       =hbar [R^2 Jp I0+J2 I2]/(6 Iq J2).                (5)

The full spatial integral of the nonminimal improvement is a
boundary term; it vanishes because these fixed-time mode
profiles are Schwartz in space. Thus (5) is also the physical
total energy for every xi in this family. It is positive and
finite at every finite R. Since 1<=Jp/J2<=4, its exact bounds are

    hbar R^2 I0/(6Iq) <=Energy_R
       <=hbar [2R^2 I0/(3Iq)+I2/(6Iq)].

Using (2), Bq2/B2<=width^2 gives the explicit rational energy
coefficients recorded by the certificate. In particular the
family does NOT have a uniform energy or UV-resource bound.

## Complete-line ANEC is positive for this family

F is smooth compact with strictly positive frequency support.
Hence v and all its derivatives are Schwartz in s, and

    integral v'(s)^2 ds=0,
    integral |v'(s)|^2 ds=2pi I2.

The first identity follows from the convolution at zero:
q+q' cannot vanish on this strictly positive support. The
second is Parseval. The improvement integrates to
-xi[w'(infinity)-w'(-infinity)]=0. In (4) the anomalous term
therefore vanishes after complete-line averaging, leaving

    integral <:T_ell,ell:>_R ds
       =2n A_R^2 integral |v'|^2
       =(4pi/3) A_R^2 I2 >0.

All integrals here are absolutely convergent for each finite R.
Taking R to infinity is NOT exchanged with a distributional
limit state. Positive energy elsewhere on the line compensates
the uniformly negative finite segment.

For a prescribed affine length tau, dilate momenta by 1/tau:
h_{R,tau}(k)=tau*h_R(tau*k) retains unit one-particle norm.
Then u_{R,tau}(gamma(t))=tau^-1 A_R v(t/tau);
the Wick square and null stress acquire tau^-2 and tau^-4,
respectively. Thus the same obstruction holds on every fixed
nonzero finite sampling length. No physical Planck or momentum
cutoff is held fixed under R->infinity.
