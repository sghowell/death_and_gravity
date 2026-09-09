# Actual central homogeneous background time jets

All quantities use coordinate clock time u=phi. At the central slice
the hat spatial metric is the identity, N>0, canonical trace momentum
p0=0, vector data vanish, and matter density is the positive root of

    P^2=(5624995 N^4-14244998 N^2+8625003)/(500000 N).

It solves the actual lapse constraint. It is not substituted before
a fixed-phase lapse derivative. The physical metric has lapse N
and spatial metric N times the hat metric.

Let R be the homogeneous hat scale and ell=C/R^3 for the conserved
matter momentum C. Its canonical scale momentum is 3R^2 p0.
With the actual scalar Hamiltonian

    H=R^3 N[(p0-b)^2/(4a)-f0+ell^2/(2U)],

the Hamilton equations give

    R'/R=N(p0-b)/(6a),
    p0'=N[-(p0-b)^2/(4a)+f0+ell^2/(2U)],
    ell'=-3(R'/R)ell.

The code independently differentiates at fixed canonical scale
momentum and C, not fixed p0 during the scale derivative.

At u=0, a=-1/(3sqrt(N)), U=N^(3/2), b=0. The original
boundary primitive has

    I_phi=-6 N^(-3/2)-18sqrt(N)+24.

Differentiating the actual full b before evaluating gives

    b_u=3[-8 N^(3/2)+5N^2+3]/N^(3/2).

Keeping the boundary primitive is necessary. The actual on-family
momentum time derivative is

    p0'=[-6000000 N^(3/2)+1124999N^4+2875001]/
         [250000 N^(3/2)].

The hat Hubble derivative is N(p0'-b_u)/(6a). At N=1
these give p0'=-8 and H_hat'=4. The canonical background
coefficient -p0/2 must not be identified with H_hat away
from the original clock.

The exact scalar coefficients a,U,e^omega,B4,f0 are even in u,
where I_phi is even. The primitive integrand and b are odd.
The defining integral has fixed zero basepoint, preserving this
parity. The local homogeneous Hamiltonian equations and implicit
lapse, whose derivative is nonzero in the certified interval,
therefore have N,R even and p0 odd through the central data.
In particular N'=R'=ell'=0 there. Equivalently, differentiating
the lapse constraint has zero explicit time derivative and zero
linear trace derivative at this slice. These are actual background
jets, not independently assigned Hubble derivatives.

Let Fp=partial_N[N(p0-b)/(2a)] at fixed canonical phase data.
Its actual time derivative is

    Fp'=(partial_N[N/(2a)]) p0' - partial_N[N b_u/(2a)]
       =-9[1124999N^4-5000000N^2+2875001]/(1000000N).

The derivative of p0'(N) along the constraint family MUST NOT
be included in the partial lapse derivative. Doing so changes
the eventual clock gradient by
-sqrt(N)(5624995N^4-8625003)/250000. That nonzero discrepancy
is retained as a negative control.
