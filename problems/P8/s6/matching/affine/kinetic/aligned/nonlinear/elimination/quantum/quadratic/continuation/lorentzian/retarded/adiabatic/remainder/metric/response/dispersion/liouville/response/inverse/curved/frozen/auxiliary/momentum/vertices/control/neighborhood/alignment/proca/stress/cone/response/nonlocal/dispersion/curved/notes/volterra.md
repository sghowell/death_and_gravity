# Full-interval causal inverse and the original physical variables

Use the actual normal form in matching.md and write z=(eta,w),

    [P+V]z=I4 g_adapted,
    P=diag(A(u),gamma F_m(partial_u^2)), A=-6delta^2.

Here gamma=1/(64*pi^2*L^2)>0 is fixed, not adjusted in response to
a source. The background functions in Q use its stated normalization.
Every remaining term is the actual tree, selected-state response
or already-fixed profile, not an unspecified replacement operator.

## Invert the two diagonal channels, not the rank-one loop matrix

S6.86 gives the scalar causal inverse K_m in L1(0,infinity),
with no instantaneous term. The missing channel is supplied by
the actual tree coefficient A. Thus

    B=P^-1=diag(A^-1 multiplication, gamma^-1 K_m convolution).

The first entry is instantaneous and nonzero; the second is purely
integral. The isolated loop matrix itself remains singular.
Let a0=sup|A^-1|<=15625/6144 and
M0=max(a0,gamma^-1 ||K_m||_L1), both finite.

The remainder kernel is bounded by v(r)=C(1+|log r|) on 0<r<=1,
with finite C, and the same is true for all fixed diagonal
derivatives with their own constants. Then BV is Volterra with
the integrable scalar majorant

    w0(r)=a0*v(r)+gamma^-1 (|K_m|*v)(r).

The first term retains the nonzero instantaneous channel of B;
there is no instantaneous term in BV because V itself is integral.
Young's integral inequality gives w0 in L1(0,1).

For the norm sup exp[-lambda(u-u_initial)]|z(u)|, the operator
norm of BV is at most integral_0^1 exp(-lambda r)w0(r)dr.
Dominated convergence makes this tend to zero as lambda grows.
Choose a finite lambda so it is below 1/2. The Neumann series
for I+BV converges on the full fixed interval, giving a unique
continuous prepared z with

    ||z||_C0 <= 2 exp(lambda) M0 ||I4 g_adapted||_C0.

Large L enlarges constants but never changes integrability.
There is no assertion that the quantum correction is perturbatively
small, that lambda is numerically reasonable, or that the inverse
is stable in a spectral sense.

## No derivative loss in returning to the physical sources

The adapted change uses n=eta', so a C0 estimate for z alone
would not bound the original lapse. Here the first row has more
structure than the weak-log estimate needed for the other row.

After separating A eta'''', the entire first raw row is local of
order at most three. Its four-primitive kernel V_eta(t,s) and
its first output derivative are bounded, including the local
fourth-coefficient commutator. Hence the first integrated equation

    A(t)eta(t)+integral_initial^t V_eta(t,s)z(s)ds=I4 g_eta(t)

can be differentiated once without differentiating z(s). For
continuous z this gives a continuous eta' and

    ||eta'|| <= a0 [||I3 g_eta||+C1 ||z||],

where a finite admissible C1 is the sum of sup|A'|,
the diagonal norm of V_eta, and the supremum of the integral
of its first output-derivative norm. These quantities are finite
because the first row is a local-primitive kernel, not because
a logarithmic inverse arbitrarily smooths C0 inputs.

For raw prepared force g=(g_N,g_Z), put h_g=-3H g_N+H g_Z.
Integration by parts in the zero-past distribution sense yields

    I4 g_eta=-I3 g_N+I4 h_g,
    I3 g_eta=-I2 g_N+I3 h_g.

Thus no derivative of a continuous original forcing is required.
On the unit-length interval, |H|<=8/5 and the joint max norm gives

    ||I4 g_adapted|| <=13/30 ||g||,
    ||I3 g_eta|| <=47/30 ||g||.

The actual clock extrema, integration signs and these rational
bounds are checked exactly. Also v=w+H eta implies
||v||<=13/5 ||z||. Combining these inequalities with the finite
adapted inverse bound proves a finite C0 bound for the original
physical pair (n,v) against the original pair of forces. No
numerical bound for C1,lambda,||K_m|| or the final inverse norm
is claimed.

## Smoothness, preparation and recovery of the original equations

All fixed diagonal derivatives of V have integrable majorants.
For B, such derivatives act only on its smooth A^-1 multiplier;
they annihilate the lag-only scalar kernel. Translating both
endpoints in the causal kernel composition therefore preserves
integrable convolution majorants for every fixed derivative order.

For prepared inputs, differentiate a Volterra integral by moving
the derivative onto the input plus the diagonal derivative of
the kernel. Lower endpoint terms vanish. Difference quotients
with the integrable majorants justify the first derivative.
Induction and the same bounded inverse I+BV then give every
higher derivative. Smooth prepared g makes B I4 g_adapted smooth
and prepared, so z and the reconstructed n,v are smooth and
prepared. A finite extension and cutoff after the observation
interval can be used whenever compact future support is needed;
causality leaves the result on I unchanged.

Differentiating the integrated equation four times in the
zero-past distribution algebra recovers the exact transformed
retained response equations. The full new source response of
S6.85 is defined for these smooth inputs, so this is also its
original mode-response equation. Conversely, every smooth
prepared original solution integrates to the same Volterra
equation and hence is the unique solution.

With E_Z-g_Z=0, the remaining original lapse residual obeys
(partial_u+3H)(E_N-g_N)=0. Prepared zero data set it to zero.
The force transformation has therefore removed no original
equation and introduced no free homogeneous mode. The original
matter displacement follows by its single zero-charge integral.

Continuous forcing may be interpreted through the zero-past
distributional extension; the C0 reconstruction estimate above
does not assert that the forward logarithmic operator maps
arbitrary C0 sources to pointwise C0 currents.

Every operation is causal. Uniqueness on an initial subinterval
preserves preparation and excludes adding arbitrary higher-jet
solutions or independently changing the state. It does not
exclude a large transient or growing physical modes.
