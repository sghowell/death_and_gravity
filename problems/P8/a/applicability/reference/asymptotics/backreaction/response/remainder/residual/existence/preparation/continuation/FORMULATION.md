# P8(a) A.13: the causal Einstein--logarithm block and continuation boundary

This checkpoint is an **exact constant-coefficient linear calculation and
an exact preconditioner decomposition of the full A.11 map**. It is not
the Frechet derivative of the rolling prepared state, a longer SEE solution,
a proof of actual nonlinear runaway, or a prescription to change the model.

Keep the original scalar, ordinary radiation normalization, transported
state, lambda=2 sqrt(2) A eta_star^2, gamma=0, epsilon=1 and delta=10^-14.
No order reduction, state reset, future boundary condition or finite
counterterm is introduced. The dimensionless actual conformal coordinate
and physical normalization remain those of A.11: x=eta/eta_star-1 and
T0=A eta_star^2, with kappa hbar=2880 pi^2 delta T0^2.

## Exact decomposition, not a replacement of the unknown history

The A.11 actual Wick functional is

    S[u,a]=R[u]+4 pi^2 D_gammaE[u]+d(a)u,
    d(a)=-19/60-log(a/2)/2.

The trace in the unforced region is exactly

    S''-2hS'+2(u+h^2)S
       =a^2 u/(60delta)-u^2/4-h^2(u+h^2)/30.

Choose a fixed scale af>0 for a preconditioner. This is a mathematical
choice, not a frozen physical metric asserted to solve SEE. The named
linear block on a zero-past variation X=u'-ub' is

    Bf X=Lf X-(c/2) I^2 X,
    Lf=4 pi^2 D_gammaE+d(af) Id,
    beta=gamma_E-19/30-log(af/2), c=af^2/(30delta).

The full auxiliary trace map has exactly `Bf X=Gf[X]`, where Gf retains
the baseline defect, rolling Einstein coefficient, anomaly history,
auxiliary q product, nonlinear actual-state response and local-coefficient
variation. [notes/decomposition.md](notes/decomposition.md) lists them all.
Neither their size nor their growing-pole projection is inferred here.

## Exact causal resolvent and the retained poles

With the principal logarithm and Re(s) larger than every pole,

    Rf(s)=2/[log(s)+beta-c/s^2].

For every real beta and c>0 there are exactly three poles on this sheet:

    w_j=W_j(2c exp(2beta)), p_j=exp(w_j/2-beta), j=0,1,-1.

The positive pole p0 is simple; the conjugate pair lies strictly in the
left half-plane. The residues are 2p_j/(1+w_j), and the cut density is

    2/{[log(r)+beta-c/r^2]^2+pi^2}, r>0.

All poles and the cut belong to the causal inverse. A positive Volterra
series proves that its full time kernel is nonnegative, even though the
damped pole pair alone need not be. Its C0 norm is the integrated kernel.
Smooth flat input has the same norm supremum as C0 input; the result is
not an artifact of incompatible nonzero initial forcing.

## Exact rational calibration and the limited conclusion

At af=5/2 and delta=10^-14, with no numerical evaluation of Euler's
constant or Lambert W used as a premise,

    10^6<p0<2*10^6, W0<33,
    ||Rf||_[0,L] >= [exp(10^6 L)-1]/17-8/3.

In particular the norm exceeds 1000 at L=10^-5 and 10^28 at L=10^-4.
Conversely, in the exponential weight exp(-sigma t), sigma=2*10^6,
the half-line convolution norm is at most 240/769. Returning from that
weight loses exp(sigma L); it is not an unweighted continuation bound.

These L values are durations of the **linear comparator**, not extensions
of A.11's proved interval. They show why resumming the stiff term cannot
by itself certify a uniformly small macroscopic causal inverse. They do
not impose a necessary step length on the nonlinear equation.

For compact forcing g, the growing term is absent exactly when
`integral exp(-p0 t) g(t) dt=0`, with any initial/history numerator included
in g. No such condition has been proved for the actual remaining Gf[X].
Deleting the positive pole is not the same causal inverse; imposing a
future cancellation or order reducing would be a separate physical/model
choice and is not done. The concrete next estimates and open constraint/
state issues are in [notes/next-step.md](notes/next-step.md). P8 remains open.
