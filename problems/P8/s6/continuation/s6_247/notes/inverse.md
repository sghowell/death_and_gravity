# Entire causal quotient inverse and exact all-transfer bound

Use physical scalar quotient amplitudes Q=2zeta I-2c Pi with q=|P|² and Pi the longitudinal unit projector. The full covariant channel map is

B(D,q)=[[D²+2q/3,-D²/3],[-q,-D²]].

In the three amplitudes(n,zeta,b), the map is
M=[[q/3,D²+2q/3,-D/3],[q,-q,-D]].
The detector replaces D by-D. With source gauge chart
(n,zeta,b)=(D eta,w,D c+q eta),
both ordered maps reduce to a zero eta column and the SAME quotient B. The full three-source block remains gauge singular; it is not called invertible. These are source-pinned covariant projector identities, not COM substitutions at finite external momentum.

The total quotient reference is
A_q=B^T diag(F_trace(D²+q),(8/3)F_2(D²+q)) B,
using the two complete factors from notes/measure.md.

For p=lambda²+q, the exact inverse coordinate transform is

B^-1=[[1/p,-1/(3p)],
       [1/p-1/lambda²,-1/(3p)-2/(3lambda²)]].

Its causal kernel uses s_q(t)=sin(sqrt(q)t)/sqrt(q), continuously s_0=t:

R_q=[[s_q,-s_q/3],[s_q-t,-s_q/3-2t/3]].

R_q(0)=0. The full derivative's squared Frobenius norm is
[20 cos²(sqrt(q)t)-14cos(sqrt(q)t)+13]/9<=47/9<25/4.
Hence both its operator norm and that of its transpose are below C=5/2. Integration gives||R_q(t)||<Ct, uniformly in q, including its continuous zero-transfer extension. This does not solve the separate literal homogeneous lapse/shift constraint problem.

## Construct and check both causal products

Let J_q=diag(J_trace,q,(3/8)J_2,q). The complete two-mass static moments give

||J_q(t)||<=45/(ell+2)=Jmax.

Define the ORDINARY causal matrix kernel
E_q=R_q' * J_q * R_q^T.
The transpose is of matrix entries, not a reversal of time. Since R_q(0)=0, its Laplace transform is exactly

L[E_q]=B^-1 diag(1/F_trace(lambda²+q),(3/8)/F_2(lambda²+q))(B^T)^-1.

The two matrix products with the complete A_q transform are the identity. Its determinant is(8/3)F_trace F_2 lambda^4(lambda²+q)². For Re lambda>0 neither lambda²+q nor either full factor can vanish: a nonreal lambda has a nonzero imaginary part of lambda²+q, and positive real lambda gives positive p. No cut, pole or new initial branch prescription is crossed.

The forward factors and E_q are causal distributions with Laplace transforms in the right half-plane and polynomially bounded time growth. Their associative convolutions and transform uniqueness give BOTH
A_q*E_q=I delta, E_q*A_q=I delta.
This includes the complete initial boundary. No identity is restricted to strictly positive times.

## The actual triangle constants

Submultiplicativity gives

||E_q(t)|| <= C² Jmax int_0^t du int_0^(t-u) (t-u-v)dv
           =375 t³/[8(ell+2)].

Thus
int_0^T||E_q(t)||dt <=375 T^4/[32(ell+2)]
                   <125 T^4/4224<.03T^4.

For the first derivative, differentiate the last R factor. Its zero initial value removes the endpoint term without differentiating J:
E_q'=R_q'*J_q*(R_q')^T.
The corresponding bound is
int_0^T||E_q'(t)||dt <=375 T³/[8(ell+2)].

Both estimates are uniform over ALL external q>=0. For each real r, applying the Fourier multiplier and time Young/Minkowski inequalities gives the stated C([0,T];H^r) and L2([0,T];H^r) bounds for the quotient amplitudes and their dual sources, with no spatial derivative loss. No derivative in q, spatial Schwartz preservation or arbitrary higher-time bound is inferred from this estimate.

For a finite window, the graph consists of continuous zero-past quotient histories whose full forward causal distribution equals an ordinary dual source at and after the initial boundary. Sources may have a nonzero right-hand value at zero; initial atoms are not discarded. E_q maps ordinary sources into this graph and solves the equation. Conversely convolution by E_q on the full graph identity proves uniqueness. Extension past T has no effect before T. The same causal-distribution argument extends the stated L2 mapping, without claiming density of an unrelated curved graph.

## Physical normalization is not suppressed

The displayed reference is in64pi²-normalized units. For physical metric-current units the inverse restores64pi²; for the normalized gravitational force it additionally restores kappa. Therefore the physical-force estimate is

750 pi² kappa T^4/(ell+2),

not a small multiple of kappa^-1. The metric embedding has Gram
[[12,-4],[-4,4]], strictly between2I and14I. If both the dual-source and reconstructed-metric legs are converted with this embedding, their combined norm factor is at most14.

This flat proper-time convolution does not justify commuting a variable a(t)^-3 through an actual curved retarded kernel. The actual output-density normalization, gauge/clock contacts, classical tree, curved states, full profile and matter blocks remain to be restored in their original order.

Independent full two-mass spectral reconstructions verify both matrix inverse products at complex Laplace frequencies and multiple transfers, including the physical heavy scale. Separate oscillator checks retain the delta supplied by R_q'(0). These diagnostics test the identities; the contour, distributional and functional-space arguments above establish their stated scope.

A flat reference inverse with a finite bound does not establish physical stability or a controlled inverse for the actual coupled nonlinear system. Original P8 remains OPEN.
