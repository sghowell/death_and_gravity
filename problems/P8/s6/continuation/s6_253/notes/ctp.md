# Same eight-mode reference and complete finite Gaussian response

The free reference is the product of the S251 coupled two-scalar state,
its two TT states, the three fixed S55/S176 Proca states and the S240 H
state. No preparation is repeated when a background is varied. The
initial physical covariance has one4-dimensional block and six2-dimensional
blocks, reordered into eight canonical Q and eight canonical P. Under a
homogeneous physical probe, the evolving Hamiltonian has one full6-dimensional
scalar/longitudinal block and five2-dimensional blocks.

The exact Proca and H state transports are part of this statement. For
the oscillator covariance in variables(v,v_dot), define

    gL=a/(1+zeta*q), gT=a, gH=a^3, d=g_dot/(2g),
    (A,Pi)=(v/sqrt(g),sqrt(g)*(v_dot-d*v)).

This map is symplectic. Its full time connection converts the physical
Hamiltonian to diag(omega^2-d_dot-d^2,1). Both Proca corrections agree
with the entire existing source-pinned WKB operator, not a different
adiabatic state. The H correction is3H_dot/2+9H^2/4. A covariance already
expressed in(v,v_dot-d*v) receives only the scaling, not that shear twice.
The TT canonical variables remain those of the S251 state. The added
positive mass in its sampling form is a preparation device, not an
action mass or a modified physical vertex.

At a common finite regulator, let S(t) be the full canonical evolution
and V0 the same initial pure covariance. Then

    W(t,s)=S(t)[V0+iOmega/2]S(s)^T,
    O_A=z^T A z/2,
    mean_A=Tr[A V(t,t)]/2,
    connected_AB=Tr[A W(t,s) B W(t,s)^T]/2.

The last transpose is ordinary, not Hermitian conjugation. The real part
is symmetric noise. Its positivity follows from the norm of the centered
Hermitian quadratic operator applied to the fixed Gaussian vector.
With A*=S(t)^T A S(t), B*=S(s)^T B S(s),

    chi_AB=-i expectation([O_A(t),O_B(s)])
          =Tr[(A* Omega B*-B* Omega A*) V0]/2.

The effective-action retarded term is -theta(t-s)chi_AB. Its entire local
second-vertex term is -Tr[H_AB V(t,t)]/2. Thus the full physical second
matrix, including its held-W embedding contact when appropriate, cannot
be omitted. All fixed-state initial derivatives vanish in the stated
canonical representation; another chart transports the same density and
boundary phase instead of selecting a new state.

The complete S252 CTP trace applies in16 dimensions. After incorporating
the chosen initial pure-state transform, form the full relative symplectic
matrix R=S0^-1 Sminus^-1 Splus S0. For its Q/P blocks,

    alpha=(A+D+i(C-B))/2,
    Z=exp[-Tr Log(conjugate(alpha))/2].

The logarithm follows the actual path from the identity. Its continued
phase, including the metaplectic sign, is not an endpoint principal
square root. C-number action phases remain separate. The block product
at the finite regulator is a single common-regulator trace; it is not
an assertion that independently regularized continuum determinants
multiply without local terms.

The whole light/vector source-square first vertex can have zero product
reference mean and nonzero connected noise. The science tests use actual
mixed source-vertex entries to verify that distinction. Therefore a
conditional source-free Proca determinant plus an independent light
determinant does not supply the full perturbed Gaussian functional.

Both fixed stress profiles retain formal loop grade1. They are held fixed
under every probe, but their coefficient dependence on the physical
background must still be differentiated. Expand the entire inverse lapse
pivot, determinant and matching countervertices at one regulator. Keeping
the current pivot unexpanded is a convenient algebraic representation;
it is not convergence of the loop series at grade1 or an omitted-loop bound.
The source-free H determinant survives its high-order source filtration.

The unchanged product reference has the existing common time orientation.
Local finite-jet Wick insertions and mixed-mode contractions inherit the
S252 no-zero-wavefront-sum and same-state smoothing argument, now with a
finite mode index. Reconstructed shift observables retain their stated
infrared and symbol domains; no unrestricted zero-mode or inhomogeneous
extension is asserted here. The finite functional does not select local
curved subtraction coefficients, restore nonlinear Ward identities,
establish spacelike microcausality of a reduced kernel, or evaluate a
physical stress. None of these claims follows from reference normal
ordering or finite-cutoff Gaussian positivity.
