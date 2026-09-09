# Actual spatial reduction and the coupled scalar phase blocks

All momenta below are canonical densities per unit hat
volume on the chosen flat central slice. Use q=|k_hat|^2.
The physical squared momentum there is q/N, a positive
finite rescaling. It cannot change a determinant sign.

For a scalar Fourier perturbation, the exact linear
momentum constraint is solved by

Z=-grad(Delta^-1)(p-3P*chi)/8.

This is the inherited full matter-sourced York solution,
now evaluated at the actual homogeneous momentum P(N).
The vector contribution to the spatial current is
quadratic about W_i=Pi_i=0; it is retained in the
derivation and therefore contributes no omitted linear
source. No zero momentum is inverted.

Writing the unit vector along k as e, the induced shear
momentum is -(p-3P*chi)*(ee^T-I/3)/4 and its squared norm
is (p-3P*chi)^2/24. The trace invariant is dp=p/3.
The two gravitational momentum-square contributions
in the actual Hamiltonian thus combine as

-N^(3/2)*(p-3*delta*j)^2/12
 +N^(3/2)*(p-3P*chi)^2/12.

The bare p^2 cancels exactly. The surviving p term is

(N^(3/2)/2)*p*(delta*j-P*chi).

The cancellation does NOT remove the p*j coupling.
The coefficient is checked by a direct full-geometry
two-leg reconstruction with the spatial constraints
imposed, not only by this restricted scalar formula.

Use W_i=partial_i sigma and P_sigma=-div(Pi).
Then Pi_i(k)=i*k_i*P_sigma(k)/q, so
integral Pi_i dot(W_i)=integral P_sigma dot(sigma)
and j=-P_sigma. The exact remaining vector terms are

H_vector=[kappa/sqrt(N)+10^6*sqrt(N)/q]*P_sigma^2/2
         +N^(3/2)*q*sigma^2/(2*gamma_s).

In particular kappa includes the joint trace/temporal
correction. Using just gamma_t would change the answer.

## Lapse variation and gamma canonical map

Construct the complete physical quadratic density from
the frozen central invariant Hamiltonian and actual
metric determinant, inverse, curvature, matter density
and vector terms. Differentiate at fixed canonical P
and momenta before imposing P^2=P(N)^2. The stationary
correction is -F_1^2/(2*h_N), where F_1 is the linear
lapse force. At this slice its coefficients of p and
P_sigma vanish. It therefore neither produces a p^2
term nor cancels the p*P_sigma coupling.

Make the actual linear gamma swap

v=P_b/(2q), p=-2q*b.

The old symplectic form is mapped exactly to the new
one for Q=(b,chi,sigma), Pi=(P_b,P_chi,P_sigma).
The transformation has no Hubble or Theta denominator.
For time-dependent q, its canonical generator is
retained in the Euler-order accounting below.

Writing H2=Pi^T A Pi/2+Pi^T B Q+Q^T C Q/2, the full
six-by-six phase calculation reconstructs all 21 independent
entries. At the central datum,

B=q*d*e_sigma*e_b^T, d=N^(3/2)*delta=(1-N^2)/sqrt(N),
C_bb=C_bsigma=C_chisigma=0,
C_chichi=N^(3/2)*q+O(1),
C_sigmasigma=N^(3/2)*q/gamma_s,
C_bchi=O(q).

The Hamiltonian momentum Hessian has the finite limit

A_infinity=[[r,-rP,0],[-rP,N^(-1/2)+rP^2,0],
            [0,0,kappa/sqrt(N)]],
r=-1/(4*N^3*h_N)>0.

The light two-by-two determinant is r/sqrt(N)>0 and
the vector pivot is positive. Thus A_infinity is
positive definite, and so is A(q) for sufficiently
large q at each fixed datum. Its exact vector entry
is kappa/sqrt(N)+10^6*sqrt(N)/q.

Subtraction of a moving homogeneous gravitational
momentum can add a finite coefficient times v^2.
Under the gamma map this changes A_bb only by O(q^-2).
The instantaneous gamma generator vanishes when the
hat metric has zero velocity, and its time derivative
has at most the lower orders recorded below. Neither
effect changes the displayed limits or the q^2 term.

The new gamma coordinate uses the full canonical p
directly. The older mixed-boundary gamma coordinate
differs by an O(chi/q) canonical change, with its
paired matter momentum shift. This changes no leading
coefficient in the high-frequency comparison.
