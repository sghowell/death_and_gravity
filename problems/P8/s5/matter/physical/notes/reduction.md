# Matter-sourced physical reduction through degree four

## 1. Canonical variables and the density source

Start after S5.5's regular auxiliary-metric transformation and its local lapse
boundary primitive. Use the linear spatial metric gauge of FORMULATION.md.
Write `v=zeta`, `s=delta chi`. Before the final mixed matter boundary, define the
metric scalar momentum `p_old`, TT momentum `PiTT`, and matter density `pi_chi`.
At the evaluation point `a=1`, the inverse boundary transformation is

`p_old=p+3*l*s`, `pi_chi=l+P+3*l*v`.

The generating boundary is `F=3*a^3*l*v*s` in original variables. Its scalar
derivatives give both shifts. Its explicit cosmic-time derivative vanishes by
the exact free-matter equation `(a^3*l)'=0`. Here that equation is applied to
the physical unscaled l before taking fixed local evaluation units; it is not
an assertion that the scaled function `y^11/10` obeys the same differential
equation. Using just one shift is not canonical. The earlier metric boundary
is also retained, not replaced by this one.

The spatial constraint is a density equation. With the contravariant metric
momentum density written in the local spatial frame,

`C^i=partial_j Pi_bar^ij+Gamma^i_jk Pi_bar^jk-(pi_chi/2)*g^ij*partial_j s=0`.

The factor `pi_chi` on its right side is not divided by `sqrt(g)`. The density
weight cancels the trace-Christoffel term in the divergence on its left side.
Matter geometry instead enters the lapse Hamiltonian by
`eta=pi_chi/sqrt(g)-l` and `z=g^ij partial_i s partial_j s`.

## 2. All three constraints and canonical York orthogonality

Decompose the metric momentum as

`Pi_bar=[-H*(1+v)+p_old/6]*I+PiTT+H*gammaTT+L W`,

`(LW)^ij=partial_i W_j+partial_j W_i-(2/3)*delta_ij*div W`.

At homogeneous perturbative degree n the dependence on W_n is exactly
`E W_n`, with `E=Delta I+grad div/3`. All extra Christoffel factors have
positive degree. The remaining degree-n source is known from lower orders.
For nonzero transfer k,

`E(k)^-1=-(I-k*k^T/(4*k^2))/k^2`.

Thus `W_n=-E^-1*source_n` is unique in this Fourier chart. The implementation
recomputes and checks every component of the residual after the recursion.
The nonexceptional proper-subset condition supplies exactly the inverses needed
for n=1,2,3 in a four-leg integrated coefficient. It is not a prescription to
drop a sourced homogeneous mode.

The York tensor is tracefree. Its contraction with any TT metric variation is
zero after spatial integration by parts: `integral LW:delta gammaTT=0`.
Consequently its nonlinear dependence on the two scalar/matter pairs does not
alter their reduced symplectic form in this gauge. The prescribed remaining
trace and TT shifts are canonical boundaries. This proves the formal local
reduction without assuming a nonlinear global TT gauge theorem.

W_1 through W_3 suffice for the quartic Hamiltonian. An omitted W_4 first pairs
with the background derivative of H with respect to Pi_bar; isotropy makes that
derivative proportional to I, so its contraction with LW_4 vanishes. Every
other factor would raise total degree to at least five. The same observation
explains why a degree-(n-1) York recursion computes an integrated n-leg kernel.

## 3. Full invariant substitution and time boundaries

Build g inverse, determinant, Christoffels and R[g] directly. Define

`M=Pi_bar*g/sqrt(g)`, `sigma=tr(M)+3*H`,

`shear2=tr(M^2)-tr(M)^2/3`, `rho=R[g]`,

`eta=pi_chi/sqrt(g)-l`, `z=g^ij partial_i s partial_j s`.

Sigma, rho and eta start in degree one; shear2 and z start in degree two.
Their nonlinear parts are not discarded. Substitute all of them in the pinned
S5.5 invariant Hamiltonian

`h=A+B*sigma+C*rho+L*eta+D*sigma^2+E*shear2+M_m*eta^2+Z*z`.

The symbol M_m in this line is a lapse-dependent matter coefficient, not the
momentum matrix M. Lapse jets are derivatives, not Taylor coefficients.
The generic stationary formula is replayed against the pinned recurrence.
H through degree four requires only n1 and n2; the unused n3 need not be
expanded again as a physical Fourier jet. Since each substituted invariant
has at least its declared weight, omitted invariant degree five cannot affect
physical degree four.

The final local Hamiltonian density is

`sqrt(g)*h_reduced -2*H*Pi_bar:g`

` -(Hdot+3*H^2)*(6*v+3*v^2-gammaTT:gammaTT/2) -l*pi_chi`.

The second line contains the original metric canonical boundary and the matter
background-velocity subtraction. For chi=chi_bg+s, omitting `-l*pi_chi` leaves
nonzero tadpoles. The density definition is linear in the chosen coordinates,
so that particular subtraction introduces no cubic/quartic term; it must still
be included to identify the correct background and canonical perturbations.

The regular phase gamma chart is the fixed-mode canonical swap

`p=-2*q*b`, `v=P_b/(2*q)`.

Its time generator is `-H*b*P_b`, since the physical coefficient
`2*a^3*k_physical^2` has logarithmic derivative H. Apply this swap to the scalar
pair **after** the mixed matter boundary; the matter pair remains canonical.
Both phase charts require only nonzero spatial q, never inverse Theta or Lambda.

## 4. Coupled velocity substitution and quartic contact

For either canonical chart write the full quadratic Hamiltonian as

`H2=P^T*A*P/2+P^T*B*Q+Q^T*C*Q/2`, `Q=(gravitational scalar,s)`.

Then `P0=A^-1*(Qdot-B*Q)`. This is the full two-by-two matrix, not a diagonal
gravitational response with an independent matter spectator. The TT response
is `PiTT=gammaTT_dot/4`. The implementation substitutes P0 in every invariant,
source and canonical boundary before extracting the labelled coefficient.

Writing the momentum correction as `delta P=-A^-1*H3_P+O(Q^3)` gives

`L2=P0*Qdot-H2(P0)`, `L3=-H3(P0)`,

`L4=-H4(P0)+(1/2)*H3_P^T*A^-1*H3_P`.

The last expression includes the momentum contractions for both scalars and
both tensors. For four labelled external legs, the two ordered orientations
of each 2+2 split cancel the factor 1/2, leaving three unordered partitions.
Each scalar partition uses the full inverse Hessian at its own transfer.
For a rational orthogonal TT basis E, let its conjugate momentum coefficient
be defined by `PiTT=E*p_E/(E:E)`. The inverse TT momentum Hessian is then
`(E:E)/4`; summing both E contributions gives the complete tensor projector.
Rescaling or orthogonally changing that basis leaves the total invariant.

This is the stationary Legendre **contact** correction. It is not a propagator
or the exchange of an internal physical excitation over time. No on-shell or
finite-time scattering conclusion is drawn from these kernels alone.

## 5. Fixed-momentum phase regularity, not interaction control

The compact background and lapse jets are polynomials in x,y with rational
coefficients; `|x|,|y|<=1`. The geometric inverse/volume expansions truncate
exactly in the labelled algebra and introduce only rational constants. York
recursion introduces only fixed proper-subset `1/k^2` factors. The finite lapse
formula introduces powers of `1/J`, and pinned S5.5 gives `J>1/10`. The phase
gamma swap introduces only external `1/q`, also fixed and nonzero.

By induction over these finite additions/products, every fixed-momentum phase
kernel belongs to the algebra of compact-continuous bounded coefficient
functions generated by those inputs. Hence it extends continuously to r=+-1
and has a finite bound there. This is an analytic closure lemma, not a sampled
sign argument or a serialized expanded quartic majorant. Bounds can diverge
as a proper-subset momentum goes soft; no uniform soft/forward limit is claimed.
Here fixed momenta mean fixed dimensionless local momenta. A fixed-comoving
mode instead has q tending to zero in the tails, outside that fixed-transfer
bound; its global evolution is not covered.

Velocity responses have additional inverse determinants. They must remain on
the stated positive chart domains at every momentum. They do not inherit the
all-time phase regularity conclusion. Time-dependent configuration/mode
normalization introduces further derivative/connection terms and is not done
here. For that future step `q=k_com^2/d^3` in CD moving local units, so
`sqrt(d)*d_u q=-6*x*q`; the D-only drift cannot be reused.
