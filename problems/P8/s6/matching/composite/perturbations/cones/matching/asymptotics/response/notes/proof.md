# S6.11: exact finite-positive-parameter TT source response

## 1. Experiment and exact physical units

Use the constant model, positive roots, reconstructed internal canonical
scalar and physical composite frame specified in FORMULATION. In
physical time T the EH tensor kinetic coefficient for a unit real TT
polarization is `M² a_g³/(8N_g)` and similarly for f. Set `u=mT`,
`q=(kcom/m)²`, and divide the entire physical action by m. This common
constant division does not alter its equations, but must also divide
the probe action. All coefficients below refer to this divided action.

An external probe is defined by

`delta Sprobe = (1/2) integral dT d³x sqrt|g_eff| T^{ij} delta g_eff,ij`.

Set `T^{ij}=Pi(T)eij/Ae²`, `T^{0mu}=0`,
`k_i eij=0`, `eii=0`, `eij eij=1`. For the physical background FLRW
connection, the time component of its divergence is proportional to
the spatial trace, and its spatial divergence is proportional to
`k_i eij`. Both vanish for arbitrary C1 Pi(T). In particular no equation
requires Pi to be the anisotropic stress of the homogeneous internal
scalar. It is a distinct externally applied conserved TT probe.

Since `delta g_eff,ij=-Ae² Gamma eij`, its T-time coefficient is
`-Ae³ Pi Gamma/2`. Write its literal u-time coefficient as `J_u Gamma`.
Then `Pi=-2m J_u/Ae³`. The divided action uses `Jbar=J_u/m`, and we fix

`Jbar=M² epsilon sigma p/4`, `p=u²(1-u)²` on `[0,1]`, zero outside.

Thus `Pi=-M²m² epsilon sigma p/(2Ae³)`. This convention specifies the
polarized sign and all clock/volume factors. The pulse extended by zero
is C1. Its temporal Fourier transform is not compactly supported, and
no temporal-frequency cutoff is asserted. Working at linear order in
sigma, TT transversality/trace and vanishing probe density/momentum
leave the scalar and vector linear constraints unsourced. This follows
also from the rotational decomposition of the isotropic background.
No scalar/vector stability conclusion is needed or drawn.

## 2. Regular background coordinates and correct initial root

Let `e=1/y=epsilon*z`,

`Ae=(1+11epsilon*u²/4)²`, `j=11u/[(1+11epsilon*u²/4)z]`,

`R=sqrt(1+eta(1+e)³/3)`, `X=sqrt(1+e(R²-1))`.

The exact S6.10 flow becomes

`lambda=(1+e)[-X(R+j)+j e² R]/(X+eR)`,

`kappa=[X-je(1+e)]/[R+j(1+e)]`,

`z'=-z lambda`,

`R'=[(R²-1)lambda(2-e)/(1+e)-j(2-e+3eR²)]/(2R)`.

The initial conditions are `z(0)=1` and
`R(0)=sqrt(1+3(1+epsilon)³)`. The latter is not 2 at finite epsilon:
its initial derivative with respect to epsilon is `9/4`, and its
initial offset is bounded by `3epsilon` on the box below. Setting it
to 2 would change the pinned initial density to
`eta(0)=9/(1+epsilon)³`, with first derivative -27 at zero.

The scalar source and lapse expressions are

`(rho+p)/(M²m² e)=(2-e+3eR²)/(1+e)³`,

`N_g=e/(e+kappa)`, `N_f=kappa/(e+kappa)`.

For positive finite epsilon on the proved domain the lapses, scale
factors, radicands and canonical null source are positive. The full
background constraints and undivided pressure branch are inherited
by the exact source-aware reconstruction identities, not by assigning
these variables independently. Positive scalar speed permits the
analytic local potential reconstruction throughout this compact
interval. Neither density nor potential need be positive. This is not
the globally free old M1 matter action.

## 3. Literal two-metric action and flux variables

Put `g=h_g/sigma`, `f=h_f/(epsilon sigma)`. The composite tensor is

`Gamma=(h_g+y h_f)/(1+y)=epsilon sigma*(zg+f)/(1+e)`.

Divide the already divided action by `M²epsilon²sigma²/8`. The exact
quadratic action, including the external source, has integrand

`Ag g'² + Af f'² - Vg g² - Vf f² - W(g-epsilon f)²`

`+ 2p(zg+f)/(1+e)`,

where

`Ag=Ae³ z²(e+kappa)/(1+e)³`,

`Af=Ae³(e+kappa)/[kappa(1+e)³]`,

`Vg=Ae z² q/[(1+e)(e+kappa)]`,

`Vf=Ae kappa q/[(1+e)(e+kappa)]`,

`W=Ae³ z²(1-e)(1-kappa)/[(1+e)^4(e+kappa)]`.

These are the complete TT coefficients, including the internal shared
matter contribution to the relative mass. No derivatives of the
vanishing canonical normalization have been dropped: they are
contained in the exact derivative of the kinetic flux.

Define `Pg=Ag g'`, `Pf=Af f'`. The first four components of `model.matrix`
are exactly

`g'=Pg/Ag`,

`Pg'=-(Vg+W)g+epsilon Wf+pz/(1+e)`,

`f'=Pf/Af`,

`Pf'=epsilon Wg-(Vf+epsilon²W)f+p/(1+e)`.

All four initial data are zero. These are ordinary regular first-order
equations for every finite member in the proved box; no quotient by
the mass or a Hubble rate is used. Both TT polarizations obey the same
equations independently. One polarization is sufficient to falsify a
prediction subject to the stated source contract.

The locked candidate is an explicit competing action, not an asserted
exact integration of the relative field. Restricting `h_g=h_f` is
`g=epsilon L`, `f=L`. Its exact equation is

`(C_lock L')'+V_lock L=p`,

`C_lock=Af+epsilon²Ag`, `V_lock=Vf+epsilon²Vg`.

The source and output weights reduce to one. Its initial value and
flux also vanish. In particular it is incorrect to retain only Af at
finite epsilon. The six-dimensional block-diagonal flux system in
model.py evolves the full and locked candidates together.

## 4. The limit retains a physical source-coupled relative response

The frozen S6.10 forward orbit is

`z=12cosh u+2sinh u-11`, `v=z'/z`, `a=11/z`,

`R=(z'-11u)/z`, `v'=1+a-v²`, `a'=-av`.

It has `z>=1`, `1<=v²<=140/19<9`, and

`D_H=N_R-m_alg²=3(a+1)²/(4v²)+(-v²+4v-2)/4>=1/4`.

Here `N_R=(z/sqrt(v))''/(z/sqrt(v))` is the complete normalization
curvature. At epsilon zero,

`Ag=z²/v`, `Af=1`, `Vg=z²v q`, `Vf=q`, `W=z²(v-1)`.

With `chi=(z/sqrt(v))g`, the exact limiting equations are

`f''+q f=p`,

`chi''+[v²q-D_H]chi=sqrt(v)p`,

`R_0=f+sqrt(v)chi`, `L_0=f`.

The substitution uses the identity
`(F²(chi/F)')'=F chi''-F''chi`, `F=z/sqrt(v)`.
It explicitly retains the normalization derivatives. Both operators
have coefficients of order epsilon-zero, or order m² in physical
time, on this fixed window; there is no epsilon-diverging decoupling
gap. An order statement does not itself assert a cutoff-safe spectrum.

The causal initial-value kernels have zero initial value and unit
initial derivative at their source time. They are

`G_L(u,s)=sin(sqrt(q)(u-s))/sqrt(q)`

and `G_H` for `partial_u²-[D_H-v²q]`. The source-normalized output
kernel is exactly

`K_0(u,s)=G_L(u,s)+sqrt(v(u))G_H(u,s)sqrt(v(s))`.

Both source and observation factors are indispensable. This is a
retarded equation calculation. A homogeneous relative solution would
add another term if its initial data were nonzero. No causal inverse
is silently substituted into a symmetric single-copy effective action.

## 5. Positive pulse comparison and finite-error margin

On `q in [1/256,1/64]`,

`D_H-v²q >= 1/4-(140/19)/64 = 41/304 > 1/8`.

The positive Volterra iteration of the retarded equation, or its
variation-of-constants formula relative to `partial²-1/8`, gives

`G_H(u,s)>=sqrt(8)sinh((u-s)/sqrt(8))>=u-s`.

The iteration converges on a compact interval because its coefficient
is continuous and bounded there. Its kernels and inhomogeneous terms
are nonnegative, proving the comparison rather than assuming a
positive frozen eigenvalue. Meanwhile for `0<=u-s<=1`,

`(383/384)(u-s)<=G_L(u,s)<=u-s`.

The lower bound follows from `sin x>=x-x³/6`; all arguments are between
0 and 1/8. Since `sqrt(v)>=1`, the limiting relative contribution is
at least the linear kernel response. At observation u=1 let

`I=integral_0^1 (1-s)p(s) ds=1/60`.

Then `(383/384)I<=L_0<=I` and the relative term `H_0>=I`, with
`R_0=L_0+H_0`. In particular, the relative term is not suppressed by
the singular source normalization.

Suppose the full and locked normalized outputs each differ from their
limits by at most `delta=1/600=I/10`. Then

`R_epsilon-L_epsilon>=I-2delta=1/75`.

Both outputs are positive by their explicit lower bounds, and

`(2/3)R_epsilon-L_epsilon`

`>= (2/3)H_0-(1/3)L_0-(5/3)delta`

`>= I/3-(5/3)delta=I/6=1/360>0`.

Hence the relative error of the locked prediction is at least 1/3.
The next two sections prove the assumed finite errors quantitatively.

## 6. Outward-rational background continuation

Use the convex box

`0<=epsilon<=1/10000`, `0<=u<=1`,

`1/2<=z<=12`, `1/10<=R<=3`, `1/256<=q<=1/64`.

intervals.py carries exact Fraction intervals and first partial
derivatives with respect to `(epsilon,z,R)`. Each arithmetic operation
rounds endpoints outward to 32 binary places; sqrt uses integer
isqrt with a checked upper endpoint. Products and reciprocals enclose
their full operand intervals; no denominator crossing zero is allowed.
The derivative chain/product/inverse/sqrt rules therefore enclose
partial derivatives of the exact algebraic functions, not derivatives
of a rounded approximation. Floating/nonfinite endpoints and invalid
jet variable indices are rejected.

The resulting value/derivative records, replayed in bounds.build(),
show positive `X`, lapse numerator and denominator, `X+eR`, `kappa`,
`Ag`, `Af`, `C_lock`, and null source. They also show `e<1/100`.
Although `N_g=0` is possible at the formal epsilon-zero endpoint,
`N_g>0` follows for every actual member since `e=epsilon*z>0`.

The background sup-norm state Lipschitz constant is at most
`K_b=26372`; the parameter partial derivative is at most `B_b=350559`.
Include `|R_epsilon(0)-R_0(0)|<=3epsilon`. Until exit from the box,
Gronwall gives

`|b_epsilon-b_0|_infinity <= (3+B_b) exp(K_b) epsilon`

`<= 2^52763 epsilon`.

Here and below `exp(K)<=4^K=2^(2K)` for integer nonnegative K. No huge
integer denominator needs to be materialized to compare these bounds.

For the limiting orbit, z is increasing. Using `8/3<exp(1)<11/4`,

`z(1)=7exp(1)+5exp(-1)-11 < 81/8 <11`.

The lower exponential bound is its positive partial sum through 1/3!;
the upper is the sum through 1/5! plus a geometric tail bounded by
`7/4320`. Both are replayed rationally. The numerator `z'-11u` starts
at 2 and has derivative z>=1, hence `R>=2/11`. Also `R'=1-Rv`, v>=1
and R(0)=2 imply R<=2 by a barrier argument. Thus the limiting orbit
lies in `1<=z<=11`, `2/11<=R<=2`. A sup-norm error at most 1/128 lies
strictly inside the larger box, with positive rational margins in the
report. For `epsilon<=2^(-52770)` this closes the bootstrap and proves
existence on all of `[0,1]`. Analyticity on the positive-root tube and
the standard continuation theorem complete the argument; there is no
assumption of a finite-epsilon solution on this interval in advance.

## 7. Quantitative causal propagation and output transfer

Let x be the six flux variables `(g,Pg,f,Pf,L,P_L)`. Write
`x'=A_epsilon x+p b_epsilon`, x(0)=0. The same box enclosures give

`||A||_infinity<=K=33572`, `||b||_infinity<=B=12`,

`||partial_(epsilon,z,R) A||_row,1 <=D_A=68152776910250`,

`||partial_(epsilon,z,R) b||_infinity,1 <=D_b=145`.

The derivative norms sum all three absolute partial bounds in each
entry, then row sums for the matrix. Since the joint parameter/state
change has sup norm at most `2^52763 epsilon`, the mean-value theorem
on the convex box bounds coefficient differences by the displayed
derivative constants times that quantity. No q derivative is needed:
the comparison fixes the same q and is uniform throughout its band.

Using `0<=p<=1` gives `||x||<=B exp(K)<=2^67148`. Duhamel and Gronwall
for the difference, with the limiting x on the forcing side, yield

`||x_epsilon-x_0|| <= (D_A+D_b) 2^52763 2^67148 exp(K) epsilon`

`<=2^187101 epsilon`.

The full output row is `c=(z/(1+e),0,1/(1+e),0,0,0)`; the locked row
selects L. Both have row norm at most 13. Their summed joint partial
bound is at most 158. Thus the sum of the phase and output-map errors
is bounded for each output by `2^187106 epsilon`. Choosing

`0<epsilon<=2^(-187116)`

gives each error at most `2^-10<1/600`, uniformly in u and q. This gate
also lies inside the parameter box and closes the background tube.
The report computes every integer from the exact enclosure records
and checks each exponent inequality; these exponents are not a
rounded decimal representation of an integration result.

## 8. Prepared zero-source data: a positive control

Set p=0. At u=0 choose a common normalized metric value l0 and
u-velocity v0 with `|l0|,|v0|<=1`. The finite initial flux vector is

`(epsilon l0, epsilon Ag v0, l0, Af v0, l0, C_lock v0)`.

Thus `h_g=h_f=epsilon sigma l0` and their u-derivatives are equal.
The relative tensor and its derivative vanish, so the full relative
canonical field H and H' also vanish, without ignoring a rotating
normalization boundary. At epsilon zero the initial vector is
`(0,0,l0,v0,l0,v0)` and has sup norm at most one. It evolves with
`g=0`, `f=L`, so both limiting physical outputs coincide.

The initial finite-minus-limiting vector is bounded by
`3039038489 epsilon`: include epsilon Ag explicitly, and use the joint
partial bounds for Af and C_lock with the initial root offset 3epsilon.
The finite initial norm is at most 1450. The same Duhamel calculation
now has an initial-error term and no forcing term. bounds.py records

`||delta x|| <=2^187098 epsilon`,

each output error relative to its common limit `<=2^187103 epsilon`.

Consequently full-minus-locked is bounded by `2^187104 epsilon`, less
than `1/600` under the same gate. This proves an actual conditional
light-only positive control. It does not replace the physical probe
experiment with differently prepared heavy data, or establish
universal matching for arbitrary states or sources.

## 9. Uniformly small metric amplitude and precise conclusion

The zero-data sourced phase norm is at most `2^67148`; the prepared
source-free finite phase norm is at most `2^67155`. Since
`h_g=sigma*g`, `h_f=epsilon*sigma*f`, taking

`0<=sigma<=2^(-67162)`

makes each metric polarization amplitude at most `1/128<1/100`.
At sigma=0 the response is the trivial zero field; use positive sigma
for the normalized ratio. All reported linear operator ratios are
independent of its choice. This is a metric-amplitude statement, not
a uniform proper-time curvature or nonlinear remainder theorem.

For finite positive epsilon the background, TT equations, source and
constraints are all regular on the certified interval. The absolute
composite response and probe scale as epsilon*sigma and vanish in the
limit. Nevertheless the specified source-normalized locked prediction
fails by at least one third. No unit canonical column with an
unbounded individual relative metric was used to assert uniform
physical smallness.

The necessary inference is only failure of the explicitly stipulated
locked-action source-response contract on this family. A general
second-order theory with another matched source map, extra derivative
operators, memory, state restrictions or error contract has not been
excluded. Neither operator in the limit has an epsilon-parametrically
large mass, but that fact alone is not a theorem against all
nonadiabatic reductions. The prepared-data positive control makes this
distinction concrete.

There is no uniformity on an interval growing with Y or on a full CD
duration. The zero-epsilon endpoint is singular as a two-metric
geometry; individual lapse/scale factors shrink and some parent
proper-time curvature scales grow. No cutoff safety, full scalar or
vector health, global reconstructed-potential bound, nonlinear
instability, quantum state, original physical-frame M1 match, or P8
closure follows. The very small sufficient parameter is a proof of a
finite-member mathematical statement, not a viable EFT window.
