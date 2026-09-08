# Actual physical-current ultraviolet kernel

This proof is on the original clock, not a finite off-clock
ordinary-Proca replacement. Fix m>=1000, L>0 and the selected
S6.55 state. All unlabelled time derivatives are in physical u.
The auxiliary common background acoustic time has sigma'=1/a.
It is used only to represent the unchanged physical modes.

## Exact current kernel before any UV expansion

At a fixed background canonical phase coordinate, the physical
Hamiltonian current has the form

J_i=(A_i*p_u^2+B_i*omega_u^2*v_u^2)/2, i=N,Z.

Use precisely the S6.67 weights, not only a mass-insertion
vertex. The retarded kernel is i*theta times the current
commutator, plus the already-derived instantaneous second
Hamiltonian and output-normalization contacts. This is the
physical Kubo representation underlying S6.64 and S6.69.
The S6.71 source/readout identities check its equivalence to
moving-clock covariance transport, including both source
shifts; those shifts are not set equal or dropped here.

For background acoustic mode F_k and P_k=F_k'-b_k*F_k, the
positive connected two-current function at t,s is one half
the sum of

A_i(t)A_j(s) [P_k(t)P_k(s)^*]^2,
B_i(t)B_j(s) D_k(t)D_k(s) [F_k(t)F_k(s)^*]^2,
A_i(t)B_j(s) D_k(s) [P_k(t)F_k(s)^*]^2,
B_i(t)A_j(s) D_k(t) [F_k(t)P_k(s)^*]^2,

multiplied by f(t)f(s), f=1/a. There are two transverse
sectors and one longitudinal sector in physical D=3.
The output normalized current has the additional a(t)^-3;
the input is an action current and has no extra a(s)^-3.
The radial measure is k^2 dk/(2*pi^2).
The real retarded part is -2*theta*Im of this positive function.
In complex dimension replace conjugates by the prescribed
analytic minus modes before taking the physical real slice.

## Controlled finite expansion, not an assumed convergent series

On every compact time triangle, the S6.70 potentials and pumps
are smooth rational functions of k^2 with a>0,m>0. Thus they
have finite inverse-k expansions with uniform remainders
after any specified number of time derivatives.

For each fixed desired expansion and derivative order, choose
a sufficiently high finite Riccati frequency W_N. Its exact
normalized variation-of-constants equation bounds the error
by its residual divided by frequency, with a k-independent
finite-interval exponential. S6.55's all-order Cauchy matching
then makes the exact physical state differ from that reference
by an arbitrarily high inverse power, also after the prescribed
finite number of time derivatives. An overall phase cancels
in F_k(t)F_k(s)^*. The nonzero unvaried initial squeezing is
included in the error at sufficiently high N, not reset.

Expand the finite reference amplitude and its phase around
exp[-ik(sigma(t)-sigma(s))]. Taylor's integral remainder for
the residual phase gives bounded smooth coefficient functions
on the compact triangle. Multiplying the finite expansions
in the four current products yields, for k>=K>0,

retarded radial kernel =
sum_(j=-1)^4 [a_j(t,s) k^j sin(2k Delta_sigma)
             +c_j(t,s) k^j cos(2k Delta_sigma)]
 + remainder(t,s,k).

The remainder is O(k^-2) uniformly. For each desired number
of time derivatives, choose a higher finite N and also expand
through correspondingly more negative powers of k, so that
the differentiated remainder is integrable. The extra finite
negative-power terms have no higher diagonal singularity.
No derivatives in momentum of the cutoff-summed state are
needed in this argument. The finite low-k interval is smooth
and integrable, using m>0 and the original constrained
canonical equations at k=0.

In the common complex-D neighborhood, only the undifferentiated
kernel limit is needed. The chosen W4(D) plus fixed Borel
correction has relative Cauchy mixing O(k^-6) against the
corresponding dimensional high-order reference. After the
radial k^4 prefactor this is O(k^-2). The same expansion
therefore has an integrable remainder
O(k^(-2+abs(Re(D)-3))) in a sufficiently small strip.
The fixed-dimensional plus/minus construction and multiplicity
are retained. At D=3, arbitrary time differentiability follows
from the stronger all-order physical comparison above.

## The highest two powers and their normalization

Writing P_k=(-ik-b_k)F_k at leading order, the exact polynomial
in the generic current product has coefficients

k^4: (A_i-B_i)(A_j-B_j),
k^3: 2i*[b_s*(A_i-B_i)*A_j-b_t*A_i*(A_j-B_j)].

The certificate verifies these before specializing either
current. Terms involving U_t,U_s occur at lower degree;
the first phase correction contains an integral along the
lag and vanishes on the diagonal.

At physical D=3 the transverse highest pair is zero. The
longitudinal pair is

b_L=(alpha+beta,4)=4v, v=(16/(81h),1).

The normalized Q=64*pi^2 times the unscaled physical current
kernel therefore has leading radial term

8*b_L(t)b_L(s)^T*k^4*sin(2k Delta_sigma)/(a(t)^4*a(s)).

For positive lag the Abel radial integral is
integral k^4*sin(2k Delta_sigma)dk=3/(4*Delta_sigma^5).
Consequently the leading off-diagonal kernel is

6*b_L(t)b_L(s)^T/[a(t)^4*a(s)*Delta_sigma^5].

These are off-diagonal formulas; multiplication by the causal
step and extension at the diagonal still require the fixed
subtraction treated in the next note. No unrenormalized
divergent time integral is used as an ordinary function.

Let tau=t-s. The geometric ratio
G(t,s)=tau^5/[a(t)^4*a(s)*Delta_sigma^5]
is smooth up to the diagonal, because
Delta_sigma/tau=integral_0^1 1/a(s+x*tau) dx>0.
It has G(s,s)=1 and first lag derivative -3H(s)/2.
The highest diagonal singularity is therefore exactly

6*b_L(t)b_L(s)^T/tau^5,

up to one lower singular order. Four derivatives of
b_L(t)b_L(s)^T/(4*tau) have this same highest term;
every derivative of its smooth multiplier is lower order.
The rank-one pole is b_L b_L^T/8, exactly S6.69's P.
The logarithmic symbol is -2P*log(s/m), with no sign change.

The connection term is lower singular order but nonzero.
Its diagonal antisymmetric matrix, after dividing out b_k,
is [[0,8+32/(27h)],[-8-32/(27h),0]].
Discarding it would change the remainder, even though it
would not change the highest logarithmic coefficient.

For diagonal derivatives, differentiating the phase inserts
k*[sigma'(t)-sigma'(s)], with the bracket O(t-s). Thus each
extra radial power is accompanied by a lag factor. In the
explicit finite radial integrals these cancel in the singular
order. Expanding to a sufficiently deeper negative power for
each fixed derivative count leaves an absolutely integrable
remainder. This supplies the diagonal-derivative kernel class
used in the smooth Volterra argument without claiming that a
fixed O(k^-2) remainder stays integrable after arbitrarily
many time derivatives.
