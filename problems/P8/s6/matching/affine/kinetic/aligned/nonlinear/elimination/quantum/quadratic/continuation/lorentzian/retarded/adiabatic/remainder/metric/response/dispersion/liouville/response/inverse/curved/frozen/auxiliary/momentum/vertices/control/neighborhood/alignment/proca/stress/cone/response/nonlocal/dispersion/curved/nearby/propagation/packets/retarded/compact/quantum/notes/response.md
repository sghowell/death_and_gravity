# State-independent Kubo response and compact commutator

By canonical.md and state.md, the exact Fourier commutator of the
linear relational field is
[Ohat(t,k),Ohat(s,k')]=-i*hbar/kappa U_{chi,pi_chi}(t,s,k)
times the Fourier delta enforcing k+k'=0.
Indeed e_2^T U Omega e_2=-U_{2,4}; this sign is checked natively.
S6.91's normalized classical retarded response for t>s is
Ghat(t,s,k)=V_s U_{chi,pi_chi}(t,s,k).

The probe action is +kappa integral V J O, so its interaction
Hamiltonian is -kappa integral V J O. Linear response is therefore
(i*kappa/hbar) theta(t-s) V_s [O(t),O(s)] = G(t,s).
The bracket is central, making the response independent of the
chosen state in this exact quadratic algebra. For an unscaled
action probe +integral V j O, retain the factor 1/kappa instead.
These two source conventions are checked separately.

The same result follows without an interacting Dyson expansion:
solve the linearly forced canonical Heisenberg equations by
variation of constants. The source adds the c-number displacement
integral U(t,s)V_s e_4 J. Commutators are unchanged, and its
observable component is the exact classical response. For compact
sources the associated Cauchy test is Schwartz, so the corresponding
linear-field and Weyl operations are available in the construction.
No UV interaction Hamiltonian or global quadratic Fock implementer
is assumed in this reasoning.

Now take the actual smooth real compact J and f of S6.92. Their
supports are strictly time-ordered and every source/detector point
pair is spacelike for the physical matter metric. Multiplication
by the smooth positive V_s preserves source support.
Define O(f)=integral dt d^3x f(t,x)O(t,x), and similarly O(VJ).
All these smearings are valid by growth.md, not point fields.

S6.92 proves C=integral f G J>=3D/4>0 for its fixed probes and
explicit D. Hence
[O(f),O(VJ)]=-i*hbar*C/kappa times the identity,
whose central coefficient has magnitude at least
3*hbar*D/(4*kappa)>0.
This is an exact statement within the reduced quadratic scalar
CCR theory on the actual nearby local bounce.

If another state of this same algebra has a suitable Hadamard
property, the commutator is unchanged. This conditional observation
does not establish existence of that state, its renormalized
stress, or a self-consistent semiclassical nearby background.
It also does not bound the difference between this quadratic
response and any interacting candidate's response.
