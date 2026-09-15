# Formulation and exact scope

## Source, variables and order

Keep the original affine parent EFT, its matched flat vacuum and field
normalizations. Let mu>0 and n>4mu be the light and heavy mass squares,
g the H Phi Phi cubic, kappa>0 the gravitational normalization, and
nu^2>0 the dimensional loop scale. The inherited source values are
unchanged; allowing symbolic positive masses proves identities on this
specified massive domain, not a different original source.

Use the inherited harmonic metric gauge, signature(+---), D=4+2e,
0<e<=1/8, the Feynman normal sheet and the raw master convention in
S6.283. A raw Euclidean C0 carries a minus sign. The inverse matter
propagator is p^2-mass+Sigma, so Z=1/(1+Sigma').

This packet computes the entire minimal proper cubic correction and all
three minimal-gravity external residue factors, deltaT/g at order1/kappa.
The graph boundary excludes, but explicitly retains as separate work,
the non-1PI H-metric response, finite cubic/curved counterterms and
higher-EFT matching. It does not call this proper-plus-LSZ subset the
full finite three-point amplitude.

## Exact scoped claims

For all outgoing pair momenta with k_i^2=m_i, k_j^2=m_j and z=k_i.k_j,
set d_i=l^2+2k_i.l and d_j=l^2-2k_j.l. The full harmonic-projector
numerator, without a soft expansion, is
N_D=N0_D+2z*l^2-2z(d_i+d_j),
N0_D=4z^2-4m_i*m_j/(D-2).
The metric-cubic contacts and the entire on-shell self-energy derivatives
reduce the complete selected graph sum to
deltaT/g=-Gamma(-e)(4pi nu^2)^(-e) B_e/(16pi^2 kappa),
where
B_e=sum_pairs integral_0^1[
 N0_D/2*(A-i0)^(e-1)+2z*(A-i0)^e]dv
 +sum_legs m_i^(1+e)/(2(1+e)),
A=m_i*v^2+m_j*(1-v)^2-2z*v*(1-v).
There is one light-light pair and two heavy-light pairs.

The complete physical-root subtraction fixes B0 and B1=d_e B_e at0,
including their imaginary branches, with no arbitrary finite constant.
Define beta=sqrt(1-4mu/n), V=n^2-4mu*n+2mu^2 and
R3=mu+n/2-4mu(1-mu/n)*atanh(beta)/beta>0.
Then B0=-R3+i*pi*V/(n*beta). The first coefficient is the complete
convergent expression proved in notes/masters.md and proper.py.

For a fixed C1 weight W on[n,n(1+L)], 0<L<=1/2, define ONLY the named
formal combination
D_W(e)=integral_n^{n(1+L)}rho_Hh,e(s)W(s)ds
       +pi*g^2*2Re(deltaT/g)*W(n).
Here rho_Hh,e is the whole S6.291 real cut. Its leading K0/(2e)
pole cancels the known virtual delta coefficient. Its finite value is
integral_0^L[
 K0(n(1+y))W(n(1+y))-K0(n)W(n)]/y dy
 +K0(n)W(n)ln L
 +K0(n)W(n)/2*[1+ln(n/4)+L_B/F0+Re(B1)/R3].
All symbols and the explicit conditional O(e) remainder are defined
in notes/pairing.md. The loop scale and Euler constant cancel. The
mass-unit dependence cancels inside the whole coefficient.

## Deliberate nonclaims

The heavy state is perturbatively unstable. This is a formal expansion
of its pole coefficient, not an exact stable-H asymptotic state or a
uniform approximation throughout a resonance. The surviving imaginary
Coulomb pole multiplies principal values and other cuts. H-metric,
local/higher-EFT and additional mass matching have not been set to zero.
No complete physical inclusive/dressed observable, near-resonance width,
fixed-transfer Regge remainder, all-loop bound, full-source b20 sign,
vacuum-to-bounce state/domain map or original V/G/B/P8 closure follows.
