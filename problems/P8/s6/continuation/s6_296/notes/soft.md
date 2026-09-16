# Full dimensional soft kernel and finite conversion

At Born take four future physical momenta p_i, common COM energy E,
masses mu and signs eta=(-1,-1,+1,+1),sum eta_i p_i=0.
The full physical polarization contraction of the conserved soft current
is sum_ij eta_i eta_j[(p_i.p_j)^2-mu^2/(D-2)]/
[(p_i.q)(p_j.q)], retaining all16 ordered pairs, including four diagonals.

In D=4+2e define K_e as its normalized angular average after removing
omega^-2. Feynman pairing uses p(x)=xp_i+(1-x)p_j, p(x)^0=E,
m_x^2=mu+2x(1-x)(p_i.p_j-mu),b_x=1-m_x^2/E^2.
The normalized angular denominator is
J_e=2F1(1,3/2;3/2+e;b_x)/E^2.
At zero J0=1/m_x^2. Substituting y=(1-z)/2 in its log-weight derivative
gives raw integral[2ln2-2atanh(sqrtb)/sqrtb]/(1-b).
The normalized weight contributes2ln2-2, hence
J1=2[1-atanh(sqrtb)/sqrtb]/m_x^2.
The removable b0 value is retained explicitly.

Therefore K1=sum_ij eta_i eta_j integral dx[
mu^2/2+2N_ij(4)(1-atanh(sqrtb_x)/sqrtb_x)]/m_x^2.
Moving diagonals contribute6mu-4mu*atanh(beta)/beta, not their rest-frame
value. Their omission changes the finite answer.

The radial/phase normalization is
p(e)=(4pi)^(-e)Gamma(3/2)/Gamma(3/2+e),p0=1,
p1=EulerGamma-2-lnpi. Thus the full leading real rate is
(resolution/nu)^(2e)*p(e)*K_e/(8pi^2*kappa*e).

On the physical S288 i0 sheet, K0=-2ReBsoft. After pairing the virtual
soft term from the exact S278 convention, the finite conversion is
Delta=[K1+(EulerGamma-2-lnpi)K0]/(8pi^2*kappa).
It is independent of resolution and reference scale; the bounded
nonsoft real remainder still depends on resolution. Both forward and
backward real kernels vanish for all D, but the virtual amplitude
Coulomb phase is not removed. The explicit physical-sheet/UV connection
is in inclusive.md.
