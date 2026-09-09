# Positive inverse density and exact moments

Let R(p)=1/F(p)=-1/H(p). With D,U from the first-sheet note,

    rho(tau)=Im R(-tau+i0)/pi=U/(D²+pi² U²)>0,
    z=1-4m²/tau, tau>4m².

The positivity is a scalar inverse spectral identity, not a claim that R itself is a fundamental particle propagator or a UV-completed metric spectrum.

The threshold denominator is nonzero. Since sqrt(z) atanh(sqrt(z)) is analytic in z near zero, the explicit density has a differentiable endpoint expansion

    rho(tau)=675/512 sqrt(z)+O(z^(3/2)).

At infinity, writing L=log(tau/m²), the cut formula gives D=2L-14/15+O((m²/tau)L) and U=2+O(m²/tau). Consequently rho=1/(2L²)+O(L^-3)+O((m²/tau) polylog(tau)). The explicit rational/logarithmic form justifies termwise differentiation, not differentiation of an unspecified remainder.

There are no pole residues by the gap proof. R=O(1/log|p|) on large circles, so the outer Cauchy contribution vanishes. The threshold arc vanishes by continuity, and the cut integral converges absolutely with its Cauchy denominator because rho/tau is integrable at infinity. Cauchy's formula, with the cut orientation retained, yields

    R(p)=-integral_(4m²)^infinity rho(tau)/(p+tau) dtau.

There is no constant instantaneous term. Evaluating at p=0 and differentiating there give the two exact positive moments

    integral rho(tau)/tau dtau = 1/4,
    integral rho(tau)/tau² dtau = 9/(560m²).

The second identity uses H'(0)/H(0)², independently determined by the original radial integral. These moment statements follow from the analytic representation and exact endpoint data, not from a sampled spectral quadrature.

