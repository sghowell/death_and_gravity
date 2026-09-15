# Complete normal-sheet masters and their finite derivatives

## Whole radial reduction

For each on-shell pair the two-mass bubble polynomial equals
A(v)=m_i*v^2+m_j*(1-v)^2-2z*v*(1-v):
v*m_i+(1-v)*m_j-v(1-v)*(m_i+m_j+2z)=A(v).
In the complete triangle simplex put x=u*v,y=u*(1-v). Its Jacobian
has magnitude u, its denominator is u^2*A, and0<=u,v<=1.
The raw C0 convention gives
C0=-Gamma(1-e)(4pi nu^2)^(-e)/(2e)*integral(A-i0)^(e-1)dv.
Using Gamma(1-e)=-e Gamma(-e), this becomes
C0=Gamma(-e)(4pi nu^2)^(-e)/2*integral(A-i0)^(e-1)dv,
while B0_pair=Gamma(-e)(4pi nu^2)^(-e)*integral(A-i0)^e dv.
These are the whole masters. Keeping the Feynman branch is essential
when A has a physical root. A singular real integral at e=0 is not
a replacement for its analytic boundary.

B0_on=Gamma(-e)(4pi nu^2)^(-e)*m_i^e/(1+2e).
Substitution in the entire proper-plus-external result gives
deltaT/g=-Gamma(-e)(4pi nu^2)^(-e)*B_e/(16pi^2 kappa),
B_e=sum_pairs[N0_D/2*J_pair(e)+2z*M_pair(e)]
 +sum_legs m_i^(1+e)/(2(1+e)),
where J_pair(e)=integral(A-i0)^(e-1) and
M_pair(e)=integral(A-i0)^e.

## Full physical light-light moment

Set beta=sqrt(1-4mu/n),0<beta<1. For A=mu-n*v*(1-v),
v=(1+w)/2 and reflection reduce the full moment to
(n/4)^alpha*integral_0^1(w^2-beta^2-i0)^alpha dw.
Split at w=beta. On the inner segment set w=beta*u; on the
outer segment use t=1-u^(-2). The exact result is
M_alpha=(n/4)^alpha*beta^(2alpha+1)/2*
 [B_(1-beta^2)(alpha+1,-alpha-1/2)
  +exp(-i*pi*alpha)*B(1/2,alpha+1)].
Initially Re(alpha)>-1 makes the split root integrable. Analytic
continuation applies to this whole expression, not each pole
separately. At alpha=e-1 the two beta-function poles cancel.

Put c=1-beta^2 and write
J_e=(n/4)^(e-1)*beta^(2e-1)/2*H_e,
H_e=B_c(e,1/2-e)-exp(-i*pi*e)*B(e,1/2).
The identity
H_e=[c^e-exp(-i*pi*e)]/e+U_e-exp(-i*pi*e)*V_e
subtracts the entire common endpoint pole, with
U_e=integral_0^c t^(e-1)[(1-t)^(-1/2-e)-1]dt,
V_e=integral_0^1 t^(e-1)[(1-t)^(-1/2)-1]dt.
At t=0 the bracket differences gain one power of t; near the
remaining square-root endpoint the integrable logarithmic powers
from e derivatives are harmless. Since0<c<1, U has a strictly
separated upper endpoint. Thus these subtracted functions and their
first two e derivatives are continuous on0<=e<=1/8. The quotient
has its unique removable value, not an arbitrary finite prescription.

For explicit convergent first derivatives use t=1-z^2:
V_e=2*integral_0^1(1-z^2)^e/(1+z)dz.
Then
V0=2ln2, V1=2ln(2)^2-pi^2/6,
U1=2*integral_beta^1[
 ln(1-z^2)/(1+z)-2ln z/(1-z^2)]dz.
The integrand has only an integrable logarithm at z=1; the apparent
ln z/(1-z^2) singularity is removable. Direct evaluation gives
H0=-2atanh(beta)+i*pi,
H1=(ln(c)^2+pi^2)/2+U1-V1+i*pi*V0.
Consequently
J_ll,0=2H0/(n*beta),
J_ll,1=2[H1+ln(n*beta^2/4)*H0]/(n*beta).
This fixes the complete complex first coefficient. It retains the
Coulomb phase and involves ordinary convergent logarithmic integrals
only after the necessary analytic subtraction.

The entire light bubble logarithm is
L_ll=integral_0^1 ln(A-i0)dv
    =ln mu-2+2beta*atanh(beta)-i*pi*beta.
It follows either by integrating the split-root logarithm directly or
differentiating M_alpha at0. The imaginary interval has length beta.

## Positive heavy-light moments

For either heavy-light pair A=n*v+mu*(1-v)^2 is strictly positive.
Under v=(1-t)/(1+t),
A=n*(1-beta^2*t^2)/(1+t)^2, -dv=2dt/(1+t)^2.
Thus
J_hl(e)=2n^(e-1)*integral_0^1
 (1-beta^2*t^2)^(e-1)*(1+t)^(-2e)dt,
J_hl,0=2atanh(beta)/(n*beta),
J_hl,1=integral_0^1 ln A/A dv,
L_hl=integral_0^1 ln A dv.
No heavy-mass series or massless limit is used. Positivity and the
compact parameter interval justify the needed derivatives.

## Full finite cubic coefficient

For each pair N0 at D4 is4z^2-2m_i*m_j and its e derivative is
2m_i*m_j. Therefore
B1=sum_pairs[m_i*m_j*J_pair,0
 +(N0_D4/2)*J_pair,1+2z*L_pair]
 +sum_legs m_i*(ln m_i-1)/2.
All three pairs and all three leg terms are included.
The zeroth coefficient is
B0=-R3+i*pi*V/(n*beta),
R3=mu+n/2-4mu(1-mu/n)*atanh(beta)/beta.
Expansion of the entire common Gamma/scale factor gives
deltaT/g=B0/(16pi^2 kappa e)
 +[B1+(EulerGamma-ln(4pi nu^2))*B0]/(16pi^2 kappa)+O(e).
No D0-only guess can supply the missing B1.

Exact residuals check the whole polynomials, radial and root Jacobians,
Cayley denominator, primitives, subtracted integrands and Gamma
coefficients. Independent high-precision tests compare the entire
beta-function moment with nonsingular physical-segment integrals and
the first derivative with the convergent expressions above. Those
finite numerical examples support, but are not a formal proof of,
the analytic all-domain argument.
