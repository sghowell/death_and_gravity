# Full-D transform and every rung

DistinguishDkin=D=sqrt(s(s-4)) fromDdim=4+2e.
The transverse dimension isd=2+2e, andeta_e=V_e/(8pi*kappa*D).
The Euclidean transverse Fourier identity with the original(2pi)^-d
measure gives

int d^dq/(2pi)^d exp(iqb)/q^2
=Gamma(e)/(4pi^(1+e))*b^(-2e).

Multiplying the Born numerator and flux gives the impact phase

chi_e(b)=eta_e*Gamma(e)*(pi*b^2)^(-e).

For each fixedb>0,
chi_e-eta/e -> eta[2/V-gamma_E-ln(pi*b^2)].
The dimensional derivativeV_e'=2 supplies the finite2/V term.
It is fixed by the existing numerator, not by a new matching choice.

The inverse Riesz transform in transverse dimensiond is

int d^db exp(-iqb)*(b^2)^(-a)
=pi^(d/2)*4^(d/2-a)*Gamma(d/2-a)/Gamma(a)
 *(q^2)^(a-d/2).

Initially it is understood in its convergence/analytic-continuation
strip and then as the associated tempered distribution. Apply it with
a=ne at fixed positive integern to the nth Gaussian exchange term.
For nonzerotau=q^2 this gives

A_n=(2D/i)(i*eta_e)^n/n!*Gamma(e)^n*pi^(-ne)
    *pi^(1+e)*4^(1-(n-1)e)
    *Gamma(1-(n-1)e)/Gamma(ne)
    *tau^((n-1)e-1).

The first rung is exactlyV_e/(kappa*tau).
WithL=n-1 and division by that full-D Born, use
Gamma(e)=Gamma(1+e)/e andGamma(ne)=Gamma(1+ne)/(ne) to obtain

R_L=(i*eta_e/e)^L/L!*H_L(e),
H_L=Gamma(1+e)^(L+1)*Gamma(1-Le)/Gamma(1+(L+1)e)
    *(tau/(4pi))^(Le).

This formula is independently calibrated before its Laurent expansion.
For each fixedn one can choosee small enough that the continued
Gamma(1-(n-1)e) is away from its ultraviolet poles.
There is no claim that one such positivee works for every rung order
inside an unregulated small-b integral.
