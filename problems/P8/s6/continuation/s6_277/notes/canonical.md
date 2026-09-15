# Full homogeneous canonical flow and normalized constraint weights

Let V=e^(3alpha) and use one explicit diagonal anisotropy beta. The raw
mixed momentum and metric are

    pi*gamma=(Palpha/6)I+diag(Pbeta/4,-Pbeta/4,0),
    gamma=e^(2alpha)diag(e^(2beta),e^(-2beta),1).

Direct differentiation gives pi:gamma_alpha=Palpha and
pi:gamma_beta=Pbeta. Its traceless squared norm gives
s=Pbeta^2/(8 kappa^2 V^2). The complete one-form, including homogeneous
M1 and heavy h, is therefore the usual absolute canonical one-form.
These identities do not identify the data with any Gaussian reference.

For the normalized Hamiltonian H(u,N,p,m,s,eta,ph), the full homogeneous
Hamiltonian is kappa V H, per unit coordinate torus volume. The rates are

    alpha_u=H_p/3,
    p_u=-H+mH_m+2sH_s+phH_ph,
    m_u=-mH_p, s_u=-2sH_p,
    eta_u=10^100 H_ph,
    ph_u=-10^100 H_eta-phH_p,
    M1_u=H_m.

For positive conserved Pbeta,
beta_u=N Pbeta/(2 kappa V R^(1/4)). Its sign can instead be retained
by using the signed conserved Pbeta. Both cyclic coordinates are integrated,
not held fixed. The full spatial Proca zero sector is invariant; diagonal
Bianchi-I geometry and momenta remain diagonal by covariance. All spatial
derivatives vanish, so the spatial momentum constraints vanish identically.
The temporal solution is the complete original T(u,N,x), not T=0 off clock.

For any density-coordinate scalar A define

    L_H A =
       A_p[-H+mH_m+2sH_s+phH_ph]
       -H_p[mA_m+2sA_s+phA_ph]
       +10^100[A_eta H_ph-A_ph H_eta].

This is {A,kappa V H}. It is NOT an antisymmetric bracket in the two
normalized arguments A,H: in particular L_C C=-C C_p.
Independent differentiation in the absolute alpha,Palpha,h,PH variables
checks the density factors and the heavy conversion, including mixed
heavy/shape fixtures. The canonical proof itself follows by the chain rule
from the entire one-form and holds off the constraint.

The fixed source profile functions are coefficients of the original local
action; they are not replaced by the classical energy or pressure of this
comparison solution. No quantum expectation or scalar-center subtraction
is made in this classical result.

The regular-side equations are equivalent to the full homogeneous
Euler equations of that action: R,N,Gamma stay positive, the original
mixed trace/temporal Legendre system is invertible, and C_N is nonzero
away from the endpoint. Unitary clock remains timelike. The original
covariant clock equation follows from the remaining equations and their
Noether identity; it is not removed as an additional physical assumption.
