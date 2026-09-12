# Continued metric second variations before the dimension limit

Use g=a(eta)^2[-deta^2+Eij dxidxj], E=exp gamma, detE1, sigma=log a. Let f=a^(d-3), h=sigma', c=sigma''+(d-1)h^2 and b=2d sigma''+d(d-1)h^2. All primes in this paragraph denote eta derivatives.

The exact conformal scalar identity is R=a^-2(Rbar+b). With K=E^-1 E'/2 and trK0, Rbar=R_sp+trK^2. Its quadratic compact spatial integral is

    Rbar2=||gamma'||^2/4-||grad gamma||^2/4+||div gamma||^2/2,

while Rbar1=partial_i partial_j gammaij. These identities give the Einstein and scalar-curvature-square Hessians without an isotropic-only extrapolation.

For the Ricci tensor,

    R00=Rbar00-d sigma'',
    Rij=Rbarij+(d-1)h Kij+c Eij,
    R0i=Rbar0i.

Here Rbar00=-trK^2 and tr(E^-1 Rbar_sp)=R_sp. Thus the spatial quadratic Ricci-square difference, before multiplying f, is

    Ricbar1_squared_spatial
    +(d-1)h[q gamma.gamma'-2(gamma P).(gamma'P)]/2
    +2c R_sp2.

The flat linearized invariants for a spatial tracefree perturbation are

    Riem1^2=||gamma''||^2-2q||gamma'||^2+2||gamma'P||^2
              +q^2||gamma||^2-2q||gamma P||^2+(P.gamma.P)^2,

    Ric1^2=||gamma''||^2/4+q gamma''.gamma/2
              -(gamma''P).(gamma P)-||gamma'P||^2/2
              +q^2||gamma||^2/4-q||gamma P||^2/2+(P.gamma.P)^2/2,

    R1^2=(P.gamma.P)^2.

The d+1 dimensional Weyl square is Riem^2-4Ric^2/(d-1)+2R^2/[d(d-1)]. Its quadratic action is f times this flat combination because the background is conformally flat.

## Proper-clock spatial operators

Set A=a^(d-2), F=a^(d-4), S=T-2V and source time jets Gamma0,Gamma1,Gamma2. The spatial Hessian differences are

    H_R=-A p^2 S Gamma0/2,

    H_R2=-A[2d H'+d(d+1)H^2]p^2 S Gamma0+2F p^4 W Gamma0,

    H_Ric=A p^2(T-V)[Gamma2+(d-2)H Gamma1]
             -2A[H'+(d-1)H^2]p^2 S Gamma0
             +F p^4(T/2-V+W)Gamma0,

    H_W=4(d-2)A p^2(T-V)[Gamma2+(d-2)H Gamma1]/(d-1)
          -2(d-3)A[H'+(d-2)H^2]p^2 S Gamma0/(d-1)
          +2(d-2)F p^4 S Gamma0/(d-1)
          +2(d-2)F p^4 W Gamma0/d,

    H_Riem=H_W+4H_Ric/(d-1)-2H_R2/[d(d-1)].

These follow by integration by parts in detector time before d is restricted, then d_eta=a d_t and division by a for the proper-clock measure. The original R_old sign is unchanged.

They imply the exact continued Euler variation

    H_Euler=-a^(d-2)(d-3)[(d-2)H^2+2H']p^2(T-2V)Gamma0.

It vanishes atd3 but its dimension derivative need not vanish. It is therefore kept in the original fixed pole variation. The volume potential3m^4 has zero unimodular second variation at every d.

An independent direct computation expands the metric and inverse metric through the mixed detector/source coefficient, builds Christoffel symbols and the full Riemann tensor, contracts all invariants and varies the resulting density. It agrees with these operators for tensor,vector,scalar and noncommuting tensor pairs in dimensions3,4,5,6. The mixed metric coefficient is(D Gamma+Gamma D)/2; it is not omitted in favor of linearized curvature products. The formulas above are derived analytically, not inferred from that finite integer sample.
