# S6.42: source-aligned curled trace

Status: exact source, quadratic constraints, canonical frequency and
nonlinear source-bound checkpoint; not a full acceptance verdict.
Original P8 remains OPEN. This is separate from frozen ancestors.

Keep the complete S6.41 source-centered mass term with mu=11/9, but
replace its unshifted curl by the curl of the specified local one-form

    W_mu=T_mu-B(phi,x)*partial_mu phi,
    B(u,x)=-d(u)*(1+x)/2,
    d(u)=15H(u)/(4h(u)), h=(1+u²)^3, H=4u/(1+u²),
    x=g^(mu,nu)*partial_mu phi*partial_nu phi (source signature -+++).

Thus the literal new density beyond S6.37 is

    -(mu/2)*g^(mu,nu)*(T-Tstar)_mu*(T-Tstar)_nu
       -(zeta/4)*F(W)_mu_nu*F(W)^mu_nu,
    mu=11/9, zeta=zeta_physical/(M² tau²).

These are normalized coefficients after the same overall factors as
S6.41. B_physical=B_normalized/tau. The physical metric, scalar clock,
free matter, original mass-centering source and counterterm are kept.
B is a local coefficient function, not an externally prescribed time
source. It has no Hessian dependence. No old action is overwritten.

## Precisely scoped claims

1. Derive the full shifted stationary source Sstar=Tstar-B*dphi and
   check its background and every first variation from the all64
   stationary connection, not only an isolated vector ansatz.
2. In the exact nonlinear spatial chart, verify

       Sstar_normal=(4p²-1)*(K_hat-3H*s)+(3/2)*s*Q_lower,
       s=sqrt(-x), 4p²-1=(s²-1)/h.

   Here Q_lower is the unchanged S6.37 coefficient-ODE solution; it
   is not spatial momentum. Q_lower=Q_lower,x=0 at x=-1. The proposed
   identity makes the source quadratic in deviations from the
   actual rolling solution while preserving the primary lapse null.
3. The complete quadratic physical action, including metric/matter
   constraints, separates into the old CD/M1 system and a positive
   three-polarization Proca block. All five scalar/vector principal
   cones are exactly the original matter cone. The original tensor
   block is unchanged. The first-order center remains regular.
4. Derive the actual time-dependent canonical vector equations at fixed
   comoving momentum, with q=k_com²/a² and all normalization derivatives.
   A positive canonical-frequency floor is not a stationary S-matrix
   gap, a nonlinear health theorem or a cutoff.
5. Quantify the nonlinear source on a stated field/derivative domain
   using the actual lower-coefficient ODE. Do not turn an algebraic
   source bound into a retarded-inverse or EFT-remainder claim.

The certified frequency range is 0<zeta<=1/2000: at fixed comoving
momentum every canonical vector frequency obeys Omega²>=q+1985 in
normalized units, uniformly in real u. See notes/modes.md for the
actual normalization derivatives and the homogeneous vector chart.
This is a quadratic rolling-background result, not a stationary gap.

The nonlinear source and electric-curl bounds are stated in
notes/source.md on the original closed x tube and with explicit
spatial-derivative norms. The shifted source is nonzero at second
order and the mass term has cubic light-vector interactions.

Zero curl and zero
spatial momentum are separate rank/chart controls. This candidate
still requires nonlinear secondary constraints, a controlled heavy
elimination with an explicit state and errors, loop/cutoff estimates,
and the adopted vacuum/finite-gravity V/G/B tests. A scalar/vector
quadratic decoupling is not UV completion or P8 closure. The scope is
exact symbolic verification with written proofs, not formalization.
