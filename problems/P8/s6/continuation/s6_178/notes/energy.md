# Actual constrained Hamiltonian and a finite-time energy theorem

Use A=sqrt(kappa*zeta)W, J=sqrt(kappa/zeta)S, m^2=1/zeta.
For spatial LOWER components and the physical +--- CD metric, the
uneliminated canonical vector density is

    a/2 (Adot_i-partial_i A0)^2-F_ij^2/(4a)
      +a^3 m^2 A0^2/2-a m^2 A_i^2/2
      -a^3 A0 J0+a A_i J_i.

The light contact +kappa a^3 S_mu S^mu/2 is present as well.
The canonical momentum is pi_i=a(Adot_i-partial_i A0).
Integrate pi_i partial_i A0 by parts in space, then vary A0:

    A0=(J0-div pi/a^3)/m^2.

No conserved-current or transverse-source condition was imposed.
The reduced driven Hamiltonian, BEFORE its separate light contact, is

    H=integral[pi^2/(2a)+F_ij^2/(4a)+a m^2 A_i^2/2
               +(div pi-a^3 J0)^2/(2a^3 m^2)-a A_i J_i].

The contact Hamiltonian is
-kappa a^3(S0^2-S_i^2/a^2)/2. Its temporal constant cancels the
a^3 J0^2/(2m^2) in H, leaving the spatial contact. This cancellation
does not remove the temporal source from A0 or from the equations.

Define a POSITIVE REFERENCE energy, not the full driven Hamiltonian,

    E0=integral[pi^2/(2a)+F_ij^2/(4a)+a m^2 A_i^2/2
                  +(div pi)^2/(2a^3 m^2)]
      =Ee+EB+Em+ED.

The equations are

    Adot=pi/a-grad(div pi)/(a^3 m^2)+grad J0/m^2,
    pidot=div F/a-a m^2 A+a J_sp.

The compact source, zero initial data and causal finite propagation
on the globally hyperbolic CD metric justify every spatial integration
by parts. The exact energy identity is

    E0'=Hubble[-Ee-EB+Em-3ED]
          + integral[pi.J_sp+a A.grad J0
                         +(div pi)(div J_sp)/(a^2 m^2)].

The curl contracted with grad J0 vanishes after integration, but the
mass and div-momentum terms do NOT vanish. Cauchy-Schwarz across
Ee,Em,ED gives

    E0'<=3|Hubble|E0+sqrt(2E0) F_J,
    F_J^2=a||J_sp||2^2+||div J_sp||2^2/(a m^2)
                         +a||grad J0||2^2/m^2.

Use sqrt(2E0+epsilon^2) and then epsilon->0 to justify the norm
inequality at zero energy. For v=sqrt(2E0) with zero initial data,

    v(t)<=integral_(t0)^t exp[(3/2)integral_s^t |Hubble|] F_J(s) ds.

On I, integral_I |Hubble|=4 log(5/4), so every retarded subinterval has
propagator weight at most B=(5/4)^6=15625/4096. There is no temporal
resonance denominator. All spatial derivatives commute with these
equations because the coefficients depend only on time; the same
argument gives the differentiated energy estimate when its source
norm is finite.

This is a finite-time causal mapping between explicit energy/source
spaces. It does not assert a bounded whole-time Fourier inverse at
the Proca pole, nor physical EFT validity for arbitrarily high momentum.
The exact finite-time resonant oscillator control illustrates the
distinction. Independent forced transverse and full longitudinal
Fourier modes on the actual CD geometry check the energy calculation,
not replace its all-momentum proof.
