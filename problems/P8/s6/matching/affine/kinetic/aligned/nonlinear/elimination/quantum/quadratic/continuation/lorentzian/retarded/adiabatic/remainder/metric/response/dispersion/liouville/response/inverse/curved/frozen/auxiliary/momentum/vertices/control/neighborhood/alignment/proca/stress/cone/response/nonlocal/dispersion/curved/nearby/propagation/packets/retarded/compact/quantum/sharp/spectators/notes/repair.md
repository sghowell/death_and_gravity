# Necessary repair and controls

Use the exact physical canonical comparison from principal.md. If a
finite Hermitian pencil with positive K_full has no real characteristic
outside either matter lightcone, P(1) and P(-1) must be positive
semidefinite. Otherwise its negative minimum at an endpoint becomes
positive at large absolute speed and forces a root beyond that endpoint.
This is a necessary condition; it is not asserted sufficient for a
healthy, hyperbolic or interacting system.

Restriction to the canonical clock direction gives

    1-c^2+DeltaK_11+2 DeltaC_11-DeltaH_11 >= 0,
    1-c^2+DeltaK_11-2 DeltaC_11-DeltaH_11 >= 0.

Together these imply

    DeltaK_11-DeltaH_11 >= (c^2-1)+2 abs(DeltaC_11).

Averaging the full light restrictions also yields the Loewner necessary
condition DeltaK-DeltaH >= diag(c^2-1,0). Consequently, with operator
norms in the specified OLD kinetic canonical comparison,

    ||DeltaK|| + ||DeltaH|| >= c^2-1 > 5*10^-6.

No smallness assumption on C_full or the heavy/cross blocks was used.
For prescribed exact nonnegative bounds epsilon_K,epsilon_H with sum
strictly below delta=5*10^-6, the positive margin
delta-epsilon_K-epsilon_H certifies failure of principal cone repair
under the theorem's other hypotheses. repair.budget validates the strict
domain. Its example uses epsilon_K=epsilon_H=delta/4. These are conditional
inputs, not actual UV error estimates or authorizations to change a model.

## Why checking one direction is insufficient

Set DeltaK=DeltaH=0 and DeltaC_11=(c^2-1)/2. The isolated clock pencil is

    s^2+(c^2-1)s-c^2=(s-1)(s+c^2).

The positive-speed root is exactly luminal, but the negative-speed root
is -c^2<-1. At s=-1 the pencil is -2(c^2-1)<0. This explicit control
prevents a one-orientation mixed-term change from being called a repair.

## Why changing the light block is a separate route

The algebraic change DeltaH=diag(-(c^2-1),0), DeltaK=DeltaC=0 gives
K_new=H_new=I in the two-mode comparison. It saturates the clock-direction
necessary norm bound. This is deliberately a permissive algebraic
control: the theorem does not exclude genuinely changed light actions.
It supplies no covariant operator, new on-shell solution, all-sector
health proof, EFT validity scale, Lorentzian vacuum or matching model.
None of the old action, background, quantum state or certificates is
overwritten by this comparison.
