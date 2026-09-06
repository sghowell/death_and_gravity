# Full actual map and exact energy-constraint propagation

All variables and physical prefactors are fixed in
[FORMULATION](../FORMULATION.md). Bars always mean the same actual A.8
metric and transported state, evaluated at the same conformal time.

## 1. An exactly conserved switch

Writing the common physical prefactor as `F=hbar/(pi^2 A^4 eta_star^8)`,
the source density and pressure are

    rho_ext=F c/a^4,
    p_ext=F[c/(3a^4)-c'/(3h a^4)],
    Theta_ext=F c'/(h a^4).

Direct differentiation gives `rho_ext'+3h(rho_ext+p_ext)=0` on any
unknown geometry with positive `a,h`. Prescribing `c`, rather than an
independent pressure, is important. Using `p_ext=rho_ext/3` while
switching `c` would leave the nonzero conservation defect `F c'/a^4`.

The baseline defect used to define `cbar` is itself conserved, because
the Einstein tensor, ordinary radiation and actual quantum stress are
all conserved. Consequently the pressure obtained from the displayed
formula equals the old pressure defect. The old geometry therefore
solves the complete *forced* equations with `cbar`, not just their
density component. It does so in the original free past as well.

On the exact radiation past `a'=1`, `u=0`, the actual quantum density
bracket is `h^4/960`; hence `cbar=-h^4/960`, not zero. This rules out
describing the source as globally compact while keeping that past.
It does not prevent extinguishing it before the new physical slab.

## 2. Exact integrated trace and constraint

For the actual mode functional from A.10, set

    S=Rmode-K/2+u/10, q=S/a^2,
    J'=u'S, ell'=h u^2,
    P_rho=h^4/960+[-hS'+(h^2-u)S+J]/8-ell/32.

The quantum density is `F P_rho/a^4`; both histories and the quantum
radiation constant are fixed by the original state. A.8's exact trace
becomes

    P_trace = -a^2(q''+2h q')/8-u^2/32-h^2(u+h^2)/240.

With the physical `delta` relation, the full traced equation is

    q''+2h q' = u/(60delta)
                -[u^2/4+h^2(u+h^2)/30]/a^2
                +8c'/(h a^2).                                    (1)

In particular the large gravitational factor `1/(60delta)` is kept.
It is not replaced by the A.10 response constant. Define

    C=3a'^2-3-2880delta(P_rho+c),
    Dtrace=-6u/a^2-2880delta(P_trace+c'/h)/a^4.

Using the actual history equations gives the exact identity

    C'=a^4 h Dtrace.                                             (2)

At the starting slice, `c=cbar` and all metric/state data are unchanged,
so `C=0` exactly. A solution of (1) therefore satisfies the complete
energy equation throughout the transition. The pressure equation then
follows from density and trace. Trace evolution alone, with a nonzero
initial `C`, would not have this implication.

## 3. Eliminate the switch derivative without discarding it

The metric satisfies `a''=-ua` and `h'=-u-h^2`, so

    (1/h)'=1+u/h^2.

Writing `P=a^2q'`, integration of (1) is exactly

    P(s)=P0+ integral_0^s [a^2u/(60delta)-u^2/4-h^2(u+h^2)/30] dr
          +8[c(s)/h(s)-c(0)/h(0)
              -integral_0^s c(r)(1+u(r)/h(r)^2) dr],              (3)
    q(s)=q0+integral_0^s P(r)/a(r)^2 dr.

No estimate for `chi'/L`, `c'`, or a higher cutoff jet is needed in
the uniform contraction norm. This is an integration by parts of the
complete conserved source, not an omission of its pressure.
The finite initial data are the actual `qbar(0),qbar'(0)`.

## 4. The complete actual-Wick fixed point

At the last free slice the potential and all of its jets vanish. In
the unchanged named prescription the exact logarithmic split is

    S=Rmode[u]+4pi^2 D_gammaE[u]+d(a)u,
    d(a)=-19/60-(1/2)log(a/2), d'=-h/2,
    D_gammaE[f]=-(1/(8pi^2)) integral log(t-r) f'(r) dr.

Its inverse is the causal positive pole-plus-cut operator from A.10,
denoted `J_gammaE`. No pole or finite curvature term is removed.
The complete original prepared history is retained in `Rmode`.

For `X` in the ball of FORMULATION, put `W(s)=integral_0^s X` and
`u=ubar+W`, solve `a''=-ua` with the common initial `a,a'`, and obtain
`P,q` from (3). The unknown and its perturbation are zero before the
new starting slice. Define

    F[X] = a^2(q'+2h q)-abar^2(qbar'+2hbar qbar)
           -(Rmode'[u]-Rmode'[ubar])
           -d(a)X-[d(a)-d(abar)]ubar'
           +(h/2)W+[(h-hbar)/2]ubar.                            (4)

The actual map, including every local and state term, is

    X = (J_gammaE/(4pi^2)) F[X].                                (5)

The singular linear history cancels only in the **difference**: its
remaining kernel in (5) starts at the new slice. The nonlinear mode
functional still includes all earlier scattering. A.10's shared-history
bound controls that interaction.

The initial actual value `S=a^2q` agrees by definition. Equation (5),
once smoothness is established, equates its derivatives. It therefore
implies `S=a^2q` identically, not a separate state ansatz. Equation (1)
is then the trace of the actual transported state's stress.

## 5. Quantitative existence and physical interpretation

[The full contraction proof](contraction.md) establishes, with exact
rational arithmetic, the domain and estimates

    ||J_gammaE/(4pi^2)|| < 21/100,
    Lip(F) < 4/7, Lip((J_gammaE/(4pi^2))F) < 3/25,
    ||(J_gammaE/(4pi^2))F[0]|| < 63*10^-9.

The closed ball of radius `10^-6` is preserved since
`63*10^-9+(3/25)10^-6 < 10^-6`. Banach's theorem gives its unique
continuous fixed point and the explicit error estimate in FORMULATION.
No unknown numerical solution is inserted to supply these constants.

[The regularity proof](regularity.md) upgrades this particular
flat-start Picard construction to a smooth potential on the **same**
slab. The actual original state, propagated on its smooth positive
metric, remains positive and Hadamard. Equations (1), (2), and (5)
then prove the full forced SEE pointwise. On the final half-slab the
source and all its derivatives vanish, giving the unforced equation
with the unchanged ordinary radiation normalization.

This is not an exact SEE on the old A.8 metric. The new metric and state
co-evolve after the common preparation. Nor does this tiny conformal
interval prove whole-window shadowing or a cosmological focusing claim.
