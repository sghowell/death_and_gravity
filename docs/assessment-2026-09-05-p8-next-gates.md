# P8 next-gate research handoff after A.7 and S5.11

Recorded 2026-09-05. **This is a research plan with candidate derivations,
not a new certificate or ledger result.** The completed checkpoint is commit
`2f55aded79b1c503fa3a8f6b5cbe41c4d70054cb`: all 826 P8 tests, both new
read-only certificate replays, source hashes and new-subtree Ruff checks pass.
P8 remains open. Neither next task currently requires a new user choice.

The immutable immediate inputs are
[A.7](../problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/FORMULATION.md),
certificate SHA256
`cc4ec02bcb23e0ef01d0fb58ad4b95a3bb602a5bbbfc6871d044b2a09a78e963`, and
[S5.11](../problems/P8/s5/matter/physical/control/tree/loops/tensor/scalars/FORMULATION.md),
certificate SHA256
`b00d124cf5c7ab712cb7f3198766f1aa542065d801741273db37855a5824d27c`.
Implement future gates in new subtrees; do not rewrite either input.

## A: quantitative actual semiclassical residual on the existing metric

The strongest immediate bounded task is to compare the actual stress to the
positive A.3 reference and then bound the actual Einstein-plus-radiation
residual on the target plateau. Keep the exact A.7 prepared metric/state,
actual conformal clock, zero cosmological term, and the already named
`lambda=2*sqrt(2)*A*eta_star^2`, `gamma=0` prescription. This need not first
compare the exact clock to an epsilon-linearized one.

The following candidate reconstruction should be independently derived,
tested and certified against A.7's full conserved trace before promotion.
Write `h=a'/a`, `U=-a''/a`, and, with zero preparation-past integrals,

    S=Rmode-K/2+U/10, J=integral U'*S d_eta, L=integral h*U^2 d_eta.

The `U/10` term accounts for the Box-R anomaly without changing the
subtraction prescription. Candidate exact physical components, after
factoring out `hbar/(pi^2*a^4)`, are

    rho=h^4/960+[-h*S'+(h^2-U)*S+J]/8-L/32,
    p=(5*h^4+4*h^2*U)/2880
        +[S''-3*h*S'+(3*h^2+U)*S+J]/24+(U^2-L)/96,
    E=(3*h^4+2*h^2*U)/960
        +[S''-4*h*S'+4*h^2*S+2*J]/16+(U^2-2*L)/64.

Required checks include the exact trace/anomaly dictionary, conservation,
the actual common-past integration constant, and negative controls for
discarding either history. A coincident-stress bound is not a QSEI bound.

Using auxiliary span `3*eta_star`, the local logarithmic factor is
`log((3/2)*y*f^(1/4))+5/6`. A sufficient proposed absolute cap is 3 on the
active history `1<=y<=3`. The candidate kernel bounds are
`|K^(j)|<=delta*L_j/eta_star^(j+2)`, where A.7 supplies every b_j:

    L0=3*b1+3*b0,
    L1=3*b2+3*b1+2*b0,
    L2=3*b3+3*b2+4*b1+(4+delta_bar*b0)*b0.

Candidate rounded constants for
`|T_actual-T_ref| <= hbar*(delta*C1+delta^2*C2)/(pi^2*A^4*eta_star^8)`
are `(C1,C2)=(14000,8*10^14)`, `(236000,9*10^14)`, and
`(360000,17*10^14)` for density, pressure and EED. They require a new exact
majorant replay; they are not A.7 certificate claims. A useful test target is
one-percent reference-scale accuracy for `delta<=10^-14`, including physical
epsilon=1 when `t_star>=2*10^7*sqrt(d)`.

For the actual SEE residual, combine this with the exact A.5 frozen defects.
On the plateau put `z=delta/y^4<=1/32`. The elementary proposed enlargement

    [(1-z)^(-2)-1]/z <= 2048/961, y^12>=4096

uses the continuous value 2 at z=0. It explains the literal denominator
`1922=2*31^2` in the candidate bound

    |D_i_actual| <= delta^2/(A^2*eta_star^4)
        *[c_i/1922+2880*(C1_i+delta*C2_i)], c_i=(3,5,9).

Verify the exact coupling restoration
`epsilon*kappa*hbar/pi^2=2880*delta*A^2*eta_star^4` and component signs.
Completion means an actual finite-amplitude residual certificate on this
plateau. It does not mean an exact/nearby SEE solution, stability, control
over the off-shell preparation interval, or a cosmological focusing theorem.
Afterwards, exact-solution control and the full two-point-function estimate
needed for perturbed QSEIs remain separate gates.

## B: a specified heavy-metric parent screen

The next recommended bounded task is a source-aware parent test under the
existing [S6 contract](../problems/P8/s6/FORMULATION.md). A physical-data or
original-equation residual bridge would consolidate S5.11, while adding R^2
would broaden local-operator coverage; neither alone supplies the missing
parent or its finite matching coefficients.

Proposed family: singly coupled Hassan--Rosen bimetric gravity with positive
Einstein coefficients Mg^2 and Mf^2, all positive-field-metric canonical
clock/M1 matter coupled only to the physical g metric, and

    beta0=-3*beta1, beta4=-beta1, beta2=beta3=0, beta1>0.

First independently verify the proportional Minkowski equations, both
spin-2 residues and the mass, then diagonalize with the physical g-source
retained. A TT Schur expansion is a useful independent check, not a complete
scalar/source-projector matching proof. Eliminating a metric via a different
field's equation must not be relabelled as integrating out the massive
eigenstate: these procedures and source terms require separate accounting.
See [Hassan--Schmidt-May--von Strauss, section 2 and Appendix A](https://arxiv.org/pdf/1303.6940)
and [Gording--Schmidt-May](https://arxiv.org/pdf/1807.05011).

In parallel derive both lapse equations before gauge fixing and test exact
CD background compatibility before any derivative truncation. A candidate
obstruction to certify is: on the regular positive-lapse, positive-scale-factor
flat-FLRW dynamical branch, a physical g bounce forces `y=a_f/a_g=1`, and
the two acceleration equations then imply

    -2*(Mg^2+Mf^2)*Hdot_g=(rho+p)_g.

If the independently derived equations establish this over the stipulated
family, positive-metric canonical matter would contradict the CD target
`Hdot_g(0)=4/tau^2`. The falsifiable completion condition is either a regular
CD branch satisfying the full parent equations or a certified obstruction
for this exact parameter/matter family. A failure excludes only that parent
ansatz, not general bimetric models, the DHOST witness, its row or P8.

No common parent, bounce heavy-spectrum bound, quantum remainder or UV
completion is assumed by this proposal. A healthy quadratic Minkowski test
would not supply any of those missing results.
