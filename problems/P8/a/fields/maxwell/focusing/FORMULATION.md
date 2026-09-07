# P8-A.17 — photon-QSEI incompleteness with geometric contraction history

This proves a conditional incompleteness theorem for an actual free-Maxwell
semiclassical model in a stated global FLRW spacetime. It assumes no initial
pointwise SEC/Ricci sign. A short contraction history and future geometric
bounds conditional on extension replace that premise; they remain explicit
hypotheses, not consequences of the QSEI or observations.

## Global model and exact quantifiers

Let M=I_t times R^3, with I_t an open interval containing 0, be the stipulated
spacetime with smooth positive scale factor and metric

    ds^2=dt^2-a(t)^2 dx^2, H=a'/a.

This is a global spacetime description, not merely a coordinate chart inside
an unspecified larger spacetime. Let one physical free Maxwell field be in
any Hadamard state on M, in the explicit A.16 prescription with real constant
beta_M. Neither state homogeneity nor a quasifree/zero-mean condition is added.
Assume the actual semiclassical equation

    G_FK+Lambda g_FK=-kappa(T_Maxwell+T_other), kappa>0,
    R_UU=Lambda-kappa(E_Maxwell+E_other).

The extra source, if present, obeys the explicitly specified E_other>=ell
on the relevant comoving segments. Lambda and Einstein/Newton normalization
are fixed separately. Extra left-hand curvature tensors must be included
with their correctly signed effective source before using this dictionary.
Setting T_other=0 and Lambda=0 is an available pure-Maxwell specialization,
not silently assumed by the generic theorem.

Fix tau>0, r>0 and h>0. The past interval [-r tau,0], including a smooth
neighborhood of its endpoints, must ACTUALLY exist in M and satisfy

    H(t)<=-h/tau,
    |H^(j)(t)|<=d_j/tau^(j+1), j=0,1,2,3.

Future caps have a different quantifier: IF a comoving normal reaches proper
time tau (equivalently tau belongs to I_t), assume on that segment

    |H^(j)(t)|<=c_j/(tau-t)^(j+1), 0<=t<tau, j=0,1,2,3.

All d_j,c_j are finite nonnegative constants. The future caps are conditional
on that extension; no extension of an incomplete curve or metric through a
missing endpoint is assumed. They are not asserted on a shorter curve that
already ends before tau. This distinction must be retained in applications.
The caps are not a uniform absolute-curvature bound: their right sides may
diverge at tau. Conversely their validity on a hypothetical long extension
does not follow from a verified initial patch or from A.16 alone.

Define

    delta=kappa hbar/(8pi^2 tau^2),
    sigma=tau^2[max(Lambda,0)+kappa max(-ell,0)].

Upper bounds on these quantities may be used. All finite beta_M dependence
is retained; zero type D is still the separately named beta_M=0 choice.

## Explicit sufficient inequality and conclusion

Let V(c)=186(c0^4+2c0^2 c1)
 +|beta_M|(18c3+90c1^2+90c0 c2+108c0^2 c1), and similarly V(d). Set

    q_d=3d0^2/4+3d1/2, q_c=3c0^2/4+3c1/2,
    C_P=[7/2+(11/5)d0 r+(2/3)q_d r^2]^2/r^3+13r V(d)/12600,
    C_F=[(7/2)(1+2c0)+(25/12)q_c]^2+13 V(c)/1080,
    M=3h+39h^2r/35-18/5-delta(C_P+C_F)-13(1+r)sigma/35.

If M>=0, no comoving normal can reach tau. In particular the upper endpoint
of I_t is at most tau, no future-directed timelike curve from t=0 has proper
length greater than tau, and the stipulated spacetime is future timelike
geodesically incomplete. For the calibration below M is strictly positive.

This is an incompleteness theorem, not a proof of curvature blow-up, absence
of any larger extension, a fundamental singularity, or EFT validity at the
endpoint. Reversing the physical time orientation gives the corresponding
past-incompleteness statement with the reversed contraction history.

## Exact macroscopic-scale conditional calibration

Take

    r=1/100, h=3/2, d=c=(2,4,16,96), |beta_M|<=1,
    delta<=10^-8.

At the worst finite-beta cap,

    V=16704,
    C_P=54968328329/4375, C_F=363631/240,
    C_P+C_F=2638797936917/210000 <13000000.

With sigma=0 the margin is
2398243151869/3000000000000 >3/4. Even sigma<=1 leaves margin >1/3.
For kappa=8piG, delta=ell_P^2/(pi tau^2), ell_P^2=G hbar (c=1).
Thus the quantum requirement is compatible with cosmological tau by an
enormous separation of scales. No observed curvature/history/source values
or actual cosmological solution are inferred from that dimensional fact.
These caps are a geometric calibration, not a radiation witness: for example
classical radiation with |H|tau=2 has |Hdot|tau^2=8, not 4.

## Nonvacuity and exclusions

An explicitly defined mollified-quintic metric is smooth, obeys these
geometric history/cap assumptions, violates the pointwise SEC, and is future
timelike and null complete. Therefore the geometric assumptions alone do
not imply the conclusion. It is deliberately NOT an allowed Maxwell SEE
solution with the stated small quantum/source constants; the independent
source-mismatch control makes this distinction quantitative.

The theorem's analytic field input is the pinned A.16 all-Hadamard Maxwell
QSEI, not an arbitrary-field energy axiom. The geometric countercontrol is
not a state/solution witness. No interacting-QED result, all-spacetime
renormalization independence, observed-universe application, theorem of
fundamental EFT validity, or completion of original P8(a)/P8 is claimed.
Machine certificates verify the identity chain and exact constants; the
Hadamard, Sobolev and geometric implications have written proofs.
