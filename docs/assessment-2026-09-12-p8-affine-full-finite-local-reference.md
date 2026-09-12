# P8 continuation: complete original finite-local curved reference

S6.227 bounds the COMPLETE original finite-local remainder on the same weighted space as the two-channel reference inverse. It closes the finite-local comparison step, not the actual curved nonlocal scalar/clock/matter inverse, stability, nonlinear control, or original P8. Original P8 remains OPEN.

## Exact local factorization and retained density boundary

In the actual S225 conformal curvature coordinates, the original finite remainder is

Lrem=alpha rDG+6alpha(TD SG+TG SD)+gamma TD TG,
alpha=-2U/3+(5m^2/3)a^2,
gamma=-2U^2+10m^2a^2U+(5m^4/2)a^4,
T=3w-c.

Its COMPLETE Euler Hessian, with D*= -D in conformal Lebesgue density, is

Rloc=Bc* A+A^T Bc+D* H D+D* J+J^T D+V0,

A=alpha[[18,-6],[10,-4]],
H=alpha[[24,-18],[-18,12]],
J=-alpha'[[0,0],[10,-4]],
V0=gamma(3,-1)^T(3,-1).

The direct action identity retains the nonzero boundary

D{alpha[cD'(10wG-4cG)+(10wD-4cD)cG']}.

Every finite-P gradient contribution is included in Bc; no q c^2 term is invented. Coefficients are not commuted through Green operators. The coefficient bounds use the unchanged mass1000 and actual scale factor, with a^2<3, a^4<6, U<18 and |U'|<112. They imply ||A||<110m^2, ||H||<185m^2, ||J||<286m^2 and ||V0||<160m^4.

## The required source-time adjoint bound

Let Y=Bc^-1_ret and Z=(Bc*)^-1_ret. Green reciprocity gives Z(t,s)=Y_adv(s,t)^T. Time reversal of the primal equation retains its coefficient bound and gives

||Z(t,s)||<=(5/2)(t-s)exp(20(t-s)),
||partial_s Z(t,s)||<=(5/2)exp(20(t-s)).

This is a SOURCE-time derivative estimate; it is not a stronger output-time adjoint estimate. The original initial delta in D* cancels the lower integration-by-parts boundary. Consequently ZD* is the ordinary kernel partial_s Z. Omitting that delta gives the wrong result on a source with nonzero initial value; the independent boundary diagnostic detects it.

Writing g=sigma-20>0, the complete ordered conjugate is

Omega=Z Rloc Y
 =A Y+Z A^T+(ZD*)H(DY)+(ZD*)JY+ZJ^T(DY)+ZV0Y.

Weighted Young bounds on L2_tH^r, uniformly in spatial momentum, give

||Omega||<=6825m^2/(4g^2)+3575m^2/g^3+1000m^4/g^4.

At the original S226 weight sigma_M=2m exp(32+8M), g>=100m and

||Omega||<=6825543/40000000<1/4.

This is a same-space bound on the full local remainder, not an inference from its differential order.

## Ordered finite-local inverse and precise graph

For an optional stated bounded channel multiplication V with ||V||<=M,

Afinite=Bc*(Fdiag+V)Bc+Rloc
       =Bc*(Fdiag+V+Omega)Bc.

Omega is a bounded causal OPERATOR, not a time-multiplication matrix. The original channel inverse gives

||Kdiag(V+Omega)||<=5(M+1/4)/(32+13M)<5/13<1/2.

The middle inverse has norm at most20/(123+32M). The composed inverse Y Kfinite Z has norm at most

125/[(123+32M)(sigma_M-20)^4].

The physical inverse retains RIGHT multiplication by64pi^2 kappa a^4, giving

48000pi^2 kappa/[(123+32M)(sigma_M-20)^4].

Removing the weight costs exp(sigma_M T). Neither this factor nor kappa is suppressed, so no physical smallness or stability follows.

The compatible composed graph is explicit: s=Yx, x belongs to the middle graph, y=(Fdiag+V+Omega)x lies in the weighted source space, and Bc*y is also an ordinary source there. Then Bs=x, so merely measurable V never multiplies an arbitrary distribution. For a source f, first y=Zf, then x=Kfinite y, then s=Yx. The reverse implication follows from prepared uniqueness in the three ordered factors. This is not a maximal metric-distribution domain or the unrestricted S222 graph. The external S225/S226 assessments use the same precise graph boundary; their frozen inputs remain unchanged.

## Independent checks and completed fresh validation

Exact independent tests cover the arbitrary-coefficient full Euler Hessian, density boundary, actual coefficient estimates, source-time Green reciprocity, initial delta, all six ordered conjugate terms and both inverse products. Reciprocity is checked by exact constant-coefficient kernels and nonconstant sixth-order jets. Finite causal block solves at q=0,7/5,10000 detect deleting J, deleting Rloc or reversing the Neumann order. These are diagnostics, not numerical continuum certificates.

An early private comparison required expansion of a gradient derivative before exact cancellation. No physical formula changed. Formatting and lint passed before freezing. Preflight verified17 source inputs and20 report fields in83.71 seconds. Final private science passed253 tests in83.31 seconds; repository science passed253 in86.80 seconds.

Fresh original-SymPy ordinary replay passed278 tests in2308.60 seconds. Standalone CLI replay passed. Complete P8 regression passed46336 tests in4503.42 seconds, exit code0. Its707-file snapshot SHA is

`59d0ad83ca01e2cfb90ee164e2c66b577542c8d4d9006bac343fc5eada7dd496`.

The full launcher retained the frozen S219 helper-directory PYTHONPATH allowance. Only full regression used the audited exact-GCD adapter:128 original tuple comparisons passed, with final counters36087 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

The39166-character native report was transferred in four checked chunks. It contains20 fields,31 named identities,33 scalar entries,24 gates,nine controls and184 rejected inputs. Its SHA is

`48c4549159564a4ece8cfd0eb274b24d5e76296f8bd4870b1ccfc2c1ef0d6347`.

The17 frozen source hashes and native report hash are rechecked before exact21-file publication:17 inputs, report, this assessment, CLAIMS and README. The homogeneous successor, private coupled research and unrelated P4/P9 changes are excluded. The analytic arguments are written proofs, not FORMALIZED.

## Remaining work

The actual homogeneous quantum-plus-classical tensor successor has passed ordinary and CLI replay and is awaiting full regression. The private finite-transfer scalar/clock/matter track is developing an actual smooth bounded-external-momentum inverse. Neither this finite-local result nor those pending results give unrestricted-space stability, nonlinear/finite-coupling common-parent control, a physical cutoff, or original V/G/B/P8 closure. No user intervention is currently required.
