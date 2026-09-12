# P8 continuation: actual curved scalar coordinates and causal reference inverse

S6.225 supplies the actual curved scalar curvature coordinates, the complete original finite local Hessian, and a two-channel causal reference inverse on the unchanged bounce slab. The middle nonlocal factors are still the specified flat trace/shear factors. This is not the actual full curved quantum-force inverse, and original P8 remains OPEN.

## Actual curvature, not a frozen flat substitution

On the unchanged proper-time slab[-1/2,1/2], a=(1+t^2)^2 and dt=a d eta. Write h=aH, U=h'+h^2=a''/a, where the primes in this paragraph are conformal-time derivatives. The actual interval bounds are

|h|<=5/2, |h'|<=175/16, U<=275/16<18, a^4<6.

The conformal duration is2/5+atan(1/2)<1. For the original exponential synchronous spatial scalar chartQ=2wI-2cPi, q=P^2, define

S=(D^2+3hD+2q/3)w-(D^2+3hD)c/3=(a^2/6)delta R_old,
W=-qw-D^2c.

Literal curvature calculation checks the scalar identity and every linearized Weyl component;24 of the256 components are nonzero. The curvature-coordinate matrix isBc=B0+C D, withC=[[3h,-h],[0,0]], norm below8 and derivative norm below35. This is a spatial synchronous quotient, not a reconstruction of all lapse/shift/clock constraints.

## The whole original finite local action is retained

Before the original factor64pi^2, the fixed density is

5m^4/2+5m^2 R_old/3-Weyl^2/30+Euler/90-R_old^2/18.

ForTi=3wi-ci and rDG=a^2 delta_D delta_G R_old, the exact mixed density is the curvature-square part-4SD SG-(4/45)WD WG plus

Rrem=-(2/3)U rDG-4U(TD SG+TG SD)-2U^2 TD TG
 +(5m^2/3)a^2[rDG+6(TD SG+TG SD)+6U TD TG]
 +(5m^4/2)a^4 TD TG.

Here

rDG=24wD'wG'-8(wD'cG'+cD'wG')+4cD'cG'
 -20q wD wG+4q(cD wG+cG wD).

The background-curvature, second-metric-variation, Einstein and volume pieces are not discarded. The remainder has total differential orders0,1,2. Differential order alone does not supply a compatible-space norm or smallness. Compact Euler and box-R variations are variational boundaries, not a claim of pointwise zero density or absent flux on arbitrary histories.

An independent nonlinear warped-ADM curvature calculation checks the same full mixed terms. Atq0, w=c/3, division by the tensor norm8/3 gives the S189 fourth coefficient-1/30 and kinetic coefficient5m^2 a^2/6-U/3, with the original prescription unchanged.

## Both causal curvature-coordinate inverses

LetY=Bc^-1 andZ=(Bc*)^-1 be the zero-past inverses. The formal adjoint is

Bc*=B0^T-C^T D-(C')^T.

The derivative-of-coefficient term cannot be omitted. The complete Green boundary isz^T L y'-z'^T L y+z^T C y. The primal velocity Volterra equation gives

||partial_t Y(t,s)||<=(5/2)exp(20(t-s)),
||Y(t,s)||<=(5/2)(t-s)exp(20(t-s)).

The separately derived adjoint equation gives the corresponding bounds with exponent108 forZ and its output-time derivative. These are uniform in every spatial momentum; no inverse-transfer division is introduced.

The reference operator is

Aref=Bc* Fdiag,flat Bc,
Fdiag=diag(Ftrace,(8/3)F2).

Its inverse has the ordered formY Kdiag Z. Using the source-pinned bounded spectral primitiveJdiag andZ(s,s)=0 gives the ordinary double kernel

E(t,s)=integral_s^t du integral_s^u dv
       Y(t,u) Jdiag(u-v) partial_v Z(v,s).

The bound is

||E(t,s)||<=375(t-s)^3 exp(108(t-s))/16,
||E||<=375T^4 exp(108T)/64<6T^4 exp(108T)

on the statedC_tH^r andL2_tH^r source spaces, allq and every realr, without spatial derivative loss. The original shear pole and complete causal initial boundary remain.

For clarity, the composition uses its compatible graph: s=Yx withx andy=Fdiag x in the stated source space andBc*y an ordinary source there, all as causal distributions. Existence follows fromy=Zf,x=Kdiag y; uniqueness follows in reverse order. It is not an assertion that an informal product is defined on every arbitrary metric distribution, or that the full S222 graph is dense or invariant.

The physical conformal force reference is output multiplication by1/(64pi^2 kappa a^4) followed byAref. Its inverse retains RIGHT multiplication by64pi^2 kappa a^4, with bound2250pi^2 kappa T^4 exp(108T). The amplitude-to-metric factors from S224 also remain. This is not a kappa-small feedback or stability estimate.

## Independent checks and completed fresh validation

The independent tests cover literal nonlinear scalar/Weyl geometry, the entire finite local action, both formal adjoints, forced zero/nonzero-transfer coordinate equations and noncommuting finite causal matrix products. An early private probe needed expansion before comparing an exact coefficient. A finite adjoint fixture atq0 was accidentally singular for its arbitrary coarse coefficient; changing only that diagnostic coefficient to a nonsingular value fixed the fixture. Neither the actual background nor the continuum formulas changed.

Private integrated science passed245 tests. Preflight checked17 inputs and20 fields in82.52 seconds. Final private science passed245 in83.57 seconds, and repository science245 in83.95 seconds. Formatting and lint passed before freezing.

Fresh original-SymPy ordinary replay passed270 tests in2304.28 seconds. The standalone CLI replay passed. Complete P8 regression passed45768 tests in4450.38 seconds, exit code0. Its703-file snapshot SHA is

`7f310cea0543a4f05526105cdbff6e2c014791b19542136f2dc37796662e887f`.

The full launcher retained the frozen S219 helper-directory PYTHONPATH allowance. Only full regression used the audited exact-GCD adapter:128 original tuple comparisons passed, with final counters36075 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

The49340-character native report was transferred in five checked chunks. Its20 fields contain38 named identities,55 scalar entries,24 gates,nine controls and169 rejected inputs. All17 frozen input hashes and the native report SHA were verified:

`5c4a3649bc435475d6bfd65901af47b63da2543faa9bc9bcd6364e00fac9a02d`.

The exact21-file publication manifest contains the17 scientific inputs, report, this assessment, CLAIMS and README. Pending weighted/finite-local successors and unrelated P4/P9 edits are excluded. Written causal and functional proofs are not FORMALIZED.

## Remaining work

Actual nonlocal curved mass/state/contact matching, incorporation of the full finite local remainder in the inverse, classical/matter coupling and the complete S222 quantum-force graph still require further work beyond this checkpoint. Weighted and finite-local successors are under fresh verification. A new actual homogeneous nonlocal shear calculation is in private development. None of these references alone supplies physical stability, finite-coupling/nonlinear common-parent control, a cutoff, original V/G/B or P8 closure. All original primitive statuses remain unchanged; there is no user-intervention blocker.
