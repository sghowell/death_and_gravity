# S6.231: original scalar unitarity and conditional dispersion matching

## Fixed input and observables

Use S177's complete anchored canonical family at the original S182 retained vacuum. The S182 nonconstant retuning has no field jets below2048 and therefore changes no quartic amplitude; its constant vacuum term is retained separately. The independent full functions, higher vertices, affine source and source contacts remain. The scalar mass squared is1, distinct from the Proca mass1000. The fixed quartic coefficients are lambda=10^-600 and gamma=1024*10^-800.

For S=s>=4 and x=cos(theta), t=-(S-4)(1-x)/2 and u=-(S-4)(1+x)/2, the original amplitude is

A=2lambda[(s-2)²+(t-2)²+(u-2)²]+3gamma stu-8gamma=a+b x².

With beta=sqrt(1-4/S), define normalized identical-channel partial waves by t_l=beta integral(A P_l)/(64pi), so S_l=1+2i t_l. This convention includes both identical two-body state normalizations. It is not the distinguishable16pi convention applied to the same labelled amplitude.

## Exact original-tree result

Only l0 and l2 occur. Their complete formulas, the elastic squared amplitude and first absorptive coefficient are proved in notes/tree.md and notes/normalization.md. For the fixed coefficients, t0>0 for S>4, t0>=5 abs(t2), and t0 is strictly increasing from zero to infinity. It crosses1/2 once at E_U in(10^133,2*10^133). At E=10^134 it exceeds50000.

An existing exact unitary normalized channel has abs(Re t0_exact)<=1/2. Its distance from the original tree therefore exceeds(1-10^-5)t0_tree at the latter energy. This bound does not suppose corrections are small. E_U is a nominal real-part ceiling, not a theorem that the real tree satisfies exact unitarity below it.

## Conditional physical-matching result

The separately explicit hypotheses are:

1. An actual nongravitational scalar S matrix with physical scalar pole mass1, canonical positive LSZ normalization and the required crossing and positive unitarity exists.
2. The specified pole-subtracted forward amplitude is analytic at crossing point s=2. No unresolved cut crosses that neighborhood; in particular an interacting massless forward cut is not retained.
3. Its relevant twice-subtracted dispersion relation has the specified pole accounting and vanishing infinity arc. The positive absorptive representation therefore applies to physical b2=B''(0)/2.
4. For every S in[E²/2,E²], the norm of the FULL angular amplitude error is at most eta times the original tree L2 norm, where0<=eta<1 and E²/2>=16.

None of these physical matching or full-amplitude hypotheses is proved here. In particular tree mass and residue do not prove their quantum values. All loops, higher operators, heavy states, thresholds and finite matching effects are included in the full error.

Under those hypotheses, b2>=45(1-eta)² gamma² E^8/(131072pi²), hence the strict scaled inequality in README. An additional tolerance abs(b2-4lambda)<=4lambda delta with delta>=0 therefore requires

1+delta > 9(1-eta)²(E/10^125)^8.

For E10^125, eta<=1/2 and delta<=1 this is impossible. It is a conditional error tradeoff, not a UV-parent exclusion.

## Boundary and verification

The ordinary Einstein pole is not discarded to manufacture a finite-gravity partial wave. Original quantum decoupling, physical cutoff, nonlinear/bounce control, unrestricted curved response and V/G/B/P8 remain open. The tensor comparison does not change the original response equation, poles or contour.

The six source modules, two tests, two root documents and seven notes form17 immutable inputs. Native, direct, ordinary and CLI replay use original SymPy. Only full P8 regression may use the audited exact-GCD adapter. Exact tests are distinct from continuum optical/dispersion reasoning and its explicit premises; these written arguments are not FORMALIZED.
