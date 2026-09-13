# S6.234: separate-model first-loop input

## Model and separate renormalization prescription

Use the unchanged V2S-T1 classical action from S233, with canonical light mass1, D=2lambda/gamma, M_H²=D+2, g=1/8192 and C=-g²(3/D-2/D²). This is not the original affine/DHOST/Proca action.

The perturbative prescription is explicit: dimensional regularization in4-2epsilon; MSbar subtraction of UV poles at mu=1; a heavy-one-point-zero condition; and finite local light mass/kinetic subtractions at s=1. Tree renormalized parameters keep their S233 values. No finite four-point matching condition is imposed in this checkpoint. Counterterms in this separately specified loop prescription do not edit a frozen original finite prescription.

## Exact regulated Gaussian identity

For a finite Euclidean lattice or equivalent specified finite-dimensional regulator, K_H=-Delta+M_H² is positive and G_H=K_H^-1. Integrating H exactly gives the field-independent Gaussian determinant and

S_eff=1/2<phi,K_phi phi>-C integral phi4/24
      -g²<phi²,G_H phi²>/8.

Since G_H<=1/M_H², this action is bounded below by the canonical light quadratic form plus the positive valley quartic from S233. This exact finite-regulator statement is distinct from removing the regulator while preserving a quantum theory.

The full light Hessian is K_phi+W, where

W=-C M_(phi²)/2-g² M_(G_H phi²)/2-g² M_phi G_H M_phi.

Here M_f means multiplication by f. Its local and bilocal terms both matter. The complete formal first-loop light-field Taylor jet through degree4 is

Gamma_one=Tr(G_phi W)/2-Tr(G_phi W G_phi W)/4+counterterms.

This is one-light-particle-irreducible, and includes mixed heavy-light loops. It is not an evaluated renormalized four-point matching amplitude.

## Complete first light two-point coefficient

The heavy one-point counterterm cancels the explicitly retained one-point-reducible heavy source tadpole. The quartic tadpole and mixed bubble remain. With iPi denoting the insertion, the inverse has sign s-1+Pi. At mu1,

Pi_MS(s)=-C/(32pi²)+g² B0_MS(s;1,M_H²)/(16pi²),
B0_MS=-integral_0^1 Log[x M_H²+(1-x)-x(1-x)s-i0]dx.

The finite on-shell coefficient is Pi_OS(s)=Pi_MS(s)-Pi_MS(1)-(s-1)Pi_MS'(1). With F=(1-x)²+xM_H² and alpha=x(1-x)/F,

Pi_OS(s)=g²/(16pi²) integral[-Log(1-alpha(s-1))-alpha(s-1)]dx.

The mixed first-sheet cut begins at(M_H+1)²; its pseudothreshold(M_H-1)² is not an additional first-sheet cut. Higher-loop light cuts are not excluded.

## Uniform finite-order comparison

For abs(s-1)<=R and r=R/M_H²<1,

abs(Pi_OS(s))<=g² abs(s-1)²/[96pi² M_H4(1-r)].

At R=10^196 the actual r<32/625. The analytic quotient q=Pi_OS/(s-1) obeys a rational bound epsilon<10^-209 and epsilon/(1-epsilon)<10^-209. Hence the algebraic reciprocal of the one-loop-truncated inverse has only its anchor simple pole in this disk, unit residue, and a propagator ratio error below10^-209. This does not bound the exact propagator or omitted loops.

The unsubtracted first coefficients separately satisfy abs(Pi_MS(1))<10^-7 and0<Pi_MS'(1)<10^-207. These are calculated finite-order coefficients; the on-shell conditions are not presented as an independent proof of exact physical mass and LSZ.

## Heavy first absorptive coefficient

The heavy bubble has the identical internal factor1/2:

Im Pi_H,first(s+i0)=g² sqrt(1-4/s)/(32pi), s>4.

This agrees with the leading identical two-body decay rate Gamma_H,first=g² sqrt(1-4/M_H²)/(32pi M_H). Exact parameter brackets put Gamma_H,first/M_H between10^-208 and10^-207. Its outgoing formal second-sheet displacement has negative imaginary part. No exact resonance pole/width or stable quantum heavy atom is asserted.

## Boundary

The17 source inputs have exact identities, independently reconstructed diagrams, complex quadrature diagnostics, rational inequalities and written all-domain arguments. They are not FORMALIZED. The full real four-point matching, all-loop error, continuum quantum construction, original covariant bounce/state dictionary and finite-gravity Regge problem remain open. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter.
