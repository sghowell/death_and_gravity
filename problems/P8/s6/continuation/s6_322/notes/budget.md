## Uniform O(W) remainder budget

Use operator inputs ||A||,||B||<=1 and B physical TT. In the stated
tube, |u'|<2u,|c'|<2W, |N_A|<8,|N_B|<4,
|X_B|<16u,|Y_A|<32|c'| and |Z|<4u|c'|.
Also |A_d|>u/8, |B_d/c'|>1/8,
|A_d+B_d+Z|>W/8. Multiplying the three ordering corrections by
|u'c'| gives at most327680W<10^6W.

Replacing u'*N_A/A_d by u'*p.A.p/(p.Q_A) costs at most8192u:
|N_A-p.A.p|<=14u and |Q_A^2/2|<=4u^2 yield the conservative
coefficient224+4096=4320<8192. Multiplication by the physical
singleton current, whose normalized bound is32, costs<10^6W.
No such replacement is made without retaining its error.

For247 regular graphs, the scaled light-propagator factor is<=16W:
the possible cut is Q_A, q_c or total Q, or there is no light cut.
The S312 arbitrary-tensor vertex budgets apply unchanged. With
R3=3*138240*30000000, R4=3*10616832*30000000, coefficients are

 regular_m=(33*4+177*8)*1024^2*16,
 regular_g=177*1024^4*16*1800000*max(1,R3,R3^2,R4).

They multiply g^2 W/(n*kappa) and W/[kappa^2*(delta+W^2)].
Divide by Am>4g^2/n^3 or AG>8/(kappa*delta), respectively.

For the double-external leading term, assign radiation to its four
hard legs and group LL,LR,RL,RR in each hard channel. The unshifted
Born paired currents are<=2048sqrt(tau), and their recoil changes
are<=16384W. Thus their product is<=2*10^10*(tau+W^2), and its
difference from the Born product is<=10^10*W*(sqrt(tau)+W).
The stress numerator bounds remain |N|<2*10^5, assigned change
<10^6W and Born change<10^7W. Hard-denominator differences remain
<200W*(sqrt(tau)+W). Telescoping numerator, product and inverse
gives the same S313 forward cancellation, with these enlarged
constants. In particular delta*(sqrt(tau)+W)/(tau+W^2)<=2.

Use140 as a safe overcount for ordering/off-shell and assigned-stress
errors. Contact currents are bounded by256 with change32768W;
two-leg heavy-channel currents by128 with change16384W.
The phase has |sqrt(rho)|<2 and |sqrt(rho)-1|<16W.
The independent arithmetic script records each resulting term and
finds final coefficients

 matter:13632146432,
 gravity:1952093286852638820978919582900252803544196317184/9.

They are strictly below10^40 and10^60. Thus the operator
remainder estimate is

 ||K2-K20|| < B2*W, B2=10^40*n^2+10^60,

where K20(A,B) is the full Born product of the two normalized soft
currents. K20 depends on Q_A/u' and n_c but not c. These are derived
conservative envelopes, not assigned unknown physical matching data.
The integer check does not establish the topology or uniform tube
by itself; those are the preceding written arguments.

## Cauchy, current product, compatible faces

For E2=K2-K20, the anisotropic radii imply
 ||partial_c E2||<=B2/eta,
 ||partial_ac E2||,||partial_bc E2||<=B2/(eta^2*u),
 ||partial_abc E2||<=B2/(eta^3*u^2).
Combine these with the four pair-current terms in the product rule:

 |partial_abc G387|
 <=B2/kappa^(3/2) *
 [3104/eta+2*1566/eta^2+530/eta^3]/u.

Exact arithmetic at unchanged n=10^200/512+2,kappa=10^800 puts
this coefficient below10^-723. The compatible-face rectangle is bounded by10^-723*c*I(a,b).

For u>0 the hard operator extends across the c face and individual
a/b faces. At u=0 the pair current is O(u), while K20 and E2 stay
bounded, giving a common zero axis limit. The origin is likewise
zero. This supplies the same trimmed-cube FTC argument as S321.

Summing the three pair choices gives10^-723*J, not three times
that bound, because each contributes the corresponding c*I(a,b).
Together with the188 class, the1349 nonsingleton temporal terms
have rectangle<2*10^-723*J. The3767 three-singleton
hard-core terms, the complete probability subtraction, real-virtual
matching, all-N summation, quantum state, Regge and original P8
closure would STILL remain open.


The regular budgets deliberately overcount13,117,117 graphs by the
older33,177,177 inventories. The source vertex estimates allow an
arbitrary complex spatial first tensor of unit Frobenius norm.
The140 factor likewise overcounts the ordering and assigned-stress
corrections. The n^2 and constant coefficients are derived from the
original positive matter and Einstein Born lower bounds; neither is
a free matching datum.

No abstract uniform estimate is inferred just from a passing integer
inequality. The topology, analytic gaps, complex Lipschitz comparison,
off-shell ordering correction and forward grouping above supply the
proof hypotheses; the exact integers audit its conservative budget.

Integrate the mixed derivative first on an epsilon-trimmed rectangle.
The compatible faces justify epsilon->0. Since the kernel1/(a+b)
is locally integrable, the result is c*I(a,b). The retained S321
bounds integral(J/(abc))<=24*x^2 and integral(J^2/(abc))<=18*x^4
on0<a,b,c<x remain available. These partial-class integrals do not
subtract the complete probability or provide its real-virtual match.

