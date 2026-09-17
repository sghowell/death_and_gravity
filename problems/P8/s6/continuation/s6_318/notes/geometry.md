# Weighted geometry and conserved temporal inversion

## Statement and domain

Use S317's complete temporal pure-soft recursion, unchanged canonical
Einstein vertices, physical future null free rays of energies w_i>0,
and unit-Frobenius spatial TT leaves. Work away from exact internal
propagator poles. For a nonempty subset S let W=sum_(i in S)w_i,
Q=(W,Q_sp), v=Q_sp/W and delta^2=1-|v|^2=Q^2/W^2.

For a temporal symmetric current H write A=H_sp and define

N(H,Q)=max( ||A||F, ||A*v||2/delta, |v^T*A*v|/delta^2 ).

At a singleton, delta=0 and transversality makes the last two numerators
zero; define the seed norm as the TT Frobenius norm. At generic larger
subsets delta>0.

The uniform bound proved with the other notes for n=|S|:

N(H_S,Q_S) <= C_n * W^n/[product_i w_i] / kappa^((n-1)/2),

C_1=1,
C_n=64*sum_(set partitions pi of S, k=|pi|>=2)
      (k+1)!*32^(k+1)*13^k*product_(A in pi) C_|A|.

An explicit envelope is

C_n <=2*(2*10^10)^(n-1)*n!   for n>=1.

This controls every angular configuration and nested collinear approach
on the generic domain, with explicit soft-energy factors. It is NOT
energy-independent boundedness and does not assign values to excluded
exactly-collinear propagator points. Leaf norms multiply the RHS if
the leaves are not normalized.


## 1. Parent geometry and transfer of a child norm

Choose a unit spatial n along Q_sp, with any unit n when Q_sp=0.
Then Q_sp=rho*W*n, 0<=rho<=1, delta^2=1-rho^2.

For any child subset A of energy W_A and velocity v_A=Q_A,sp/W_A,
future-null summation gives

W_A-n.Q_A,sp
 =sum_(i in A)w_i*(1-n.n_i)
 <=(1-rho)*W <=delta^2*W.

Cauchy-Schwarz gives
|Q_A,perp|^2<=2*W_A*(W_A-n.Q_A,sp)<=2*delta^2*W_A*W.
Also

delta_A^2=1-|v_A|^2 <=2*(1-n.v_A)<=2*delta^2*W/W_A,
|n-v_A|^2=2*(1-n.v_A)-delta_A^2<=2*delta^2*W/W_A.

Let R=W/W_A>=1 and M=N(H_A,Q_A). By the norm definition,

||H_A*n|| <=||H_A*v_A||+||H_A||*|n-v_A|
           <=2*sqrt(2)*delta*sqrt(R)*M,

|n^T H_A n|
 <=|v_A^T H_A v_A|+2|n-v_A|*||H_A*v_A||
   +|n-v_A|^2*||H_A||
 <=8*delta^2*R*M.

These statements also hold for TT singleton seeds, where delta_A=0.
They avoid dividing by |Q_A,sp|, so no direction singularity is hidden
when a child's spatial momentum vanishes.

With P=I-n*n^T, decompose the child's spatial tensor exactly as

H_A = A_A+delta*B_A+delta^2*C_A,

A_A=P H_A P,
delta*B_A=n*(P H_A n)^T+(P H_A n)*n^T,
delta^2*C_A=(n^T H_A n)*n*n^T.

The coefficient norms satisfy

||A_A||F<=M,
||B_A||F<=4*sqrt(R)*M,
||C_A||F<=8*R*M,

hence their coefficient-l1 norm is <=13*R*M.
The first coefficient is purely transverse, the second has one
longitudinal spatial index, and the third has two.

## 2. Graded momentum decomposition

Relative to the same parent n and null vector k=(1,n), every child
momentum has the exact form

Q_A = W_A*k+delta*b_A+delta^2*c_A,

where b_A is purely transverse spatial and c_A purely longitudinal
spatial. The bounds above imply

||W_A*k||<=sqrt(2)*W_A,
||b_A||<=sqrt(2*W_A*W),
||c_A||<=W.

Their coefficient-l1 sum is <=(2*sqrt(2)+1)*W<4W.
The negative root momentum has the same bound. Momentum conservation
holds COEFFICIENT BY COEFFICIENT in the formal grading variable:
transverse child coefficients sum to zero and the root longitudinal
coefficient is the negative sum of the children's coefficients.


## 4. Conserved temporal propagator closes the weighted norm

S317's finite Noether induction supplies conservation of the COMPLETE
amputated source. Align n with z. The exact conserved-root temporal
propagator identity gives

Hxx=-(Rxx-Ryy)/(2W^2 delta^2)+Rzz/(2W^2),
Hyy= +(Rxx-Ryy)/(2W^2 delta^2)+Rzz/(2W^2),
Hxy=-Rxy/(W^2 delta^2),
Hxz=-Rxz/W^2, Hyz=-Ryz/W^2,
Hzz=(Rxx+Ryy-delta^2*Rzz)/(2W^2).

After extracting B/W^2, the component bounds are
|Hxx|,|Hyy|<=3/2, |Hxy|<=1,
|Hxz|,|Hyz|<=delta, |Hzz|<=3*delta^2/2.

Thus ||H||F^2 <=(51/4)*(B/W^2)^2 <16*(B/W^2)^2.
Also ||H*v||^2/delta^2 <=(17/4)*(B/W^2)^2,
and |v^T H v|/delta^2 <=3B/(2W^2).
The common weighted norm is therefore <=4B/W^2.

This cancels the root angular denominator with explicit uniform
constants. No separate off-shell descendant is assumed conserved:
only the completed root uses the Noether identity.
