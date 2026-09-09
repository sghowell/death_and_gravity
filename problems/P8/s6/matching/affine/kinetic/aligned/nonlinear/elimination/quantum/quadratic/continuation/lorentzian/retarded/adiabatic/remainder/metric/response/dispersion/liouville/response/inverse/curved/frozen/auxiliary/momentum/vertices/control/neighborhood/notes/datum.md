# A constraint-compatible same-action off-clock family

Use the positive lapse N of the actual clock gauge u=phi.
At phi=0 the spatial hat metric is I, gravitational
momentum vanishes and all spatial Proca data vanish.
The physical spatial metric is N times this hat metric:
e^omega=sqrt(N) here. Its lapse is still N. This point
chart does not change which metric couples to matter.

The original fixed-basepoint primitives have I=Q=0
on this central slice, but

I_phi=-6*N^(-3/2)-18*sqrt(N)+24

must be retained. It is not zero away from N=1.
The full actual central Hamiltonian has already been
matched invariant by invariant in the parent. Writing
its old background part as A(N), the homogeneous
Hamiltonian with canonical matter momentum P is

H_hom=A(N)+(P^2-1/100)/(2*sqrt(N)).

The lapse constraint at FIXED canonical P is solved by

P(N)^2=(69/4+6*epsilon)/N
       -(2849/100-4*epsilon)*N
       +(45/4-10*epsilon)*N^3,
epsilon=10^-6.

Take its positive root. This is an actual constraint
family, not a lapse chosen independently of the field
data. All homogeneous momentum constraints vanish.
At this slice the reconstructed trace and temporal
vector also vanish, since b=c=p=j=0.

The fixed-phase lapse Hessian, evaluated only AFTER
differentiation, is

h_N=[-3*a1-a2*N^2+15*a3*N^4]/(2*N^(7/2)),
(a1,a2,a3)=(23/8+epsilon,2849/200-2*epsilon,9/8-epsilon).

Differentiating H along P(N) instead would compute a
different quantity and would not establish the lapse
pivot. The code checks the fixed-phase expression.

Let gamma_t,gamma_s be the actual retuned mass functions,
delta=N^-2-1, and

kappa=gamma_t-(3/2)*N^2*delta^2.

The joint temporal pivot is -N^(5/2)/kappa; the mixed
lapse/temporal pivot vanishes at this datum. These are
derived with the trace Legendre map retained, not from
an isolated Proca mass sign.

## Continuous small interval

Use |N-1|<=10^-28. For each rational function, expand
its numerator and denominator about one. Bound the
nonconstant terms by their absolute coefficient sums
times the corresponding powers of 10^-28. A positive
denominator lower bound gives a continuous rational
interval, not a grid of sample values.

The explicit intervals give P^2>1/200,
N^(7/2)*h_N<-2, gamma_t>9/10, gamma_s>9/10 and
kappa>9/10. Furthermore,

|P-1/10|=|P^2-1/100|/(P+1/10)
        <=10*|P^2-1/100|<10^-24.

Thus these data lie in the already certified nine-
invariant auxiliary polydisc and its lapse disc.
Only dc=P-1/10 is nonzero among its nine invariants.
The family is punctured at N=1 for the new diagnostic.

The homogeneous reduced Hamiltonian is smooth, the
auxiliary pivots and the trace/matter Legendre map are
nonzero, and N stays positive locally. Its ordinary
finite-dimensional evolution therefore has a local
smooth solution through each datum; the background
time jets used below are finite. This elementary local
existence statement asserts neither a global new bounce
nor a uniform nonlinear spatial Cauchy theorem.
