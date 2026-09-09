# Vary the actual Euler equations before freezing

The quadratic Hamiltonian has real matrices A=A^T,
C=C^T and B in the notation of spatial.md. Its Legendre
transform is

L=hat_a^3*[Q'^T alpha Q'/2+Q'^T beta Q
           -Q^T gamma Q/2],
alpha=A^-1, beta=-A^-1 B, gamma=C-B^T A^-1 B.

Here gamma names the sign convention, not a positive
matrix. The ACTUAL Euler equation is

alpha Q''+(alpha'+3*H_hat*alpha+beta-beta^T)Q'
 +(gamma+beta'+3*H_hat*beta)Q=0,

where H_hat=hat_a'/hat_a in the clock time coordinate.
These derivatives follow the actual local homogeneous
solution. Freezing H2 before taking them would give a
different diagnostic.

## Why the time-jet orders are controlled

About a nearby homogeneous isotropic state with zero
spatial vector background, the full canonical ADM
reduction has the following momentum orders. Linear
spatial curvature is O(q*v); the lapse Schur term is
the square of one linear force. The York solution has
one inverse gradient. Longitudinal vector coordinates
enter through q*sigma^2 and momenta through P_sigma,
with their electric correction O(P_sigma^2/q).

Consequently the original v-v Hessian is at most O(q^2),
v-p and v-P_chi/v-P_sigma entries are at most O(q),
and the momentum-only entries are O(1), except for the
explicit decaying electric correction. After the gamma
swap, A is O(1). B's b column is at most O(q), and its
chi/sigma columns are O(1). All these coefficients are
rational in q with smooth local background coefficients.
Their first time derivatives have the SAME momentum
orders: the background jets are finite and differentiating
q gives -2*H_hat*q, not an extra power of q.

At the present datum A_infinity is invertible. Thus
alpha and alpha' are O(1). In particular beta' and
3*H_hat*beta are O(q) only in their first column and
O(1) in their other columns. This retains arbitrary
finite actual jets; it does not bound them by zero or
assert a numerical uniform time-evolution estimate.

At the central datum B has only B_sigma,b=d*q and A
has no scalar/vector kinetic off-diagonal entry.
Therefore the exact Legendre correction is

gamma_bb=-d^2*q^2/A_sigmasigma(q).

No time derivative above can cancel this q^2 term.
The lower two-by-two coordinate block has positive
leading diagonal q*(N^(3/2),N^(3/2)/gamma_s); its
off-diagonal entries have lower order. The b-chi
coordinate coupling is at most O(q).

## An exact asymptotic sign

Let E0(q)=gamma+beta'+3*H_hat*beta be the constant
coefficient of the Euler-first frozen operator.
The determinant therefore satisfies

lim_(q->infinity) det(E0(q))/q^4
 =-d^2*N^3/[gamma_s*(kappa/sqrt(N))]
 =-(1-N^2)^2*N^(5/2)/(gamma_s*kappa).

The code independently checks this determinant identity
with arbitrary real coefficients in EVERY allowed O(q)
and O(1) remainder slot. It also matches d, the vector
pivot, both gradient coefficients and the absent bare
b^2 term to the fully spatially reduced action.

For each fixed real 0<|N-1|<=10^-28 the displayed limit
is strictly negative. Thus det(E0(q))<0 for sufficiently
large q. For such q, alpha is positive definite. For
real lambda define the frozen matrix polynomial

D_q(lambda)=alpha*lambda^2
 +(alpha'+3*H_hat*alpha+beta-beta^T)*lambda+E0(q).

It has real coefficients, det(D_q(0))<0 and positive
leading coefficient det(alpha)*lambda^6. The intermediate
value theorem gives at least one positive real root.
This is a genuine growing solution of the FROZEN local
three-scalar system. No assertion of simplicity or of
a projected response residue is required.

The physical momentum rescaling q_physical=q/N and
the positive physical lapse do not alter the sign
argument. It is not necessary to interpret the
isolated vector mass ratio as a physical coupled cone.

## Negative control at the original clock trajectory

At N=1 the new q^4 coefficient vanishes. Lower orders
then matter and cannot be discarded. If one falsely
freezes the Hamiltonian before varying, the constant
matrix would have

lim det(gamma)/q^3=-1/100.

The actual original Euler-first scalar block from
S6.78 instead gives, including the vector's unit
leading coefficient,

lim det(E0_clock)/q^3=1199/100>0.

Both are exact checks. The second value retains the
actual original time jets and agrees with the positive
coupled scalar energy construction. The negative
off-clock conclusion does not come from reusing the
false on-clock freeze.
