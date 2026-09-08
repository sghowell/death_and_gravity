# All-equation and physical-constraint interfaces

The source-hashed wrapper rebuilds S6.38 and, through it, the full S6.37
action, the original CD/M1 input and the adopted P8 target/contract.
No ancestor source or report is edited. Independent geometry and
dynamics modules reconstruct both trace maps and every original
quotient block; these interfaces compare the resulting formulas.

## Eight retained traces, without an omitted connection equation

At quadratic order on the actual rolling solution, let M be the
invertible 60-dimensional projective-quotient Hessian, E its embedding
in the unrestricted 64 components, and N the eight-trace map on that
quotient. The literal all-block calculation gives W=M^-1 N^T and
D=N W, with D=(3/8)[[1,-7],[-7,1]] tensor diag(1,-1,-1,-1).
The wrapper checks both M W=N^T and M_full E W=T_full^T, the latter
retaining all 64 Euler equations.

D is invertible. Define R=W D^-1 and P=I-R N. Exact checks give

    N R=I8, N P=0, P²=P, trace(P)=52,
    R^T M R=D^-1, P^T M R=0.

Thus every quotient perturbation decomposes into its eight fixed trace
components and a 52-dimensional complementary kernel. This kernel has
nondegenerate restricted M: a vector in ker N orthogonal to ker N is
W l, and N W l=D l=0 forces l=0. An independent test constructs a
full nullspace basis and verifies its 52-by-52 restricted determinant.
The exact quadratic reduced mass action is consequently
(C-Cstar)^T D^-1(C-Cstar)/2, with a plus sign. This proves the rank-two
action premise without a selected connection ansatz. It does not assert
a nonlinear inverse on an open field tube for the new kinetic family.

For rank-one null-Schur choices this eight-trace inverse is not used to
discard a remaining multiplier: the separate four-trace multiplier
action is retained until it forces T=Tstar. The actual all-64 Euler
identity is in geometry.action_identity(); the independent small
matrix in dynamics.null_schur() is explicitly only a control.

## Full scalar embedding and time-dependent null constraints

The wrapper places a general rank-two curl form, with arbitrary time
and spatial lapse-source shifts, next to the complete original CD/M1
scalar action before eliminating either temporal vector, lapse or
shift. Exact algebra shows that the shift equation is unchanged and
the full auxiliary determinant is

    -4 q² Theta² det(D2^-1+q Z).

On the stated regular punctured chart all four auxiliaries can be
reconstructed. The pure two-vector velocity submatrix of the full
reduced action is exactly (D2+Z^-1/q)^-1/2. This establishes the physical
embedding required by the rank-two inertia proof; a negative value in
this restriction is a negative value of the full velocity form.

There is also a direct time-dependent constraint check for the null
branch. Write the unreduced action as

    L0(v,v_dot,s,s_dot,n,n_dot;u)
      +2 a³ q b(Theta*n-v_dot-ell*s/2).

Let E_n0 denote the lapse Euler derivative of L0. The shift fixes
n=(v_dot+ell*s/2)/Theta, while the full lapse equation reconstructs
b=-E_n0/(2 a³ q Theta). The remaining Euler equations equal those of
the substituted L0, by the exact variational chain rule:

    E_v,reduced=E_v0-d_u(E_n0/Theta),
    E_s,reduced=E_s0+ell*E_n0/(2 Theta).

All time-dependent factors are retained. An independent nonconstant
coefficient example checks these identities, including lapse dynamics
and the induced second derivative of v. Therefore the reconstructed
lapse equation is not an unaccounted constraint on the Ostrogradsky
momentum. This supplements the actual full-M1 highest-Hessian and
Legendre identities; it is not a replacement toy for the actual action.

## Background and normalization

The independently reconstructed Vstar, Ustar, gradient shifts and
curl lapse terms agree before center specialization. All 64 background
sources vanish, so the new terms preserve the original solution.
The original J numerator has only even nonnegative powers with a
strictly positive constant, and its denominator is positive on all
real time; its proof is replayed, not replaced by samples.

The A=1,B=0 regular rank-one action agrees exactly with S6.38 before
the temporal-vector solve. Nonunit M²=3,tau=2,zeta_physical=5 gives
normalized zeta=5/12 and the isolated mass control 8/5 physically,
32/5 in normalized units. That control is not the coupled health
test, a heavy-gap estimate or a cutoff claim. The physical metric,
clock and free chi, and all previously fixed proof scopes, are unchanged.
