# Positive subset gaps and the complex weighted current

## Positive real center and relative complex tube

Fix real positive wi, fixed real unit directions ni, fixed complex
spatial TT leaf tensors, and a generic real tree configuration.
For each leaf put zi=wi*(1+ei), |ei|<=epsilon. The analytic action
continues bilinear dot products and the same square-root recoil
branches. Absolute values occur only in estimates, not in the
definition of a continued amplitude.

For every subset S set W=sum wi, Qreal=(W,sum wi ni),
v=Qreal_sp/W, delta^2=1-|v|^2. The original weighted current norm
uses this REAL center v and delta throughout the complex tube.

For any subset,
 |sum zi-W|<=epsilon W,
 |Qcomplex^2-Qreal^2|<=(2epsilon+epsilon^2)Qreal^2.
The second inequality follows termwise from
Q^2=2 sum_(i<j)wi wj(1-ni.nj), whose real coefficients are nonnegative.
Hence every generic pure-soft inverse remains nonzero, uniformly
as real angular invariants tend to zero. Its relative gap has no
multiplicity factor. The temporal energy inverse remains nonzero too.

## Conserved complex temporal inverse

Let Z=Qcomplex^0, u=Qcomplex_sp/Z. If R is the spatial block of a
complete conserved source, its remaining components are
J0i=(R*u)i and J00=u^T R u. Put T=I-u*u^T; transpose is bilinear.
The exact spatial temporal inverse is

 Hsp=[-T*R*T+(1/2)*T*tr(T*R)]/Qcomplex^2.

This follows directly by applying the unchanged trace-reversed
inverse and the temporal projection. A symbolic nine-entry check
verifies it for arbitrary u and arbitrary symmetric R. It does not
assume real u, a conserved isolated diagram, or a projected free
off-shell external block.

Choose the real parent axis n along v; for v=0 choose any real unit
axis. For epsilon<=1/100, real future geometry gives
 ||u_T||<=sqrt(2)*epsilon*delta/(1-epsilon)<delta/50,
 |1-u_z|<2delta^2, |u_z|<2.
For the first inequality use cancellation of the real transverse
sum, followed by Cauchy-Schwarz:
sum wi|ni_T|<=sqrt(2W sum wi(1-n.ni)).
The longitudinal bound uses
sum wi(1-n.ni)=(1-|v|)W<=delta^2 W.

Set D=diag(1,1,delta). Then T=D*B*D, where the blocks of B satisfy
 ||B_TT||F<2, ||B_TL||<1, |B_LL|<6.
Thus ||B||F^2<4+2+36<64.
If each source TT entry is bounded by delta^2*L, each TL entry
by delta*L, and LL by L, write
 R=delta^2*D^(-1)*Rbar*D^(-1), ||Rbar||F<=3L.
The numerator factors EXACTLY as
 delta^2*D*[-B*Rbar*B+(1/2)*B*tr(B*Rbar)]*D.
Its middle norm is at most(3/2)*8^2*3L=288L.
Since |Qcomplex^2|>delta^2 W^2/2, every component of the real-center
weighted norm is at most576L/W^2<1024L/W^2.
The D factors give the velocity and double-velocity weights; |v|<=1.
No complex rotation or Hermitian replacement of the action is used.

## Complex extension of the all-order grading

Use the S318 decomposition relative to the real parent center.
The child tensor transfer remains valid for complex fields: the
geometric vectors are real and all field estimates use Hermitian
Frobenius norms. The child coefficient-norm sum remains
13*(W/W_child)*N_child.

Continuing energies relatively multiplies each momentum coefficient
budget by at most1+epsilon. The real budget below4W is therefore
below5W. Momentum conservation still holds coefficientwise in the
auxiliary angular grading variable.

The constant TT vertex vanishes for the same arbitrary Rosen metric.
That identity is algebraic and remains valid for complex transverse
metric coefficients. The real proper transverse rotation still
implies TT even, TL odd and LL even. Thus source TT begins at
delta^2, TL at delta, LL at delta^0, at every valence.

The S315 coefficient-l1 vertex norm estimate uses triangle and
submultiplicative inequalities valid for complex coefficients.
Each root partition with k children now has source budget
 25W^2*(k+1)!*32^(k+1)*13^k
 *product[(W/W_child)*N_child]*kappa^(-(k-1)/2).

Combining it with1024/W^2 gives the recurrence
 D1=1,
 Dn=25600 sum_(labeled partitions,k>=2)
       (k+1)!32^(k+1)13^k product D_child.

Energy and coupling exponents are exactly those in S318, so
 N_center(H_S(z),Qreal_S)
 <=Dn W_S^n/[product_i(wi)*kappa^((n-1)/2)].

The nonnegative EGF is
 d=z+819200*((1-416d)^(-2)-1-832d).
An exact barrier at cap1e-13 and radius5e-14 gives
 Dn<=2*(2e13)^(n-1)*n!.
This is an auxiliary majorant series, not quantum convergence.
