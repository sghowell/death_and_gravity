# Literal new action and all-connection lift

This is S6.80, a new local source/contact candidate under the S6 matching
contract. Nothing in the frozen previous model is overwritten. Work in
the same timelike-clock tube, physical matter metric, hat chart and curled
one-form W=T_trace-B*dphi. The exact old normal source is

S = delta*K_hat+c = delta*(K_hat-3*H/N)+(3/2)*Q/N,
delta=(N^-2-1)/h.

The original affine source is covariantly parallel to dphi, not merely
normal on homogeneous data. `degeneracy.source()` reconstructs it from
the full Hessian of phi and checks all four components. The aligned shift
has the same property. S has no lapse or temporal-vector velocity in the
full spatial hat chart. Replace only

N*sqrt(h_hat)*U*(T-S)^2/(2*gamma_t)

by N*sqrt(h_hat)*U*T^2/(2*gamma_t). Thus the literal added density is
N*sqrt(h_hat)*U*(T*S-S^2/2)/gamma_t. All scalar/matter terms, original
boundary primitive, Q, margin, Maxwell curl and mass functions are kept.
This is not a field redefinition of the old theory.

For the all-connection lift let L be the original 4-by-64 retained trace
map, E the 64-by-60 quotient embedding, M the unchanged full Hessian and
W_lift the old quotient inverse lift. The frozen identities give

M*E*W_lift=L^T, D=L*E*W_lift,
D=diag(gamma_t,-gamma_s,-gamma_s,-gamma_s).

This (+---) convention is the actual retained mass-square convention;
it is not an instruction to change the physical (-+++) metric. Write
Svec=(S,0,0,0), Cstar for the original stationary connection, and
shift=L*Cstar-Svec. The added retained density is

(L*C-shift)^T D^-1 Svec - Svec^T D^-1 Svec/2.

At fixed light fields it is linear in C, with added source
L^T D^-1 Svec and zero Hessian. Its new stationary connection is
Cstar-E*W_lift*D^-1*Svec. Direct multiplication checks all 64 new Euler
equations and L*Cnew-shift=0. The added source annihilates all projective
gauge columns and the original 56-dimensional complement. The old
invertibility and complementary Gaussian elimination therefore persist.
No unspecified affine operator is being assumed to supply the change.

Along the old clock S0=S1=0 in all four components. With N=1+e*n and
K_hat=3H+e*k, the actual fixed-basepoint Q equation gives
Q_NN=-3*h'/h^3 and

S2=-2*n*(k+3*H*n)/h-(9/4)*h'*n^2/h^3.

For C=C0+e*C1, T=e*T1+e^2*T2 and S=e^2*S2+e^3*S3, the changed
density starts at degree three, C0*T1*S2; degree four is
C0*(T1*S3+T2*S2-S2^2/2)+C1*T1*S2. Both terms are retained.
The old and new background and full quadratic action agree exactly.

The new action is even in the complete W, not just its temporal
component: curl and mass terms are quadratic and the source is gone.
Mass and metric coefficients still depend on the light fields. Even
light-vector interactions remain; this is not a decoupled free theory.
