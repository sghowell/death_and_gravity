# Literal full-target clock tree, boundary and matter charge

Use the full analytic S6.109 scalar/tensor functions inside the
current S6.174 parent. Their literal clock jets are checked through
R_XXX and F_XX. These suffice for this quadratic calculation;
the global functions are not replaced by their clock polynomials.

In the regular hat metric chart write the homogeneous action per
kappa and per reference coordinate volume as

    exp(3vhat) [9a(N)(H+vhat')^2/N
                 +3b(N)(H+vhat')+N f(N)
                 +U(N)(ell+sigma')^2/(2N)].

Here U=R^(-3/4), a=-R^(1/4)/3. The exact transformed coefficients
b and f include the S6.174 primitive I(u,1)=0 with its full
weighted integration-by-parts boundary. Its relevant lapse jets are

    I_N=0, I_NN=-3h'/h^3, h=(1+u^2)^3.

The actual b_NN contains -I_NN, and f_NN contains
-derivative_u I_NN. The latter is nonzero at the bounce even
though I_NN itself vanishes there. No placeholder primitive
value is held fixed while differentiating.

Expanding the literal density to degree two, retaining the
weighted vhat*vhat' boundary, gives

    L2=-3vhat'^2+(J+w0^2/2-3Theta^2)n^2+6Theta n vhat'
         +sigma'^2/2+w0 n sigma'-3ell vhat' sigma,

    ell=1/[10(1+u^2)^6], delta=1/(2h), w0=ell(3delta-1).

The full calculated J agrees with the independently derived
S6.174 Hamiltonian lapse pivot. Both scalar/matter kinetic
normalizations, the vanishing n*vhat and vhat^2 terms, and
ell'+3H ell=0 are checked. Omitting the primitive shifts J(0)
by -9: the actual J(0)=243/160 would acquire the wrong sign.

The M1 equation fixes a^3(sigma'+w0 n+3ell vhat). Prepared
zero charge perturbation therefore selects
sigma'=-w0 n-3ell vhat. Varying the original equations and
performing the fixed-charge Routh reduction agree exactly.
Direct on-shell velocity substitution is not the variational rule.

The reduced density is

    L_eff=-3vhat'^2+6Theta n vhat'+(J-3Theta^2)n^2
           -3ell w0 n vhat-(9/2)ell^2 vhat^2.

Now set n=eta' and vhat=w+H eta-delta eta'. Weighted Euler
variation of this actual density agrees with an independent
physical-current transformation. The complete derivative-order
matrix is [[4,3],[3,2]], with top coefficients
(-6delta^2,6delta,-6delta,6). In particular

    A=-6delta^2, |A|>=6144/15625, |A^-1|<=15625/6144.

These are continuous compact-slab bounds. A is a coefficient in
the prepared constraint reduction, not an asserted particle
kinetic sign. No division by H, Theta or spatial momentum occurs.
