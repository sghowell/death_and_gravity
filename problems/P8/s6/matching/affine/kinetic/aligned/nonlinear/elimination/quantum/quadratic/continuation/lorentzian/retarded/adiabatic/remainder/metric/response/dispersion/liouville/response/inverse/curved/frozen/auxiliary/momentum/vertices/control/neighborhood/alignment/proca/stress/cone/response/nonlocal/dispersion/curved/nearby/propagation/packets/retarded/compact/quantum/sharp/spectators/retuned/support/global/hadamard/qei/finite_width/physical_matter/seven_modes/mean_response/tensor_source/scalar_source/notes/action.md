# A scalar quadratic action before homogeneous mean restriction

Use the unchanged retuned Hamiltonian and spatial-only frame.
Keep the homogeneous lapse N, scale A, trace density p and
matter charge density l independent before differentiating.
Write the scalar spatial metric A^2 exp(2zeta) I and natural
canonical density momenta Pzeta,Pchi. Set r=Pzeta/A^3,
c=Pchi/A^3, q=k^2/A^2. The exact trace and matter densities are

    p_full=exp(-3zeta)(p+r/3),
    l_full=exp(-3zeta)(l+c).

The linear momentum constraint fixes the tracefree mixed
density (divided by A^3) to

    sigma=(k_i k_j/k^2-I/3)(3l chi-r+9p zeta)/4.

Its divergence and trace are checked in all three Cartesian
directions. Its square is (r-9p zeta-3l chi)^2/24.
The background shear vanishes, so higher-order tracefree
York terms cannot enter this quadratic Hamiltonian.

Denote the literal scalar Hamiltonian coefficients by
a_j,b_j,c_j,d_j,e_j,f_j,g_j, where j is a lapse derivative:

    Hcal=a p^2+b p+c+d l^2+e S+f R3+g |grad chi|^2.

They are reconstructed from the ACTUAL S6.102 invariant
rows and their original retuned lapse jets, not fitted at
the bounce. The tensor/Proca invariants are zero here.

The integrated quadratic Hamiltonian per A^3 before lapse
elimination is Q+nu L+nu^2 D/2, with

    Q=a_0(9p^2 zeta^2/2-2pr zeta+r^2/9)
       +9c_0 zeta^2/2
       +d_0(9l^2 zeta^2/2-6l zeta c+c^2)
       +e_0(r-9p zeta-3l chi)^2/24
       +g_0 q chi^2+2f_0 q zeta^2,
    L=(-3a_1 p^2+3c_1-3d_1 l^2+4f_1 q)zeta
       +(2a_1 p+b_1)r/3+2d_1 l c,
    D=a_2 p^2+b_2 p+c_2+d_2 l^2.

The curvature term follows directly from
A^3 exp(3zeta) R3=A exp(zeta)(-4 Delta zeta-2|grad zeta|^2).
After integration by parts its constant-lapse quadratic
coefficient is 2A k^2 zeta^2; the mixed lapse coefficient
retains the unintegrated linear curvature 4q zeta.

The independent source differentiates the canonical-volume
expression twice and checks every term above. Eliminate
only the NONZERO-MOMENTUM lapse: nu=-L/D. This gives

    H2=A^3[Q-L^2/(2D)].

The homogeneous lapse has NOT been solved or frozen in this
step. On the actual clock D=-2J with the positive retuned J.
At every finite time continuity supplies a regular nearby
homogeneous parameter domain. Differentiating H2 with respect
to N retains the third lapse jet in D_N. The direct shifted-
coefficient calculation independently checks this derivative.
