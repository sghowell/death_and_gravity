# Actual lapse variation before the dimensional limit

Use the exact original alpha=4/(9h), beta=28/(81h),
h=(1+u^2)^3, and Yt=m^2 alpha n,Ys=m^2 beta n. Their full time
derivatives and k'=-Hk are retained. Let delta=D-3.

Write the compact action as integral a^D[A_D n''^2+B_D n'^2+C_D n^2].
The zeroth dimensional jets are exactly S6.61. For r=1+u^2,q=k^2,
the first derivatives with respect to D at D=3 are

    A2,1=0,
    B2,1=-574m^2/(19683r^6),
    C2,1=-2m^2[259q r^2+32994u^2+2106]/(19683r^8),

    A4,1=-2896/(98415r^6),
    B4,1=-4[1389q r^2+94370u^2+8070]/(98415r^8),
    C4,1=-4[1995q^2 r^4+q r^2(149764u^2+49984)
            +6795480u^4+3319008u^2-107352]/(295245r^10).

All six coefficients are checked against the direct dimensional
expansion, not inserted as the defining calculation.

Let W_D=D_u+D H and W=W_3. Compact reduction must use W_D.
For a raw density with coefficients c_ij of n^(i)n^(j),

    A_D=c22,
    B_D=c11-c02-W_D(c12)/2,
    C_D=c00-W_D(c01)/2+W_D^2(c02)/2.

Thus differentiating the compact coefficients includes
partial_D W_D=H. The normalized Euler operator is

    E_D=2C_D n-2W_D(B_D n')+2W_D^2(A_D n'').

Its first dimensional jet is not obtained just by replacing
A,B,C with their first jets in the four-dimensional operator:
one must additionally include

    -2H B0 n'+2W(H A0 n'')+2H W(A0 n'').

Direct variation of the unreduced density agrees exactly with
the compact result at both dimensional orders. The omitted
measure term has the nonzero order-four bounce fixture
-1024/6561 for n''=1, all other n jets zero, k=0.

## Counterterm sign and interpretation

With epsilon=(3-D)/2, the pole counterterm is -Q_D/(32 pi^2 epsilon).
Its normalized lapse variation expands as

    -E0/(32 pi^2 epsilon)+2E1/(32 pi^2)+O(epsilon).

The displayed finite quantity is its evanescent component only.
It must be combined with the dimensionally continued bare
contact/retarded kernel and common normalization before limiting.
It is not a standalone physical finite effective action or
self-adjoint four-dimensional response: differentiating the
dimension-dependent adjoint weight is part of this component.
No loop norm, finite stability, cutoff or Wilson-coefficient
boundary value follows from these coefficient formulas alone.
