# Literal mass update, not an unchanged-Hessian shortcut

Use the S6.80 source-free action and its same one-form
W=L*C-B*dphi in the physical metric. Its retained mass
matrix is Dold^-1=diag(1/gamma_t,-1/gamma_s,-1/gamma_s,-1/gamma_s).
Replace it by eta=diag(1,-1,-1,-1). This is the actual
retained mass-square convention; the physical metric
remains the original (-+++) matter metric.

The added local density per physical volume is
W^T*(eta-Dold^-1)*W/2. Thus, writing shift=B*dphi,

Delta M=L^T*(eta-Dold^-1)*L,
Delta J=-L^T*(eta-Dold^-1)*shift,
Delta c=shift^T*(eta-Dold^-1)*shift/2.

The connection Hessian CHANGES. Both its source and
contact change are necessary. The source-free stationary
connection from S6.80 remains stationary, since its
retained trace is shift; the code checks all 64 Euler
equations directly with the new Hessian and source.
The generic shift symbol is the existing fixed B(phi,x),
not permission to choose another function.

Let R be the original full retained right inverse,
L*R=I and M_old*R=L^T*Dold^-1. The update gives
M_new*R=L^T*eta, hence M_new*(R*eta)=L^T and
L*(R*eta)=eta. These all-component identities reconstruct
the new inverse retained response and its Lorentz mass
matrix. The original 56-dimensional complementary
projector is annihilated by L; the Hessian update has
no complementary or complementary-retained block.
All projective columns and source annihilations are
independently checked, rather than inferred by deleting
four arbitrary connection coordinates.

On the sixty-dimensional quotient, the matrix determinant
lemma gives det(M_new)/det(M_old)=det(I+(eta-Dold^-1)*Dold)
=gamma_t*gamma_s^3. The frozen domain has both gammas
strictly positive, so this ratio never vanishes there.
The four retained and 56 complementary directions
remain a nondegenerate quotient. This is a local
classical affine elimination result, not unitarity of
an unspecified full affine quantum parent.

Since Dold=eta on the entire clock trajectory, the
added W^2 density starts at cubic physical field
degree about W=0. The clock and full classical
quadratic theory are unchanged. Its nonlinear
vertices and quantum insertions are not unchanged.
