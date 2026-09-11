# Actual-state leading source-only tensor noise

Use exactly the S6.186 full local centered Proca stress
in the actual unchanged all-order CD Gaussian state,
with mass1000, kappa=1e800 and its pinned covariant
renormalization prescription. Its complete constrained
three-mode, internal/external momentum proof gives

    Var T[f] <1e50 N(f)^2,
    N(f)^2=sum_(j=0..3)||partial_t^j f||L2^2+||grad f||L2^2.

A spatial TT test is embedded as the spatial block of
the physical four-tensor; the Frobenius norm introduces
no additional polarization or off-diagonal factor.
The covariance is that of the local quadratic stress,
not a Gaussian assumption on its higher moments or a
product of separately smoothed fields. Fixed c-number
stress counterterms cancel only in centering. The mean
prescription and scalar retuning are not changed.

For the compatible class, v=Gadv q is a valid compact
TT stress test and the source-only observable identity
h[q]=T[v]/sqrt(kappa) holds. The preceding energy bound gives

    Var h[q] <5000^2*1e50/kappa *D(q)^2,
    STD h[q] <5e-372 D(q).

For gamma=2h/sqrt(kappa),

    STD gamma[q] <1e-771 D(q).

The machine packet also uses the smaller exact rational
S6.186 variance majorant before comparing to these displays.
All spatial momenta and the actual wave operator are
retained. No finite radial grid or physical UV cutoff
stands in for the continuous proof.

This is the source-only tree-propagated component.
The homogeneous part cancels on this detector class,
so no initial-state factorization is needed for this
identity. The calculation does not establish the
marginal or joint state of a full interacting tensor
field. It does not bound the full self-energy inverse,
intrinsic/cross correlations at subsequent order, or
a finite-coupling error. Smallness of this leading
observable is not a quantum-background verdict.
