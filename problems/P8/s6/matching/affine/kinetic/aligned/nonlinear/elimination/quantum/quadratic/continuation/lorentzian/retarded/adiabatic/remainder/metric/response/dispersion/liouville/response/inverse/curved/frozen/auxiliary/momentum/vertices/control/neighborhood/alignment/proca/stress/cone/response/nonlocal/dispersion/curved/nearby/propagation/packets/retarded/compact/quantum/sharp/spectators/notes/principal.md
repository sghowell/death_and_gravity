# Fixed-frame light block and the all-dimension principal theorem

## The actual comparison

Write kc,km>0 for the inherited clock and matter kinetic coefficients,
ell for the homogeneous matter velocity, and c for physical clock speed.
The inherited action gives

    K = [[kc+km*ell^2, km*ell], [km*ell, km]],
    H = [[kc*c^2+km*ell^2, km*ell], [km*ell, km]].

Here H=e^2 G/N^2 uses ONE physical matter frame. Equivalently divide the
coordinate-frequency gradient by omega_m^2; do not normalize each mode
by its own speed. With

    T = [[1/sqrt(kc), 0], [-ell/sqrt(kc), 1/sqrt(km)]],

T^T K T=I and T^T H T=diag(c^2,1). The exact whole-interval parent
enclosure implies c^2-1>delta=5*10^-6, with both kinetic coefficients
strictly positive. This is an actual on-shell nearby background, not an
independent off-shell choice of lapse and momenta.

## Finite Hermitian principal extensions, including mixed terms

Fix one spatial direction and one background point. Assume the reduced
physical principal pencil has the finite-dimensional Hermitian form

    P(s)=s^2 K_full+2s C_full-H_full,

where K_full is positive definite. All coefficients are finite at this
point. C_full can be any Hermitian mixed time-space matrix. Allow arbitrary
heavy blocks and arbitrary cross-block entries in all three matrices.
The canonical light restriction is I+DeltaK, DeltaC and
diag(c^2,1)+DeltaH. The real symmetric four-mode calculation in pencil.py
is an exact algebraic anchor; the following proof has no dimension bound.

Embed the canonical light clock unit vector as x=(1,0,...,0). Its Rayleigh
value is

    x*P(s)x=s^2(1+DeltaK_11)+2s DeltaC_11-c^2-DeltaH_11.

All heavy and cross blocks disappear from this restriction. If DeltaK and
DeltaH vanish, x*P(c)x+x*P(-c)x=0. Thus for at least one sign sigma,
lambda_min(P(sigma*c))<=0. On that ray, the bound

    lambda_min(P(sigma*r)) >= r^2 lambda_min(K_full)
                              -2r ||C_full||-||H_full||

is positive for sufficiently large finite r. Eigenvalues of finite
Hermitian matrices are continuous in their entries. The intermediate
value theorem therefore gives a zero eigenvalue at some r>=c>1.
If the initial minimum is already zero, that point itself is a root.
There is a real characteristic with |s|>=c. No assumption that all other
roots are real is required; a separately ill-posed branch would not make
a healthy causal completion either.

When C_full=0 this reduces to the familiar Rayleigh statement that the
largest eigenvalue of K_full^-1 H_full is at least c^2. That special case
is not used to silently exclude mixed time-space couplings.

This theorem assumes an actual reduced physical polynomial pencil and an
explicit light embedding. A rational frequency-dependent constraint or
heavy-field Schur complement is not automatically such a pencil. Nor is
an EFT Wilson coefficient automatically a bare parent light block. The
Proca application below is proved separately using its actual constraints.
