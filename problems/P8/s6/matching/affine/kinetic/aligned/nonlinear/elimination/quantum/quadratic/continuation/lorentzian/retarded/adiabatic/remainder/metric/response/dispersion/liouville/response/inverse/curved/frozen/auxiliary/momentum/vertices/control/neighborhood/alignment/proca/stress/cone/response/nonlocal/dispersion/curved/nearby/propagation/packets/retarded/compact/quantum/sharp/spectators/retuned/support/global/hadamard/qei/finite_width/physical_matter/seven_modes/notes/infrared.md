# Global distributional existence and the zero-frequency end

Fix an exact finite R>=3 and consider any times in [-R,R].
The actual and preparation scales lie between 1 and
A_R=(1+R^2)^2. The previous energy argument gives

    Proca floor=min(1,10^6/A_R),
    Proca ceiling=max(10^6,A_R)(1+q),
    tensor floor=2/A_R^3,
    tensor ceiling=A_R(1+q)/2.

All floors are strictly positive. The common relative energy
growth rate 697 bounds both sectors on this whole strip.
Gronwall and the energy eigenvalue comparison consequently
bound every transfer by C_R sqrt(1+q). Differentiating
U'=M U gives polynomial frequency bounds for every fixed
number of time derivatives, because M and all its time
derivatives have finite polynomial frequency degree on the strip.

At m=1000, the nonzero Proca covariance eigenvalues are
those of (Q+P)/2. In the longitudinal direction Q=omega,
P=1/omega; in the two transverse directions Q=m^2/omega,
P=omega/m^2. Using m<=omega<=m sqrt(1+q) gives

    ||C_P,flat|| <=501 sqrt(1+q).

It is smooth at k=0. The tensor covariance norm is exactly
1/k+k/4, at most (1+k^2)/k. After reference preparation,
its anchor norm is bounded by a constant times
(1+k^2)^2/k. Near zero the radial three-dimensional integral is

    integral_0^1 k(1+k^2)^2 dk=7/6.

Actual compact-time propagation is bounded near k=0 and adds
no new negative power. The Cartesian TT projector has norm one
for k!=0, so it likewise preserves local integrability. Its
value at the single origin is irrelevant to the Lebesgue
distribution. No delta measure or classical zero-mode atom
is being attached there.

The scalar factor keeps the already established S6.99/S6.101
infrared completion and degree-nine global transfer bound.
Its high-frequency two-point majorant has degree 21; the
spectator bounds above require a smaller polynomial degree,
so a common degree-21 high-frequency majorant is available
for the combined phase on each compact time strip. Any
specified time or spatial derivatives increase this by a
finite degree. Spatial Schwartz smearing then defines a
continuous two-point distribution.

An infrared-cutoff tensor kernel is smooth: all finite time
derivatives obey integrable compact-frequency bounds, and
each spatial derivative only inserts another bounded power
of k. Dominated differentiation therefore proves smoothness
in both spacetime legs. The angular projector need not be
differentiable at zero for this argument. It is the kernel's
spacetime variables, not the momentum origin, that are being
differentiated.
