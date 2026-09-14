# All twelve actual invariant images and full physical volume

For every real w in ||w||<=2*10^20, reconstruct the
complete fields as in notes/geometry.md. These are
deterministic classical phase-space points, not twelve
independent Gaussian oscillator coordinates.

Use EXACTLY the source definitions in notes/source.md.
Write f_X for the full physical A_s field bound on
this support ball. With V=exp(3v) and ||V^+-1||<2,
the following outward majorants retain all products:

    |p| <=2 f_Pi_v/3,
    |G| <=6(1+P) f_Pi_W,
    |pm-1/10| <=2 f_delta_Pi_M+10^-329,
    |ph| <=2 f_Pi_H,
    |eta| <=f_eta,
    |sh| <=48 f_pi_TF^2,
    |el| <=8 f_Pi_W^2/zeta,
    |ma| <=48 zeta [(1+P)f_W]^2,
    |wm| <=2 f_W^2,
    |gm| <=2 f_M1_A1^2,
    |gh| <=2 f_gradH^2,
    |curv| <=18[6(9e)+18(9e)^2].

The Gauss factor6 conservatively covers all derivative
components and V^-1. The matter density shift keeps
the background:
pm-1/10=(delta Pi_M)/V+(1/10)(V^-1-1).
The latter is at most(6/10)||v||<10^-329.
The magnetic bound includes all antisymmetric F=dW
components, both inverse metrics and their complete
contraction; no polarization is suppressed.

Every displayed quantity is strictly below10^-260.
This is checked by exact rational comparison for ALL
twelve source-pinned entries, uniformly in space.
Since10^-260<10^-250/2, the entire image is strictly
inside the original S266 invariant box. The actual
heavy ratio is the unchanged one in that theorem.
The full positive lapse and temporal roots, pivots
and coefficients therefore apply. This does not
claim off-bounce time control or solve the remaining
global translation charges classically.

Let N0=Nstar(0) at the original reference invariants.
The unchanged full source gives |N0-1|<10^-375.
The full derivative bound |N_i|<2 throughout the
convex twelve-invariant box permits integration along
the entire straight invariant segment. Consequently

    |Nstar(z)-N0|<=2 sum_i |z_i|<=24*10^-260.

This does not treat the lapse as its linear Gaussian
reference. The actual invariant curve in canonical
phase can be nonlinear; the uniform box and full
implicit derivative bound are enough.

The physical volume is the original complete
exp(3v)U(Nstar), U=R_full(0,Nstar)^(-3/4), spatially
averaged with normalized torus measure. The S265
whole-strip derivative ceiling is |U_N|<=12.
Using only the center value3/2 would be invalid.
Use also U<2, exp(3v)<2 and
|exp(3v)-1|<=6|v|. Then

    |F(z)-F0| <=12*10^-329+24*(24*10^-260)
                  <10^-256,
    |F0-1| <=12*10^-375.

The factor24 in the first line permits the extra
exp(3v) factor in the lapse variation. These are
bounds for the complete source function on all
reconstructed points, not Taylor-polynomial
substitutions. S261's incorrect physical slope
is neither used nor repaired inside its archive.

Any smooth 0<=chi<=1 with core radius10^20 and
support strictly inside the allowed radius2*10^20
gives the SAME named extension

    Fext=F0+chi(F-F0).

It remains positive everywhere, and
||Fext-1||infinity<10^-255. The background F0 is
not reset to1. This uniform pointwise estimate,
together with the positive unital S268 coherent
map, gives ||Q_V(Fext)-I||<10^-255 on the full
infinite Hilbert space and on its reducing
zero-translation-charge subspace. It is a finite
BOUNCE-SLICE coherent observable estimate, not
an evaluated original interacting volume mean.
