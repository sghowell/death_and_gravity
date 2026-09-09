# Independent center calculation without labelled geometry

Set ell=1/10, h=a=1, H=Theta=0, Lambda=-1/2, w=1/20,
Jnew=243/160. For wave direction x the independently solved
linear momentum is

    pi_xx=ell chi/2,
    pi_yy=(p-ell chi)/4+PT/2,
    pi_zz=(p-ell chi)/4-PT/2,
    pi_yz=pi_zy=PS/2.

Its trace is p/2. Thus dp1=p/3 and the shear invariant is
tr(pi^2)-(tr pi)^2/3. Expanding the inverse determinant of the
actual LINEAR metric coordinate gives

    dc1=Pchi-3ell v,
    dc2=(15/2)ell v^2+(ell/2)(T^2+S^2)-3v Pchi,
    R1=4k^2 v,
    R2=-10k^2 v^2-(k^2/2)(T^2+S^2) at zero output.

For the pure scalar metric the exact curvature is
-4 Delta v/(1+2v)^2+6|grad v|^2/(1+2v)^3.
Its quadratic coordinate average is -10|grad v|^2 by integration
by parts. The TT quadratic average is -|grad t|^2/4; mixed
scalar-TT terms vanish by transversality and integration by
parts. The TT first volume variation vanishes, so replacing
sqrt(g)R2 by R2 in that TT calculation changes nothing at this
order. It WOULD change the pure scalar coefficient.

The other quadratic invariants are

    G2=k^2 chi^2, j^2=k^2 Px^2,
    electric=10^6(Px^2+Py^2+Pz^2),
    magnetic=(2/10^6)k^2(Wy^2+Wz^2),
    vector=Wx^2+Wy^2+Wz^2.

The independent force is

    F1=-dc1/20+R1/4,
    F2=-dc2/20+R2/4-9dp1^2/8-dc1^2/4-j^2/4
        +3 shear+electric/4+magnetic/8+3 vector/4+3G2/4,
    L2=3dc1/40-3R1/8,
    n1=F1/(2Jnew),
    n2=[F2+(13821/800)n1^2+L2 n1]/(2Jnew).

Differentiating the independent rho2, p2 and omitted-n2
polynomials in all fourteen fields reproduces all 588 entries
of the three separately reconstructed matrices. In particular

    M_rho[PT,PT]=M_rho[PS,PS]=-2/135,
    M_rho[Py,Py]=M_rho[Pz,Pz]=-200000/81.

The large numerical Proca entry uses the original unrescaled
vector momentum and its 10^-6 kinetic coefficient. It must not
be confused with a coefficient in unit-normalized Proca fields.
The pressure matrix differs from the density matrix only in the
matter-gradient entry, as the physical formula requires.
