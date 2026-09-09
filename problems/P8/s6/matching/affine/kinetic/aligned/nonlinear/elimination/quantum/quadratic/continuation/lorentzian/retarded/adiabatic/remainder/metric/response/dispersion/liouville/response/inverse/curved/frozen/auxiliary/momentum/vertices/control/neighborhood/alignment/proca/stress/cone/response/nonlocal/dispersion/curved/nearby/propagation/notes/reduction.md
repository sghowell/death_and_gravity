# Full moving-background reduction and Euler-first principal

Let Hcal=a p0^2+c p0+f+m ell^2 be the homogeneous Hamiltonian per
hat-frame volume, with p0 the trace density and ell the matter
density. The symbols a,m here are coefficient functions, not the
physical scale or Proca mass. Then

    Hhat=(2a p0+c)/3,
    p0'=-a p0^2-c p0-f+m ell^2.

The exact geometric/York calculation in geometry.py expands the
full volume, metric, inverse, scalar curvature, matter gradient,
trace and shear invariants for every unordered scalar phase pair.
Its local hat scale is normalized to one, with the scale dependence
restored by q=k_comoving^2/R^2 and the measure R^3 du.

The moving Hamiltonian includes -2 Hhat pi:g, the background
matter translation, and

    (p0'+3 Hhat p0)(6v+3v^2)/2.

The 3 Hhat p0 term cannot be omitted off center. It comes from
differentiating the volume-weighted background trace in the
canonical one-form. The conserved matter charge, not a time-varying
density held incorrectly constant, supplies the canonical matter
boundary. On Hcal_N=0 the full ten phase-pair coefficients give

    H2=m(P-3 ell v)^2+(2/3)a ell p chi
       +(g q-a ell^2)chi^2+2r q v^2-F1^2/(2h),
    F1=alpha p+beta(P-3 ell v)+4r_N qv,
    alpha=Hcal_Np/3, beta=2m_N ell, h=Hcal_NN.

All lapse derivatives here hold canonical phase fixed. The shear
coefficient is -8a/3 for this literal scalar-tensor action.
The last term is actual linear lapse elimination, not a new
constraint or an ansatz for the spatial quadratic.

## Crossing-regular momentum chart

Set v=P_b/(2q), p=-2qb. Since q'=-2 Hhat q, integrating the
R^3-weighted canonical one-form by parts adds -Hhat b P_b
to the transformed Hamiltonian. Its momentum Hessian has determinant

    -(beta^2 r-2h m r+8m q r_N^2)/(h q).

This can vanish at some finite positive q; it is not declared an
all-frequency chart. For sufficiently large q with m,h,r_N nonzero,
the ordinary Legendre transform gives the exact leading matrices

    K=[(beta^2-2mh)/(8m r_N^2), -beta/(4m r_N);
       -beta/(4m r_N),          1/(2m)],
    mixed/q=diag(B,0), B=alpha/r_N,
    frozen_potential/q=[-2 Hhat B+r B^2, -4a ell/3;
                        -4a ell/3,       2g].

Compute the Euler equation before the principal limit. The symmetric
leading mixed term contributes B'+Hhat B after the time derivative
of both q and R^3. Consequently

    G=[B'-Hhat B+r B^2, -4a ell/3;
       -4a ell/3,       2g].

Using only the frozen potential gives a different, incorrect
gradient matrix. The remaining antisymmetric mixed terms are lower
order in q. There is no division by Hhat, alpha or a crossing
coefficient anywhere in this construction.

## Literal physical-metric specialization

For this action, writing e=eomega,

    a=-3Ne/4, m=N/(2e^3), g=Ne/2, r=-m, r_N=-m_N.

The change chi=sigma-ell b has unit determinant and diagonalizes
both principal quadratic forms:

    K_clock=-h/(4m_N^2), K_matter=e^3/N,
    G_clock=B'-Hhat B-m B^2-N e ell^2, G_matter=N e.

Time variation of ell changes lower-order terms, not these principal
forms. Since dt=Ndu and physical q=qhat/e^2, physical speed squared
is e^2/N^2 times a generalized eigenvalue of (G,K). The matter
value is exactly one, and

    c_clock^2=-4e^2 m_N^2 G_clock/(N^2 h).

These expressions agree with the independent full S6.83 central
calculation, including its Euler time derivatives. Positive K and
positive characteristic values imply positive G by this same
real invertible basis.
