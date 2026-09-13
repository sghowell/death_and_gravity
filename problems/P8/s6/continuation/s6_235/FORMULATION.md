# S6.235: complete first-loop representation and separate symmetric value matching

## Frozen input and new prescription

Keep the S233 V2S-T1 classical model and the S234 base first-loop scheme. In particular the light mass and residue are fixed on shell THROUGH THIS ORDER, the heavy one-point condition remains, and the base vertex/heavy-mass UV poles use MSbar at mu1. The original affine/DHOST/Proca parent is a different model.

Define a new matching prescription, V2S-T1-OS4, by adding one finite four-light contact at first loop order:

deltaC_fin=-A1_base(s0,s0,s0), s0=4/3.

The point is on shell with s+t+u=4, but is subthreshold and is not a real2-to2 scattering configuration. The new condition is explicitly recorded; no frozen finite prescription is edited. It fixes a value, not derivatives or the full angular amplitude.

## Complete first-loop formula

Let n=M_H²=D+2 and A_s=C+g²/(n-s). With the common loop factor outside the scalar parameter integrals,

A1_base=(16pi²)^-1 sum_channels[
 A_s² B_MS(s)/2+2A_s g² Cbar(s)
 +g4(Dbar(s,t)+Dbar(s,u))].

All six ORDERED boxes remain. Their first argument crosses two light lines and their second two heavy lines, so Dbar(s,t) is not interchangeable with Dbar(t,s).

At equal external virtuality e, define

Q_C=(1-z)+nz-ez(1-z)-s(1-z)²xi(1-xi),
Q_D=Q_C-tz²eta(1-eta).

Cbar integrates(1-z)/Q_C over z,xi in[0,1]; Dbar integrates z(1-z)/Q_D² over z,xi,eta in[0,1]. Set e=1 for on-shell light fields. The Feynman boundary prescription is retained. The off-shell all-zero external jet instead has e=s=t=0. B_MS(s)=-integral Log[1-sx(1-x)-i0]dx at mu1.

## UV and first elastic checks

The complete pole is Delta sum A_s²/(32pi²). It cancels against the combined tree insertions of

deltaC_UV=-3C²Delta/(32pi²),
deltag_UV=-CgDelta/(32pi²),
deltaM_H²_UV=g²Delta/(32pi²).

The finite quadratic/one-point conditions remain those of S234. The complete forward first cut for4<s<M_H² agrees with beta integral A_tree²/(64pi), including the full rational angular dependence. These are first-order statements, not an exact S matrix.

## Subthreshold sign and value correction

At s=t=u=4/3 the parameter denominators are uniformly positive and both exceed2/3+D z. Complete integral bounds give

0<B_symmetric<1/2,
0<C_symmetric<462/D,
0<D_symmetric<462/D².

The actual parameters additionally give C_symmetric>375/D and abs(A_symmetric)>(19/10)g²/D, with A_symmetric negative. Consequently

A1_base(s0,s0,s0)<-5985g4/(64pi²D²)<0.

The complete tree there is4g²/[D²(3D+2)]>0, and their magnitude ratio exceeds10^190. The large ratio concerns this extremely small cancelled value; it is not a proof of strong coupling or a UV no-go.

An independent upper bound gives

abs(A1_base)<14793g4/(192D²),
0<deltaC_fin/(24q)<14793g²(D+2)/[768(D-1)]<10^-6.

The contact-only classical comparison C to C+deltaC_fin thus leaves remaining heavy-square quartic above(1-10^-6)q. It remains nonnegative and coercive. The quantum effective potential and the bare minimum after ALL finite counterterms are separate matters.

## Independent off-shell control and remaining work

The full constant-background two-field Hessian yields an independently evaluated off-shell zero-momentum jet with its complete bubble/triangle/box multiplicities. That jet is not used in place of the on-shell matching point.

Exact parameter inequalities and analytic primitives supplement independent diagram, cut and high-precision checks. Written arguments are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter.

This checkpoint does not prove the remaining physical-angle real remainder, higher quantum b20/b21/b40 matching, all-loop errors, continuum UV completion, original covariant bounce/state domain or finite-gravity Regge bound. Original V/G/B/P8 remain open.
