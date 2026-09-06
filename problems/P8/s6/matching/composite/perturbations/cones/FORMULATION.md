# P8-S6.8.COMPOSITE: a scoped all-mode principal-cone no-bounce theorem

Pin S6.7.COMPOSITE certificate
`018295c17af163cc80801e3dc2a9d7309184e4d1f5eca319b2a59ea85f064090`,
including its full S6.6 source-aware background lineage. The action and
physical metric are unchanged: constant positive Einstein coefficients G,F,
constant real HR beta_n, constant positive m4, and a single canonical scalar coupled to
`g_eff=g(alpha I+beta sqrt(g^-1 f))^2`, alpha,beta>0.
The positive affine couplings alpha,beta are constants, not functions of
the scalar or time.

Consider a smooth solution on a connected time interval with regular
positive finite common-flat FLRW scales and lapses on the positive root.
Write `y=b/a`, `c=Nf/Ng`, `r=alpha+beta*y`, `s=alpha+beta*c`;
physical time is `dT=Ng*s*dt`, physical scale `Ae=a*r`.
Use the full source-aware equations, and `rho+p>=0` for the canonical
source. Do not impose separate conservation of the induced g/f sources.

The additional, explicit **all-mode formal principal contract** includes
the following necessary TT/vector subset (no scalar health is assumed proved):

1. Both tensor principal speeds relative to the composite matter cone are
   no greater than one, at every time in the interval.
2. The relative vector has nondegenerate positive inertia and strictly
   positive formal principal speed at every time. Equivalently in the
   regular all-k vector chart used here, Xi>0 and cV,eff²>0.

Then `c=y`, `Q>0`, `B=0` and y is constant throughout the interval, where
`Q=P-alpha*beta*r^2*p`, `B=Ng*bdot-Nf*adot`. Moreover,

    dH_eff/dT = -alpha*r*(rho+p)/(2G) <= 0.

Consequently no contraction-to-expansion transition of any degeneracy
can occur on this interval. There is no Hubble division, simple-root
assumption, fixed parent mass/duration identification, or prescribed
scalar potential in this implication.

A separate pointwise consequence on any pressure-branch slice Q=0 is:
either a tensor principal cone exceeds the physical matter cone, or the
vector speed vanishes (for regular Xi>0). Scalar calculations cannot repair
this simultaneous formal-principal requirement for that slice.

This is **not** a theorem excluding every composite or bimetric EFT,
arbitrary matter couplings, finite-band controlled matching, or UV
completion. Massive tensor principal cones need not govern a light-only
EFT below a verified heavy threshold. No such threshold, background cutoff,
subcutoff amplification, quantum front velocity or positivity amplitude is
computed here. Relaxing the explicit speed requirement changes the class;
it is not a proved pathology or a proved solution. P8 remains open.

Proof and countercontrols: [notes/proof.md](notes/proof.md).
Source boundary: [notes/sources.md](notes/sources.md).
Report: [certificates/composite-cones.json](certificates/composite-cones.json).
