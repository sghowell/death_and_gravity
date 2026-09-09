# Cosmological proof with the nonminimal field cost retained

## Source and geometry, with no photon field substitution

This branch reuses only the pinned A.17/A.18 geometric identities,
history tube, sampler norms and the direct global-FLRW square argument.
Its scalar field inequality and physical reference are freshly derived
in field.md and anomaly.md. For an actual admitted SEE source,

    R_UU=3(H'+H^2)=Lambda-kappa(E_scalar+E_other).

If E_other>=ell and Sigma=Lambda_+ + kappa max(-ell,0), the field
inequality and the one-sided physical Wick-square cap give

    integral R_UU f^2 <= Q[f],
    Q[f]=7 kappa hbar/(144pi^2)||L_H f||^2
          +(kappa Phi_*^2/3)||G_H f||^2
          -kappa integral E_conf f^2+Sigma||f||^2.         (1)

All source coarsenings are upward. The beta_S reference is not removed.
An additional gravitational R^2 term needs its own effective source
budget; it is not identified with or silently absorbed into the
reference-state amplitude cap.

Let r=1/100 and p(v)=3v^2-2v^3. On the actual past use
u(s)=p((s+r tau)/(r tau)); on the hypothetical future use
g(s)=p((tau-s)/tau). Their value and first derivative match at zero;
at the two outer endpoints both traces vanish. The joined f is C1
and H2_0, not C-infinity. It is used only if the smooth positive
global metric actually reaches tau, so there is a compact smooth
strip around its support. The field inequality extends by H2 density.

The exact geometric identity is

    R_UU f^2=3(f'-Hf)^2-3f'^2+(3Hf^2)'.

Write K=3H(0) and J[g]=integral_future(3g'^2+R_UU g^2).
Since u,u'>=0 and H<=-h/tau on the actual past,

    3 integral_past[(u'-Hu)^2-u'^2]
        >=[3h+39h^2r/35]/tau.

Thus (1) implies

    tau(J[g]+K)
        <=18/5 +tau Q[f] -[3h+39h^2r/35].                 (2)

For the actual A.18 tube, h=19/10, d=(21/10,9,70,800).
The proof of those bounds is a geometric comparison of one actual
function and its first three derivatives, not independent choices
of jets. No initial pointwise SEC estimate is inserted into (2).

## Exact scalar and anomaly costs

Set q_v=3v0^2/4+3v1/2. The proper sampler moments give

    tau^3||L_H u||^2
       <=[7/2+(11/5)d0*r+(2/3)q_d*r^2]^2/r^3 = P,
    tau^3||L_H g||^2
       <=[(7/2)(1+2c0)+(25/12)q_c]^2 = F.

For the future estimate Hubble jets have relative powers
(tau-s)^(-j-1); the cubic's quadratic endpoint zero makes the
weighted integrals finite. The relevant exact unit moments are

    integral p^2=13/35, integral p'^2=6/5,
    integral p''^2=12, integral (p'/v)^2=12,
    integral (p/v^2)^2=13/3.

The rational upper radicals used above are 7/2,11/10,2/3,25/12;
their squared comparisons are strict. These norms are purely geometric.
The scalar coefficient relative to delta is 7/18, NOT the photon value.
The anomaly factor remains 8/2880=1/360, with the fresh V_S from
anomaly.md. Therefore

    C_S(beta_S)=(7/18)(P+F)
        +(13r/12600)V_S(d)+(13/1080)V_S(c)
        =C0+Cbeta |beta_S|.

For the stated distinct d,c these are exactly

    C0=373108140471263/75000000,
    Cbeta=3345309847373/10500000,

both less than 5,000,000. Hence the quantum/reference part of
tau Q is at most 5,000,000 delta(1+|beta_S|).
No bound on |beta_S| independent of tau is needed.

## Fresh field-strength cost

The second proper-clock term is G_H f=f'-3Hf/2.
On the past, triangle inequality yields

    tau ||G_H u||^2
       <=[sqrt(6/5)+(3d0*r/2)sqrt(13/35)]^2/r
       <[11/10+d0*r]^2/r = Cphi_P.

On the future, integral (p(v)/v)^2=4/5, so

    tau ||G_H g||^2
       <=[sqrt(6/5)+(3c0/2)sqrt(4/5)]^2
       <[11/10+(27/20)c0]^2 = Cphi_F.

The physical a^2 conversion in field.md is what makes this a
proper-field-amplitude bound with these constants. Replacing it
by a flat field amplitude or dropping the first derivative of a
would change the theorem.

The exact values are

    Cphi_P=1256641/10000,
    Cphi_F=169/4,
    Cphi=1679141/10000.

As a separate polynomial check, expanding the squared absolute
bounds before using coarse radicals gives

    past <=6/(5r)+3d0/2+117d0^2*r/140,
    future <=6/5+27c0/10+9c0^2/5.

The coarse expressions exceed these by, respectively,

    1/(100r)+7d0/10+23d0^2*r/140,
    1/100+27c0/100+9c0^2/400,

both positive. This independently verifies that the first-derivative
state cost has not been underestimated.

The source weight is ||f||^2/tau=13(1+r)/35=1313/3500.
Inserting all three costs into (2), the margin under the exact gates is

    M >=2009079/350000-18/5-1/20
         -(1679141/10000)/5000 -(1313/3500)*5
       =63325013/350000000
       =9/50+325013/350000000 >9/50.                       (3)

cosmology.py checks the symbolic identities; independent.py rebuilds
every rational cost, gate, thermal constant and comparison using
Fraction with no SymPy or imports of the production expressions.

## Global conclusion and exact conditional quantifiers

If the comoving normal reached tau, the future identity would give

    J[g]+K=3 integral_0^tau(g'-Hg)^2 >0.

Equality would require g(s)=a(s)/a(0), impossible at g(tau)=0
because a(tau)>0. But (2),(3) make this quantity negative. Thus
the upper endpoint of the stipulated global I_s is at most tau.
The inextendible comoving normal in M has finite proper length.
All timelike curves from s=0 satisfy d proper_length<=ds,
hence have total length at most tau in this global spacetime.

There is no compact-Cauchy theorem import for the noncompact R^3.
No sampler is continued through an absent endpoint. The separate
future geometry and field bounds are only required IF tau is reached.
This distinction is essential for the actual thermal control, whose
amplitude grows near its own earlier endpoint.

With kappa=8piG, delta=ell_P^2/(pi tau^2); the scalar quantum gate
is compatible with cosmological tau for every fixed finite beta_S.
The field gate separately requires kappa Phi_*^2<=3/5000 and the
source gate requires Sigma tau^2<=5. These are explicit conditional
physical budgets, not measured properties of the universe. The
theorem does not establish endpoint EFT control, a fundamental
singularity, or inextendibility in every larger spacetime.
