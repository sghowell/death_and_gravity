# Full current source, not an off-clock jet truncation

background.data binds R and F to parent.fixed_functions from S253, with
the full S238 additions, both S240/earlier fixed profiles, and all three
vacuum constants. These are literal full functions. The finite source is
the unchanged S250 physical source

    JH=G kappa u^2 (1-X)^8 exp(-10^420 X^2)/2 - J1 V(X),
    J1=-G/(32 pi^2), V=(1-X)^1024/[X^1024+(1-X)^1024], X=N^-2.

It is not divided by kappa in the physical KG equation. No source
subtraction is changed. The affine-clock map leaves the lapse and shift
unchanged and gives a_phys=R^-1/4 a_hat. The full coefficient families are

    M=R^1/4, U=R^-3/4, C3=R^3/4/2, Cchi=R^-1/4,
    B=-U R_u/(2N)-I,
    Fhat=U[F+9 R_u^2/(16 R N^2)]-I_u/N,
    I_N=3 U R_u R_N/(4 R N), I(u,1)=0.

The primitive is the exact integral with this initial value. Every mixed
derivative in the stored full parameter substitution acts on it. It is
not independently varied or replaced by a finite lapse polynomial.
The same action-to-ADM map proved in S253 applies before restriction.

At Hbar=0 write the physical heavy fluctuation H=sqrt(kappa) h. Dividing
the action by kappa a_hat^3 gives its entire quadratic addition

    Z hdot^2/2 - Y q h^2/2 - N U nH h^2/2
       + [partial_N(N U JH/sqrt(kappa)) n
          +3 N U JH v/sqrt(kappa)] h,
    Z=U/N, Y=N Cchi, q=P^2/a_hat^2.

The nonzero linear density N U JH h/sqrt(kappa) remains. Nothing sets the
off-shell heavy equation to zero. The independent lapse/volume variation
re-enters the full kinetic, gradient, mass and source density and checks
both its first and second field derivatives. The source has no spatial
derivatives in unitary clock gauge. At Hbar0 its matter charge and rate
mixing vanish, but its lapse-H and volume-H contacts do not.

Consequently the exact four-mode block in coupled.py covers this entire
source by its second matter field: c_H=w_H=0, d_H=partial_N(N U JH/sqrt(kappa)),
v_H=3 N U JH/sqrt(kappa), mass_HH=-N U nH. These and the full M1 coefficients
are listed explicitly. The more general symbolic two-matter block also
retains arbitrary such contacts; they are only shown lower order after
the whole reduction. The order-eight/1024 source filtration from S250
licenses S253's finite clock vertices, not deleting these terms at fixed
nonzero displacement. The huge but fixed mass and localizer affect
remainder constants and thresholds; no smallness is inferred from them.
