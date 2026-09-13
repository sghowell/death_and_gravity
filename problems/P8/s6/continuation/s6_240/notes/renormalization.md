# Complete scalar dimensional counteraction and finite matching

## Full-dimensional mode coefficients

Use spatial D=3-2epsilon, scale a, omega²=n+p²/a², z=(p²/a²)/omega², lambda=omega'/omega=-H z and U_D=D H'/2+D²H²/4. Here R=2D H'+D(D+1)H²; this explicitly fixes the curvature convention. For a polynomial in H,H',... and z, the total derivative is

dt=sum_j H^(j+1) partial_(H^j)-2H z(1-z)partial_z.

The complete normalized WKB expansion begins W=omega+P2/omega+P4/omega³, with

P2=-U_D/2-dt(lambda)/4+lambda²/8,

P4=-P2²/2-dt²(P2)/4+5lambda dt(P2)/4+[dt(lambda)/2-3lambda²/2]P2.

Put B2=dt(P2)-2lambda P2 and d0=(D H+lambda)/2. After removing the common a^-D, the physical rho and pressure coefficients multiplying omega^(1-2j), for j0,1,2, are

e0=1/2, e2=d0²/4,
e4=(P2²+d0 B2-d0²P2)/4,

p0=z/(2D),
p2=[(2-2z/D)P2+d0²]/4,
p4=[(2-2z/D)P4-(1-2z/D)P2²+d0 B2-d0²P2]/4.

All D dependence is retained. Independent tests expand the full phase-cancelled energy and pressure readouts, not just the Riccati equation, and recover every displayed coefficient.

## Full radial integral and evanescent variation

Changing to physical momentum removes a^-D against the measure. For each polynomial coefficient at order j, define

F(epsilon)=sum_r c_r(3-2epsilon) rising(3/2-epsilon,r)/Gamma(r+j-1/2).

Let q=2-j and pref=8sqrt(pi)(-1)^q/q!. In units n^(2-j)/(64pi²), at mu1 the exact radial pole and finite component are

pole=pref F(0),

finite=pole[H_q-log n]+pref F'(0).

This comes from the full D-dimensional gamma-function integral and its Laurent expansion. It includes the D dependence in the physical pressure and in the tensor contractions.

The full covariant scalar pole density before64pi²epsilon is

n²-nR/3+2a2,

a2=R²/72+(Riemann²-Ricci²)/180,

Ricci²=D²(H'+H²)²+D(H'+D H²)²,

Riemann²=4D(H'+H²)²+2D(D-1)H^4.

The compact divergence in a2 is a variational boundary term; there is no assumption about an infinite-time boundary flux. For any of these local densities f(H,H'), complete lapse variation gives

rho_D=-f+H f_H+(H'-D H²)f_H'-H dt(f_H'),

P_D=-rho_D-dt(rho_D)/(D H).

The last expression simplifies to a polynomial and extends through H=0. Independent lapse AND scale-factor Euler-Lagrange calculations, retaining lapse derivatives and scale-factor second derivatives, verify these expressions for1,R,R²,Ricci²,Riemann² in spatial dimensions3,4,5.

Crucially, the pole counteraction is varied in D+1 dimensions BEFORE taking D=3. Its expansion contains the finite counterstress -2 partial_D T_pole. Subtracting it adds2 partial_D T_pole to the component finite part. All six component/order identities then match the stress of the finite scalar local density

[(3/2-log n)n²+(log n-1)nR/3-2log n a2]/(64pi²).

Subtracting only the component pole at D3 leaves extra finite curvature terms and is not this covariant prescription. The test retains that explicit negative control. This scalar calculation does not borrow the different finite vector evanescence from S176.

## Complete finite stress and Ward identity

At D3 write ell=log n. The coefficient at order j multiplies n^(2-j)/(64pi²):

j0: rho=ell-3/2, P=3/2-ell;

j1: rho=-2H²(ell-1), P=2(3H²+2H')(ell-1)/3;

j2: rho=ell(6H²H'+2HH''-H'^2),
P=-ell(18H²H'+12HH''+9H'^2+2H''')/3.

Each order satisfies the full exact Ward identity. The complete adiabatic MODE subtractions also obey dt(e_j)+(1-2j)lambda e_j+D H p_j=0 at arbitrary spatial D. Restoring a^-D omega^(1-2j) proves the pointwise mode Ward identity at each order before integration; no unexamined radial boundary term is needed. The four-dimensional identity a2=Weyl²/120-Euler/360+R²/72 is checked; the dimensional variation was not replaced prematurely by this four-dimensional identity.

The renormalized rho_H_ref and P_H_ref are defined as the full exact-SLE radial integral minus these complete fourth-order adiabatic mode subtractions, PLUS the displayed covariant finite local stress. The bounds prove absolute convergence with five time derivatives. The full finite density and the mode Ward identity give the renormalized Ward identity. This is a fixed state and prescription, not state-dependent normal ordering.

In particular the flat heavy Gaussian pressure is p_v_H=(3/2-log n)n²/(64pi²), with rho_v_H=-p_v_H. It is not zero. A mass-one minimally coupled light scalar in the formal flat limiting-action vacuum has p_v_Phi=3/(128pi²) at the same mu1. Both are retained in the new profile. This says nothing about quantizing the interacting light sector on the curved clock reference.

The scalar adiabatic/DeWitt-Schwinger equivalence provides primary context, while the full D-dependent finite matching above is derived here; see [del Rio and Navarro-Salas](https://arxiv.org/abs/1412.7570).
