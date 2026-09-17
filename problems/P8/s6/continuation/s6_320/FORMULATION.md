# All-finite relative-energy complex tubes and Cauchy bound

Retain the original selected covariant action, canonical normalization,
original couplings and the complete untruncated tree source. Use the
same recoil map at fixed real E,u and emission directions, with
5/4<=E<=2 and positive associated elastic transfers t,u.
The leaf polarizations are fixed complex spatial TT tensors of unit
Frobenius norm. No matching coefficient or real point API changes.

Choose any generic positive center energies wi with sum wi<=1/8.
Continue only the energies zi and the same analytic recoil branches.
The closed polydisc is

 |zi-wi|<=epsilon*wi, epsilon=1e-12.

The complete finite-N tree is holomorphic on a neighborhood of this
closed polydisc, and for every N>=1 obeys

 |M_N(z)|/A0
 <=3*n^2*N!*(2e17)^N/[kappa^(N/2)*product_i(wi)],

where n=HEAVY_MASS2 and A0=Am+AG is the same positive, unexpanded,
REAL elastic Born normalization. The wi on the right are center
energies, not complex zi. Pure-soft angular and hard-angle approaches
are uniform on the generic domain; exact internal pole points remain
excluded.

For the analytic scaled amplitude
G_N(z)=product_i(zi)*sqrt(rho(z))*M_N(z)/A0, put
K_N=6*n^2*N!*(4e17)^N/kappa^(N/2).
Every nonnegative integer multiindex alpha then satisfies at the center

 |partial^alpha G_N|
 <=product_i(alpha_i!)*K_N/product_i[(epsilon*wi)^alpha_i].

The proof extends the real-center weighted current norm to complex
energies. Positive pair-invariant sums keep every pure-soft inverse
away from zero. The complex conserved temporal inverse, arbitrary
Rosen identity, reflection grading and all-valence vertex estimates
close a nonnegative coefficient recurrence. Unchanged hard-core
ownership and complex recoil/cut estimates give the full bound.
Ordinary multivariable Cauchy supplies the derivative estimate.

The current recurrence is
D1=1,
Dn=25600*sum_(labeled root partitions,k>=2)
 (k+1)!32^(k+1)13^k*product D_child,
with Dn<=2*(2e13)^(n-1)*n!. Its pure-current tube allows relative
shifts up to1/100; the full hard source uses the smaller1e-12 tube.

A complete three-ray/four-tree source has a nonzero transverse root
residue47089/2601 at a nearby signed-energy pole. With a harder
spectator this pole lies in a naive global-total-W disc. The control
concerns that invalid block-holomorphy shortcut, not a positive-real
physical divergence or every full hard-amplitude residue.

The valid relative tube does NOT reach wi=0. The derivative bound
therefore does not establish an integrable all-N overlap subtraction,
infrared cancellation, inclusive probability or physical-series
convergence. Hard/evanescent matching, quantum state, unitarity,
absolute complex Regge and common-parent bounce remain open.
Original V/G/B/P8 remain OPEN; scoped P8(a) is unchanged.
