# Complete virtual pole, physical branch and inclusive reference

Use the S288 physical continuation of the frozen S278 pair integral:
M_s=(-atanh(beta)+i*pi/2)/(s*beta), beta=sqrt(1-4mu/s).
For a spacelike channel a, M_a=one quarter of
integral_0^1 dx/[mu-a*x(1-x)].
Do NOT substitute real s>4mu separately into the principal square roots
of S278's compact analytic radical formula: taking the cut boundary after
that substitution can reverse the real branch. S288 already certifies
the correct i0 sheet and the retained imaginary Coulomb coefficient.

The Feynman-current pair integral at physical incoming s has
L_s=integral dx/[mu+(s-4mu)x(1-x)]=4atanh(beta)/(s*beta).
Writing N_s=((s-2mu)/2)^2-mu^2/2 and analogous crossed pair N_t,N_u,
the complete ordered real kernel is
K0=2mu+4[N_s L_s-N_t L_t-N_u L_u]=-2Re Bsoft.

Under x=(1+z)/2 and reflection, the pointwise K0 density equals twice
the S288 positive attenuation density. This connects the current and
analytic descriptions without dropping any self terms or imaginary phase.
S288 bounds therefore also imply0<=K0<=Q^2/(5mu), Q=s-4mu.
At mu1,E<=2 this is<=144/5, stronger than the coarse800 used in S296.
Use the existing coarse bound unless a sharper constant is needed.

## Local UV subtraction before IR pairing

The S293 full raw quartic graph sum is
C Gamma(-e)(4pi nu^2)^(-e)[-2T_e+E_e]/(16pi^2 kappa).
At zero T0=Bsoft,E0=5mu/3. The local UV pole is
-5Cmu/(48pi^2 kappa e); its covariant renormalizing counterterm removes
that UV pole before the IR comparison. Its independent finite part is
NOT chosen. The remaining infrared pole is C Bsoft/(8pi^2 kappa e).

The S294 whole cusp/box/noncusp sum has pole
g^2 Htree[sum V J0-2mu]/(16pi^2 kappa e)
=g^2 Htree Bsoft/(8pi^2 kappa e).
Its full metric endpoint and RH UV counterterm are handled covariantly;
independent finite RH/heavy-residue/higher matching stays in the hard
amplitude, not in the soft conversion. The off-shell vertices, contact
bubble and residues are not individually interpreted as the full pole.

Thus the complete selected renormalized virtual pole is
A0 Bsoft/(8pi^2 kappa e), A0=C+g^2 Htree.
The real Born-rate interference is2Re of this divided by A0,
so it cancels K0/(8pi^2 kappa e).
The imaginary Coulomb pole remains an amplitude phase, not a real-rate
divergence. Nothing here removes it from fixed-transfer gravity/Regge.

The hard reference is the finite coefficient after the EXACT stated
S278 soft division, including any independent finite hard matching.
The selected inclusive formula must be truncated consistently:
1+2Re(deltaA_hard/A0)+Delta_soft+R_real,
not an unqualified |A_hard|^2 which includes an uncomputed loop square.
A bound on Delta_soft+R_real does not fix deltaA_hard or establish
positivity, complete physical unitarity, full-source matching or P8 closure.

The fixed-domain continuity proof removes pure metric terms only AFTER
the full Ward identity. After q projection all anisotropic vectors span
at most two dimensions, where tr(AB)-tr(A)tr(B)/(2+2e) is a positive
bilinear bounded by Frobenius norms. The normalized sphere pushforward is
rho_e=(1/2+e)/pi*(1-r^2)^(e-1/2)<=5rho0/4 on a fixed unit disk.
Equivalently its folded perpendicular coordinate v has normalized weight
(1+2e)v^(2e)dv dphi/(2pi). No noninteger helicity count is introduced.

With J_phase,e=J_phase,0*(rprime/r0)^(2e),1-J_phase,e<3omega,
the exact-minus-soft kernel is bounded by
[B^2+(128B+12288)/omega]/kappa,B330000.
The radial omega^(1+2e) and the fixed-domain weight provide an integrable
majorant. Dominated convergence therefore yields exactly the D4 S295
physical remainder, not an assumed finite prescription.
