# Whole scalar ADM Hamiltonian and its metric contacts

The full minimally coupled scalar phase-space Hamiltonian is

H=N/2[pi²/sqrt(h)+sqrt(h)(h^ij partial_i H partial_j H+m²H²)]
+beta^i pi partial_i H,

with the shift term Weyl ordered, h=a² exp(Q), N=1+n_lapse and pi held canonical. An independent full Legendre transform derives it from the ADM scalar Lagrangian; no scalar or TT restriction is imposed.

At the reference use the five physical features

F=(pi/a³,partial_x H/a,partial_y H/a,partial_z H/a,mH).

For direction D=(n_D,beta_D,Q_D), put tau_D=tr Q_D and B_D=tau_D I/2-Q_D. The first HAMILTONIAN feature has diagonal blocks

M_D=(n_D-tau_D/2, n_D I+B_D, n_D+tau_D/2)

and the energy/gradient off-diagonal blocks a beta_D and their transpose. The current is J_D=-a³ F^t M_D F/2.

The entire mixed second Hamiltonian feature C_DG has energy coefficient

tau_D tau_G/4-(n_D tau_G+n_G tau_D)/2,

gradient block

(B_D B_G+B_G B_D)/2+n_D B_G+n_G B_D,

and mass-field coefficient

tau_D tau_G/4+(n_D tau_G+n_G tau_D)/2.

There is no second shift or lapse-lapse term in this canonical chart. In particular the noncommuting matrix anticommutator cannot be replaced by a trace-only expression. The instantaneous current contact has the NEGATIVE sign -a³<F^t C_DG F>/2.

The literal full nilpotent matrix exponential, determinant and inverse calculation checks every matrix entry. Independent high-precision finite mixed derivatives of the complete noncommuting matrix-exponential Hamiltonian agree as well.

For Fourier canonical variables z=(H,pi), let F_p map z to( pi/a³,ipH/a,mH). The complete first canonical kernel at incoming/outgoing momenta k,l is

K_HH=a³[(n_D+tau_D/2)m²+a^-2 k^t(n_D I+B_D)l],
K_Hpi=-i k dot beta_D,
K_piH=+i l dot beta_D,
K_pipi=a^-3(n_D-tau_D/2).

The exact bridge a³ F_(-k)^t M_D F_l=K_D(k,l) retains both shift phases and the mass/volume term.

Finally |tr Q|<=sqrt(3)||Q||F. Since sqrt(3)/2+1<2, full block norm bounds give

||M_D||op<=v(D)=|n_D|+a|beta_D|+2||Q_D||F,
||C_DG||op<=v(D)v(G).

These estimates apply to the full complexified component matrices using the Hermitian Frobenius norm; no separate physical-state positivity is assigned to complex source directions.
