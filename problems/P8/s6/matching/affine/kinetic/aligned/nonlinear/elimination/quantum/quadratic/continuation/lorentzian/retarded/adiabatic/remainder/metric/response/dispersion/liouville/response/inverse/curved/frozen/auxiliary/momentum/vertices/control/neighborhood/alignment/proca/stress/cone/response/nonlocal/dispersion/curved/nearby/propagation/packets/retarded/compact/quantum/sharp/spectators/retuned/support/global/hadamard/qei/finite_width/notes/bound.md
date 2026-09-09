# The finite-width inequality and its constants

Let g be real smooth with compact support in (-1/100,1/100).
Use the actual comoving observable and
rho_test=(O_u^2+a^-2 sum_i O_i^2)/2.
For the reference family in notes/state.md and every state in
the same generalized Hadamard class, S6.100's difference argument
and the estimates here prove

    integral g(u)^2 [<rho_test>_omega-<rho_test>_ref](u) du
      >= -(hbar/kappa) [(C_L+C_R)||g||_2^2
                       +C_A (T_g^2+S_g^2)].

The finite constants, with K=32 and C_E from notes/remainder.md, are

    C_L=153269043200 exp(196/25)/(3 pi^2),
    C_A=exp(1/25)/(160 pi^3),
    C_R=257 C_E^2/(16 pi^2 32^8).

T_g,S_g are the explicitly tabulated derivative-norm combinations
in notes/fourier.md. These constants contain no unevaluated
reference-state integral or unknown high-frequency threshold.

For clarity, the normalization is as follows. The state modes
carry hbar/(2 kappa); the test energy carries 1/2; the radial
spatial measure is k^2 dk/(2 pi^2); and the reference functional
has d_alpha/pi. The high-frequency approximate term therefore
starts with hbar/(8 pi^3 kappa), before its initial-mode norm,
Fourier decay, or exact/error split is bounded.

Below K, extend the alpha integral to the whole real line and
use Parseval. The reference's initial regular covariance norm
is <16, its evolved norm is <=16*182 exp(196/25), and its
time/spatial observation constants are 16 and 1. Also
integral_0^K k^2 sqrt(1+k^2) dk
<=K^4/4+K^3/3. These give C_L.

Above K, apply |exact|^2<=2|approximate|^2+2|error|^2 to each
derivative mode BEFORE integrating. The factor two is retained
in BOTH terms. Three integrations by parts give C_A.
For the error, full-line Parseval and
integral_K^infinity k^(-9) dk=1/(8K^8) give C_R.
The two spatial directions orthogonal to momentum do not create
extra copies: sum_i k_i^2=k^2 is kept exactly.

All norms are in the normalized proper-clock variable u, on the
same action/source normalization as S6.99--100. This is an explicit
finite-width majorant, deliberately much weaker than an optimal
bound. It is not a claimed small remainder relative to S6.100's
leading short-sampling coefficient.
