# Entire initial state mismatch, not a reset

The input is the actual S240 SLE, prepared once before t0=-1/2. It is not the full W6 approximation or the exact KG comparison initialized from W6 at time-1. Their complete Bogoliubov errors are used only to bound the actual initial covariance.

Write Omega_star²=n+k²/16. The unchanged S240 proof gives |beta_SLE|<=10^83 Omega_star^-11 and the full exact-comparison/W6 error <=10^62 Omega_star^-13. Both include all internal momenta. A normalized Bogoliubov change acts on a pure graph by

r_new=(alpha r+beta phase)/(alpha+beta phase rsharp),
r_new-r=beta phase(1-r rsharp)/(alpha+beta phase rsharp).

Here |phase|=1, |alpha|>=1; a phase choice makes alpha positive. For |r|<1/50, |beta|<1/2, the difference is below2|beta|. This is the exact transformation, including alpha-minus-one and all denominators.

The FULL W6 frequency has the marker graph

r6(zeta)=(omega-W6+i zeta d6)/(omega+W6-i zeta d6),
d6=(theta+W6'/W6)/2, theta=3H at preparation.

For an arbitrary W, not only its truncated expansion, define

F_W=W²-omega²+zeta²[theta'/2+theta²/4+W''/(2W)-3W'^2/(4W²)].

Direct exact algebra proves

zeta r6'-2iomega r6-zeta(lambda+theta)(1-r6²)/2
=2iomega F_W/(omega+W-i zeta d6)².

S240's whole sixth-frequency residual starts at marker14. The unique coefficient recursion therefore agrees with rhat through order10 (indeed further); no physical WKB evolution is inferred.

On the inherited entire marker disk |zeta|<=Omega_star/10^8, the full normalized iterate obeys |W6/omega-1|<=10^-8, |d6|<14 and |omega|>=Omega_star/2. Consequently the normalized numerator is less than3*10^-7 and the denominator exceeds2-3*10^-7, giving |r6|<1/10 and a nonzero denominator on the FULL disk. At physical marker1 it is also below1/50.

Since the disk radius exceeds2, Cauchy bounds the ENTIRE coefficient tail from order11 by

(1/5)*(10^8)^11 Omega_star^-11 <10^88 Omega_star^-11.

Combine the full graph transformations, comparison residual and coefficient tail. With nu=sqrt(99/100)Omega_star,

|r_actual(t0)-rhat(t0)| <(2*10^83+4*10^62+10^88)nu^-11 <10^90 nu^-11.

The initial error is nonzero. Independent full six-iterate computations at1300 and1500 decimal digits for k/sqrt(n)=0,1,4 resolve a nonzero tail after multiplication by Omega_star^11 and agree beyond150 relative digits. They are diagnostics of the exact argument, not its uniform proof.

All admitted source histories have the same entire zero germ at t0. The exact initial covariance, frame and each local coefficient jet are consequently source-independent there: e1(t0)=e2(t0)=0. This does NOT set e0(t0) to zero or reminimize the state on a varied history.
