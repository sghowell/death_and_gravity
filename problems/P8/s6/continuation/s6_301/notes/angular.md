# Full fixed-ball dimensional continuation

Use the S300 exact physical external momenta. Four massive scalars have
mass squared one; any finite number of outgoing null legs has total
energy R<=1/8. Incoming COM scalar energy is in [5/4,2]. All hard momenta
are kept in four physical dimensions. The additional soft graviton
uses the D=4+2e projector and angular phase convention of S296.

Let n be a unit vector in R3 and t in [0,1]. Embed each spatial external
momentum into R4 with last component zero and define

v_t=(sqrt(t)*n,sqrt(1-t)), P_t=I-v_t*v_t^T,
T_t=sum_A (P_t*pvec_A)*(P_t*pvec_A)^T /
                 (p_A^0-sqrt(t)*pvec_A.n).

All momenta in this formula are outgoing. Equivalently, use future
momenta and the incoming/outgoing signs. Whole momentum conservation
makes the covariant current transverse. Its D-dimensional contraction is

F_e(t,n)=tr(T_t*T_t)-tr(T_t)^2/(2+2e).

The equality follows by eliminating the temporal components of a
conserved symmetric tensor; the certificate checks it for generic
conserved tensors and an actual recoiled radiative state. Physical
integer-dimensional transverse projection and the explicit rational
D dependence fix the continuation. This does not assign a noninteger
number of independent polarizations. For noninteger D, F_e inside the
ball is not assumed nonnegative.

The normalized projection of a spatial unit sphere in 3+2e dimensions
onto its first three coordinates gives t~Beta(3/2,e), n uniform on S2.
For positive integer 2e this follows by ordinary spherical coordinates;
the beta density defines the dimensional continuation for e>0. It is
a positive normalized measure on a fixed domain, with normalization

h(e)=Gamma(3/2+e)/[Gamma(3/2)*Gamma(e)]=e*A(e),
A(e)=Gamma(3/2+e)/[Gamma(3/2)*Gamma(1+e)].

Set H(t)=mean_S2 F_0(t,n), U(t)=mean_S2 tr(T_t)^2, K0=H(1).
The exact continuation, including the evanescent trace term, is

K_e = K0+c(e)*U(1)
 + e*A(e)*integral_0^1 sqrt(t)*(1-t)^(e-1)
       *[H(t)-K0+c(e)*(U(t)-U(1))] dt,
c(e)=e/[2(1+e)].

The subtracted form uses the EXACT normalization of the beta measure.
Once the integrability in notes/bounds.md is proved, differentiation
at e=0 is justified and gives

K1 = U(1)/2
   + integral_0^1 sqrt(t)/(1-t)*[H(t)-H(1)] dt.

Neither the first term nor the radial difference may be omitted.
K0 is the full S300 radiative soft index, including every null/null pair.
When there is no radiation this K1 agrees with S296's independently
derived massive Feynman-parameter expression.

A null leg exactly parallel to n at t=1 gives a bounded but potentially
direction-dependent representative. Choosing its value on that
zero-area set does not change any angular integral. No pointwise
continuous collinear extension is asserted. The tests exercise the
bounded representative without treating it as a physical directional
limit.

For stable numerical calibration, put t=(1-u^2)^2. The radial weight
sqrt(t)*(-dt)/(1-t) is 4*(1-u^2)^2/[u*(2-u^2)] du.
This removes the artificial square-root endpoint in a preliminary
quadrature. It is a numerical chart, not the proof of convergence.
