# Actual bounce covariance symbols from the complete regular chart

The original crossing uses the S251 physical canonical swap

Qb=-Pv/(2a^3 q), Pb=2a^3 q Qv,

where Qv=sqrt(kappa)v. Retain the full symmetric boundary shear before
counting ultraviolet orders. At the ACTUAL reference bounce H=Theta=0,
the complete central B matrix is exactly zero, even with Jc,A,Tcorr kept.
Thus its cleaned momenta equal the unsheared central momenta at this
time. This fact is checked against the full finite-q S221 matrix, not
an asymptotic truncation.

The whole current values are

a0=1, E0=-1/2, ell0=1/10, F0=1199/800,
J0=243/160+rho_fixed(0)-7 pressure_fixed(0)/8.

The two fixed normalized profiles are retained as their original functions.
The clock germ used for these finite jets agrees with the entire parent
to order eight; no off-clock function has been replaced. The code checks
J against the existing full reference coefficient, with both profiles.

Let cs=sqrt(F/J)>0. The complete central principal pair is

K0=[[2J/E^2+ell^2,ell],[ell,1]],
G0=[[2F/E^2+ell^2,ell],[ell,1]].

The existing matrix
R=[[E/sqrt(2J),0],[-ell E/sqrt(2J),1]]
satisfies R^T K0 R=I and R^T G0 R=diag(cs^2,1).
Literal equality with the S251 normalization is checked.

In normalized cleaned coordinates Q=R a^(-3/2)x and
p=R^(-T)a^(3/2)pi. The previously written S251 full two-cone
classical-symbol proof gives, for its ACTUAL fixed reference state,

Cov(pi,pi)=diag(k cs/a,k/a)/2+O(1).

This is not a newly minimized instantaneous state. It follows from the
arbitrary-depth positive-frequency comparison, full symbol remainders,
regular chart overlaps and smooth-sampling anomalous-block estimate
already supplied there. At the bounce the exact coordinate rows are

v=Pb/(2sqrt(kappa)a^3 q),
n=[ell E Ps-(E+3Tcorr/(2q))Pb]/(2J sqrt(kappa)a^3).

No Tcorr term is dropped from these exact rows. Substituting R gives
the leading metric momentum row
[ sqrt(2J)/E, ell ]/(2sqrt(kappa)a^(3/2)q)
and lapse row
[-1/sqrt(2J),0]/[sqrt(kappa)a^(3/2)].
The remaining lapse row term is O(q^-1) times the same complete metric
momentum combination. With q=k^2/a^2, therefore,

Cvv(k)=(2J cs/E^2+ell^2)/(8kappa k^3)+O(k^-4/kappa),
Cvn(k)=-cs/(4E kappa a^2 k)+O(k^-2/kappa),
Cnn(k)=cs k/(4J kappa a^4)+O(1/kappa).

The big-O constants depend on the fixed complete reference and are not
claimed explicitly bounded here. Both scalar branches contribute to Cvv;
the matter term ell^2 is not deleted. At the actual E0=-1/2, Cvn has a
positive leading coefficient. This is consistent with the exactly
vanishing LINEAR v-n commutator at Theta0=0.

## Radial coincidence asymptotics

Use the original R3 Fourier convention d^3k/(2pi)^3 and a fixed positive
lower endpoint k0. The radial measure is k^2 dk/(2pi^2). For a sharp
upper diagnostic endpoint Lambda,

Cvv_[k0,Lambda] =
 (2J cs/E^2+ell^2)/(16pi^2 kappa) log(Lambda/k0)+O(1/kappa),
Cvn_[k0,Lambda] =
 -cs Lambda^2/(16pi^2 E kappa a^2)+O(Lambda/kappa),
Cnn_[k0,Lambda] =
 cs Lambda^4/(32pi^2 J kappa a^4)+O(Lambda^3/kappa).

The finite lower-endpoint terms are included in the exact principal
integrals. At the bounce these give
Cvv log coefficient=(8sqrt(JF)+1/100)/(16pi^2 kappa),
Cvn quadratic coefficient=cs/(8pi^2 kappa).
Multiplying by -9/4 gives the displayed gauge-mean jet coefficients.

These are divergences of unrenormalized COORDINATE contractions.
The cancellation in the physical-volume identity is retained at every
finite band. It does not by itself define an interacting renormalized
physical observable. In particular one may not identify Lambda with the
heavy mass, infer a Wilsonian matching scale, turn asymptotic remainders
into numerical finite errors, or equate these coordinate jets with the
actual interacting metric/lapse mean.

Independent full noncommuting positive 4x4 comparison matrices check the
leading powers and coefficients. They are diagnostics of the symbol
algebra, explicitly not replacements for the actual S251 preparation.
