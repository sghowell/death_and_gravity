# Exact repartition and no double-counting

At the same finite two-leg mask, both decompositions represent the complete endpoint row:

E5 = Q203 + Pfin + PUV = Qnew + U5.

Therefore

Qnew = Q203 + Pfin + PUV - U5,

and

Jactual,K = (R_actual-unit,K+F_unit,K+Qnew,K)
           + (C_unit,K+U5,K).

No finite endpoint term, actual-state contribution or contact is discarded. The mathematical subtraction is changed only by an exact regrouping; the physical finite covariant prescription is not altered.

The bound must be rebased from S6.201's independently controlled R_actual-unit+F_unit, which is below 2e48 M[D]N61[Gamma] with tail 2e52 M[D]N61[Gamma]/K. Adding Qnew to the already combined S6.205 known piece would count the old finite endpoint terms twice, and is not done.

Using ||D||L2 <= M[D] and Y^2 = N61^2+X46^2 gives

|known_new| < 5e48 M[D]Y[Gamma],
|known_new-known_new,K| < 5e53 M[D]Y[Gamma]/K.

The unchanged norms are

M[D]^2 = ||D||L2^2+||grad D||L2^2,
N61[Gamma]^2 = sum_(r=0)^6 ||partial_t^r Gamma||L2^2+||grad Gamma||L2^2,
X46[Gamma]^2 = sum_(r=0)^4 ||(1-Delta)^3 partial_t^r Gamma||L2^2.

Both external metric factors supply 4/kappa. The endpoint displays are 8e-752 and 16e-747/K, while the new known actual piece has 2e-751 and 2e-746/K.

The fixed local action Hessian from S6.204 is still only a target. It is not added as an independently established quantum term. The exact UV-symbol integral and complete contact must still be matched in the original prescription, including all finite/divergent coefficients and subleading regulator artifacts.
