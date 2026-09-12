# Complete local Euler-Hessian factorization and boundary

Use the source-pinned S225 complete finite local density, not a reduced or retuned action. Let y=(w,c), T=(3,-1)y and B y=(S,W). Its remainder is

Lrem=alpha rDG+6alpha(TD SG+TG SD)+gamma TD TG,

with alpha andgamma in FORMULATION.md. Production verifies this exact bridge against the entire frozen S225 remainder, and independently checks alpha andgamma against its original finite R_old^2, R_old and volume coefficients.

The only spatial-transfer terms in rDG are
-20q wD wG+4q(cD wG+cG wD).
Since W=-qw-c'', these are precisely the spatial terms in

WD(10wG-4cG)+(10wD-4cD)WG.

Integrating its two second-time-derivative terms by parts with alpha retained generates
alpha[10(wD'cG'+cD'wG')-8cD'cG']
+alpha'[cD'(10wG-4cG)+(10wD-4cD)cG'].

Subtracting those from the actual temporal part of alpha rDG gives the statedH andJ matrices. The6alpha T S terms supply the first row ofA, and gamma TD TG suppliesV0.

More precisely, if Lfact is the bilinear expression

(B yD)^T A yG+yD^T A^T B yG
+yD'^T H yG'+yD'^T J yG+yD^T J^T yG'+yD^T V0 yG,

then the exact POINTWISE identity is

Lrem-Lfact
=D[alpha(cD'(10wG-4cG)+(10wD-4cD)cG')].

The boundary is exhibited, not silently called zero. The resulting Euler-Hessian identity

Rloc=B* A+A^T B+D* H D+D* J+J^T D+V0

holds as a local differential-operator identity for arbitrary smooth alpha,gamma andh. The independent test computes the Euler derivative directly from the full original bilinear density and detects deletion ofJ. A second residual verifies the whole spatial-transfer coefficient.

On compact variations the displayed total derivative is the usual variational boundary. The local Euler operator then extends uniquely to causal distributions, including derivative atoms at the initial time. It does not introduce an alternative boundary prescription, discard a one-current contact or use a background equation to delete terms.
