# Causal curvature-coordinate inverse without an inverse-transfer pole

The channel matrix is

B=[[D^2+2q/3,-D^2/3],[-q,-D^2]],
det B=-D^2(D^2+q).

For Re(lambda)>0 neither lambda nor lambda^2+q vanishes forq>=0. On zero-past histories use I^2, the causal inverse ofD^2, and Lq^-1, the causal inverse ofD^2+q. The latter has kernel s_q(t)=sin(sqrt(q)t)/sqrt(q), with s_0=t.

Subtract one third of W from S to obtain (D^2+q)zeta=S-W/3. Then D^2 c=-W-qzeta. The exact resolvent identity

q I^2 Lq^-1=I^2-Lq^-1

eliminates the apparent inverse-transfer issue and gives

zeta=Lq^-1(S-W/3),
c=(Lq^-1-I^2)S-(Lq^-1/3+2I^2/3)W.

Both inverse products are checked as rational matrix identities. Their causal kernels give

R_q(t)=[[s_q,-s_q/3],[s_q-t,-s_q/3-2t/3]].

There is no singularity atq0, and R_q(0)=0. With x=cos(sqrt(q)t), the first derivative matrix is[[x,-x/3],[x-1,-x/3-2/3]]. Its squared Frobenius norm is(20x^2-14x+13)/9. The exact maximum on[-1,1] is47/9. The simpler entrywise enclosure55/9 is also below25/4 and is what the uniform bound uses:

||R'_q(t)||<5/2, ||R_q(t)||<5t/2.

This is valid for everyq>=0. No derivative of the cosine, which would grow withq, is used. The nonzero-q kernel has a removable square-root notation, not an unbounded inverse-q multiplier.

Independent forced time-domain wave solves atq0,1e-12,1 and10000 reconstruct both curvature channels and the alternative double-primitive cancellation. Further exact matrix solves at small and largeq check both gauge signs and detect deletion of the shear channel.

These are constant-coefficient flat reference coordinate identities. In a curved time-dependent background, coefficient multiplication does not automatically commute with time primitives or wave resolvents. The same algebra cannot be copied there without a separate ordered derivation.
