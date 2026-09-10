# Actual one-loop forward matching target

Use S6.110's same canonical polynomial potential with light mass one:
V=Phi^2/2 + M H^2/2 + G H Phi^2/2 + lambda4 Phi^4/24.
Here M denotes the heavy mass **squared**, g=G^2, D=M-2,
g=2 lambda D^3, and lambda4=g(3D-2)/D^2.
All actual rational parameters are imported, not independently retuned.

The exact tree amplitude is A0=-lambda4+g sum_z(M-z)^(-1).
At forward transfer t=0, u=4-s and b2=[(s-2)^2]A(s,0).
The full tree coefficient, not only a derivative truncation, is b2_tree=4lambda.

Bound the complete renormalized one-loop b2 for this model in the fixed
scheme stated in notes/subtraction.md. Work on the complex disc |s-2|<=1,
strictly below every relevant subthreshold singularity. Include the
t-channel angular dependence and all mixed heavy/light contributions.
Do not impose a radial cutoff or expand the heavy inverse in y/M.

The once-fixed light pole and residue are from S6.112. The target is a
quantitative one-loop error relative to the tiny tree coefficient, not
an inference based solely on small dimensionless couplings. The latter
would fail: the naive lambda4^2/tree coefficient ratio is enormous.

Higher-loop errors, full absorptive/contour control, finite gravity and
a common rolling-bounce parent are not established by a finite-order
forward coefficient. Original V/G/B and P8 remain open.
