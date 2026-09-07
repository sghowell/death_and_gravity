# P8-S6.15.AUXILIARY: actual-background no-bounce theorem

This is a new screen of the specified auxiliary parent, not a change to
S6.14 or a closure of P8. Its main conclusion is independent of tensor
cones, a vacuum, large-time assumptions and the auxiliary tensor inverse.

## Action and physical dictionary

The four-dimensional action is

\[
 S=-\sum_{i=g,f}{G_i\over2}\int\det(e_i)R[g_i]
 -2\int\det(u)\left[B+\operatorname{tr}\left(u^{-1}\sum_i p_i e_i\right)\right]
 -2\sum_i b_i\int\det(e_i)+\epsilon S_m[h,\psi],\qquad
 h=u^T\eta u,\quad g_i=e_i^T\eta e_i.
\]

All coefficients are constant. The positive Einstein coefficients are
\(G_i>0\); \(p_i,b_i,B\) are real, of either sign. The index set
\(I=\{i:p_i\ne0\}\) is nonempty. There is no Einstein or other derivative
term for \(u\), no direct matter on the two Einstein metrics, and no
additional interaction operator. This includes the separately named
endpoint-cosmological extension of S6.13, not arbitrary trimetric actions.

Use the original P8(b) convention \(+---\),
\(R_B=-6(DH+2H^2)\), \(G_{B00}=3H^2\), and EH density
\(-G_iR_B/2\). This is not the opposite A/FK curvature convention.
\(h\) is the actual physical matter metric, not its vacuum approximation.

On a connected regular common spatially flat homogeneous interval, take
\(e_i=\operatorname{diag}(n_i,a_i,a_i,a_i)\) and
\(u=\operatorname{diag}(n_u,a_u,a_u,a_u)\), with all entries positive and
smooth, in the parent coframe/root domain. Define

\[
 dT=n_u dt,\quad d\tau_i=n_i dt,\quad
 R_i={a_u\over a_i},\quad c_i={R_i n_i\over n_u},\quad
 H_u={d\log a_u\over dT},\quad H_i={d\log a_i\over d\tau_i}.
\]

The actual isotropic \(u\)-only source obeys
\(n_h=\epsilon(\rho+p)\ge0\). Homogeneous positive-field-metric canonical
scalars with \(\epsilon>0\) are sufficient, including zero rolling
velocity. The main theorem only needs this isotropic physical-source NEC,
not a full matter perturbation model. Matter conservation is part of an
actual solution, not an extra interaction-NEC premise.

## Main theorem and controlled-background corollary

The full Einstein equations and their Bianchi identities imply, without
division by a Hubble rate,

\[
 H_i=R_iH_u\ (i\in I),\qquad R_i'=R_i(1-c_i)H_u,
 \qquad K=\sum_{i\in I}{G_i\over R_i^2}>0,
\]

\[
 2K H_u'-H_uK'=-n_h,\qquad
 \left({H_u\over\sqrt K}\right)'=-{n_h\over2K^{3/2}}\le0.
\]

A contraction-to-expansion sign change of the actual physical Hubble rate
is therefore impossible anywhere on the connected regular interval,
including arbitrarily degenerate bounces. No sign of an individual link,
cone bound, flat vacuum, tail completeness or auxiliary tensor denominator
is used. A disconnected Einstein metric is not added to \(K\).

If a proposed matching map identifies the actual \(h\) and its proper
clock with CD, then the target
\(H_{CD}(T)=4T/(\tau^2+T^2)\) at \(T=\pm L\tau\) cannot be approximated
at both endpoints with normalized error strictly below
\(4L/(1+L^2)\). At \(L=1/2\) the necessary maximum error is at least
\(8/(5\tau)\). This concerns an actual full-parent solution and explicit
physical metric/clock dictionary; it can constrain even an asymmetric
light-only proposal that claims that background match. It does not infer
such a dictionary from a reduced action or identify the parent's matter
content with the old DHOST clock plus free \(\chi\).

## Separate supplementary results

The proof also records an independent, weaker cone/three-slice argument;
an actual same-action flat-vacuum calibration; and the mixed-sign,
positive-vacuum-FP denominator bound. These are separate corollaries with
extra hypotheses, not hidden premises of the main theorem.

The background no-bounce result includes auxiliary tensor singularities
\(\sum_i p_i/R_i=0\). A tensor inverse or physical full-cone inference
does not. The identically disconnected case \(I=\varnothing\) is excluded:
it has genuine arbitrary-\(u\) controls and no invertible auxiliary
construction. Additional operators, an auxiliary kinetic term,
non-NEC/quantum sources, anisotropy/spatial curvature, singular coframes or
a different matter frame change the hypotheses. No cutoff, universal UV
exclusion, scalar/vector health, all-parent exclusion, or P8 closure follows.

Read the [proof](notes/proof.md) and [source audit](notes/sources.md).
