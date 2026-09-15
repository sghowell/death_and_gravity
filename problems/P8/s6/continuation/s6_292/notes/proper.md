# Entire pair, contact and external-factor reduction

## Full tensor algebra before integration

Let all three matter momenta be outgoing. For any pair write
k_i^2=m_i, k_j^2=m_j, z=k_i.k_j,
a=k_i.l, b=k_j.l and d_i=l^2+2a,d_j=l^2-2b.
The complete on-shell/off-shell stress tensors are
T_i=2k_i k_i+k_i l+l k_i-eta*a,
T_j=2k_j k_j-k_j l-l k_j+eta*b.
The harmonic projector contracts as
T_i:T_j-tr(T_i)tr(T_j)/(D-2).
Contraction using the full Gram matrix of(k_i,k_j,l) gives
N_D=4z^2-4m_i*m_j/(D-2)+4z(b-a)-2z*l^2
   =N0_D+2z*l^2-2z(d_i+d_j).
This is a complete polynomial identity, not an asymptotic soft rule.
Independent literal component contractions in D4,5,6 verify the
pair numerator, its entire constant soft term and the contact trace.

Division by l^2*d_i*d_j therefore reduces the whole pair to four
masters. With the inherited raw negative-Euclidean C0 convention,
deltaGamma_pair/g=-
 [N0_D*C0_ij+2z*B0_pair-2z*B0_on_i-2z*B0_on_j]
 /(16pi^2 kappa).
The sign follows from the literal product(ig)(-iT)^2 i^3, the
i*C0/(16pi^2) loop integral and division by the tree ig; it is not
selected by demanding a later real/virtual cancellation.

The metric-cubic contraction is
eta.P.T_i=-4m_i/(D-2)+d_i-l^2.
Its complete graph contributes
+[4m_i/(D-2)*B0_on_i+A0_i]/(16pi^2 kappa).
The leftover massless tadpole is scaleless and IR integrable. No
finite bubble polynomial or evanescent tensor numerator is deleted.

## All contacts and external poles

Momentum conservation implies sum_(j!=i)z_ij=-m_i.
The exact on-shell masters and derivative in D=4+2e obey
A0_i=m_i*(1+2e)/(1+e)*B0_on_i,
Sigma'_i=m_i*(3+2e)/(1+e)*B0_on_i/(16pi^2 kappa).
Since the inverse propagator is p^2-m_i+Sigma, each external factor
is1-Sigma'_i/2. Summing the three pairs, three contacts and three
external factors gives the entire selected result
deltaT/g=-[
 sum_pairs N0_D*C0_ij+sum_pairs2z*B0_pair
 +sum_legs m_i*(1+2e)/(2(1+e))*B0_on_i
 ]/(16pi^2 kappa).

The pair data are
(mu,mu,(n-2mu)/2),
(n,mu,-n/2),
(n,mu,-n/2).
The last two entries are separate diagrams. No expansion in mu/n or
replacement by massless external momenta enters this formula.

## Separate UV and IR origins

In the inherited raw convention the B0 UV pole is-1/e. Pair triangles
are UV convergent. The full pair bubble terms have combined UV
coefficient+sum m_i, the contact bubbles/tadpoles have-3sum m_i,
and external normalization has+2sum m_i, in units1/(16pi^2 kappa e).
These cancel for the selected proper-plus-LSZ sum.

The derivative part of the on-shell self energy instead has separate
IR pole+m_i/(16pi^2 kappa e). Thus its external sqrtZ has IR pole
-m_i/(32pi^2 kappa e). The raw combined Sigma' coefficient-3m_i/e
contains UV and IR pieces and is not the IR pole. Keeping their
origins separate prevents an incorrect factor or sign in the pair.

The complete proper real IR coefficient is
+4mu(1-mu/n)*atanh(beta)/beta,
and the three external factors give-mu-n/2.
Their sum is-R3. The light-light pair additionally leaves
+i*pi*V/(n*beta), in the same units1/(16pi^2 kappa e).

UV cancellation here is not a full-amplitude cancellation or a
physical cubic beta function: the separate H-metric and local matching
terms of notes/source.md remain. As a scope cross-check only, the
massive quartic calculation of Gonzalez-Martin and Martin,
https://arxiv.org/pdf/1711.08009 section3.1 equations3.3--3.7,
also distinguishes proper-plus-external cancellation from non-1PI
metric terms. That paper is not a derivation of this unequal-mass
cubic formula; all tensor algebra above is independent.
