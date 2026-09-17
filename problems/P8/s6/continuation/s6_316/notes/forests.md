# Pair branches, matchings and exact finite inclusion-exclusion

Choose a labeled external Phi root. A pair branch is an Einstein cubic
vertex incident to two prescribed external h leaves and one internal
h edge. A given h leaf has one adjacent vertex. It cannot belong to
two different pair branches, so shared-label pair intersections are
empty. Every set of pair branches present in a graph is a matching.

For a fixed matching M of k pairs, collapse each cubic branch and its
internal edge to a single distinguished h leaf. This is a bijection
from trees containing M to the original tree inventory at N-k labels.
The distinguished leaves can be off shell: graph counting does not
depend on their momenta or wave tensors. Restoring the branches is
unique. The extra cubic vertices add no contact/heavy marker.

There are m(N,k)=N!/[2^k k! (N-2k)!] matchings of size k. This follows
by selecting2k labels, ordering them, and dividing by the pair
orientations and by permutations of the k pairs.

For a tree with p actual pair branches, its weight in the alternating
sum over nonempty contained matchings is
sum_(k=1)^p (-1)^(k+1)*binomial(p,k)=1 if p>0, and0 if p=0.
This proves the amplitude union identity graph by graph. It is valid
with the original signed/complex graph weights and kinematics; no
positivity or convergence is needed for this finite sum.

Taking formal graph weights gives
R_N=sum_k (-1)^k m(N,k) T_(N-k).
To obtain the graph EGF, set m=N-k. The coefficient of y^(m+k) in
T_m/m!*(y-y^2/2)^m is
(-1)^k T_m/[2^k k! (m-k)!], exactly the required R_N/N!.
Thus R(y)=T(y-y^2/2) as a formal graph series.
This scalar composition is NOT a composition rule for amplitudes,
which still depend on the collapsed momenta and off-shell tensors.

The first no-pair inventories, with C,g formal markers, are:
N0: C+3g^2+3;
N1: 5C+21g^2+21;
N2: 33C+177g^2+177;
N3: 274C+1770g^2+1770;
N4: 2758C+20646g^2+20646.
At C=g=1 these are7,47,387,3814,44050.

An independent labeled root-cut recursion suppresses exactly two-h-leaf
currents, without using the closed formula, and recovers these five
inventories. Explicit matching enumeration and the EGF identity are
also checked to higher finite orders. The combinatorial bijection and
Boolean identity, not those calibrations, prove the arbitrary-N result.
