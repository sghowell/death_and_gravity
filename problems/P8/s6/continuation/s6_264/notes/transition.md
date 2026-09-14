# Complete canonical chart transport, not an endpoint reset

The original S251 map is
T_clean=shear_central * T_central * shear_outer^(-1),
where each shear changes its original canonical momentum by
-a^3(B+B^T)y/2. Both time-boundary terms are necessary.
T_central exchanges the first position and momentum with coefficient
2a^3q; its time derivative is retained in the original dynamics.

In balanced phase use D=diag(sqrt(P),sqrt(P),1/sqrt(P),1/sqrt(P)).
The ENTIRE matrix M=D T_clean D^(-1), with q=P^2/a^2, has leading
diagonal L=diag(-E/Theta,1,-Theta/E,1). This leading matrix alone
is not the finite-P estimate.

The exact symplectic identity M Omega M^T=Omega gives
M^(-1)=-Omega M^T Omega. The code checks both this identity and
the full inverse multiplication. It avoids an expensive generic
symbolic inversion without approximating any matrix entry.

Put w=1/P and multiply each entry of(M-L)/w, and each entry
of(M^(-1)-L^(-1))/w, by
a^6 E^4 Theta^4 Delta(a^2 w^2)^2.
All32 results are exactly reconstructed polynomials and given
outward coefficient majorants in the inherited overlap box.
At the switch |E|,|Theta|>=1/4, |E|<=1/2, |Theta|<=2,
1<=a<=2 and Delta>=1/8, so every remainder-entry coefficient
is below10^20. The all-entry sum of either full remainder is
less than16*10^20/P<1.

The all-entry sum of each principal diagonal is at most12.
Consequently both full Euclidean balanced transition norms are
less than16. Energy metrics on both sides have eigenvalues in
[10^-4,100]; conversion therefore costs at most
16 sqrt(100/10^-4)=16000 in root energy.
The exact same state passes through the switch. There is no
instantaneous ground-state reset, new preparation or discarded
symmetric time boundary.

Independent finite fixtures reconstruct both shears from the original
generator and Legendre momentum relation, then compare the complete
matrix and its inverse norm. These fixtures are checks on the algebra,
not a replacement for the evaluated all-P>=10^64 estimate.
