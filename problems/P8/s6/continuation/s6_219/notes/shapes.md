# Complete scalar cutoff shapes

Let p=|P|, u=n.Phat, q=3u^2-1, a=(1+t^2)^2. Use unnormalized I/I,I/S,S/I with S=diag(-1,-1,2). Unit Frobenius channels multiply the first by1/3 and each cross by1/sqrt18.

The complete corrected S216 product is first convolved with its geometry, before angular averaging or endpoint deletion. The raw current coefficient of r^(1-q) is the sum over all j+degree=q and every source jet r<=j. There are45 radial/source-jet rows in the three scalar channels through total grade4. In particular the LC,CL,CC j1 nonlog products are nonzero; the tracefree S217 special cancellation is not extended to them.

The leading azimuthal densities are f0_tt=1/(2a) and f0_tS=f0_St=-(3u^2-1)/(4a). The complete centered second coefficients of source jets0,1,2 are
\[
\begin{aligned}
 f^{tt}_{20}&=-m^2a/4-p^2(1+u^2)/(16a)+3(t^2-1)/2,\\
 f^{tt}_{21}&=-a'/8,\qquad f^{tt}_{22}=-a/8;\\
 f^{tS}_{20}&=-m^2aq/4-p^2(u^2-2)/(8a)-3(3t^2+1)q/4,\\
 f^{tS}_{21}&=-a'q/8,\qquad f^{tS}_{22}=aq/16;\\
 f^{St}_{20}&=-m^2aq/4-p^2(u^2-2)/(8a),\\
 f^{St}_{21}&=a'q/4,\qquad f^{St}_{22}=aq/16.
\end{aligned}
\]
These include j0/degree2,j1/degree1,j2/degree0. The uncentered first grade obeys the exact translation identity. Centering changes only source jet0 by
\[
 -p^2(1-u^2)[f_0-u f'_0+(1-u^2)f''_0]/8.
\]

The original intersection shell is integrated before expanding: its angular endpoint is p/(2K), not zero. S209's exact uncentered hemisphere and positive grazing-strip functionals are applied to every q0..3 row. All30 possible cutoff-power/source-jet outputs are checked. The complete K^2 and K^0 terms cancel. Independently, centered integration gives
\[
 A_3=-p(8\pi^2)^{-1}\int_{-1}^1|u|f_0du,\quad
 A_{1r}=(4\pi^2)^{-1}\int_{-1}^1|u|
 [\delta_{r0}p^3(3-5u^2)f_0/16-p f_{2r}/2]du .
\]
Thus A3_tt=-p/(16pi^2 a), A3_tS=A3_St=p/(64pi^2 a), and
\[
\begin{aligned}
 A^{tt}_{10}&=\frac p{\pi^2}[m^2a/32+p^2/(64a)+3(1-t^2)/16],\\
 A^{tt}_{11}&=pa'/(64\pi^2),\quad A^{tt}_{12}=pa/(64\pi^2);\\
 A^{tS}_{10}&=\frac p{\pi^2}[m^2a/64-5p^2/(256a)+3(3t^2+1)/64],\\
 A^{tS}_{11}&=pa'/(128\pi^2),\quad A^{tS}_{12}=-pa/(256\pi^2);\\
 A^{St}_{10}&=\frac p{\pi^2}[m^2a/64-5p^2/(256a)],\\
 A^{St}_{11}&=-pa'/(64\pi^2),\quad A^{St}_{12}=-pa/(256\pi^2).
\end{aligned}
\]
The two ordered operators obey the formal proper-time transpose, including first and second derivatives of coefficients, but are not equal.

For the complete tensor form, each degree-d azimuthal coefficient has polynomial degree at most d+4<=8. Trace shifts add no momentum, and the new constraint geometries have no higher degree. Therefore the original per-row Cauchy bound B_j controls their angular C2 norm by6121 B_j rho^-d. The exact small-transfer shell/grazing bound8192 and full shape-functional bound64 give S209's original sums unchanged. For large transfer |P|>K/4, bound the whole original lost UV integral and each shape separately; q4 retains log(K/m)<=K/m. All resulting external degrees are at most6. The complete conversion remainder is below1e40 M X46/K for K>=2m.

The q4 shell remains nonzero and decaying; its bulk logarithm is not removed. Independent high-precision ray integration includes q0..4, all three source jets in each channel, and the exact original positive strip. Subtracting the displayed K^3 and K terms gives the expected1/K asymptotics. These fixtures supplement the uniform proof; they do not extrapolate a fixed-transfer expansion to all momenta.
