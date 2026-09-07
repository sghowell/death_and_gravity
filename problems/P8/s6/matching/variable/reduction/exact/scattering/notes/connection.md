# S6.35: coefficient-16 own-f connection and fixed-window comparison

This note concerns the classical **homogeneous own-f tensor equation**, with
the physical metric prescribed, in the same off-shell probe and positive exact
branch as S6.33. It is not the coupled-relative coefficient-80 equation of
S6.21. The canonical remainder bound supplied by `potential.py` is an explicit
premise below. A free exterior is introduced only to name finite-endpoint
coordinates; it is not a physical continuation or vacuum. No source or state
has been integrated out implicitly.

## 1. Exact reference and two independent fundamental bases

The dimensionless inner comparison is

\[
y_{xx}+\frac{16}{1+8x^2}y=0.
\]

With \(z=\sqrt8x=\sinh t\) and \(y=\sqrt{\cosh t}\,\psi\), its literal
pullback is

\[
\psi_{tt}+[\rho^2+V(t)]\psi=0,
\qquad \rho=\frac{\sqrt7}{2},\quad V=\frac34\operatorname{sech}^2t.
\tag{1}
\]

Indeed the equation before removing the first derivative is
\(y_{tt}-\tanh(t)y_t+2y=0\). The old coefficient 80 would instead give
\(\rho^2=39/4\); those constants are not imported.

Put \(a=-1/4+i\rho/2\), \(b=\bar a\). An exactly normalized real basis is

\[
y_e={}_2F_1(a,b;1/2;-8x^2),\qquad
y_o=x\,{}_2F_1(a+1/2,b+1/2;3/2;-8x^2).
\tag{2}
\]

At zero its Cauchy matrix is the identity, hence \(W_x(y_e,y_o)=1\).
The argument is nonpositive, so it does not cross the ordinary Gauss cut
\([1,\infty)\). Conjugate numerator parameters make (2) real. Direct
substitution gives Gauss equations with \((a+b,ab)=(-1/2,1/2)\) and
\((a+b+1,(a+1/2)(b+1/2))=(1/2,1/2)\), respectively.

For positive large x, write
\(y_e=E_+x^{1/2+i\rho}+\bar E_+x^{1/2-i\rho}+O(x^{-3/2})\), with the
corresponding differentiated expansion, and similarly with \(O_+\) for
\(y_o\). The oscillatory remainders have magnitude \(O(x^{-3/2})\), not a
claimed finite-x numerical bound. The exact connection coefficients are

\[
E_+=\frac{\sqrt\pi\,\Gamma(i\rho)8^{1/4+i\rho/2}}
 {\Gamma(-1/4+i\rho/2)\Gamma(3/4+i\rho/2)},\qquad
O_+=\frac{\Gamma(3/2)\Gamma(i\rho)8^{-1/4+i\rho/2}}
 {\Gamma(1/4+i\rho/2)\Gamma(5/4+i\rho/2)}.
\tag{3}
\]

Their Wronskian obeys \(4\rho\operatorname{Im}(E_+\bar O_+)=1\).
At negative x parity leaves E unchanged and reverses O. Thus, independently
of the Jost construction below, the radial connection is
\(S\operatorname{diag}(1,-1)S^{-1}\), where
\(S=\left(\begin{smallmatrix}E_+&O_+\\\bar E_+&\bar O_+\end{smallmatrix}\right)\).
The connection formula behind (3) is
[DLMF 15.8.2](https://dlmf.nist.gov/15.8.E2), converted from its **normalized
Olver** \(\mathbf F={}_2F_1/\Gamma(c)\) to the ordinary function used here.
The extra denominator Gamma factors in the transformed functions must also
be converted; simply replacing bold F by ordinary F would be wrong.

## 2. Jost data and genuinely nonzero mixing

The right Jost function is

\[
\psi_R(t)=e^{i\rho t}{}_2F_1\left(-\tfrac12,\tfrac32;
 1-i\rho;\frac{1-\tanh t}{2}\right).
\tag{4}
\]

It is asymptotic to \(e^{i\rho t}\) at plus infinity. At minus infinity,
the ordinary Gauss connection in
[DLMF 15.10.21](https://dlmf.nist.gov/15.10.E21) gives

\[
\psi_R\sim Ae^{i\rho t}+Be^{-i\rho t},\qquad
A=\frac{\Gamma(1-i\rho)\Gamma(-i\rho)}
 {\Gamma(3/2-i\rho)\Gamma(-1/2-i\rho)},\qquad
B=i\,\operatorname{csch}(\pi\rho).
\tag{5}
\]

Here \(w_t=-2w(1-w)\), so the direct Jost pullback gives numerator
parameters \((-1/2,3/2)\), not a reflectionless integer-index potential.
Gamma reflection and recurrence yield

\[
|A|^2=\coth^2(\pi\rho),\quad |B|^2=\operatorname{csch}^2(\pi\rho),
\quad |A|^2-|B|^2=1.
\tag{6}
\]

For the mathematical convention of unit left incoming wave and no right
incoming wave, transmission is \(1/A\), reflection is \(B/A\), and reflected
Wronskian flux is \(\operatorname{sech}^2(\pi\rho)>0\). This flux ratio is
not \(|B|^2\). Numerically, only as corroboration, \(|B|\simeq0.0313499\)
and the reflected flux is about \(0.000981852\). No particle-production or
physical spatial-scattering interpretation is needed.

There is a simple exact lower margin. Polynomial division shows

\[
0<\int_0^1\frac{x^4(1-x)^4}{1+x^2}\,dx=\frac{22}{7}-\pi.
\]

Together with \(\rho<4/3\), this gives \(\pi\rho<17/4\).
The exponential series gives
\(e<49/18<11/4\) and \(e^{1/4}<4/3\), hence
\(e^{\pi\rho}<e^{17/4}<(11/4)^4(4/3)=14641/192<80\).
Consequently

\[
|B|=\frac1{\sinh(\pi\rho)}>\frac1{40}.
\tag{7}
\]

No numerical special-function enclosure is needed for (7).

## 3. Time frequency, radial log waves, and orientation

The left-to-right coefficients of \(e^{\pm i\rho t}\) transform by

\[
T_t=\begin{pmatrix}\bar A&-\bar B\\-B&A\end{pmatrix},\qquad\det T_t=1.
\tag{8}
\]

With \(x=u/\sqrt\delta\), the coefficients of
\(|u|^{1/2\pm i\rho}\) instead transform by

\[
T_u=\begin{pmatrix}
-\bar B&\bar A e^{i\varphi_\delta}\\
A e^{-i\varphi_\delta}&-B
\end{pmatrix},\qquad
\varphi_\delta=2\rho\log\frac{4\sqrt2}{\sqrt\delta},\qquad\det T_u=-1.
\tag{9}
\]

The left positive time-frequency wave has the *negative radial* exponent.
That reversal explains the determinant in (9); canonical Cauchy evolution
still has determinant +1. In real cosine/sine amplitudes, put
\(A e^{-i\varphi_\delta}=a_\delta+i b_\delta\) and
\(\beta=\operatorname{csch}(\pi\rho)\). The real radial map is

\[
\begin{pmatrix}a_\delta&b_\delta+\beta\\
b_\delta-\beta&-a_\delta\end{pmatrix}.
\tag{10}
\]

Its singular values are \(|A|\pm\beta=\coth(\pi\rho/2),
\tanh(\pi\rho/2)\). The logarithmic delta phase cannot be discarded, even
if B were artificially set to zero. Formula (9) is the reference's
asymptotic radial connection, not the exact finite-endpoint frame used next.

## 4. Actual finite-window endpoint coordinates

Let \(u=T/\tau\), \(k=b^3/N>0\), and \(Y=\sqrt{k}\,Q\), using the actual
S6.33 branch. With zero prescribed tensor q, its homogeneous equation is

\[
Y_{uu}+\left[\frac{16}{\delta+8u^2}+R_\delta(u)\right]Y=0.
\tag{11}
\]

The **entire** canonical pump is in \(R_\delta\). Suppose the separately
verified bound \(|R_\delta|\le44\) holds on \(|u|\le1/100\),
\(0<\delta\le1/100\). We use only the narrower
\(0<\delta\le10^{-6}\). Define

\[
r=\sqrt{u^2+\delta/8},\quad
t=\operatorname{arsinh}\frac{\sqrt8u}{\sqrt\delta},\quad
Y=\sqrt r\,\psi,\quad Z=(\psi,\psi_t/\rho)^T.
\]

At either endpoint the literal map from \((Q,Q_u)=(Q,\tau Q_T)\) is

\[
Z=\begin{pmatrix}r^{-1/2}&0\\
-u/(2\rho r^{3/2})&\sqrt r/\rho\end{pmatrix}
\begin{pmatrix}\sqrt{k}&0\\k_u/(2\sqrt{k})&\sqrt{k}\end{pmatrix}
\binom Q{Q_u}.
\tag{12}
\]

In particular, omitting \(k_u\) changes the actual endpoint data. Constant
M factors cancel in the transfer, but the physical time and kinetic factors
in (12) do not. The transformed potential is

\[
\psi_{tt}+[\rho^2+V(t)+E_\delta(t)]\psi=0,
\qquad E_\delta=r^2R_\delta.
\tag{13}
\]

Set \(L=1/100\), \(t_*=\operatorname{arsinh}(\sqrt8L/\sqrt\delta)\).
Extend (13) **mathematically** by the free equation
\(\psi_{tt}+\rho^2\psi=0\) outside \([-t_*,t_*]\), matching both Cauchy
data continuously. No spacetime metric, potential profile, source, or state
is asserted outside the actual interval.

The normalized plane-wave frame

\[
V_\rho(t)=\frac1{\sqrt2}
\begin{pmatrix}e^{i\rho t}&e^{-i\rho t}\\
i e^{i\rho t}&-i e^{-i\rho t}\end{pmatrix}
\]

is unitary. If \(F_\delta(t_*,-t_*)\) is actual Z evolution, then its
finite-window coefficient matrix is exactly

\[
\mathcal S_\delta=V_\rho(t_*)^*F_\delta(t_*,-t_*)V_\rho(-t_*).
\tag{14}
\]

There is no hidden factor two in the matrix norm: the same constant
normalization multiplies coefficients at both ends. This is also exactly
the scattering matrix of the fictitious free-exterior equation.

## 5. Energy bound, complete tails, and an explicit surviving margin

For any real potential W, the generator of Z is
\(\left(\begin{smallmatrix}0&\rho\\-\rho-W/\rho&0\end{smallmatrix}\right)\).
Its symmetric part has eigenvalues \(\pm |W|/(2\rho)\), so

\[
\left|\frac{d}{dt}\log\|Z\|\right|\le\frac{|W|}{2\rho}.
\tag{15}
\]

For the reference, \(\int_{\mathbb R}V=3/2\), hence every forward or
backward interval evolution has norm at most
\(B_0=e^{3/(4\rho)}<2\). A fully rational proof uses \(\rho>5/4\) and
the series bound
\(e^{3/5}<1+3/5+(3/5)^2/2+[(3/5)^3/6]/(1-3/20)=1549/850<2\).
The factor one half in (15) matters; the old large-rho shortcut is not used.

Compare the free-exterior equation with the full reference (1). The
potential difference is E on the central interval and \(-V\) on both
tails. Since \(dt=du/r\), its integrated central absolute value is bounded by

\[
\int_{-t_*}^{t_*}|E|dt\le44\int_{-L}^{L}r\,du
 <44\left[L^2+2L\left(\frac1{2000}\right)\right]
 =\frac{121}{25000}.
\tag{16}
\]

Here \(\sqrt{\delta/8}<1/2000\) for the chosen delta range. Both reference
tails together give

\[
\int_{|t|>t_*}Vdt=\frac32\left(1-\frac L{\sqrt{L^2+\delta/8}}\right)
 \le\frac{3\delta}{32L^2}\le\frac3{3200}.
\tag{17}
\]

The first inequality follows by rationalizing the numerator; no tail
asymptotic is substituted at a finite point. Therefore the entire
perturbation integral divided by rho is

\[
J<\frac45\left(\frac{121}{25000}+\frac3{3200}\right)
 =\frac{2311}{500000}<\frac1{200}.
\tag{18}
\]

Duhamel's identity and (15) give the **full-line coefficient** comparison

\[
\|\mathcal S_\delta-T_t\|\le2B_0(e^{J/2}-1)
 <4(e^{1/400}-1)<\frac4{399}<\frac1{80}.
\tag{19}
\]

For completeness, at an integration time s the later reference factor
and earlier perturbed factor together are bounded by
\(B_0\exp[\tfrac12\int_{-\infty}^{s}|\Delta W|/\rho]\).
Integration against \(|\Delta W|/\rho\) yields the factor
\(2(e^{J/2}-1)\), not an extra square of \(B_0\). Integrability justifies
passing finite Cauchy endpoints to the free-wave coefficient limits. At
the last step, the geometric-series bound is
\(e^{1/400}<1/(1-1/400)\).

Define \(B_{\rm eff}=-[\mathcal S_\delta]_{21}\). Combining (7) and (19),

\[
|B_{\rm eff}|>\frac1{40}-\frac4{399}
 =\frac{239}{15960}>\frac1{80}.
\tag{20}
\]

This is a fixed physical interval of duration \(2\tau/100\), not an inner
interval whose duration shrinks with delta. The error bound is deliberately
nonzero as delta tends to zero: direct comparison still retains the bounded
outer-region remainder. A stronger vanishing-error theorem would require
the actual outer dressing and is not asserted here.

## 6. What is and is not established

Equations (12), (14), and (20) give a specified, uniformly nonzero
homogeneous transfer mixing in exact canonical endpoint coordinates, under
the separately checked 44 remainder premise. Incoming data are prescribed.
With zero homogeneous initial data and zero prescribed q the solution is
identically zero; nothing here creates a response without data or forcing.

The original q is an off-shell prescribed physical-metric tensor, not the
free chi matter field or a conserved external matter source. Extending a
canonical ODE by free waves is not proof of a parent vacuum, an asymptotic
physical S matrix, a cutoff, a fixed spatial-frequency EFT band, or a
quantum state. No retarded inverse has been inserted into a single-copy
action. The action's other sectors, nonlinear terms, source matching, and
original P8 completion remain separate obligations.
