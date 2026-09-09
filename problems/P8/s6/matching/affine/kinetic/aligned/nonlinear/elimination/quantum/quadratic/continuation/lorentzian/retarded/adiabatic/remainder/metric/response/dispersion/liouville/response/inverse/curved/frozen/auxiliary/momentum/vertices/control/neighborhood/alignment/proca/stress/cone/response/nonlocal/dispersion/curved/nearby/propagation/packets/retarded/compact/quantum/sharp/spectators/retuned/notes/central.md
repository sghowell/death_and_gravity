# Central positive-matter family and why the larger margin was selected

Keep the total margin epsilon symbolic. At u=z=0 let y=N^2, with
N>0 and 9/10<=X=N^-2<=11/10, hence 10/11<=y<=10/9. Derive the
constraint and time jets at fixed canonical phase BEFORE substituting
a solved matter charge. The exact central expressions are

    W=w^2=-[ (1000 epsilon-1125)y^2
              +(2849-400 epsilon)y-600 epsilon-1725 ]/100,
    D=(3000 epsilon-3375)y^2+(2849-400 epsilon)y+600 epsilon+1725,
    B=(2200 epsilon-2475)y^2+(8849-400 epsilon)y-1800 epsilon-5175,
    A=(8 epsilon-9)y^2-60y+24 epsilon+69,
    PIV=-D/(400N^3), c_clock^2=B/D, 1-c_clock^2=100A/D.

W is a rescaled matter square; the physical conserved charge square
at central R=1 is W/N. A negative W is not assigned a real field state.
The lapse velocity and physical Hubble value vanish at the symmetric
center, but their actual derivatives are retained in the acceleration.

## Explicit weaker-margin negative control

For epsilon=1/10000 and N=10001/10000, exact rational substitution gives
W>0, PIV<0, c_clock^2>1 and dHphysical/dtau>0. The approximate values
are 0.00880207, -2.99545, 1.00103 and 3.99570 respectively; the proof
uses exact signs, not the decimals. The rational flow is analytic near
this point, with M, PIV and W nonzero. The reconstructed real analytic
solution is consequently a strict local classical bounce and retains a
fast clock on some open interval. No quantitative interval or old
compact response witness is transferred to this different datum.

Thus increasing the margin to 10^-4 repaired a small box but did not
remove the central-tube obstruction. This calculation motivated the
selected total epsilon=1/200, rather than silently declaring the
first favorable small-box result a general repair.

## Exact selected-parameter continuum proof

For epsilon=1/200 define

    p(y)=-100W=-1120y^2+2847y-1728,
    D(y)=-3360y^2+2847y+1728,
    B(y)=-2464y^2+8847y-5184,
    A(y)=-(224/25)y^2-60y+1728/25.

p is strictly increasing on the FULL tube because
p'(y)>=p'(10/9)=3223/9>0. At y_cap=501/500, p=2619/12500>0.
Consequently W>0 forces y<y_cap; the apparent remainder of the tube
does not describe a positive-matter central datum.

On [10/11,y_cap], D decreases and remains positive, B increases and
remains positive, and A decreases and remains positive. These claims
follow from their displayed linear derivatives and exact endpoint
values; none is inferred from sampled signs. In particular,

    c_clock^2 >= B(10/11)/D(10/11)=49753/93129 > 1/2,
    1-c_clock^2 >=100 A(y_cap)/D(10/11)
                 =129954/485046875 > 1/4000.

The positive D gives PIV<0. Together with W>0, N>0 and
M=-1/(4N^2), the physical scalar kinetic and gradient forms are
positive, as are the unchanged tensor and ordinary-Proca principal
blocks. This proof covers every positive-matter central datum in the
declared tube, not just the later chosen nearby initial lapse.

## Actual physical bounce acceleration on the continuum

The complete derivative calculation gives, for this selected parameter,

    dHphysical/dtau = P_acc(y)/Q_acc(y),
    P_acc=-(75264y^4-845488y^3+2124921y^2-1708562y+313365),
    Q_acc=25y*(-1120y^2+949y+576).

central.py derives the lapse second derivative from the unreplaced
time jets and independently checks agreement with the full new
physical-Hubble chain rule. This is not a derivative taken along the
family of constraint-parametrized charges.

Set y=10/11+(501/500-10/11)t. For each of these two polynomials,
write its degree-n power coefficients a_j and compute exact Bernstein
coefficients b_k=sum_{j<=k} a_j*binomial(k,j)/binomial(n,j).
Native reconstruction checks the polynomial equals
sum b_k binomial(n,k)t^k(1-t)^(n-k). For 0<=t<=1 these basis functions
are nonnegative and sum to one, so the polynomial is between its
smallest and largest b_k. Every actual coefficient is positive.
Their ratio gives the uniform acceleration lower

    207386958411035601/60630859375000000 > 3.

At each admitted central datum the literal coefficients and the
fixed-basepoint primitive have analytic neighborhoods: N>0, D>0 and
the primitive's integration path stays at positive lapse. The nonzero
PIV and positive W give a local analytic constraint/flow solution,
and this positive acceleration makes it a strict bounce. This is
qualitative local existence for the continuum; the quantitative
|u|<=10^-7 construction is only for the separately specified N0=1+10^-6.
No common time radius for the whole central family, whole off-center
phase tube, global nearby completeness or UV admissibility follows.
