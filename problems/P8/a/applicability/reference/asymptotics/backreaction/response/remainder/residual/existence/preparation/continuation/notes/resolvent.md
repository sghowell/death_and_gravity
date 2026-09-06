# Principal-sheet poles, causal inverse and exact duration bounds

This is a theorem for the isolated constant-coefficient block identified
in the full-map decomposition. Beta is fixed by the named prescription;
it is not an adjustable counterterm used to remove a pole.

## 1. Normalization and analytic symbol

For zero-past data, A.10 gives the Laplace multiplier
`4pi^2 D_gammaE=(log(s)+gamma_E)/2`. Adding d_f gives

    Lf(s)=(log(s)+beta)/2,
    Bf(s)=[log(s)+beta-c/s^2]/2,
    Rf(s)=1/Bf(s)=2s^2/[s^2(log(s)+beta)-c].

The principal logarithm has -pi<arg(s)<pi. The double primitive, not a
single primitive, produces c/s^2. The minus sign is inherited from the
positive u/(60delta) in the full trace. Changing this sign or omitting
the local finite coefficient changes the operator. The formulas are
initially Laplace formulas to the right of all poles, not a prescription
to Fourier transform a growing response along the imaginary axis.

## 2. Complete pole count

Let s=r exp(i theta) and A=log(r)+beta. A nonreal zero obeys

    A sin(2theta)+theta cos(2theta)=0,
    c/r^2=A cos(2theta)-theta sin(2theta).

There is no zero on the imaginary axis. For 0<theta<pi/2, elimination
gives c/r^2=-theta/sin(2theta)<0, impossible. On the positive real axis
`log(s)+beta-c/s^2` is strictly increasing from minus to plus infinity;
there is exactly one zero p0>exp(-beta).

For pi/2<theta<pi the two real equations reduce to

    c exp(2beta)=[-theta/sin(2theta)] exp[-2theta cot(2theta)].

The right side goes from zero to infinity. Its logarithmic derivative is

    1/theta-4cot(2theta)+4theta csc(2theta)^2
     =[(sin(2theta)-2theta cos(2theta))^2
          +4theta^2 sin(2theta)^2]/[theta sin(2theta)^2]>0.

Hence exactly one zero is in the upper left half-plane, and its conjugate
is the only lower-half zero. The negative real axis is the excluded cut;
its banks cannot be zeros because their imaginary parts are +/-pi.
There are precisely three poles on the sheet, not an unbounded list of
Lambert branches.

Put w=2(log(s)+beta). The root equation becomes

    w exp(w)=z, z=2c exp(2beta)>0.

The sheet restriction is -2pi<Im(w)<2pi. Thus the three roots are W0(z),
W1(z),W-1(z), with Im(W1) in (pi,2pi). Write p_j=exp(w_j/2-beta).
All are simple: the derivative of the denominator is
`p_j(1+w_j)`, which cannot vanish here. Their residues in Rf are

    b_j=2p_j/(1+w_j).

In particular p0,b0>0. No claim is made that exp(p0 t) alone is a
nonzero zero-past homogeneous solution: the causal zero-data solution
is unique and zero. The pole describes the forced response and the
response to any compatible initial/history numerator.

## 3. Full causal kernel, including the damped pair

The upper/lower banks of the negative cut give the density

    rho_c(r)=2/{[log(r)+beta-c/r^2]^2+pi^2}>0.

A keyhole contour applied to Rf(z)/(z-s) gives, for real s>p0,

    Rf(s)=sum_j b_j/(s-p_j)+integral_0^infinity rho_c(r)/(s+r) dr.

The outer-circle integral is O(1/log R), and the inner-circle integral
is O(epsilon^3), since Rf(z)=O(z^2) near zero. The cut integral converges
at both ends. The resulting retarded kernel is

    Kf(t)=b0 exp(p0 t)+2Re[b1 exp(p1 t)]
                +integral_0^infinity rho_c(r) exp(-rt) dr, t>0.

Its Laplace transform exists for Re(s)>p0. The cut is locally integrable
in t because its exponentially weighted integral is finite; the three
residue terms are locally integrable as well. This defines a causal
convolution on C[0,L] for every finite L.

Positivity of the **full** kernel should not be guessed from the complex
pair. It follows instead from A.10's positive normalized logarithmic
kernel j_beta, whose transform is 2/(log(s)+beta). Set

    H=(c/2)(j_beta*t),
    Kf=j_beta+H*j_beta+H*H*j_beta+... .

For sigma>p0, its weighted L1 ratio is
`Hhat(sigma)=c/[sigma^2(log(sigma)+beta)]<1`. Therefore the series
converges in weighted L1, is nonnegative and has exactly transform Rf.
Uniqueness identifies it with the contour kernel. This also proves the
unique causal inverse property without deleting a pole. Smooth input
flat near the initial endpoint has smooth flat output by differentiating
onto the input; for generic continuous data no automatic C1 gain is used.

## 4. Norms and explicit rational growth

Positivity gives the exact unweighted C0 operator norm

    C_f(L)=integral_0^L Kf(t)dt.

The supremum is the same on smooth inputs flat near zero: approximate
the constant input 1 from below away from zero and use local L1
domination. It is therefore not caused by a forbidden endpoint jump.

Integrating the residue formula, the positive cut can be dropped for a
lower bound. Since Re(p1)<0 and Im(w1)>pi>3,

    |2Re[(b1/p1)(exp(p1 L)-1)]|
        <=8/|1+w1|<8/3.

Thus `C_f(L)>[2/(1+w0)](exp(p0 L)-1)-8/3`.

Take delta=10^-14, af=5/2. Elementary inequalities suffice:
`0<gamma_E<1`, `0<log(5/4)<1/4`, `2<log10<5/2`, `1/2<log2<1`.
The first is the standard integral comparison for the harmonic-number
limit. The logarithm bounds follow from integrating 1/x, e<3, and the
positive exponential series: the sum through degree four at 5/2 already
exceeds 10. Consequently

    -53/60<beta<11/30,
    [s^2(log(s)+beta)-c] at s=10^6 < -16400000000000/3,
    [s^2(log(s)+beta)-c] at s=2*10^6 > 76900000000000/3.

The monotone positive-root equation proves `10^6<p0<2*10^6`.
Also `w0<491/15<33`, so

    C_f(L)>[exp(10^6 L)-1]/17-8/3.

The finite sum for e through degree four is 65/24>8/3. Hence rational
powers of 8/3 show C_f(10^-5)>1000 and C_f(10^-4)>10^28, without treating
a decimal Lambert evaluation as evidence. The code records the strict
margins and independently replays them with Fraction.

The weighted half-line C0 norm, with ||g||_sigma=sup exp(-sigma t)|g(t)|,
is exactly `Rf(sigma)` for sigma>p0. At sigma=2*10^6 the same logarithm
bounds give `Rf(sigma)<240/769`. On a finite interval this is an upper
bound as well. Returning to the unweighted norm loses exp(sigma L), so
this useful weighted estimate is not a macroscopic unweighted bound.

Finally the named scale dependence cancels from the Lambert argument:

    z=4 exp(2gamma_E-19/15)/(15delta),
    [p0/(af T0)]^2=192pi^2/[kappa hbar W0(z)].

This only restores the proper-time units of the frozen comparator. It
does not establish validity of semiclassical physics at that timescale.

## 5. The exact no-growing-response condition

For compactly supported forcing g, the growing part after its support is

    b0 exp(p0 t) integral exp(-p0 s)g(s)ds.

The cut decays and the conjugate pair is exponentially damped. Thus the
growing component vanishes **if and only if** the displayed weighted
moment is zero. A nonnegative nonzero smooth compact forcing fails it.
More generally, for bounded forcing on the whole half-line, the cut
kernel is L1: its integrated density behaves as r^3 near zero and as
1/[r log(r)^2] at infinity. The damped pair is L1 too. Subtracting the
tail of the exponentially weighted forcing integral leaves a bounded
term. Hence the same full half-line moment is necessary and sufficient
for a bounded causal response in this **linear** problem. It depends on
the future of the forcing; its vanishing is not freely available as a
condition in the original prepared-state initial-value problem.
Initial/history contributions must be included in the total numerator
before evaluating that condition; a density constraint alone is not
shown to impose it. In the full map g=Gf[X] is still unknown and nonlinear.

Deleting b0/(s-p0) leaves a nonzero Laplace inverse defect. Canceling it
using a future boundary condition or projecting the forcing changes the
specified causal problem. Order reduction is likewise a distinct model
approximation, not a rearrangement used by this checkpoint.
