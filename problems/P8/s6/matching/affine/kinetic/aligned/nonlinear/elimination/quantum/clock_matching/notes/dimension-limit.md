# Uniform regulator limit for the prepared subtracted mode family

This proof concerns the formal dimensional continuation used by S6.52
and the [stated counterterm prescription](prescription.md). It does
not posit a Hilbert space in a noninteger dimension. The real physical
state remains exactly the fourth-order Gaussian preparation of S6.50.

## Analytic continuation at fixed momentum

Take real u in I=[-1/2,1/2], real comoving k>=0, m>=1000 and the
complex disk |D-3|<=1/4. The physical a=(1+u^2)^2 is between one
and Amax=25/16. Set omega^2=m^2+k^2/a^2, z=(k^2/a^2)/omega^2,
and continue the S6.52 canonical mode equations and physical readout
coefficients in D. In particular, use D-1 transverse polarizations,
not two until the limit. Omega,z are independent of D.

For each polarization use W=omega*S with
S=1+P2/omega^2+P4/omega^4. The P2,P4 are the actual S6.52
dimensional reference coefficients with every Hubble jet replaced
by derivatives of H=4u/(1+u^2). Differentiation includes
z'=-2Hz(1-z); derivatives of fourth-order coefficients are not
truncated at the Hubble jets needed to define the coefficients.

The exact complex polynomial-box bound is coefficientwise: after
D=3+eta, write P(u,z,eta)/(1+u^2)^n with rational coefficients;
sum abs(c_ijk)*2^-i*4^-k bounds it throughout the disk and I.
The denominator is at least one because u is real. Every rational
reconstruction is checked. This bounds complex eta, not just its
real interval. The report retains all six coefficient envelopes
P2,P4,B2,B4,A,B, using the notation of S6.50/51.

Exact lower-order cancellations and the independent abstract WKB
residual identity imply, uniformly on this domain,

    |S-1|<1/2, |W'/W|<4, |rho_res|<=C/omega^4,
    C=10,000,000.

The two actual residual envelopes are below 1,646,583 and 1,747,521,
respectively, so C is a stated conservative bound, not an unspecified
constant. Also |Im integral W du| is at most P/m+Q/m^3<1/2, where
P,Q are the reported absolute P2,P4 envelopes; the interval has
length one. Since Re W>omega/2, its positive-real square-root
branch is analytic throughout the disk.

Use the analytic pair

    f_+=(2W)^(-1/2) exp(-i integral W du),
    f_-=(2W)^(-1/2) exp(+i integral W du),

with phase zero at u=-1/2. Do not conjugate the complex dimension.
For real D these are conjugates. Their bilinear Wronskian is i.
The phase estimate and exp(x)<=1/(1-x) for 0<=x<1 give
|f_+|,|f_-|<=2/sqrt(omega). The canonical physical rates obey
|d_T|,|d_L|<=13/4, whence

    |f_+'-d f_+|, |f_-'-d f_-|<4 sqrt(omega).

Indeed |W|<3omega/2 and |d+W'/(2W)|<=21/4;
3/2+21/(4m)<2. These bounds hold for either sign of the phase.

Prepare exact solutions v_+,v_- of the continued mode ODE to match
the respective f,f' at the left endpoint. The coefficients, initial
data and their chosen square-root branch are analytic in D. Uniform
Picard iteration on the finite time interval gives analytic exact
solutions at each fixed momentum. At real D the pair is conjugate;
at D=3 it is precisely the S6.50 prepared physical mode and conjugate.

## Uniform evolution and subtraction envelopes

In the reference basis the variation matrix is

    [ -i*rho_res*f_+*f_-    -i*rho_res*f_-^2 ]
    [  i*rho_res*f_+^2       i*rho_res*f_+*f_-].

The Wronskian, variation constraint and forced equation are checked
independently as exact bilinear identities. Because
|f_+ f_-|=1/(2|W|)<=1/omega and |f_+^2|,|f_-^2|<=4/omega,
its induced one-norm is at most 5C/omega^5. Put
nu^2=m^2+k^2/Amax^2. Each evolution column has norm <=exp(J),
and differs from its initial unit column by at most exp(J)-1, with
J=5C/nu^5<1/4. This is a finite-interval Gronwall estimate.

For both energy and pressure write the analytic bilinear readout as

    [A(v_+'-d v_+)(v_-'-d v_-)+B omega^2 v_+ v_-]/(2a^D).

The actual clock alpha,beta remain real and unchanged. Uniformly
|A|<=5/4, |B|<=3/2 and |a^-D|<=1. For example the largest
pressure kinetic ratio is at most (3+1/4)/(3-1/4)=13/11<5/4;
the pressure potential ratios are at most 5/11. Actual lapse
weights satisfy the stronger physical S6.50 bounds.

The reference product estimates above bound this bilinear form
by 13omega times the product of coefficient one-norms. Thus the
exact/reference difference is at most
13omega[exp(2J)-1]<=52omega*J per polarization. Summing with
|D-1|+1<=13/4 and using omega<=Amax*nu gives the uniform bound

    |combined exact/reference difference|<=1000*C*Amax/nu^4.

This bound is deliberately looser than the sharper physical D=3
bound used for the final numerical estimate.

The S6.51 generic reference-tail identity is algebraic and remains
valid for complex coefficients. Use |S^-1|<=2, the complex box
envelopes, the preceding |A|,|B| bounds and |c1|<=9/4. The computed
per-polarization tail constants are below 481,483 and 481,856.
Consequently Ctail=1,000,000 bounds both, and

    |combined reference minus orders 0,2,4|<=Ctail/omega^5.

The factor (|D-1|+1)/4 is at most 13/16<1. The code retains the
full rational tail bound, not just a large-momentum series.

## Passing the regulator limit through the finite integral

The radial measure coefficient is
2/[2^D*pi^(D/2)*Gamma(D/2)]. It is analytic and bounded on the
compact disk: the reciprocal Gamma function is entire
([NIST DLMF, section 5.2](https://dlmf.nist.gov/5.2)). The usual
fixed-scale dimensional normalization is also analytic there.
For k>0 define k^(D-1) by the real logarithm. Its modulus is
k^(Re D-1). The combined preceding envelope is therefore bounded
at k>=1 by constants times k^(-7/4)+k^(-11/4), and at 0<k<=1
by a constant times k^(7/4). Both ends are integrable, uniformly
over the disk and the full physical time interval. No numerical
bound on the Gamma prefactor is needed for this compactness argument;
the final physical error estimate uses the sharper explicit D=3
integrals, not an unknown prefactor.

The subtracted exact integrand is analytic in D at each momentum.
The uniform integrable envelope permits dominated convergence at
D=3; it also gives holomorphy in the disk interior by integrating
around triangles and applying Morera's theorem. Denote this finite
integral F(D). The regulated prescription is exactly F(D) plus
the S6.52 meromorphic local Gamma terms, with the stated local
counterterm density subtracted. F(D) has no pole. Hence its finite
part is F(3), and the full finite first variation is

    F(3) + radial local finite coefficient
         - epsilon coefficient of the continued counterterm variation.

This establishes the limiting step needed to add S6.51's physical
finite integral to the present matched local coefficients. It does
not supply derivatives of F(3), all-order Hadamard admissibility,
noncompact time control, other loops, higher metric variations,
unknown UV matching data, a corrected bounce or original P8 closure.
