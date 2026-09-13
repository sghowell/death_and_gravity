# Physical ordered endpoints and complete radial grades

For real spatial directions the correctly ordered retarded commutator is minus the imaginary part of the detector annihilation pair times the source creation pair. Equivalently it is the positive imaginary part of detector SHARP times source annihilation, with detector phase exp(+i theta_t) and source phase exp(-i theta_s). The algebraic sharp changes i to-i in mode coefficients; it does NOT conjugate the external source direction. Complex Fourier sources are inserted bilinearly after the real tensor components are formed.

For a single complete source pair U and g=1/(W_k+W_l), set L U=dt(gU). The exact identity is

integral U exp(-i theta)
=exp(-i theta) sum_(j=0)^4 i(-i)^j g L^j U
-i integral L^5 U exp(-i theta).

Differentiating the full right endpoint sum gives U exp(-i theta)+i L^5 U exp(-i theta); the remainder cancels the latter. Every lower endpoint is zero only by the source's common prepared germ. There is no final source cutoff. This identity fixes all odd signs and leaves the entire bulk, whose all-transfer estimate is not supplied here.

To extract all radial coefficients, use the order-four polynomial ring in x=1/r. Put phase=(Wbar_k+Wbar_l)^-1, gbar=a phase. For each of00,01,10,11, the left coefficient is the algebraic sharp of A or B; the initial source is the corresponding A/a or B/a. Repeatedly apply dt(gbar row), and add gbar row to the next independent source-derivative slot. At endpoint j only degrees0..4-j are required. Its coefficient is

Im[i(-i)^j phase sharp(amplitude_D) source_row_(j,r)].

The unintegrated prefactor is r^(1-j). Before the radial measure, the full metric geometry must also be multiplied and its y=Px degree convolved with these coefficients. All60 ordered endpoint/source slots and140 scalar coefficients are retained. After the complete angular average there are35 rows indexed(j,r,degree), with j0..4,r0..j,degree0..4-j.

In d spatial dimensions the radial power is r^(d-j-degree). Thus at d=3 the logarithmic coefficient has j+degree=4. The complete spatial difference has one quadratic power row, (j,r,degree)=(0,0,2). Every other row of total grade below4 vanishes identically after the full angular calculation; all19 such identities are checked. Logarithmic source orders3 and4 vanish as full dimension-dependent expressions, so no evanescent derivative of an omitted high source jet is lost.

The trace/gradient01 and10 directions are distinct. Odd endpoints are not discarded, even if some leading slots happen to vanish. Independent full numerical time-Taylor evaluation of all six unexpanded Riccati iterates checks all60 endpoint slots at two large radii for each of three distinct fixtures, including noninteger dimension. These are cross-checks of the complete coefficients, not bounds on the full causal bulk.
