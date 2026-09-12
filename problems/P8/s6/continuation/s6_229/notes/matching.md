# Full scalar leading matrix, original finite contacts and proper normalization

Use the actual S216 complete two-form contractions, with their passive-dimensional terms and both mixed constraint contractions. At zero transfer their leading combination is

(TT00+TT01+TT10+TT11)/4+LL/4+(LC+CL)/8+CC/16.

For arbitrary symmetric directions D,G the angular average is

C(d)tr(DG)+B(d)tr(D)tr(G),

C(d)=(2d^2+d-8)/[2d(d+2)],
B(d)=(d^3-4d^2-3d+16)/[4d(d+2)].

The literal trace pair D=G=I_d is T(d)=(d-1)(d^2-5d+8)/4. Both separately ordered trace/tracefree cross averages vanish. These facts follow from the complete pair, not deletion of the longitudinal mode.

Reconstruct C(d),B(d) FIRST, then contract with fixed PHYSICAL scalar directions Q=2wI_3-2cPi. The dimensionally continued physical matrix is

M(d)=C(d)[[12,-4],[-4,4]]+B(d)[[36,-12],[-12,4]].

At d3,

M=[[4,-4/3],[-4/3,8/5]]
  =L0^T diag(4,52/45)L0,
L0=[[1,-1/3],[0,-1]], det M=208/45.

The fixed-physical first dimension jet is

M'(3)=[[46/15,-46/45],[-46/45,22/25]].

In particular C'(3)=91/450 and B'(3)=4/225. A dimension-varying test family I_d gives a DIFFERENT matrix jet, [[4,-8/9],[-8/9,22/25]]. The code checks this difference and does not substitute the latter into a fixed-physical finite part.

With normalized output 64pi^2/a(t)^3, the corrected positive Kubo phase gives leading radial matrix

32 M k^4 sin(2k Delta_sigma)/[a(t)^4 a(s)].

Its Abel integral is 24M/[a(t)^4 a(s)Delta_sigma^5]. Four derivatives of M/tau are24M/tau^5. This matches the original flat scalar block

L0^T diag(Ftrace,(8/3)F2)L0 D^4.

Ftrace and F2 are the COMPLETE original factors, with original finite constants -4 and -1/30, original massive cuts and the retained inverse shear pole. In particular the fixed fourth local matrix is L0^T diag(-4,-4/45)L0. No logarithm-only or massless reference is substituted.

As in S228, at frozen scale a0 the original continued measure and leading amplitude are a0^(-d-2)k^(d+1)dk. Proper momentum k=a0 p removes this factor for EVERY nearby dimension. The scalar invariant coefficients agree throughout that neighborhood, so both their physical value and their first jet agree before fixed-physical contraction. This fixes the highest finite contact in the original prescription. Equality only at d3 would not suffice.

The resulting Gaussian scalar Euler block after normalization by kappa a^3 has coefficient

gamma=1/(64pi^2 kappa).

Neither this output normalization nor its inverse is moved through a retarded time kernel. All finite-P gradient terms, lower curved contacts and the actual state remain in the remainder described next.
