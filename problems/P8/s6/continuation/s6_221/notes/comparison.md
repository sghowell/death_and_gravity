# Two-switch propagator, low band and Sobolev extension

For |P|>=100 cover the unit interval by the left outer chart, central chart and right outer chart, switching at t=-3/16 and3/16. Both charts are valid at each switch, and every ordered subinterval crosses at most two switches. The original phase variable Z is continuous there because its smooth first-order ODE is globally regular on I.

Each chart's root energy has Gronwall factor exp(1e28 times its segment length). A switch costs at most(A lambda)^2, by converting its outgoing energy to Z and then to the incoming energy. Initial and final conversions add two further factors. Since A=1e12 and there are at most two switches,
||U(t,s;P)||<=1e72 lambda^6 exp(1e28|t-s|).
Segments of zero length and endpoints at the switches cause no problem; the same bound follows by continuity. No matching condition is inserted into the physical equations.

For |P|<=100, the preceding source-pinned S220 first-order generator bound is1000(1+|P|²)^2. Its direct matrix Gronwall estimate is at most exp[1000(10001)^2]<exp(1e12) on I. This includes the finite-q complementary-chart degeneracy. The low/high threshold is mathematical and does not alter the parent theory or supply a physical cutoff.

The total time length is at most1. Since log10<3 (already e³>1+3+9/2+27/6>10), 72log10+1e28<1e29. Both momentum regions therefore satisfy the single bound
||U(t,s;P)||<=exp(1e29)(1+|P|²)^6.

For compact momentum support, ordinary finite-dimensional ODE existence and variation of constants give the unique phase solution. Apply the uniform bound pointwise and use Plancherel and Minkowski:
||Z||Linf H^r <= exp(1e29)[||Z(s)||H^(r+12)+||f||L1 H^(r+12)].
The twelve, not six, spatial derivatives follow because the amplitude weight is(1+|P|²)^6. The H^r norm uses the square root of the weighted Fourier integral.

Approximate arbitrary indicated data and forcing by compact-momentum smooth inputs. The bound gives convergence in C_t H^r; pointwise ODE continuity and domination justify time continuity for both the free and Duhamel terms. For uniqueness in this comparison class, Fourier localization gives a finite-momentum linear ODE with zero initial data and forcing. Each localized difference vanishes, hence so does the full difference. This is the unique solution of the declared coefficient-sector phase system in the resulting distributional class.

At P=0 the coefficient system has a finite measure-zero Fourier extension. In R3 Sobolev spaces that assignment changes no norm or solution. It must not be identified with the separate literal global homogeneous constraint problem. Similarly f is mathematical additive phase forcing; no physical quantum-source-to-f map is established by this argument.

The propagator is thus uniform in momentum only as a map with twelve lost spatial derivatives. It is not uniformly bounded on H^r to itself, does not close the higher-time/spatial-derivative Gaussian feedback estimate, and its enormous constant does not demonstrate smallness relative to kappa. A different or sharper compatible-space argument is needed for the full quantum inverse.
