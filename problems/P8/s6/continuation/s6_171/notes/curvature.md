# Explicit curvature continuation and finite Euler stress

This is a separately named curved extension GY14-SAT8-MR/EC-N0.
The inherited flat calculations did not determine an evanescent
gravitational prescription. Select fixed trace4 and pole counterterms

    [a E_D-c C_D^2]/(2Q epsilon),
    a=11/360, c=1/20, D=4-2epsilon.

No additional finite R^2 is added in THIS E_D,C_D basis. Total box R
terms are boundaries under compact metric variations. This convention
is a choice, not a claim that curved stress has no finite ambiguity.

The local heat-kernel input is independently checked against
[del Rio et al., v2](https://arxiv.org/pdf/1703.00908), equations106,
109 and110. Its rank-four spin trace yields, without total derivatives,

    (5R^2-8Ricci^2-7Riemann^2)/360.

The sign/continuation is calibrated to the inherited +M^4/(Qepsilon)
counterterm and R_E=-R_P8. The finite tensor below is derived directly
from the metric variation; no unconverted anomaly-sign convention is
imported from the paper.

In general dimension the literal polynomial differs from aE_D-cC_D^2:

    literal-EC=(D-4)[Ricci^2/(10(D-2))
                    -(D+1)R^2/(60(D-1)(D-2))].

Consequently the finite counterterm difference literal minus EC is

    [-Ricci^2/20+R^2/72]/Q
      =[-C4^2/40+E4/40-R^2/360]/Q.

On FLRW the finite stress difference is the variation of -R^2/(360Q),
not zero. For example, at the actual bounce H=0,Hdot=4, variation of
R^2 gives rho=576, so this difference has rho=-8/(5Q) per copy.
No earlier frozen flat result is overwritten by selecting the new basis.

For spatially flat FLRW in D=d+1 dimensions,

    R=-2d Hdot-d(d+1)H^2,
    Ricci^2=d^2(Hdot+H^2)^2+d(Hdot+dH^2)^2,
    Riemann^2=4d(Hdot+H^2)^2+2d(d-1)H^4,
    C_D^2=0,
    E_D=d(d-1)(d-2)[4H^2 Hdot+(d+1)H^4].

Although C_D^2 itself vanishes, its metric variation also vanishes
on a conformally flat metric: it is quadratic in the Weyl tensor.
The Euler variation must be kept before the regulator limit.
Introduce lapse N and v=partial_t log a, H=v/N. Integrating its
action once by parts gives the exact reduced density

    -d(d-1)(d-2)(d-3) a^d v^4/(3N^3).

Both the unreduced lapse/scale variation and this reduced action yield

    rho_E=-d(d-1)(d-2)(d-3)H^4,
    P_E=d(d-1)(d-2)(d-3)[H^4+4H^2 Hdot/d].

Multiplying by a/(2Qepsilon), then setting d=3-2epsilon and taking
the finite limit, gives per Dirac copy

    rho_EC=11H^4/(60Q),
    P_EC=-11[H^4+4H^2 Hdot/3]/(60Q).

This polynomial is conserved through H=0 without division by H.
Its trace is aE4/Q. Setting d=3 before multiplying the pole would
incorrectly erase it. For all42 copies, Q>144, |H|<=2, |Hdot|<=4
give rho_EC<1 and |P_EC|<2. All remaining finite geometric vacuum
polarization is retained in the complete state/projector remainder;
the Euler term is not asserted to be the entire curved tensor.
