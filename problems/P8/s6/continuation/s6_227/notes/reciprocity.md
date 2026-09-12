# Retarded/advanced reciprocity and the SOURCE-time derivative

For B=L D^2+Q+C(t)D, the time-reversed operator in r=-t is

Brev=L D_r^2+Q-C(-r)D_r.

It has the same leading and spatial matrices and the same bound8 for its first-order coefficient. The original primal velocity Volterra proof therefore applies on every reversed subslab, with exponent20. Its retarded inverse constructs the ADVANCED inverse Yadv ofB in the original time direction.

Formal-adjoint Green reciprocity gives

Zret(t,s)=Yadv(s,t)^T,

whereZret is the retarded inverse ofB*. To prove the identity, pair compact sources with the advanced and retarded solutions and use the complete bilinear integration-by-parts boundary from S225. The retarded response vanishes before the source and the advanced response after the detector, so the outer boundaries vanish without deleting an initial atom. Both distributional inverse identities and uniqueness identify the kernels. This argument uses the actual formal adjoint, including-Cprime^T.

The reversed primal kernel and its first output derivative then give

||Zret(t,s)||<=(5/2)(t-s)exp(20(t-s)),
||partial_s Zret(t,s)||<=(5/2)exp(20(t-s)).

The second line is a SOURCE-time derivative. It is not a replacement of the S225 forward-time derivative bound with exponent108. At identical ordered times, Zret is also not simplyYret^T.

An independent exactly solvable constant-h, q0 matrix kernel verifies both ordinary differential equations, both initial derivative normalizations and reciprocity. A second nonconstant-h recurrence verifies the exchanged-source kernel identity through sixth order at nonzeroq. The finite Taylor jet is a diagnostic, not the proof of global reciprocity.

## The full initial-delta cancellation

Extend an ordinary source f by zero to the past. Distributionally D(theta f)=theta f'+delta f(0) when f is smooth on the right. Since D*=-D,

ZD*f=-integral_0^t Z(t,s)f'(s)ds-Z(t,0)f(0).

Integration by parts has an interior lower-boundary term+Z(t,0)f(0), which cancels the retained initial delta. The upper term vanishes becauseZ(t,t)=0. Thus

ZD*f=integral_0^t partial_s Z(t,s)f(s)ds.

The same identity extends toL2 sources by the bounded kernel formula and distributional continuity; a right-hand trace is not required of the limiting L2 function. The independent test uses nonzero f(0), and omission of its initial delta produces a nonzero defect.

With g=sigma-20>0, the time-weighted scalar majorants give

||Y||,||Z||<=(5/2)/g^2,
||DY||,||ZD*||<=(5/2)/g

on weightedL2_tH^r, all spatial momenta. The ordinaryY derivative also has no missing upper term becauseY(t,t)=0. These are the four operator bounds needed for the complete local conjugation.
