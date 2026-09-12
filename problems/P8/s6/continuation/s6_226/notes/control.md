# A bounded correction can retain a growing reference pole

The inverse theorem is deliberately an existence theorem in an exponentially weighted space, not a stability theorem. An explicit SEPARATE comparison proves the distinction.

At p=4m^2, define the constant channel correction

Vcmp=diag(0,(8/3)A2(4m^2)).

A2 is the unchanged original radial shear factor. Its positive-axis monotonicity gives A2(4m^2)>0. Its radial denominator at this point is1+(1-y^2)>=1 after scaling, so

A2(4m^2)
<=1/30+integral_0^1 y^2(30-20y^2+3y^4)/30 dy
=1/30+3/14=26/105.

Thus ||Vcmp||<208/315<1. The weighted theorem applies withM=1.

Nevertheless the corrected shear factor atp=4m^2 is exactly

(8/3)F2(4m^2)+(8/3)A2(4m^2)=0.

The actual radial derivative A2'(p)>0 for p>0 shows that this zero is simple. Atq=0 it is a simple positive-real-Laplace pole atlambda=2m, so the retarded inverse has a growing reference contribution. The chosen weight sigma_1=2m exp(40) lies strictly to its right and retains it. There is no contradiction between that pole and a bounded weighted inverse.

This comparison is not imposed on the actual physical action or its finite renormalization prescription. It does not diagnose an actual physical instability. It is an independent counterexample to the false inference "bounded weighted inverse implies no growing pole."

Other controls retain the original shear pole, detect the failure of the positive exterior argument at zero time weight, and use a noncommuting finite two-time/two-channel model. In that model, the direct full inverse equals both correctly ordered Neumann expressions, but not a reversed expression that assumesKV=VK. An exact eight-step partial sum has the predicted geometric tail. These finite matrices check order and error accounting, not convergence of a continuum discretization.
