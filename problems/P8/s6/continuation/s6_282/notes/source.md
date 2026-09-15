# Whole source and normalization

S281 pins the entire current R,F, all vacuum constants, formal loop marker
and heavy source. S280 derives the Einstein contact, both scalar exchanges
and cubic graviton exchange directly from the original action; all four
physical helicities and the Ward identities fix the normalization.

This successor imports that exact source and sewing function unchanged.
For x=a.n,y=b.n,z=a.b, h=s/(s-4mu), let

    G=1-z^2-x^2-y^2+2zxy,
    N=(h-1)^4+(z-xy)^4-6(z-xy)^2 G+G^2.

Its exact integrand check is

    sewn_helicities/(32pi)
      =(s-4mu)^2 N/[256pi kappa^2(h-x^2)(h-y^2)].

The normalized spherical average supplies the whole normal-channel cut.
Neither a high-energy approximation nor a massless external limit is used.
The inherited check dictionary is copied, and warm-parent-order tests check
that building the new source cannot mutate the parent's contract. S279
remains rejected and its frozen bytes are not repaired.
