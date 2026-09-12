# Complete mixed coefficients and the finite continuum polynomial

The endpoint index is j = 0,...,4, the spatial Taylor degree is n = 0,...,4, and the source time-jet index is r = 0,...,j. The detector is held fixed while source-time derivatives are taken. For each complete source-jet row, use the S6.202 joint time/spatial domain and the actual unit-W8 modes, including the longitudinal constraint and all nine physical pairs.

The analytic estimate is made before taking the final current imaginary part. The detector factor is continued with the analytic Schwarz operation; ordinary conjugation at the same complex argument is not holomorphic. Cauchy's coefficient formula gives delta^-n nu^-n, without the factor two needed for a geometric-series remainder. No test Fourier transform, Borel-state coefficient or moving regulator mask is differentiated.

The row coefficients c_jr are

(1),
(2,1),
(8,5,1),
(48,33,9,1),
(384,279,87,14,1).

Their sums are 1, 3, 14, 91 and 765. The factor nine already sums the physical polarization pairs; no transverse-only restriction is made.

For fixed real external direction, the homogeneous coefficient of spatial degree n satisfies the complete bound in FORMULATION.md. For each of the ten finite cells, the massive radial moment is integrable on the fixed high band |k| >= m. Dominated convergence removes the original two-leg regulator. Only after this limit can the coefficient tensors be integrated independently of P, giving a genuine finite homogeneous spatial polynomial of degree n. The finite-cutoff integral is not asserted polynomial because its intersection domain depends on P.

Independent tests calculate actual mixed source-time/spatial coefficients at internal scales 1000, 10^6 and 10^12. They use noncommuting source/detector matrices, the joint polarization frame, every source time jet and all nine physical pairs. Source-time Cauchy jets retain the complete recurrence; spatial coefficient extraction uses both twelve and sixteen nodes at 135-digit precision. Agreement is checked relative to the actual coefficient scale, with a tiny absolute floor, before the individual-jet and complete-row bounds are checked.

These numerical checks test the actual coefficient implementation at representative points. The continuous bounds rely on the written S6.202 analytic estimates, not finite sampling.
