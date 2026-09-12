# Exact original lost shell, including grazing angles

Fix P != 0 and K > |P|. Let u = n.phat. The original second-leg condition is

r^2 - 2r|P|u + |P|^2 <= K^2.

Its positive boundary is

r_* = |P|u + sqrt(K^2-|P|^2(1-u^2)).

The pair mask removes r in (r_*,K] exactly for u < |P|/(2K). This includes a small positive-u grazing strip, which must not be discarded at finite K. The branch is the positive root; at the grazing point it equals K because 2K^2-|P|^2 > 0 under the stated K > |P| assumption. The exact algebra checks the squared identity and this sign separately.

The triangle inequality bounds the lost-shell thickness by |P|. For the leading radial integrand r^3, its exact lost moment is

(K^4-r_*^4)/4.

Writing epsilon = |P|/K, its normalization by K^3|P| tends to -u at fixed u < 0. The normalized lost moment is uniformly bounded on its support; the positive-u strip has width epsilon/2 and contributes zero in the limit. Dominated convergence therefore yields the weight (-u)_+ on the entire sphere, including the grazing region.

The original pair-band minus the one-k-ball comparator is minus this lost-shell integral. This geometric sign is distinct from the positive j = 0 current sign.

Independent tests integrate the exact finite-epsilon lost shell against the actual leading angular polynomial in all five tensor channels, down to epsilon = 10^-12. Boundary identities, the finite grazing threshold and triangle-inequality bounds are also checked at multiple angles and epsilon values. These tests do not replace the written uniform shell argument.
