# Global spacelike sign and a decay-preserving envelope

Set s=-t, t>=0, d=t+1, K=4m^2-1, b=A/(m^2-A), and h=d b.
The same complete finite on-shell subtraction is exactly

    g_x,R(-t) = -(K+d) log(1+h) + K h.

Its second derivative with respect to d is the second derivative of the
original kernel at s=1-d, and its value and first derivative vanish at
d=0. These independent identities fix the expression without an unjustified
complex logarithm branch manipulation. Here all logarithms are real with
positive arguments.

Concavity of the original g_x on s<=1 makes g_x,R(-t) nonpositive, strictly
negative for interior x and positive coupling. Also 0<=b<=1/K. Since
h-log(1+h)>=0, which follows from its derivative h/(1+h) and zero anchor,

    -g_x,R(-t) <= d log(1+d b) <= d log(1+d/K).

After parameter integration,

    0 <= -f_R(-t) <= C d log(1+d/K),   C=2 N Y/Q,
    |f_R(-t)|/d^2 <= C log(1+d/K)/d <= C/K.

This covers all t, including momenta above the fermion mass. It retains
the logarithmic-over-t decay of the inserted free propagator. The last,
coarser constant loses that useful behavior. Neither bound alone is a
bound on a complete outer diagram; local forest subtractions and other
propagators must still be treated.

The ratio |f_R(-t)|/d can grow logarithmically in this fixed-order
functional. No claim of inverse positivity at arbitrary energy follows.
A hypothetical far-ultraviolet zero of a fixed-order approximation would
not by itself exclude the full candidate.
