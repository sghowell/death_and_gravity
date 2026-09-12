# Exact full two-time projector reduction

Use the literal tracefree spatial stress matrix diag(-D,-D,0,D) in the ten-field order (E,B,mA0,mAsp). Let PT_k = I-k k^T/(k.k), and let C_k v = k cross v. Projectors here are bilinear analytic matrices; a complex vector with zero bilinear square is outside their domain even if its Euclidean norm is nonzero.

For each source or detector time, the TT amplitude matrix is

Q_D = A D + B C_k^T D C_l,
A = m^2 f_kT f_lT - p_kT p_lT,
B = f_kT f_lT/a^2.

The complete two-time TT pair product is

tr[PT_k Q_G(s) PT_l Q_D#(t)^T].

Expanding it gives the four base/magnetic trace contractions used in the implementation. No magnetic momentum is replaced by the other one.

The mixed scalar amplitudes are

c_TL = m omega_l f_kT f_lL - (m/omega_l) p_kT p_lL,
c_LT = m omega_k f_kL f_lT - (m/omega_k) p_kL p_lT.

Their two-time products multiply, respectively,

n_l^T D# PT_k G n_l,
n_k^T D# PT_l G n_k.

The LL amplitude is

c_LL = omega_k omega_l f_kL f_lL
       - m^2 p_kL p_lL/(omega_k omega_l),

with geometric product (n_k^T D# n_l)(n_k^T G n_l). The n normalizations cancel into rational bilinear projector denominators in the complete expressions.

These four sectors account for all nine physical pairs. The absence of a direct mA0 coefficient in a tracefree spatial stress contraction is not permission to remove the constrained longitudinal electric or mass readouts.

Source and detector scale factors and amplitudes are independent in the exact algebra. The detector is analytically Schwarz-continued, not conjugated at the same complex argument. Each sector keeps its own inverse frequency sum before the source derivatives are taken. A numerical control confirms that replacing them with a common transverse frequency changes the actual curved result.
