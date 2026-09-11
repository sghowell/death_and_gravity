# Regulated field identity, including evanescent coordinate terms

Write the complete first slope as
k_D=r_D-fp_D=k0+epsilon k1+..., and let
K_D=1+h k_D+h^2 t_D with t_D finite near epsilon zero.
The MS field factor has only poles:

    Z_MS=1+h z11/epsilon+h^2(pure poles),
    z11=-2NY/Q.

The physical-OS field satisfies Z_H(q_H,D)=Z_MS(q)/K_D(q).
Here q denotes fixed MS inputs. The FIRST bare-coordinate
identity requires

    q_H,D=q+h q1_D+...,
    q1_D=-k_D E_Phi q,

with E_Phi=2L partial_L+G partial_G+2g partial_g+
Y partial_Y (use fundamental G or g in each expression).
M, mF and a have zero first shift. The physical light mass
is fixed independently by its prescribed mass reference.

It is incorrect to replace q1_D by its epsilon-zero value
inside a pole-containing expression. Since
E_Phi r_D=2r_D, E_Phi fp_D=fp_D and E_Phi z11=z11,

    q1_D partial z11/epsilon=-k_D z11/epsilon.

Its finite term is -k1 z11. This is exactly the same finite
term as the product -k_D z11/epsilon on the right-hand side
Z_MS/K_D. They cancel only after BOTH legs are retained.

At fixed hybrid inputs the finite two-loop OS field
counterterm is -t_H. Including the first parameter
re-expansion gives

    Fin[Z_H]_(h^2)
      =-t_H-k1 z11+k0(2r0-fp0),
    Fin[Z_MS/K_D]_(h^2)
      =-t_MS+k0^2-k1 z11.

Therefore the complete finite MS second normalization is

    t_MS=t_H-k0 r0.

The single and double kinetic pole differences also cancel
under the same first coordinate shift. This identity does
not rely on fitting t_MS or erasing the pre-OS hybrid slope.

An epsilon-independent coordinate extension is possible
only if the compensating evanescent vertices are retained.
The displayed full regulated coordinate map is the convention
used here; omitting both would leave a spurious finite term.
A private preliminary finite-only calculation was corrected
before final verification and freeze; no earlier frozen
scientific source or report is changed.
