# Prepared leading vector response

[S6.43](FORMULATION.md) keeps the S6.42 physical action unchanged and
derives its actual quadratic nonlinear source from the original ODE.
The constrained canonical force is J=(g²*S2)'/g, not the raw source.

On the stated compact window, output momenta and prepared zero-data
class, exact continuous bounds give ||v-zeta*J||<=zeta*A/80. The
physical spatial-vector error is <=6*zeta*E and the temporal error
against S2 is <=500*zeta*E. The source envelope is supplied explicitly
by C³ lapse/trace jets, with Fourier convolutions retained.

See the [retarded proof](notes/retarded.md) and [light-jet proof](notes/light-jets.md).
This is a leading perturbative conditional estimate, not a full
nonlinear solution, quantum effective action, loop/cutoff bound,
stationary gap, V/G/B verdict or original P8 closure.

Read-only replay with all local source roots on PYTHONPATH:

    python -m p8_aligned_response.verify --check

The report pins all local sources and fully rebuilds the frozen S6.42
ancestry. Tests include original temporal constraints, independent
canonical jets, a manufactured exact forced solution, preparation,
Fourier convolution, physical units, strict domains and mutation rejection.
