# Actual nearby classical cones

S6.89 extends the S6.83 central cone calculation along an actual
S6.88 nearby classical bounce on |u|<=10^-7. The clock cone stays
strictly outside the physical matter cone; both scalar principal
kinetic and gradient matrices stay positive.

The proof consists of:

- [Full geometric reduction and Euler-first principal](notes/reduction.md).
- [Primitive-free rational background evolution](notes/evolution.md).
- [Whole-box bounds and the actual-solution bridge](notes/enclosures.md).
- [Classical characteristic comparison](notes/separation.md).
- [Scope, limitations and source comparison](notes/scope.md).

The report pins every local source and proof, fully rebuilds the
unchanged S6.88 parent, certifies exact residuals and continuous
bounds, and records all exact rational-system polynomial fingerprints.
Native science uses unmodified SymPy; rational fields are ordinary
public arithmetic operations, not a library patch. The written
geometric and characteristic arguments are not proof-assistant
formalization.

With all P8 source roots on PYTHONPATH, run:

    python -m p8_proca_nearby_cones.verify --check

Omit --check to print a fresh report. Neither command overwrites a
certificate. The project-level audited snapshot runner replays the
complete P8 test set; scoped tests use the same static namespace
collection setup without the optional full-suite arithmetic adapter.
