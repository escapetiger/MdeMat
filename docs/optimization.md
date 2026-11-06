# Numerical and performance review

The final numerical checks are indexed in [status.csv](status.csv). Correctness
repairs include Hermitian test functions, physical boundary measures, valid
Gauss interpolation nodes, conserved-variable time integration, and complete
accepted-step restart histories. They take precedence over tuning arithmetic.

The [low-rank study](evidence/lowrank-study.json) measures the same complete
solver-call scope for dense and factor evolution, including operator creation.
It records a warmup, repeated timings, MATLAB/processor information, thread
policy, ranks, weighted errors, mass defects and factor payload. Its independent
reference is the dense semidiscrete exponential; this does not measure continuum
spatial or velocity error. Timings exclude reference construction, assessment
and I/O. The native memory value is a process snapshot containing the reference
and earlier session allocations, rather than the solver's allocation or peak.

For this small problem the dense solver is faster. Compression reduces factor
payload at loose tolerances, but QR/SVD, metric validation and temporary factors
cost more than the saved dense arithmetic. Tight tolerances approach the dense
Euler time-error floor while retaining more modes. These observations support
an explicit optional factor API; they do not justify making low rank the default.

Useful next optimizations, each requiring a larger controlled benchmark, are:

* Apply diagonal velocity operators and the rank-one scattering average
  directly to factors, avoiding dense velocity matrices. Preserve the current
  weighted operator-action and mass tests as references.
* Reuse factorizations of immutable physical metrics within one evolution,
  with explicit ownership and invalidation, after profiling shows repeated
  metric work dominates. Avoid a global cache keyed only by array dimensions.
* Profile assembly, linear solves and observer/I/O phases separately on
  representative multidimensional radiation/plasma studies. Keep rendering
  density independent of numerical resolution and retain numeric-only entry
  points for batch studies.

## Tools refinement

The current [tools study](evidence/tools-study.json) measures collection
construction, original-component selection, explicit copy/update, FEM plans and
sampling, physical error reduction, primitive IO, optional diagnostic hashing,
timing intervals, environment metadata and edge rendering. Every case checks
its answer before and after timing. The inputs use a private random stream;
temporary files and the owned figure are removed. The report records the native
process, warmup, repetitions, timing scope and input payloads. Memory values are
before/after process snapshots, not allocation costs or peaks.

Batch edge rendering reduced both artists and elapsed drawing cost for the
measured graph. Each endpoint was checked against the same independent edge
list. MAT was faster for the measured dense array, while the small-leaf case
favoured JSON; there is no universal format winner. Both formats retained exact
types, shapes and values. The solve timings with and without a diagnostic hash
include construction, factorization and solve, so their small-case ordering is
not evidence that hashing speeds up a solve. Hashing remains explicitly optional.

No old-class performance baseline was fabricated after removing Dataset and
Database. Struct values simplify ownership, but that design change alone does
not establish a speedup. The implementation also avoids resampling accepted
numeric frames when only a time-dependent exact reference needs refreshing;
paired-frame tests preserve the numerical/time correspondence.

Further changes should follow a representative application profile. In
particular, isolate hash cost, measure repeated large-field selection, and assess
high-dimensional/low-rank display memory before adding caches or materializing
full tensors. The current measurements do not justify a global cache, a new IO
schema or another graphics abstraction.

No current result establishes a general speedup, positivity after truncation,
AP behavior of the new low-rank consumer, or a memory peak. Existing tests verify
bounded methods and the stated physical regimes; extending those claims needs
new mathematical and numerical evidence.
