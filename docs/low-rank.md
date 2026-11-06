# Low-rank methods

The bounded matrix-factor implementation is present; its task-level acceptance
is tracked in [tasks.csv](tasks.csv). The dense DG/Hermite/radiation paths remain
available independently. Mathematical equations and component-to-code mapping
are in [mathematics/low-rank.md](mathematics/low-rank.md).

| Responsibility | Implemented entry |
|---|---|
| Optional value representation, metrics, axes, sampling and numeric snapshots | `mdemat.data.tensor.MatrixFactors` |
| Physical weighted QR, small-core SVD and controlled truncation | `mdemat.common.lowrank.weightedQr`, `compress`, `add` |
| Separated action without a full Kronecker matrix | `mdemat.common.operator.applySeparated` |
| Fixed-rank K-S-L projector splitting | `mdemat.common.ode.evolveFixedRank` |
| Rank-adaptive Euler/compression with optional scalar moment correction | `mdemat.common.ode.evolveAdaptive` |
| Selected periodic discrete-velocity transport/scattering consumer | `mdemat.application.radiation.solveLowRankTransport` |
| Independent references, time/rank studies and total evolution costs | `tests/numerical/TestLowRankDynamics.m`, `benchmarks/lowRankTransportStudy.m` |

The representation is `F = U*S*V'`, with explicit spatial and velocity axes.
Its physical metric is `trace(F'*Mx*F*Mv)` for fixed Hermitian positive definite
product metrics. Native dense/sparse arrays need no wrapper. Factor payload is
reported separately from metrics, operators, workspaces and reference arrays;
small-factor storage does not establish whole-program memory or runtime savings.

`Rank` denotes representation width. Materialization is explicit. Compression
reports discarded physical norm and rank-cap saturation; that estimate is not a
PDE error bound. Unitarly changed factors are compared through represented
fields and physical moments. The fixed-rank integrator avoids core inversion;
its global splitting order is one. Adaptive Euler/compression is a different
first-order method and does not claim positivity, AP behavior or stiff stability.

Tests use independent dense actions/exponentials and a separately projected
`ode45` reference, plus weighted identities, small singular values, rank changes,
physical mass and completed-step restart. Benchmark records distinguish time
error from compression effects and process snapshots from measured memory peaks.

The old `Separable`, `SeparableFunction`, `SumSpace`, `SeparableGraph` and
`SeparableOdeIntegrator` keep their original meanings; they are not aliases for
this optional numeric factor representation. Existing macro-micro/Hermite paths
retain their own scaling, affine components, moments and discretizations.

CP/Tucker/TT/HT formats and nonseparable joint metrics require a concrete future
consumer and their own contracts. They are not empty framework packages or
unverified capabilities of the present matrix implementation.
