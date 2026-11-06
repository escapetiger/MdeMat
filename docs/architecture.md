# Architecture

This records the implemented MATLAB architecture, including the tools redesign
specified in [tools-plan.md](tools-plan.md). Graph/mesh plotting now belongs to
its data types; configuration, samples, styles and analysis results use values;
Timer and Visualizer own their respective resources.

**Current acceptance: PASS.** Native numerical and symbolic tests, production
coverage, source analysis, packaging, fresh installation and project checks passed.
Current [test evidence](evidence/tools-required-tests.json),
[package identity](evidence/tools-package-report.json) and
[installation evidence](evidence/tools-installed-package.json) describe this release.
Task and evidence status is maintained in [tasks.csv](tasks.csv) and [status.csv](status.csv).
Project decisions are recorded in [decisions.md](decisions.md).

## Responsibilities

| Layer | Namespace | Owns |
|---|---|---|
| Application | mdemat.application | Plasma and radiation models, physical scaling, coupling and scenarios |
| Specialized solvers | mdemat.solver | Elliptic, hyperbolic and parabolic formulations, boundaries, constraints and execution |
| Common algorithms | mdemat.common | ODE, linear/nonlinear algebra, geometry, quadrature, interpolation, approximation, symbolic, FEM and spectral methods |
| Supporting data | mdemat.data | Variables, function descriptions, geometry, mesh/topology and representation metadata |
| Tools | mdemat.tools | Small predicates, value configuration, primitive I/O, timing, numerical analysis and graphics resource management |

Data-owned plot methods are optional local display interfaces; they do not
attach graphics state to data or construct upper-layer numerical methods.
Tools consume primitive samples/weights/coordinates and do not import
State/FEM/mesh classes or interpret execution history. Domain sampling
belongs to common/application/execution. Native MATLAB error/warning/MException diagnostics replace the former
custom exception wrappers; IDs identify their semantic owner explicitly.

ODE owns generic stages, steps and history. PDE solvers and coupled applications
compose it. JFNK solves nonlinear equations; fitting/projection is distinct from
interpolation. Mesh descriptions own topology; geometric algorithms own tracing,
intersection and refinement. FDM and low-rank methods have bounded initial
consumers; they do not introduce universal wrappers or empty package trees.

~~~mermaid
flowchart TD
    A[Application / experiment] --> S[Specialized solver]
    A --> C[Common algorithms]
    S --> C
    C --> D[Supporting data]
    D --> T[Foundational tools]
    C --> T
    S --> D
    S --> T
    A --> D
    A --> O[Optional analysis / visualization / I/O]
    O --> T
~~~

Arrows are source dependencies. Common modules compose through small capabilities
without cycles. Supporting data never constructs a solver or numerical method.
Numerical analysis receives sampled arrays and physical weights. Visualizer
receives plain layout specifications and sampled datasets; custom callbacks draw
only in session-owned axes. It does not sample spaces or retain a State, timeline,
reference function or solver factory. Execution/common/application functions own
observation scheduling, FEM sampling, reference evaluation and physical reconstruction.

Config is a small value class. Scheme composes it and exposes a read-only numerical
snapshot; observation settings and Timer/Visualizer resources are separate. Changing
display components does not change numerical checkpoint identity or rebuild an ODE
workspace. The accepted-step history remains owned by the ODE/execution layer.

## Mathematical and execution boundaries

Separate problem data, discrete scheme and execution. Prescribed forcing and
initial/boundary functions are inputs; an independent solution reference is an
optional assessment input supplied by the experiment. Solving requires no plot,
file output or error study. Configuration loading is a service, not a superclass.

A mapping describes input/output dimensions and optional derivative capabilities.
A space describes admissible objects and finite representation where applicable.
An ansatz reconstructs from coefficients; approximation determines coefficients.
Operators expose actions, with matrix/adjoint/derivative capabilities when useful.
These distinctions do not require a universal inheritance hierarchy.

Keep generic volume/face/source/weighted-adjoint assembly in common FEM/operator
code. Equation-specific flux, compatibility, gauge and boundary policy belongs
with its solver. Physical macro-micro/Hermite reconstruction, closure, filtering,
limiting and scaling policy belongs with the relevant application.

## Representation contracts

| Quantity | Convention |
|---|---|
| Coordinates | dimension by points |
| Field values | components by points |
| Basis values | basis functions by points |
| DG coefficients | local DOFs times elements by components; local index fastest |
| Kinetic coefficients | Explicit spatial, velocity/mode and component axes and packing |
| Time series | Explicit observation axis; no implicit squeeze |
| DG face samples | Explicit element or sided-trace request |

Current [Function](../toolbox/+mdemat/+data/+function/Function.m) and
[FiniteElementSpace](../toolbox/+mdemat/+common/+fem/FiniteElementSpace.m) expose different
layouts. Convert at named, tested boundaries. Preserve all auxiliary-gradient
components, mode order, normalization, mesh maps, physical weights and units.
Dense/sparse MATLAB arrays remain the default; complex Fourier paths retain
declared conjugation semantics. Optional factors have explicit expansion costs.

Use value records for descriptions/snapshots and handles where identity is
needed for workspaces/resources. Published solutions must not alias mutable
stage buffers or caches. Cache keys identify mesh, basis/order/normalization,
coefficients, boundaries, time-step and split-phase changes.

Preserve d(Mu)/dt versus M du/dt, accepted/stage/multistep state, and SLDG's
two-time mass/transfer maps. Constant-mass MOL is not a substitute for those
contracts. Restart initially occurs at completed-step boundaries and includes
required method history; derivable caches are rebuilt.

| Record | Content |
|---|---|
| Solution | Fields, coefficients, representation descriptions, observation times and provenance |
| SolveReport | Execution status/reason, iterations, residual, warnings and timing scopes |
| ValidationReport | Criterion, measurements, reference identity, bounds and PASS/FAIL/UNVERIFIED |
| RunRecord | Resolved inputs, code/environment provenance, artifact identities and separate reports |

Successful execution does not establish accuracy. Missing assessment is
UNVERIFIED. Persist versioned descriptions and numeric data, with explicit
caller conversions, rather than opaque mutable handle graphs.

## Repository layout

~~~text
MdeMat.prj, buildfile.m, setupMdeMat.m
toolbox/
  +mdemat/
    +application/    plasma, radiation
    +solver/         elliptic, hyperbolic, parabolic
    +common/         ode, linalg, nonlinear, approximation, operator,
                     geometry, quadrature, interpolation, symbolic, fem, spectral,
                     fdm, lowrank
    +data/           variable, function, geometry, mesh, index, tensor
    +tools/          base, environment, analyze, io, visual
    +internal/       implementation details when justified
  resources/basis/   portable numeric runtime metadata
  doc/              native introductory help
  examples/         callable user demonstrations
tests/, benchmarks/, buildUtilities/, third_party/
docs/               developer material, never toolbox release content
outputs/, release/  generated artifacts
~~~

Only toolbox is the normal source path root. Resolve resources relative to the
package, not the caller's directory; do not recursively add the repository.
Optional Symbolic generation/export belongs in buildUtilities, with mathematical
implementation in common.symbolic and tested numeric exports for runtime use.

Data mesh classes now use `@Graph`, `@SeparableGraph` and `@PolytopalMesh`
folders with external `plot.m` methods; their fully qualified class names are
unchanged, and Triangulation inherits its plot method. Config/readConfig and the
reduced visual/analyze functions are current source APIs. Timer shares the
analyze package with error, convergence and performance measurements; empty
except and timing packages are removed. plotConvergence takes
explicit scale/error arrays; constructing layouts and choosing empirical or
prescribed reference orders is separate from collecting convergence observations. The removed Dataset,
Database, Strategy, FigureManager, TiledPlotter, Analyzer and Metric classes have
no permanent compatibility proxies in tools. See [tools-plan.md](tools-plan.md)
and the [file-level migration map](tools-file-review.csv).
