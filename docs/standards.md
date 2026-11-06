# MATLAB toolbox standard

Applies to new and migrated code. Target design is in [architecture.md](architecture.md),
scope and sequencing in [plan.md](plan.md), and status in [README.md](README.md).
This is project policy adapted to MATLAB R2025b on Windows. Other releases and
platforms require separate evidence.

## Source and public APIs

* Use English source/help/errors/configuration, UTF-8, LF, four-space indentation
  and a target line length of 100 columns. Do not obscure equations to wrap them.
* Use lowercase semantic packages, UpperCamelCase classes/public properties,
  lowerCamelCase functions/methods/variables and UPPER_SNAKE_CASE constants.
  Conventional mathematical names remain valid. Filenames match public symbols.
* Document calls, shapes/axes, units/normalization, numeric types, sparsity,
  real/complex support, ownership, outside-domain behavior, errors and See also.
* Validate public inputs using arguments blocks where suitable. Validation may
  convert/reshape inputs; explicitly check exact type/axes when required.
  Avoid redundant validation in already-validated private hot loops.
* Use explicit name-value options and validated structs; reject unknown keys.
  Name-value keys use lowerCamelCase; value-record fields use UpperCamelCase.
  Configuration text must not execute arbitrary MATLAB expressions.
  Pure readConfig parses OPTION values without knowing a model schema.
  Config accepts only registered keys; an extension explicitly defines its
  value/validator pair. Scheme owners retain their validators and reject
  numeric input changes after construction. Derived workspace and observation
  options are separate from the immutable model inputs; see [tools-plan.md](tools-plan.md).
* Prefer small value records and limited capabilities to universal solver,
  Array/Tensor or configuration inheritance frameworks.

## Numerical contracts

State equation/formulation, spaces, mass/weight convention, flux orientation,
boundaries/gauge, requested fields and supported regimes. Follow the representation
contracts in the architecture. Keep all auxiliary components and normalization.
Distinguish interpolation, projection, general fitting and future rank truncation.

Prefer solves/factorizations to explicit inverses. Preserve convergence flags,
residuals, iteration counts and preconditioner identity. Reuse only valid caches;
declare failures and fallbacks. Preserve temporal/mass/transfer and splitting
semantics unless a separately assessed numerical change is intended.

Use independent continuous or demonstrably better references and physical
quadrature for assessment. Declare norm, refinement sequence, expected behavior
and acceptance bounds before evaluation. Separate algebraic, spatial, temporal,
quadrature and velocity errors; future rank/operator errors are separate too.
Near-zero reference norms and non-finite values need explicit policies.

## State, resources and side effects

Value/handle ownership, copies and cache invalidation are documented and tested.
Library calls do not clear the caller workspace, close unrelated figures, change
working directory, persist global paths/warnings, install dependencies or create
pools. Use fixtures/onCleanup for temporary resources and state restoration.
Accept a stream or preserve RNG state on success and exceptions.

Examples are callable or controlled demonstrations with explicit parameters,
seeds and output destinations. Source-relative resource loading must work from
a foreign directory. The repository current directory can otherwise hide path
dependencies, so installed-artifact checks use a fresh foreign-directory process.

Use stable project error identifiers and retain diagnostic context. Keep solution,
execution, assessment and run records distinct. Persist numeric data and versioned
schemas, not class-name-driven handle deserialization. Checkpoints include the
method history needed for supported restart boundaries. Output creation is
explicit, refuses overwrite by default, and reports partial failures.

## Verification

Use base MATLAB unit testing with explicit selection outside toolbox. Test
contracts, invalid inputs, ownership, independent identities and numerical
failure paths. Update tests/help with each capability; do not defer them until
delivery. Preserve existing requirements or document an explicit replacement.
Do not weaken tolerances or blanket-skip failures to make migration green.

Run codeIssues on new/migrated source: no errors, with warnings fixed or narrowly
justified. Keep vendor/baseline findings visible. No whole-project automatic fixes.
The project target is 80 percent migrated-production statement coverage with
justified exclusions; public/failure paths remain mandatory. Coverage is not an
accuracy certificate. Advanced MATLAB Test facilities are optional for core users.

Each numerical family has its own predeclared refinement/reference/invariant
checks. Optional dependencies may be skipped by a minimal suite, but their
dedicated capability job fails if the required environment is missing.

Measure startup, assembly, solve, assessment, rendering and I/O separately.
Performance claims require controlled warmups, repeated runs and stated hardware,
settings/seeds. Memory reports name the native process, units and scope, and
distinguish snapshots, sampled maxima and actual peaks. WSL launcher RSS is not
Windows MATLAB memory. Include instrumentation/operator/workspace costs; keep
hardware-sensitive timing thresholds out of ordinary unit tests.

## Help, build and distribution

Native help and GettingStarted are primary. Doxygen/Perl is optional development
tooling. Build commands are documented as usable only after implemented.
Create a real MATLAB project and explicit build tasks; do not hand-write fake
.prj or .mlx files or rely on newer-release-only packaging behavior.

Package only an explicit staged toolbox tree using R2025b project packaging or
ToolboxOptions/packageToolbox. Keep toolbox identity stable and versions numeric.
Review requiredFilesAndProducts plus dynamic/resource/optional dependencies;
installed products do not automatically become mandatory dependencies or prove
license availability.

T073 must explicitly exclude **all repository docs/**, including
**docs/references/*.pdf**, along with Git metadata, tests, benchmarks, build
scripts, old namespaces and generated output. Developer references are not
toolbox redistribution assets. Preserve project license and verify external-code
provenance/obligations before bundling any runtime vendor implementation.

Release acceptance requires required tests/studies, working help/resources,
dependencies and a fresh installed-artifact check without checkout paths. Local
packaging does not authorize publication. Record exact checks and limits.

## Official facilities

* [MathWorks toolbox recommendations](https://github.com/mathworks/toolboxdesign).
* [R2025b toolbox packaging](https://www.mathworks.com/help/releases/R2025b/matlab/matlab_prog/create-and-share-custom-matlab-toolboxes.html).
* [ToolboxOptions](https://www.mathworks.com/help/releases/R2025b/matlab/ref/matlab.addons.toolbox.toolboxoptions.html).
* [Argument validation](https://www.mathworks.com/help/releases/R2025b/matlab/matlab_prog/function-argument-validation-1.html).
* [Code Analyzer](https://www.mathworks.com/help/releases/R2025b/matlab/ref/codeissues.html).
* [Build tool](https://www.mathworks.com/help/releases/R2025b/matlab/ref/buildtool.html).
* [Official MATLAB setup action](https://github.com/matlab-actions/setup-matlab)
  and [command action](https://github.com/matlab-actions/run-command) used by the workflow.

Naming, coverage targets, architectural boundaries and scientific gates above
are MdeMat policies, not claims of universal MathWorks requirements.
