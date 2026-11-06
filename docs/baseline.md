# Historical baseline and repair checklist

Recorded measurements and their scopes are in [status.csv](status.csv).
The original [baseline evidence](evidence/baseline.json) was executed against
main at be9c153, whose source tree matched the pinned v4 snapshot. It remains
unchanged historical evidence, including findings for files later removed.
It is not a fresh measurement of the cleaned workspace.

[Cleanup evidence](evidence/cleanup.json) records removed paths and original
Git/blob identities, plus the Gaussian shared-body comparison and native MATLAB
checks. The canonical example now calls
[gaussianSweep](../toolbox/examples/radiation/gaussianSweep.m); varyNv preserves
the original default and varyNu preserves the former desktop variant.
Both changed files passed their targeted analyzer error check. No full Gaussian
simulation or full unit-suite rerun is implied.

The recorded unit baseline is FAILED. T011 tracks subsequent repairs and
verification; this historical result does not describe the current suite. Incomplete overlaps Failed
and is not a separate population. Analyzer cleanliness does not prove runtime
binding, physical accuracy or convergence.

## Findings that established the required repairs

The table records original findings, not a list of defects still present.
Current task status and newer native evidence are linked from the development
entry; independent tests accompany the implemented corrections.

| Source/observation | Owner and check |
|---|---|
| AffineSpace inherits missing FunctionSpace | T010/T032: establish the intended representation and test it; no placeholder hierarchy |
| FiniteElementSpace and MacroMicroState call missing error helper | T010: exercise stable defined exception paths |
| MacroMicroState.copy names the wrong Vlasov type | T010/T060: correct type and verify independent copy/sharing |
| CompiledOrthogonalBasis validates normal_spherical_harmonic but switches differently | T010/T032: every family has matching selection and normalization checks |
| FiniteElementSpace.fit leaves X undefined for explicit options.X | T010/T035: explicit-point shape/order regression |
| Poisson LDG Dirichlet saves only W(:,2) for Q | T010/T040: preserve all multidimensional auxiliary components/signs |
| SeparableFunction cell derivative loops use undefined k | T010/T032: dense/sparse-pattern derivative regressions |
| MultiLevelIntegrator.copy omits the required rule | T010/T033: copy, rule and cache independence |
| MolScheme interval endpoints both read Timeline.Next | T010/T050: observation interpolation/history interval |
| ODE convergence example has identical conditional branches | T010/T050: intended study behavior and analyzer regression |
| Existing grid/config/scheme/dataset/plot test expectations fail | T011: separate code, intended contract and environment causes; retain valid requirements |
| ConfigParser evaluates text | T021/T070: validated explicit parsing with no arbitrary expression execution |
| Numeric resource loading depends on working directory | T032: source-relative portable exports and foreign-directory checks |
| RNG/path/cache state is not uniformly exception-safe | T021/T030/T033/T034: restoration and invalidation tests |
| Matrix hashing uses Java/undocumented serialization | T021/T030: supported identity/cache scope and numeric-shape/type tests |
| Fitting fallback uses normal equations/nonconjugating transpose | T030/T032: conditioning and declared complex contracts |
| Doxygen filter names a personal Windows path | T072: portable optional tooling; native help stays primary |
| Missing names only in See also text | T072: fix links/names without creating absent methods |
| Ganapol cylinder explicitly declares failure | T075: investigate against an independent reference; retain unresolved status until corrected |

WSL UNC change-notification/CMD warnings and Krylov tolerance warnings occurred
in the original execution. The real MATLAB project has since passed a native
Windows drive-path check in [native-project.json](evidence/native-project.json).
Installed optional products do not prove license availability or mandatory
dependencies. PDE/application studies, migrated coverage and final installed
verification have separate task gates and recorded scopes.
