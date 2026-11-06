# Native MATLAB development checks

Run from a checkout on MATLAB R2025b for Windows:

```matlab
buildtool analyzer       % all toolbox .m files; findings saved as JSON
buildtool test           % migrated + unit + regression + numerical; coverage gate
buildtool numerical      % numerical studies and regression/solvers
buildtool symbolic       % requires Symbolic Math Toolbox; never blanket-skips
buildtool package        % required prerequisites and package creation
% Exit MATLAB, then start a new batch process for the installed-artifact check:
buildtool verifyPackage
```

The test task uses the documented CodeCoveragePlugin at statement metric level
and writes Cobertura XML. Its measured aggregate must reach the project 80%
target. Third-party Lebedev source and optional symbolic source are explicitly
excluded and listed in the result; all migrated v4 test classes are collected.
A measured shortfall fails the task and retains evidence. Timing thresholds are
not unit-test requirements. Analyzer errors fail; every warning remains visible
and requires review before release acceptance.

`packageToolboxArtifact` builds a native temporary staging tree from an explicit
allowlist. Its stable identity/version are in `toolbox/resources/toolbox.json`; its
manifest records all payload files and static dependency-analysis results.
After tests load external class definitions, MATLAB may report their checkout
paths during dependency analysis. The packager accepts such a path only when
it names an explicit payload member and its SHA-256 matches the staged bytes.
The manifest records these relative aliases; unknown or changed external files
still fail. The caller's class cache is not cleared.
Repository `docs/` (including reference PDFs), tests, benchmark/build scripts,
Git data, old namespaces and generated outputs never enter the tree. Native
introductory help lives under `toolbox/doc/`. Static product analysis describes
dependencies it can resolve; installed products do not prove license access.
Symbolic input/generation is an optional capability, not a requirement to load
the portable numeric basis resources.

The DistMesh v1.2 snapshot is pinned in `third_party/distmesh/provenance.json`.
Its unmodified `src/*.m` files and MIT license/notices are packaged explicitly;
upstream documentation, tests, examples and startup remain checkout material.
Vendor source is analyzed with vendor findings visible and excluded from MdeMat
production coverage. The geometry example temporarily adds only its `src` path.

`verifyPackageArtifact` must be invoked from a new external MATLAB batch process
after packaging exits; it refuses cached MdeMat production code. It enters a
native temporary folder containing only the `.mltbx` and standalone verifier.
It uses documented
temporary add-on installation settings, restores settings, uninstalls its test
installation, and checks native help, executed GettingStarted, resources,
Lebedev area, an independent quadratic FDM solution, pure Poisson potential and
gradient, numeric state/space persistence, all evolution/application families,
and actual bundled DistMesh meshing/refinement. It refuses to replace a pre-existing
MdeMat installation. Generated reports/artifacts live in `outputs/` and
`release/`; this workflow does not publish them.

The GitHub workflow uses official MATLAB actions and Windows runners. Private
repositories must configure the `MLM_LICENSE_TOKEN` secret for batch licensing,
as documented by the [official setup action](https://github.com/matlab-actions/setup-matlab).
Hosted
execution and license availability still require an actual CI run; local native
evidence must not be described as a completed hosted workflow.
