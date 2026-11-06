# MdeMat

MATLAB DG/spectral research code for elliptic, transport and diffusion problems,
radiation transport and plasma simulation.

The `mdemat` namespace in `toolbox/` contains applications, specialized solvers,
common algorithms, supporting data and tools. This working tree includes the
2.0 tools redesign: value configuration and samples, explicit observations,
data-owned graph plotting, numerical analysis functions and one graphics owner.

Tools are organized into `base`, `environment`, `analyze`, `io` and `visual`.
`mdemat.tools.analyze.Timer` shares the analysis package with error, convergence,
benchmark and memory measurements. Empty exception and timing packages are removed.

**Current acceptance: PASS.** Native numerical and symbolic tests, measured
production coverage, source analysis, packaging, fresh installation and project
checks passed for the 2.0.0 release. [Task status and evidence](docs/README.md),
[package report](docs/evidence/tools-package-report.json) and
[installation checks](docs/evidence/tools-installed-package.json) identify this snapshot.

From the checkout in MATLAB R2025b:

```matlab
setupMdeMat
runTests
runTests(Suite="symbolic") % explicitly requires Symbolic Math Toolbox
```

`buildtool` runs analysis and numerical/tests with measured coverage.
`buildtool package` enforces its prerequisites before staging the toolbox;
`buildtool verifyPackage` checks installation from a fresh external process.
Developer docs, reference PDFs, tests, benchmarks and build tools are excluded
from the installed artifact.

## Values, observations and graphics

Config updates return values and require assignment. Parsing is independent:

```matlab
cfg = mdemat.tools.environment.Config();
cfg = cfg.define("iterations", 20, @(n) isnumeric(n) && isreal(n) && ...
    isscalar(n) && isfinite(n) && n > 0 && n == fix(n));
cfg = cfg.set("iterations", 40);
resolved = cfg.toStruct();
% [patch, info] = mdemat.tools.environment.readConfig(file);
```

Given an existing scheme and matching state, configure observation separately
from its frozen numerical inputs:

```matlab
plotOptions = scheme.Observation.Plot;
plotOptions.PointsPerDimension = 2;
plotOptions.Components = struct('U', 1);
plotOptions.DatasetTypes = struct('NUMERIC', 'numeric');
plotOptions.Modes = {'1d'};
scheme.setObservation('plot', plotOptions);
state = scheme.run(state);
```

Field names must belong to that model. For example, Hermite uses `D/P/E`,
AP Hermite uses `D/Psi/Q`, and reconstructed `F` uses an explicit physical
sampling adapter. No plot or exact reference is required to solve.

For independent graphics, submit already sampled values:

```matlab
field = struct('Values', [0, 1, 4], 'ComponentIds', 2);
field.Axes = {[0, 1, 2]};
samples.RUN = struct('Type', "numeric", 'Data', struct('u', field));
panel = struct('Mode', "1d", 'Field', "u", 'ComponentId', 2, 'DatasetIds', "RUN");
layout = struct('Name', "profile", 'Size', [1, 1], 'Panels', panel);
session = mdemat.tools.visual.Visualizer(layout);
session.render(samples, time=0, title="Sampled field");
fig = session.getFigure("profile");
session.finalize();
```

`Values` is components-by-points; `ComponentIds` identifies its rows. Tensor
`Axes` follows `ndgrid` column-major order; point `Coordinates` preserves its
paired values and duplicate traces. Graphs and meshes instead expose
`plot(graph, axes=ax)` and `plot(mesh, axes=ax)` without storing graphics state.

The old Dataset/Database, Strategy, FigureManager/TiledPlotter and mesh Plotter
classes have been removed. Analyzer/Metric have become explicit sampling and numerical
functions. For a scalar convergence series, use
`[h, fit] = mdemat.tools.visual.plotConvergence(ax, scales, errors, fitted=true)`;
`referenceOrder=2` instead draws a reference normalized to the first valid error. See [graphics usage and lifetime](toolbox/doc/Graphics.m),
[architecture](docs/architecture.md) and [implementation/migration details](docs/tools-plan.md).
Persisted v1 records retain their own schema and are not renamed by the new source API.

[Mathematical formulations](docs/mathematics/README.md) map equations to code.
[Optimization review](docs/optimization.md) separates measured costs from future
ideas. Read [Contributing](CONTRIBUTING.md) before edits.

DistMesh v1.2 is bundled under [third_party/distmesh](third_party/distmesh/INTEGRATION.md)
with its MIT license and pinned upstream identity. The meshing example loads its
source for the duration of the call; the default toolbox path stays explicit.

**Contributor:** Yi Cai (yicaim@stu.xmu.edu.cn).
