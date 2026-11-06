# Decision log

Accepted decisions, recorded 2026-11-05. Implementation progress belongs in
[tasks.csv](tasks.csv), not this log.

| ID | Decision | Consequence |
|---|---|---|
| D01 | Retain MATLAB and self-developed DG/spectral methods | MATLAB R2025b on Windows is the initial verified environment; do not adopt another project's runtime/language |
| D02 | Use application, specialized solver, common algorithm, supporting data and tool layers | ODE is common.ode; PDE families are elliptic/hyperbolic/parabolic; observers consume supplied views |
| D03 | Use v4 commit ee0d143074c4e426e96086082b142dd31c5309d2 only; implement on main | No older-version recovery work; preserve Git history and the v4 snapshot |
| D04 | Exclude every legacy directory, then remove excluded material under the cleanup instruction | Cleanup is recorded in evidence/cleanup.json; deleted paths are provenance records, not pending work |
| D05 | Consolidate Gaussian copies using named varyNv/varyNu sweeps | T012 is complete; numeric-body equivalence and native checks are recorded; toolbox migration remains separate |
| D06 | Publish only an explicitly staged toolbox payload | Runtime source/resources and reviewed vendor source/licenses are eligible; tests, research output, developer documentation and reference PDFs are excluded |
| D07 | Keep low rank/FDM as bounded follow-on implementations | Reserve useful operator/data contracts; avoid speculative frameworks and empty packages |
| D08 | Use one architecture, plan, task register and source ledger | Progress is shown at docs/README.md; statistics only in CSV; immutable machine evidence retains original facts |
| D09 | Complete all planned contents under the current implementation instruction | Bounded FDM and low-rank consumers are implemented alongside the complete v4 migration; their tests do not replace existing solver verification |
| D10 | Index-coordinate conversion belongs to data.index; FEM samplers belong to common.fem | Supporting mesh descriptions avoid a common-algorithm dependency; analysis/visual tools consume caller-supplied capabilities |
| D11 | Execution composes a value Config; readConfig is a pure parser | Supersedes the initial ConfigParser service: candidate updates validate before commit; model inputs freeze after construction; derived values and observation options have separate owners |
| D12 | Make volume and intrinsic face measures explicit | A collapsed ambient interval has zero volume; DG face quadrature uses a separate boundary API and retains intrinsic measure |
| D13 | Rendering styles are ordinary per-session values | Palette defaults are shared data; plotStyle produces an independent style value and curveStyle selects native graphics properties |
| D15 | Graph, SeparableGraph and PolytopalMesh plotting belongs to the corresponding data types | External plot methods replace tools Plotter classes without attaching graphics state to data; Triangulation inherits its mesh interface |
| D16 | Cut over repository source APIs in version 2.0.0 while preserving valid v1 numeric files | No permanent old-class proxies; explicit capability mappings record replacements. State reads State snapshots; execution reads checkpoint envelopes and projects v1 configuration identities |
| D17 | Numerical owners supply samples; tools own only generic values and resources | FEM owns quadrature/nested sampling, applications own physical F reconstruction, execution owns time pairing and assessment loops; weightedError and graphics never infer State or a spatial basis |
| D18 | Name sequence diagnostics by every selected component | Preserve the historical spatial-then-component two-norm; aggregate curves must not masquerade as component 1. Physical F uses dx dv on an explicitly stated truncated velocity interval |
| D19 | Consolidate analysis and timing under tools.analyze | Move Timer beside numerical/performance analysis functions; remove empty except/timing packages and migrate all live callers without old-namespace proxies |

The old source lists and review reports were merged into these owning records.
Source-defect findings are retained in [baseline.md](baseline.md); distinct low
rank design is retained in [low-rank.md](low-rank.md). A verified backup of the
uncommitted documents was made before consolidation; no Git history was rewritten.

## D14 — Bundle DistMesh

The owner requested a vendored DistMesh on 2026-11-05. Use the official pure
MATLAB v1.2 snapshot pinned by `third_party/distmesh/provenance.json`, retaining
all upstream bytes and the MIT license. Package `src/*.m` and license/provenance
records; retain upstream examples/tests/documentation as checkout material.
The geometry study scopes its vendor path, seed and progress figure, preserving
caller state. Default toolbox paths remain unchanged. Native circle generation,
refinement, quality and installed-artifact checks cover the integration; other
upstream geometries and arbitrary convergence are not certified.
